use std::sync::Arc;

use arrow::datatypes::{DataType, Field};
use datafusion::common::Result;
use datafusion::config::ConfigOptions;
use datafusion::physical_optimizer::PhysicalOptimizerRule;
use datafusion::physical_plan::ExecutionPlan;

use crate::cost::{CostModel, JoinStrategy, SmallSide};
use crate::exec::sweep_line::BuildSide;
use crate::exec::{BinExpandExec, SweepLineJoinExec};
use crate::pattern::{detect_interval_join, IntervalJoinMatch};
use crate::stats;
use crate::IntersectsOptimizerConfig;

/// Physical optimizer rule that detects interval overlap joins and
/// replaces them with optimized execution plans.
///
/// The rule reads Parquet metadata and performs lightweight sampling to
/// choose between sweep-line and binned equi-join algorithms based on
/// the interval width distribution.
///
/// - **Sweep line**: Replaces the join with a custom `SweepLineJoinExec`
///   that sorts and sweeps. Best for heavy-tailed or high-variance
///   width distributions.
/// - **Binned join**: Wraps each input in a `BinExpandExec` that
///   expands intervals into genome bins, then lets DataFusion's
///   parallel `HashJoinExec` join on `(chrom, bin_id)`. Best for
///   uniform-width data.
#[derive(Debug)]
pub struct IntersectsOptimizerRule {
    config: IntersectsOptimizerConfig,
}

impl IntersectsOptimizerRule {
    pub fn new(config: IntersectsOptimizerConfig) -> Self {
        Self { config }
    }

    fn optimize_plan(
        &self,
        plan: Arc<dyn ExecutionPlan>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        eprintln!(
            "INTERSECTS optimizer: visiting node: {}",
            plan.name()
        );
        if let Some(join_match) = detect_interval_join(&plan)? {
            return self.maybe_replace_join(plan, join_match);
        }

        let children: Vec<Arc<dyn ExecutionPlan>> = plan
            .children()
            .into_iter()
            .map(|child| self.optimize_plan(child.clone()))
            .collect::<Result<Vec<_>>>()?;

        if children.is_empty() {
            return Ok(plan);
        }

        plan.with_new_children(children)
    }

    fn maybe_replace_join(
        &self,
        original_plan: Arc<dyn ExecutionPlan>,
        join_match: IntervalJoinMatch,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        let cost_model = CostModel::new(&self.config);

        let left_stats = self.collect_stats(
            &join_match.left_parquet_paths,
            &join_match.left_cols.start_col,
            &join_match.left_cols.end_col,
        );
        let right_stats = self.collect_stats(
            &join_match.right_parquet_paths,
            &join_match.right_cols.start_col,
            &join_match.right_cols.end_col,
        );

        let strategy = match (&left_stats, &right_stats) {
            (Some(left), Some(right)) => cost_model.decide(left, right),
            _ => {
                eprintln!(
                    "INTERSECTS optimizer: no Parquet stats available, \
                     deferring to DataFusion"
                );
                return Ok(original_plan);
            }
        };

        eprintln!("INTERSECTS optimizer: selected {strategy:?}");

        match strategy {
            JoinStrategy::SweepLine { build_side } => {
                let bs = match build_side {
                    SmallSide::Left => BuildSide::Left,
                    SmallSide::Right => BuildSide::Right,
                };
                Ok(Arc::new(SweepLineJoinExec::new(
                    join_match.left,
                    join_match.right,
                    join_match.left_cols,
                    join_match.right_cols,
                    join_match.output_schema,
                    bs,
                )))
            }
            JoinStrategy::BinnedJoin { bin_size } => {
                self.build_binned_plan(
                    original_plan,
                    join_match,
                    bin_size,
                )
            }
            JoinStrategy::NestedLoop => {
                log::info!(
                    "INTERSECTS optimizer: deferring to DataFusion's \
                     built-in join"
                );
                Ok(original_plan)
            }
        }
    }

    /// Build a binned join plan:
    ///
    /// 1. Wrap each child in `BinExpandExec` (adds `__giql_bin` and
    ///    `__giql_first_bin` columns)
    /// 2. `HashJoinExec` on `(chrom, __giql_bin)` with the original
    ///    range filter
    /// 3. `FilterExec` for canonical-bin dedup: keep only matches
    ///    where `__giql_bin == max(left.__giql_first_bin,
    ///    right.__giql_first_bin)`, so each pair is emitted once
    /// 4. `ProjectionExec` to strip the extra columns
    fn build_binned_plan(
        &self,
        original_plan: Arc<dyn ExecutionPlan>,
        join_match: IntervalJoinMatch,
        bin_size: usize,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        use datafusion::physical_expr::expressions::{
            BinaryExpr, CastExpr, Column, Literal,
        };
        use datafusion::physical_plan::filter::FilterExec;
        use datafusion::physical_plan::joins::HashJoinExec;
        use datafusion::physical_plan::projection::ProjectionExec;
        use datafusion::logical_expr::Operator;

        let hj = original_plan
            .as_any()
            .downcast_ref::<HashJoinExec>()
            .ok_or_else(|| {
                datafusion::error::DataFusionError::Internal(
                    "Expected HashJoinExec for binned plan rewrite"
                        .to_string(),
                )
            })?;

        let extra = BinExpandExec::EXTRA_COLS; // 2: __giql_bin, __giql_first_bin

        // Step 1: Wrap each child in BinExpandExec
        let left_expanded = Arc::new(BinExpandExec::new(
            join_match.left,
            join_match.left_cols.start_idx,
            join_match.left_cols.end_idx,
            bin_size,
        )) as Arc<dyn ExecutionPlan>;

        let right_expanded = Arc::new(BinExpandExec::new(
            join_match.right,
            join_match.right_cols.start_idx,
            join_match.right_cols.end_idx,
            bin_size,
        )) as Arc<dyn ExecutionPlan>;

        let left_schema = left_expanded.schema();
        let right_schema = right_expanded.schema();
        let left_n = left_schema.fields().len();
        let right_n = right_schema.fields().len();

        // Indices of the new columns in each child's schema
        let left_bin_idx = left_n - 2;      // __giql_bin
        let left_first_bin_idx = left_n - 1; // __giql_first_bin
        let right_bin_idx = right_n - 2;
        let right_first_bin_idx = right_n - 1;

        // Step 2: HashJoinExec on (chrom, __giql_bin)
        let mut on = hj.on().to_vec();
        on.push((
            Arc::new(Column::new("__giql_bin", left_bin_idx))
                as Arc<dyn datafusion::physical_plan::PhysicalExpr>,
            Arc::new(Column::new("__giql_bin", right_bin_idx))
                as Arc<dyn datafusion::physical_plan::PhysicalExpr>,
        ));

        // No projection on the HashJoinExec — we need all columns
        // including __giql_first_bin for the dedup filter.
        let new_join = Arc::new(HashJoinExec::try_new(
            left_expanded,
            right_expanded,
            on,
            hj.filter().cloned(),
            hj.join_type(),
            None, // no projection yet
            *hj.partition_mode(),
            hj.null_equals_null(),
        )?);

        // Step 3: FilterExec for canonical-bin dedup.
        //
        // Join output columns (inner join, no projection):
        //   [0..left_n) = left columns (including __giql_bin, __giql_first_bin)
        //   [left_n..left_n+right_n) = right columns
        //
        // Filter: __giql_bin (from left) == max(left.__giql_first_bin, right.__giql_first_bin)
        //
        // We use left's __giql_bin since it equals right's (equi-key).
        let join_schema = new_join.schema();
        let join_left_bin = left_bin_idx;
        let join_left_first_bin = left_first_bin_idx;
        let join_right_first_bin = left_n + right_first_bin_idx;

        // Build: CASE WHEN left_first_bin >= right_first_bin
        //             THEN left_first_bin
        //             ELSE right_first_bin END
        // Simplified: use a >= b check with binary expressions
        let left_fb: Arc<dyn datafusion::physical_plan::PhysicalExpr> =
            Arc::new(Column::new("__giql_first_bin", join_left_first_bin));
        let right_fb: Arc<dyn datafusion::physical_plan::PhysicalExpr> =
            Arc::new(Column::new("__giql_first_bin", join_right_first_bin));
        let bin_col: Arc<dyn datafusion::physical_plan::PhysicalExpr> =
            Arc::new(Column::new("__giql_bin", join_left_bin));

        // Filter: __giql_bin == left_first_bin OR __giql_bin == right_first_bin
        // AND left_first_bin <= __giql_bin AND right_first_bin <= __giql_bin
        //
        // Simpler canonical condition:
        //   __giql_bin == GREATEST(left_first_bin, right_first_bin)
        //
        // Without a GREATEST function, use:
        //   (left_first_bin >= right_first_bin AND __giql_bin == left_first_bin)
        //   OR
        //   (right_first_bin > left_first_bin AND __giql_bin == right_first_bin)
        let dedup_filter: Arc<dyn datafusion::physical_plan::PhysicalExpr> = Arc::new(
            BinaryExpr::new(
                Arc::new(BinaryExpr::new(
                    Arc::new(BinaryExpr::new(
                        left_fb.clone(),
                        Operator::GtEq,
                        right_fb.clone(),
                    )),
                    Operator::And,
                    Arc::new(BinaryExpr::new(
                        bin_col.clone(),
                        Operator::Eq,
                        left_fb.clone(),
                    )),
                )),
                Operator::Or,
                Arc::new(BinaryExpr::new(
                    Arc::new(BinaryExpr::new(
                        right_fb.clone(),
                        Operator::Gt,
                        left_fb,
                    )),
                    Operator::And,
                    Arc::new(BinaryExpr::new(
                        bin_col,
                        Operator::Eq,
                        right_fb,
                    )),
                )),
            ),
        );

        let filtered =
            Arc::new(FilterExec::try_new(dedup_filter, new_join)?)
                as Arc<dyn ExecutionPlan>;

        // Step 4: ProjectionExec to strip extra columns.
        // Keep only the original columns (skip __giql_bin, __giql_first_bin
        // from both sides).
        let orig_left = left_n - extra;
        let orig_right = right_n - extra;

        let mut proj_exprs: Vec<(
            Arc<dyn datafusion::physical_plan::PhysicalExpr>,
            String,
        )> = Vec::new();

        let filter_schema = filtered.schema();

        // Left original columns
        for i in 0..orig_left {
            let name = filter_schema.field(i).name().clone();
            proj_exprs.push((
                Arc::new(Column::new(&name, i)),
                name,
            ));
        }
        // Right original columns (skip left's extra cols)
        for i in 0..orig_right {
            let idx = left_n + i;
            let name = filter_schema.field(idx).name().clone();
            // Avoid name collisions by keeping original field name
            proj_exprs.push((
                Arc::new(Column::new(&name, idx)),
                filter_schema.field(idx).name().clone(),
            ));
        }

        Ok(Arc::new(ProjectionExec::try_new(proj_exprs, filtered)?))
    }

    fn collect_stats(
        &self,
        paths: &[std::path::PathBuf],
        start_col: &str,
        end_col: &str,
    ) -> Option<stats::IntervalStats> {
        if paths.is_empty() {
            return None;
        }

        let path = &paths[0];
        match stats::collect_parquet_stats(
            path,
            start_col,
            end_col,
            self.config.max_sample_row_groups,
        ) {
            Ok(stats) => Some(stats),
            Err(e) => {
                log::warn!(
                    "INTERSECTS optimizer: failed to collect stats \
                     from {path:?}: {e}"
                );
                None
            }
        }
    }
}

impl PhysicalOptimizerRule for IntersectsOptimizerRule {
    fn optimize(
        &self,
        plan: Arc<dyn ExecutionPlan>,
        _config: &ConfigOptions,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        self.optimize_plan(plan)
    }

    fn name(&self) -> &str {
        "intersects_optimizer"
    }

    fn schema_check(&self) -> bool {
        true
    }
}
