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
    /// 2. `HashJoinExec` on `(chrom, __giql_bin)` with a combined
    ///    filter that includes both the range overlap check AND the
    ///    canonical-bin dedup, plus a projection that strips extra
    ///    columns
    ///
    /// The dedup is folded into the JoinFilter so it runs inside the
    /// hash probe — no separate FilterExec or ProjectionExec needed.
    fn build_binned_plan(
        &self,
        original_plan: Arc<dyn ExecutionPlan>,
        join_match: IntervalJoinMatch,
        bin_size: usize,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        use datafusion::common::JoinSide;
        use datafusion::logical_expr::Operator;
        use datafusion::physical_expr::expressions::{
            BinaryExpr, Column,
        };
        use datafusion::physical_plan::joins::utils::{
            ColumnIndex, JoinFilter,
        };
        use datafusion::physical_plan::joins::HashJoinExec;

        let hj = original_plan
            .as_any()
            .downcast_ref::<HashJoinExec>()
            .ok_or_else(|| {
                datafusion::error::DataFusionError::Internal(
                    "Expected HashJoinExec for binned plan rewrite"
                        .to_string(),
                )
            })?;

        let extra = BinExpandExec::EXTRA_COLS;

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

        let left_bin_idx = left_n - 2;
        let left_first_bin_idx = left_n - 1;
        let right_bin_idx = right_n - 2;
        let right_first_bin_idx = right_n - 1;

        // Equi-keys: original chrom + __giql_bin
        let mut on = hj.on().to_vec();
        on.push((
            Arc::new(Column::new("__giql_bin", left_bin_idx))
                as Arc<dyn datafusion::physical_plan::PhysicalExpr>,
            Arc::new(Column::new("__giql_bin", right_bin_idx))
                as Arc<dyn datafusion::physical_plan::PhysicalExpr>,
        ));

        // Step 2: Build extended JoinFilter with dedup folded in.
        //
        // Start from the original filter and append:
        //   - column indices for __giql_bin (left) and
        //     __giql_first_bin (left + right)
        //   - the canonical-bin dedup expression ANDed with the
        //     original range-overlap expression
        let extended_filter = if let Some(orig_filter) = hj.filter()
        {
            let mut col_indices =
                orig_filter.column_indices().to_vec();
            let orig_len = col_indices.len();

            // Append 3 new column references into the filter schema
            // [orig_len + 0] → left.__giql_bin
            col_indices.push(ColumnIndex {
                index: left_bin_idx,
                side: JoinSide::Left,
            });
            // [orig_len + 1] → left.__giql_first_bin
            col_indices.push(ColumnIndex {
                index: left_first_bin_idx,
                side: JoinSide::Left,
            });
            // [orig_len + 2] → right.__giql_first_bin
            col_indices.push(ColumnIndex {
                index: right_first_bin_idx,
                side: JoinSide::Right,
            });

            // Build filter-local column refs
            let filt_bin: Arc<
                dyn datafusion::physical_plan::PhysicalExpr,
            > = Arc::new(Column::new(
                "__giql_bin",
                orig_len,
            ));
            let filt_lfb: Arc<
                dyn datafusion::physical_plan::PhysicalExpr,
            > = Arc::new(Column::new(
                "__giql_first_bin",
                orig_len + 1,
            ));
            let filt_rfb: Arc<
                dyn datafusion::physical_plan::PhysicalExpr,
            > = Arc::new(Column::new(
                "__giql_first_bin",
                orig_len + 2,
            ));

            // Canonical-bin condition:
            //   (lfb >= rfb AND bin == lfb)
            //   OR (rfb > lfb AND bin == rfb)
            let dedup_expr: Arc<
                dyn datafusion::physical_plan::PhysicalExpr,
            > = Arc::new(BinaryExpr::new(
                Arc::new(BinaryExpr::new(
                    Arc::new(BinaryExpr::new(
                        filt_lfb.clone(),
                        Operator::GtEq,
                        filt_rfb.clone(),
                    )),
                    Operator::And,
                    Arc::new(BinaryExpr::new(
                        filt_bin.clone(),
                        Operator::Eq,
                        filt_lfb,
                    )),
                )),
                Operator::Or,
                Arc::new(BinaryExpr::new(
                    Arc::new(BinaryExpr::new(
                        filt_rfb.clone(),
                        Operator::Gt,
                        Arc::new(Column::new(
                            "__giql_first_bin",
                            orig_len + 1,
                        )),
                    )),
                    Operator::And,
                    Arc::new(BinaryExpr::new(
                        filt_bin,
                        Operator::Eq,
                        filt_rfb,
                    )),
                )),
            ));

            // Combine: original_expr AND dedup_expr
            let combined: Arc<
                dyn datafusion::physical_plan::PhysicalExpr,
            > = Arc::new(BinaryExpr::new(
                orig_filter.expression().clone(),
                Operator::And,
                dedup_expr,
            ));

            // Build extended filter schema: original fields + 3 new
            let mut filter_fields: Vec<Arc<arrow::datatypes::Field>> =
                orig_filter
                    .schema()
                    .fields()
                    .iter()
                    .cloned()
                    .collect();
            filter_fields.push(Arc::new(
                arrow::datatypes::Field::new(
                    "__giql_bin",
                    arrow::datatypes::DataType::Int64,
                    false,
                ),
            ));
            filter_fields.push(Arc::new(
                arrow::datatypes::Field::new(
                    "__giql_first_bin_l",
                    arrow::datatypes::DataType::Int64,
                    false,
                ),
            ));
            filter_fields.push(Arc::new(
                arrow::datatypes::Field::new(
                    "__giql_first_bin_r",
                    arrow::datatypes::DataType::Int64,
                    false,
                ),
            ));
            let filter_schema = Arc::new(
                arrow::datatypes::Schema::new(filter_fields),
            );

            Some(JoinFilter::new(combined, col_indices, filter_schema))
        } else {
            None
        };

        // Projection: keep only original columns from both sides,
        // strip __giql_bin and __giql_first_bin.
        let orig_left = left_n - extra;
        let orig_right = right_n - extra;
        let mut projection: Vec<usize> =
            (0..orig_left).collect();
        projection.extend(left_n..left_n + orig_right);

        let new_join = HashJoinExec::try_new(
            left_expanded,
            right_expanded,
            on,
            extended_filter,
            hj.join_type(),
            Some(projection),
            *hj.partition_mode(),
            hj.null_equality(),
            false,
        )?;

        Ok(Arc::new(new_join))
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
