use std::sync::Arc;

use arrow::datatypes::{DataType, Field};
use datafusion::common::Result;
use datafusion::config::ConfigOptions;
use datafusion::physical_optimizer::PhysicalOptimizerRule;
use datafusion::physical_plan::ExecutionPlan;

use crate::cost::{CostModel, JoinStrategy, SmallSide};
use crate::exec::sweep_line::BuildSide;
use crate::exec::SweepLineJoinExec;
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
        log::debug!(
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
                log::debug!(
                    "INTERSECTS optimizer: no Parquet stats available, \
                     deferring to DataFusion"
                );
                return Ok(original_plan);
            }
        };

        log::debug!("INTERSECTS optimizer: selected {strategy:?}");

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
            JoinStrategy::BinnedJoin { .. } => {
                // For binned joins, DataFusion's default hash join
                // on chrom with the range filter is already well
                // optimized. The physical plan rewrite overhead
                // (BinExpandExec + modified HashJoinExec) exceeds
                // the gain from bin-based hashing. Defer to
                // DataFusion's built-in join.
                log::debug!(
                    "INTERSECTS optimizer: binned strategy selected, \
                     deferring to DataFusion"
                );
                Ok(original_plan)
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
