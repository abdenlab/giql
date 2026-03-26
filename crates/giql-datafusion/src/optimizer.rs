use std::sync::Arc;

use datafusion::common::Result;
use datafusion::config::ConfigOptions;
use datafusion::physical_optimizer::PhysicalOptimizerRule;
use datafusion::physical_plan::ExecutionPlan;

use crate::cost::{CostModel, JoinStrategy};
use crate::exec::{BinnedJoinExec, SweepLineJoinExec};
use crate::pattern::{detect_interval_join, IntervalJoinMatch};
use crate::stats;
use crate::IntersectsOptimizerConfig;

/// Physical optimizer rule that detects interval overlap joins and
/// replaces them with optimized execution plans.
///
/// The rule reads Parquet metadata and performs lightweight sampling to
/// choose between sweep-line and binned equi-join algorithms based on
/// the interval width distribution.
#[derive(Debug)]
pub struct IntersectsOptimizerRule {
    config: IntersectsOptimizerConfig,
}

impl IntersectsOptimizerRule {
    pub fn new(config: IntersectsOptimizerConfig) -> Self {
        Self { config }
    }

    /// Recursively optimize a plan tree, replacing interval overlap
    /// joins with custom execution plans.
    fn optimize_plan(
        &self,
        plan: Arc<dyn ExecutionPlan>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        // First, try to match this node
        if let Some(join_match) = detect_interval_join(&plan)? {
            return self.replace_join(join_match);
        }

        // Recurse into children
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

    /// Replace an interval overlap join with an optimized execution
    /// plan based on cost model analysis.
    fn replace_join(
        &self,
        join_match: IntervalJoinMatch,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        let cost_model = CostModel::new(&self.config);

        // Collect statistics from Parquet sources
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

        // Decide on strategy
        let strategy = match (&left_stats, &right_stats) {
            (Some(left), Some(right)) => cost_model.decide(left, right),
            _ => {
                log::info!(
                    "INTERSECTS optimizer: no Parquet stats available, \
                     using default sweep-line"
                );
                JoinStrategy::SweepLine { skip_sort: false }
            }
        };

        log::info!("INTERSECTS optimizer: selected {strategy:?}");

        match strategy {
            JoinStrategy::NestedLoop => {
                // Return original plan unchanged — reconstruct from
                // the match components. This shouldn't normally happen
                // since we default to SweepLine, but handle it.
                Err(datafusion::error::DataFusionError::Internal(
                    "NestedLoop strategy should not be returned by \
                     cost model when Parquet stats are unavailable"
                        .to_string(),
                ))
            }
            JoinStrategy::SweepLine { skip_sort } => {
                Ok(Arc::new(SweepLineJoinExec::new(
                    join_match.left,
                    join_match.right,
                    join_match.left_cols,
                    join_match.right_cols,
                    join_match.output_schema,
                    skip_sort,
                )))
            }
            JoinStrategy::BinnedJoin { bin_size } => {
                Ok(Arc::new(BinnedJoinExec::new(
                    join_match.left,
                    join_match.right,
                    join_match.left_cols,
                    join_match.right_cols,
                    join_match.output_schema,
                    bin_size,
                )))
            }
        }
    }

    /// Collect statistics from the first available Parquet file.
    fn collect_stats(
        &self,
        paths: &[std::path::PathBuf],
        start_col: &str,
        end_col: &str,
    ) -> Option<stats::IntervalStats> {
        if paths.is_empty() {
            return None;
        }

        // Use the first file for statistics. For multi-file tables,
        // a more sophisticated approach would aggregate stats across
        // files.
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
