//! DataFusion optimizer for genomic interval (INTERSECTS) joins.
//!
//! This crate provides a [`PhysicalOptimizerRule`] that reads Parquet
//! metadata and lightweight sampling to choose between sweep-line and
//! binned equi-join algorithms for interval overlap joins.
//!
//! # Usage
//!
//! ```rust,no_run
//! use datafusion::execution::SessionStateBuilder;
//! use datafusion::prelude::*;
//! use giql_datafusion::{IntersectsOptimizerConfig, register_optimizer};
//!
//! let config = IntersectsOptimizerConfig::default();
//! let state = SessionStateBuilder::new()
//!     .with_default_features()
//!     .build();
//! let state = register_optimizer(state, config);
//! let ctx = SessionContext::from(state);
//! ```

pub mod cost;
pub mod exec;
pub mod logical_rule;
pub mod optimizer;
pub mod pattern;
pub mod pruning;
pub mod stats;

pub use cost::JoinStrategy;
pub use logical_rule::IntersectsLogicalRule;
pub use optimizer::IntersectsOptimizerRule;

use datafusion::execution::SessionState;
use datafusion::optimizer::OptimizerRule;
use datafusion::physical_optimizer::PhysicalOptimizerRule;
use std::sync::Arc;

/// Configuration for the INTERSECTS join optimizer.
#[derive(Debug, Clone)]
pub struct IntersectsOptimizerConfig {
    /// Threshold for p99/median width ratio. Above this, sweep line is
    /// chosen to avoid binning replication blowup on wide intervals.
    pub p99_median_threshold: f64,

    /// Threshold for coefficient of variation. Above this, sweep line
    /// is chosen because no single bin size works well.
    pub cv_threshold: f64,

    /// Maximum number of row groups to sample for width distribution.
    pub max_sample_row_groups: usize,

    /// Enable the experimental logical optimizer rule that rewrites
    /// interval overlap joins to UNNEST-based binned equi-joins.
    /// When false (default), only the physical sweep-line optimizer
    /// is active.
    pub enable_logical_rule: bool,
}

impl Default for IntersectsOptimizerConfig {
    fn default() -> Self {
        Self {
            p99_median_threshold: 10.0,
            cv_threshold: 1.5,
            max_sample_row_groups: 3,
            enable_logical_rule: false,
        }
    }
}

/// Build a [`SessionState`] with the INTERSECTS optimizer rules.
///
/// The physical rule detects interval overlap joins and replaces them
/// with sweep-line execution plans for heavy-tailed distributions,
/// deferring to DataFusion's default join for uniform data.
///
/// The logical rule (experimental, disabled by default) rewrites
/// interval overlap joins to UNNEST-based binned equi-joins at the
/// logical level, enabling DataFusion's native parallel execution.
/// Enable by setting `enable_logical_rule = true` in the config.
pub fn register_optimizer(
    state: SessionState,
    config: IntersectsOptimizerConfig,
) -> SessionState {
    use datafusion::execution::SessionStateBuilder;

    // Physical rule: sweep-line for heavy-tailed distributions
    let physical_rule: Arc<dyn PhysicalOptimizerRule + Send + Sync> =
        Arc::new(IntersectsOptimizerRule::new(config.clone()));

    let mut physical_rules: Vec<
        Arc<dyn PhysicalOptimizerRule + Send + Sync>,
    > = state.physical_optimizers().to_vec();
    physical_rules.push(physical_rule);

    let builder = if config.enable_logical_rule {
        let logical_rule: Arc<dyn OptimizerRule + Send + Sync> =
            Arc::new(IntersectsLogicalRule::new(config));

        let mut logical_rules: Vec<
            Arc<dyn OptimizerRule + Send + Sync>,
        > = state.optimizers().to_vec();
        logical_rules.push(logical_rule);

        SessionStateBuilder::new_from_existing(state)
            .with_optimizer_rules(logical_rules)
            .with_physical_optimizer_rules(physical_rules)
    } else {
        SessionStateBuilder::new_from_existing(state)
            .with_physical_optimizer_rules(physical_rules)
    };

    builder.build()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = IntersectsOptimizerConfig::default();
        assert!((config.p99_median_threshold - 10.0).abs() < f64::EPSILON);
        assert!((config.cv_threshold - 1.5).abs() < f64::EPSILON);
        assert_eq!(config.max_sample_row_groups, 3);
    }

    #[test]
    fn test_custom_config_used_by_cost_model() {
        let config = IntersectsOptimizerConfig {
            p99_median_threshold: 5.0,
            cv_threshold: 1.0,
            max_sample_row_groups: 1,
            enable_logical_rule: false,
        };
        let model = cost::CostModel::new(&config);

        // With p99/median = 6.0 > 5.0 (custom threshold), should
        // short-circuit to sweep line even though default threshold
        // would not trigger.
        let stats = stats::IntervalStats {
            row_count: 100_000,
            domain_min: 0,
            domain_max: 1_000_000,
            is_sorted_by_start: false,
            row_group_bounds: vec![],
            width: stats::WidthStats {
                median: 100.0,
                mean: 120.0,
                p95: 500.0,
                p99: 600.0,
                cv: 0.5,
                p99_median_ratio: 6.0,
            },
        };

        match model.decide(&stats, &stats) {
            JoinStrategy::SweepLine { .. } => {}
            other => panic!(
                "Expected SweepLine with custom threshold, got {:?}",
                other
            ),
        }
    }

    #[test]
    fn test_register_optimizer_adds_rule() {
        use datafusion::execution::SessionStateBuilder;

        let state = SessionStateBuilder::new()
            .with_default_features()
            .build();
        let n_before = state.physical_optimizers().len();

        let config = IntersectsOptimizerConfig::default();
        let state = register_optimizer(state, config);
        let n_after = state.physical_optimizers().len();

        assert_eq!(n_after, n_before + 1);

        let last_rule = state.physical_optimizers().last().unwrap();
        assert_eq!(last_rule.name(), "intersects_optimizer");
    }
}
