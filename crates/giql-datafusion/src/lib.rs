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
pub mod optimizer;
pub mod pattern;
pub mod pruning;
pub mod stats;

pub use cost::JoinStrategy;
pub use optimizer::IntersectsOptimizerRule;

use datafusion::execution::SessionState;
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
}

impl Default for IntersectsOptimizerConfig {
    fn default() -> Self {
        Self {
            p99_median_threshold: 10.0,
            cv_threshold: 1.5,
            max_sample_row_groups: 3,
        }
    }
}

/// Build a [`SessionState`] with the INTERSECTS optimizer rule
/// appended to the default physical optimizer rules.
pub fn register_optimizer(
    state: SessionState,
    config: IntersectsOptimizerConfig,
) -> SessionState {
    use datafusion::execution::SessionStateBuilder;

    let rule: Arc<dyn PhysicalOptimizerRule + Send + Sync> =
        Arc::new(IntersectsOptimizerRule::new(config));

    let mut rules: Vec<Arc<dyn PhysicalOptimizerRule + Send + Sync>> =
        state.physical_optimizers().to_vec();
    rules.push(rule);

    SessionStateBuilder::new_from_existing(state)
        .with_physical_optimizer_rules(rules)
        .build()
}
