//! DataFusion optimizer for genomic interval (INTERSECTS) joins.
//!
//! This crate provides a logical [`OptimizerRule`] that detects
//! `giql_intersects()` function calls in join filters and rewrites
//! them into binned equi-joins using UNNEST. Bin size is chosen
//! adaptively from table statistics when available.
//!
//! The `giql_intersects` function is a placeholder UDF emitted by the
//! GIQL transpiler's `"datafusion"` dialect. It preserves INTERSECTS
//! semantics through the SQL layer so the optimizer can match on it
//! directly, without heuristic pattern detection.
//!
//! # Usage
//!
//! ```rust,no_run
//! use datafusion::execution::SessionStateBuilder;
//! use datafusion::prelude::*;
//! use giql_datafusion::register_optimizer;
//!
//! let state = SessionStateBuilder::new()
//!     .with_default_features()
//!     .build();
//! let state = register_optimizer(state);
//! let ctx = SessionContext::from(state);
//! ```

pub mod coitree;
pub mod logical_rule;

pub use logical_rule::{IntersectsConfig, IntersectsLogicalRule};

use std::fmt::Debug;
use std::sync::Arc;

use datafusion::common::Result;
use datafusion::execution::SessionState;
use datafusion::logical_expr::{
    ColumnarValue, LogicalPlan, ScalarFunctionArgs, ScalarUDF,
    ScalarUDFImpl, Signature, TypeSignature, Volatility,
};
use datafusion::optimizer::OptimizerRule;
use datafusion::physical_plan::ExecutionPlan;

// ── Placeholder UDF ─────────────────────────────────────────────

/// Placeholder `giql_intersects(start_a, end_a, start_b, end_b)` UDF.
///
/// Exists only so DataFusion's SQL parser accepts the function call.
/// The logical optimizer rule rewrites it away before execution.
#[derive(Debug, Hash, PartialEq, Eq)]
struct GiqlIntersectsUdf {
    signature: Signature,
}

impl GiqlIntersectsUdf {
    fn new() -> Self {
        Self {
            signature: Signature::new(
                TypeSignature::Any(4),
                Volatility::Immutable,
            ),
        }
    }
}

impl ScalarUDFImpl for GiqlIntersectsUdf {
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn name(&self) -> &str {
        "giql_intersects"
    }

    fn signature(&self) -> &Signature {
        &self.signature
    }

    fn return_type(
        &self,
        _arg_types: &[arrow::datatypes::DataType],
    ) -> Result<arrow::datatypes::DataType> {
        Ok(arrow::datatypes::DataType::Boolean)
    }

    fn invoke_with_args(
        &self,
        _args: ScalarFunctionArgs,
    ) -> Result<ColumnarValue> {
        Err(datafusion::error::DataFusionError::Internal(
            "giql_intersects should be rewritten by the logical \
             optimizer rule — was the IntersectsLogicalRule registered?"
                .into(),
        ))
    }
}

/// Create the placeholder `giql_intersects` scalar UDF.
pub fn giql_intersects_udf() -> ScalarUDF {
    ScalarUDF::from(GiqlIntersectsUdf::new())
}

// ── Registration ────────────────────────────────────────────────

// ── Custom query planner ─────────────────────────────────────────

/// Query planner that includes the COI tree extension planner.
#[derive(Debug)]
struct GiqlQueryPlanner;

#[async_trait::async_trait]
impl datafusion::execution::context::QueryPlanner
    for GiqlQueryPlanner
{
    async fn create_physical_plan(
        &self,
        logical_plan: &LogicalPlan,
        session_state: &SessionState,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        use datafusion::physical_planner::{
            DefaultPhysicalPlanner, PhysicalPlanner,
        };

        let planner =
            DefaultPhysicalPlanner::with_extension_planners(vec![
                Arc::new(coitree::COITreePlanner),
            ]);
        planner
            .create_physical_plan(logical_plan, session_state)
            .await
    }
}

// ── Registration ────────────────────────────────────────────────

/// Build a [`SessionState`] with the INTERSECTS logical optimizer
/// rule, the `giql_intersects` placeholder UDF, and the COI tree
/// extension planner.
///
/// The logical rule detects `giql_intersects()` calls in join
/// filters and rewrites them into either:
/// - A binned equi-join (uniform width distributions)
/// - A COI tree join (non-uniform distributions)
pub fn register_optimizer(state: SessionState) -> SessionState {
    use datafusion::execution::SessionStateBuilder;

    let logical_rule: Arc<dyn OptimizerRule + Send + Sync> =
        Arc::new(IntersectsLogicalRule::new());

    let mut logical_rules: Vec<Arc<dyn OptimizerRule + Send + Sync>> =
        state.optimizers().to_vec();
    logical_rules.push(logical_rule);

    let udf = Arc::new(giql_intersects_udf());

    let mut scalar_fns: Vec<Arc<ScalarUDF>> =
        state.scalar_functions().values().cloned().collect();
    scalar_fns.push(udf);

    SessionStateBuilder::new_from_existing(state)
        .with_optimizer_rules(logical_rules)
        .with_scalar_functions(scalar_fns)
        .with_query_planner(Arc::new(GiqlQueryPlanner))
        .build()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_register_optimizer_adds_rule_and_udf() {
        use datafusion::execution::SessionStateBuilder;

        let state = SessionStateBuilder::new()
            .with_default_features()
            .build();
        let n_before = state.optimizers().len();

        let state = register_optimizer(state);

        // Logical rule was added
        assert_eq!(state.optimizers().len(), n_before + 1);
        let last_rule = state.optimizers().last().unwrap();
        assert_eq!(last_rule.name(), "intersects_logical_binned");

        // UDF was registered
        assert!(
            state.scalar_functions().contains_key("giql_intersects")
        );
    }
}
