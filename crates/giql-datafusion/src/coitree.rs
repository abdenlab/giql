//! COI tree interval join — build/probe execution using cache-oblivious
//! interval trees from the `coitrees` crate.
//!
//! Used for non-uniform width distributions where binning would cause
//! excessive replication. Each interval is stored exactly once in the
//! tree; queries are O(log N + k) per probe interval.

use std::any::Any;
use std::collections::HashMap;
use std::fmt;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use arrow::array::{
    Array, AsArray, Int64Array, RecordBatch, UInt32Array,
};
use arrow::compute;
use arrow::datatypes::SchemaRef;
use coitrees::{COITree, Interval, IntervalTree};
use datafusion::common::{Column, DFSchemaRef, Result};
use datafusion::execution::SendableRecordBatchStream;
use datafusion::logical_expr::{Expr, LogicalPlan, UserDefinedLogicalNode};
use datafusion::physical_expr::{EquivalenceProperties, Partitioning};
use datafusion::physical_plan::execution_plan::{
    Boundedness, EmissionType,
};
use datafusion::physical_plan::stream::RecordBatchStreamAdapter;
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, PlanProperties,
};

// ── Logical node ────────────────────────────────────────────────

/// Logical plan node representing a COI tree interval join.
///
/// The logical optimizer rule emits this when the sampled width
/// distribution is non-uniform (cost_optimal_bin > 2 * median).
/// The extension planner converts it to a [`COITreeExec`].
#[derive(Debug, Clone)]
pub struct COITreeJoinNode {
    pub left: Arc<LogicalPlan>,
    pub right: Arc<LogicalPlan>,
    /// Equi-keys from the original join (e.g., chrom = chrom).
    pub on: Vec<(Column, Column)>,
    /// Interval column names from giql_intersects() args.
    pub start_a: Column,
    pub end_a: Column,
    pub start_b: Column,
    pub end_b: Column,
    pub schema: DFSchemaRef,
}

impl Hash for COITreeJoinNode {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.on.hash(state);
        self.start_a.hash(state);
        self.end_a.hash(state);
        self.start_b.hash(state);
        self.end_b.hash(state);
    }
}

impl PartialEq for COITreeJoinNode {
    fn eq(&self, other: &Self) -> bool {
        self.on == other.on
            && self.start_a == other.start_a
            && self.end_a == other.end_a
            && self.start_b == other.start_b
            && self.end_b == other.end_b
    }
}

impl Eq for COITreeJoinNode {}

impl PartialOrd for COITreeJoinNode {
    fn partial_cmp(
        &self,
        _other: &Self,
    ) -> Option<std::cmp::Ordering> {
        None
    }
}

impl UserDefinedLogicalNode for COITreeJoinNode {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn name(&self) -> &str {
        "COITreeJoin"
    }

    fn inputs(&self) -> Vec<&LogicalPlan> {
        vec![&self.left, &self.right]
    }

    fn schema(&self) -> &DFSchemaRef {
        &self.schema
    }

    fn check_invariants(
        &self,
        _check: datafusion::logical_expr::InvariantLevel,
    ) -> Result<()> {
        Ok(())
    }

    fn expressions(&self) -> Vec<Expr> {
        vec![]
    }

    fn fmt_for_explain(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "COITreeJoin: on=[{}]",
            self.on
                .iter()
                .map(|(l, r)| format!("{l} = {r}"))
                .collect::<Vec<_>>()
                .join(", ")
        )
    }

    fn with_exprs_and_inputs(
        &self,
        _exprs: Vec<Expr>,
        inputs: Vec<LogicalPlan>,
    ) -> Result<Arc<dyn UserDefinedLogicalNode>> {
        Ok(Arc::new(COITreeJoinNode {
            left: Arc::new(inputs[0].clone()),
            right: Arc::new(inputs[1].clone()),
            on: self.on.clone(),
            start_a: self.start_a.clone(),
            end_a: self.end_a.clone(),
            start_b: self.start_b.clone(),
            end_b: self.end_b.clone(),
            schema: self.schema.clone(),
        }))
    }

    fn dyn_hash(&self, state: &mut dyn Hasher) {
        let mut s = state;
        self.hash(&mut s);
    }

    fn dyn_eq(&self, other: &dyn UserDefinedLogicalNode) -> bool {
        other
            .as_any()
            .downcast_ref::<Self>()
            .map_or(false, |o| self == o)
    }

    fn dyn_ord(
        &self,
        _other: &dyn UserDefinedLogicalNode,
    ) -> Option<std::cmp::Ordering> {
        None
    }
}

// ── Physical execution plan ─────────────────────────────────────

/// Physical execution plan that uses COI trees for interval joins.
///
/// Build phase: collect the left (build) side, group by chromosome,
/// construct a `COITree` per chromosome with row indices as metadata.
///
/// Probe phase: stream the right (probe) side batch by batch, query
/// the per-chromosome tree for each interval, emit joined output.
#[derive(Debug)]
pub struct COITreeExec {
    left: Arc<dyn ExecutionPlan>,
    right: Arc<dyn ExecutionPlan>,
    /// Equi-key column names (e.g., chrom).
    chrom_l: String,
    chrom_r: String,
    start_l: String,
    end_l: String,
    start_r: String,
    end_r: String,
    schema: SchemaRef,
    properties: Arc<PlanProperties>,
}

impl COITreeExec {
    pub fn new(
        left: Arc<dyn ExecutionPlan>,
        right: Arc<dyn ExecutionPlan>,
        on: &[(Column, Column)],
        start_a: &Column,
        end_a: &Column,
        start_b: &Column,
        end_b: &Column,
        schema: SchemaRef,
    ) -> Self {
        let properties = Arc::new(PlanProperties::new(
            EquivalenceProperties::new(schema.clone()),
            Partitioning::UnknownPartitioning(1),
            EmissionType::Final,
            Boundedness::Bounded,
        ));

        // Use the first equi-key as the chromosome column.
        let (chrom_l, chrom_r) = if let Some((l, r)) = on.first() {
            (l.name.clone(), r.name.clone())
        } else {
            ("chrom".to_string(), "chrom".to_string())
        };

        Self {
            left,
            right,
            chrom_l,
            chrom_r,
            start_l: start_a.name.clone(),
            end_l: end_a.name.clone(),
            start_r: start_b.name.clone(),
            end_r: end_b.name.clone(),
            schema,
            properties,
        }
    }
}

impl DisplayAs for COITreeExec {
    fn fmt_as(
        &self,
        _t: DisplayFormatType,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(
            f,
            "COITreeExec: on=[{} = {}]",
            self.chrom_l, self.chrom_r
        )
    }
}

impl ExecutionPlan for COITreeExec {
    fn name(&self) -> &str {
        "COITreeExec"
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn properties(&self) -> &Arc<PlanProperties> {
        &self.properties
    }

    fn children(&self) -> Vec<&Arc<dyn ExecutionPlan>> {
        vec![&self.left, &self.right]
    }

    fn with_new_children(
        self: Arc<Self>,
        children: Vec<Arc<dyn ExecutionPlan>>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        Ok(Arc::new(COITreeExec {
            left: children[0].clone(),
            right: children[1].clone(),
            chrom_l: self.chrom_l.clone(),
            chrom_r: self.chrom_r.clone(),
            start_l: self.start_l.clone(),
            end_l: self.end_l.clone(),
            start_r: self.start_r.clone(),
            end_r: self.end_r.clone(),
            schema: self.schema.clone(),
            properties: self.properties.clone(),
        }))
    }

    fn execute(
        &self,
        _partition: usize,
        context: Arc<datafusion::execution::TaskContext>,
    ) -> Result<SendableRecordBatchStream> {
        let left_plan = self.left.clone();
        let right_plan = self.right.clone();
        let schema = self.schema.clone();
        let chrom_l = self.chrom_l.clone();
        let chrom_r = self.chrom_r.clone();
        let start_l = self.start_l.clone();
        let end_l = self.end_l.clone();
        let start_r = self.start_r.clone();
        let end_r = self.end_r.clone();

        let stream = futures::stream::once(async move {
            // ── Build phase: collect left side, build COI trees ──
            let left_batches =
                datafusion::physical_plan::collect(
                    left_plan,
                    context.clone(),
                )
                .await?;

            if left_batches.is_empty() {
                return Ok(RecordBatch::new_empty(schema));
            }

            let left_concat = compute::concat_batches(
                &left_batches[0].schema(),
                &left_batches,
            )?;

            let left_schema = left_concat.schema();
            let l_chrom_idx =
                left_schema.index_of(&chrom_l)?;
            let l_start_idx =
                left_schema.index_of(&start_l)?;
            let l_end_idx = left_schema.index_of(&end_l)?;

            let l_chrom_col = left_concat.column(l_chrom_idx);
            let l_starts = left_concat
                .column(l_start_idx)
                .as_any()
                .downcast_ref::<Int64Array>()
                .ok_or_else(|| {
                    datafusion::error::DataFusionError::Internal(
                        "start column is not Int64".into(),
                    )
                })?;
            let l_ends = left_concat
                .column(l_end_idx)
                .as_any()
                .downcast_ref::<Int64Array>()
                .ok_or_else(|| {
                    datafusion::error::DataFusionError::Internal(
                        "end column is not Int64".into(),
                    )
                })?;

            // Group intervals by chromosome and build COI trees.
            // Metadata = row index in left_concat.
            let mut chrom_intervals: HashMap<
                String,
                Vec<Interval<u32>>,
            > = HashMap::new();
            for i in 0..left_concat.num_rows() {
                let chrom: &str = l_chrom_col.as_string_view().value(i);
                let start = l_starts.value(i) as i32;
                // Half-open [start, end) → end-inclusive [start, end-1]
                let end = (l_ends.value(i) - 1) as i32;
                chrom_intervals
                    .entry(chrom.to_string())
                    .or_default()
                    .push(Interval::new(start, end, i as u32));
            }

            let trees: HashMap<String, COITree<u32, u32>> =
                chrom_intervals
                    .iter()
                    .map(|(chrom, intervals)| {
                        (chrom.clone(), COITree::new(intervals))
                    })
                    .collect();

            // ── Probe phase: stream right side, query trees ─────
            let right_batches =
                datafusion::physical_plan::collect(
                    right_plan, context,
                )
                .await?;

            if right_batches.is_empty() {
                return Ok(RecordBatch::new_empty(schema));
            }

            let right_concat = compute::concat_batches(
                &right_batches[0].schema(),
                &right_batches,
            )?;

            let right_schema = right_concat.schema();
            let r_chrom_idx =
                right_schema.index_of(&chrom_r)?;
            let r_start_idx =
                right_schema.index_of(&start_r)?;
            let r_end_idx =
                right_schema.index_of(&end_r)?;

            let r_chrom_col = right_concat.column(r_chrom_idx);
            let r_starts = right_concat
                .column(r_start_idx)
                .as_any()
                .downcast_ref::<Int64Array>()
                .ok_or_else(|| {
                    datafusion::error::DataFusionError::Internal(
                        "start column is not Int64".into(),
                    )
                })?;
            let r_ends = right_concat
                .column(r_end_idx)
                .as_any()
                .downcast_ref::<Int64Array>()
                .ok_or_else(|| {
                    datafusion::error::DataFusionError::Internal(
                        "end column is not Int64".into(),
                    )
                })?;

            // Collect join pairs as (left_idx, right_idx).
            let mut left_indices: Vec<u32> = Vec::new();
            let mut right_indices: Vec<u32> = Vec::new();

            for i in 0..right_concat.num_rows() {
                let chrom: &str = r_chrom_col.as_string_view().value(i);
                if let Some(tree) = trees.get(chrom) {
                    let start = r_starts.value(i) as i32;
                    let end = (r_ends.value(i) - 1) as i32;
                    tree.query(start, end, |hit| {
                        left_indices
                            .push(*hit.metadata);
                        right_indices.push(i as u32);
                    });
                }
            }

            if left_indices.is_empty() {
                return Ok(RecordBatch::new_empty(schema));
            }

            // Build output batch using take() on both sides.
            let left_idx_arr = UInt32Array::from(left_indices);
            let right_idx_arr = UInt32Array::from(right_indices);

            let mut output_columns: Vec<Arc<dyn Array>> =
                Vec::with_capacity(schema.fields().len());

            // Left side columns
            for col in left_concat.columns() {
                output_columns.push(compute::take(
                    col.as_ref(),
                    &left_idx_arr,
                    None,
                )?);
            }
            // Right side columns
            for col in right_concat.columns() {
                output_columns.push(compute::take(
                    col.as_ref(),
                    &right_idx_arr,
                    None,
                )?);
            }

            RecordBatch::try_new(schema, output_columns)
                .map_err(datafusion::error::DataFusionError::from)
        });

        Ok(Box::pin(RecordBatchStreamAdapter::new(
            self.schema.clone(),
            stream,
        )))
    }
}

// ── Extension planner ───────────────────────────────────────────

/// Converts [`COITreeJoinNode`] logical nodes into [`COITreeExec`]
/// physical plans.
#[derive(Debug)]
pub struct COITreePlanner;

#[async_trait::async_trait]
impl datafusion::physical_planner::ExtensionPlanner
    for COITreePlanner
{
    async fn plan_extension(
        &self,
        _planner: &dyn datafusion::physical_planner::PhysicalPlanner,
        node: &dyn UserDefinedLogicalNode,
        _logical_inputs: &[&LogicalPlan],
        physical_inputs: &[Arc<dyn ExecutionPlan>],
        _session_state: &datafusion::execution::SessionState,
    ) -> Result<Option<Arc<dyn ExecutionPlan>>> {
        let Some(join_node) =
            node.as_any().downcast_ref::<COITreeJoinNode>()
        else {
            return Ok(None);
        };

        // Build the output Arrow schema from the logical schema.
        let arrow_schema: SchemaRef =
            Arc::new(join_node.schema.as_arrow().clone());

        Ok(Some(Arc::new(COITreeExec::new(
            physical_inputs[0].clone(),
            physical_inputs[1].clone(),
            &join_node
                .on
                .iter()
                .map(|(l, r)| (l.clone(), r.clone()))
                .collect::<Vec<_>>(),
            &join_node.start_a,
            &join_node.end_a,
            &join_node.start_b,
            &join_node.end_b,
            arrow_schema,
        ))))
    }
}
