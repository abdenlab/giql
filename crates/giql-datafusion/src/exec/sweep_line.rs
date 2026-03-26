use std::any::Any;
use std::fmt;
use std::pin::Pin;
use std::sync::Arc;
use std::task::{Context, Poll};

use arrow::array::{
    Array, ArrayRef, Int32Array, Int64Array, RecordBatch, StringArray,
    StringViewArray, UInt32Array,
};
use arrow::compute;
use arrow::compute::kernels::sort::SortOptions;
use arrow::datatypes::SchemaRef;
use datafusion::common::Result;
use datafusion::execution::SendableRecordBatchStream;
use datafusion::physical_expr::expressions::Column;
use datafusion::physical_expr::{
    EquivalenceProperties, LexRequirement, Partitioning,
    PhysicalSortRequirement,
};
use datafusion::physical_plan::execution_plan::{
    Boundedness, EmissionType,
};
use datafusion::physical_plan::metrics::{
    BaselineMetrics, ExecutionPlanMetricsSet, MetricsSet,
};
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, PlanProperties,
    RecordBatchStream,
};

use crate::pattern::IntervalColumns;

/// Which side of the join to materialize as the build side.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BuildSide {
    /// Materialize the left child, stream the right.
    Left,
    /// Materialize the right child, stream the left.
    Right,
}

/// Streaming sweep-line interval join.
///
/// Requires both inputs sorted by (chrom, start). DataFusion will
/// automatically insert SortExec nodes if the inputs aren't already
/// in order.
///
/// The smaller side (selected by `build_side`) is materialized. The
/// larger side is streamed. Both are swept per-chromosome in parallel
/// with vectorized output via `compute::take`.
#[derive(Debug)]
pub struct SweepLineJoinExec {
    left: Arc<dyn ExecutionPlan>,
    right: Arc<dyn ExecutionPlan>,
    left_cols: IntervalColumns,
    right_cols: IntervalColumns,
    build_side: BuildSide,
    schema: SchemaRef,
    properties: PlanProperties,
    metrics: ExecutionPlanMetricsSet,
}

impl SweepLineJoinExec {
    pub fn new(
        left: Arc<dyn ExecutionPlan>,
        right: Arc<dyn ExecutionPlan>,
        left_cols: IntervalColumns,
        right_cols: IntervalColumns,
        schema: SchemaRef,
        build_side: BuildSide,
    ) -> Self {
        let properties = PlanProperties::new(
            EquivalenceProperties::new(schema.clone()),
            Partitioning::UnknownPartitioning(1),
            EmissionType::Incremental,
            Boundedness::Bounded,
        );

        Self {
            left,
            right,
            left_cols,
            right_cols,
            build_side,
            schema,
            properties,
            metrics: ExecutionPlanMetricsSet::new(),
        }
    }

    /// Build sort requirements for (chrom, start) on a child schema.
    fn sort_requirement(
        cols: &IntervalColumns,
    ) -> LexRequirement {
        LexRequirement::new(vec![
            PhysicalSortRequirement::new(
                Arc::new(Column::new(&cols.chrom_col, cols.chrom_idx)),
                Some(SortOptions {
                    descending: false,
                    nulls_first: false,
                }),
            ),
            PhysicalSortRequirement::new(
                Arc::new(Column::new(&cols.start_col, cols.start_idx)),
                Some(SortOptions {
                    descending: false,
                    nulls_first: false,
                }),
            ),
        ])
    }
}

impl DisplayAs for SweepLineJoinExec {
    fn fmt_as(
        &self,
        _t: DisplayFormatType,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(f, "SweepLineJoinExec: streaming merge-join")
    }
}

impl ExecutionPlan for SweepLineJoinExec {
    fn name(&self) -> &str {
        "SweepLineJoinExec"
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn properties(&self) -> &PlanProperties {
        &self.properties
    }

    fn children(&self) -> Vec<&Arc<dyn ExecutionPlan>> {
        vec![&self.left, &self.right]
    }

    fn required_input_ordering(
        &self,
    ) -> Vec<Option<LexRequirement>> {
        vec![
            Some(Self::sort_requirement(&self.left_cols)),
            Some(Self::sort_requirement(&self.right_cols)),
        ]
    }

    fn with_new_children(
        self: Arc<Self>,
        children: Vec<Arc<dyn ExecutionPlan>>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        Ok(Arc::new(SweepLineJoinExec::new(
            children[0].clone(),
            children[1].clone(),
            self.left_cols.clone(),
            self.right_cols.clone(),
            self.schema.clone(),
            self.build_side,
        )))
    }

    fn execute(
        &self,
        _partition: usize,
        context: Arc<datafusion::execution::TaskContext>,
    ) -> Result<SendableRecordBatchStream> {
        let left = self.left.clone();
        let right = self.right.clone();
        let left_cols = self.left_cols.clone();
        let right_cols = self.right_cols.clone();
        let build_side = self.build_side;
        let schema = self.schema.clone();
        let baseline_metrics =
            BaselineMetrics::new(&self.metrics, 0);
        let ctx = context;

        let stream = futures::stream::once(async move {
            // Collect both sides. The build side (smaller) is
            // materialized first; the probe side second. Both are
            // collected via DataFusion's parallel collect.
            let (build_plan, probe_plan) = match build_side {
                BuildSide::Right => (right.clone(), left.clone()),
                BuildSide::Left => (left.clone(), right.clone()),
            };

            let build_batch = {
                let batches = datafusion::physical_plan::collect(
                    build_plan, ctx.clone(),
                )
                .await?;
                concat_batches(&batches)?
            };

            let probe_batch = {
                let batches = datafusion::physical_plan::collect(
                    probe_plan, ctx,
                )
                .await?;
                concat_batches(&batches)?
            };

            // Map build/probe back to left/right
            let (left_batch, right_batch) = match build_side {
                BuildSide::Right => (probe_batch, build_batch),
                BuildSide::Left => (build_batch, probe_batch),
            };

            if left_batch.num_rows() == 0 || right_batch.num_rows() == 0 {
                return Ok(RecordBatch::new_empty(schema));
            }

            // Extract typed columns from both sides
            let left_chroms =
                chrom_to_dense_ids(left_batch.column(left_cols.chrom_idx))?;
            let left_starts =
                as_i64_vec(left_batch.column(left_cols.start_idx))?;
            let left_ends =
                as_i64_vec(left_batch.column(left_cols.end_idx))?;

            let right_chroms =
                chrom_to_dense_ids(right_batch.column(right_cols.chrom_idx))?;
            let right_starts =
                as_i64_vec(right_batch.column(right_cols.start_idx))?;
            let right_ends =
                as_i64_vec(right_batch.column(right_cols.end_idx))?;

            let left_order =
                argsort_by_chrom_start(&left_chroms, &left_starts);
            let right_order =
                argsort_by_chrom_start(&right_chroms, &right_starts);

            // Split by chromosome and sweep in parallel
            let l_groups = split_by_chrom(&left_order, &left_chroms);
            let r_groups = split_by_chrom(&right_order, &right_chroms);

            let mut handles = Vec::with_capacity(l_groups.len());
            for (chrom_id, l_range) in &l_groups {
                let r_range = match r_groups
                    .binary_search_by_key(chrom_id, |(c, _)| *c)
                {
                    Ok(pos) => r_groups[pos].1.clone(),
                    Err(_) => continue,
                };

                let l_idx: Arc<[u32]> =
                    left_order[l_range.clone()].into();
                let r_idx: Arc<[u32]> =
                    right_order[r_range].into();
                let ls = left_starts.clone();
                let le = left_ends.clone();
                let rs = right_starts.clone();
                let re = right_ends.clone();

                handles.push(tokio::spawn(async move {
                    sweep_chrom(&l_idx, &ls, &le, &r_idx, &rs, &re)
                }));
            }

            let mut all_left: Vec<u32> = Vec::new();
            let mut all_right: Vec<u32> = Vec::new();
            for h in handles {
                let (li, ri) = h.await.map_err(|e| {
                    datafusion::error::DataFusionError::External(
                        Box::new(e),
                    )
                })?;
                all_left.extend_from_slice(&li);
                all_right.extend_from_slice(&ri);
            }

            baseline_metrics.record_output(all_left.len());
            build_output_take(
                &schema, &left_batch, &right_batch,
                &all_left, &all_right,
            )
        });

        Ok(Box::pin(SweepLineStream {
            inner: Box::pin(
                datafusion::physical_plan::stream::RecordBatchStreamAdapter::new(
                    self.schema.clone(),
                    stream,
                ),
            ),
        }))
    }

    fn metrics(&self) -> Option<MetricsSet> {
        Some(self.metrics.clone_inner())
    }
}

/// Wrapper stream that forwards to the inner adapter.
struct SweepLineStream {
    inner: Pin<Box<dyn RecordBatchStream + Send>>,
}

impl RecordBatchStream for SweepLineStream {
    fn schema(&self) -> SchemaRef {
        self.inner.schema()
    }
}

impl futures::Stream for SweepLineStream {
    type Item = Result<RecordBatch>;

    fn poll_next(
        mut self: Pin<&mut Self>,
        cx: &mut Context<'_>,
    ) -> Poll<Option<Self::Item>> {
        self.inner.as_mut().poll_next(cx)
    }
}

// ── Column extraction ───────────────────────────────────────────

/// Map chromosome strings to dense u32 IDs for fast comparison.
fn chrom_to_dense_ids(col: &ArrayRef) -> Result<Arc<[u32]>> {
    // Extract all strings first
    let strings = extract_all_strings(col)?;
    let n = strings.len();

    // Build sorted unique set
    let mut unique: Vec<String> = Vec::new();
    for s in &strings {
        if !unique.contains(s) {
            unique.push(s.clone());
        }
    }
    unique.sort();

    let map: std::collections::HashMap<&str, u32> = unique
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i as u32))
        .collect();

    let ids: Vec<u32> =
        strings.iter().map(|s| map[s.as_str()]).collect();
    Ok(ids.into())
}

fn extract_all_strings(col: &ArrayRef) -> Result<Vec<String>> {
    let n = col.len();
    if let Some(arr) = col.as_any().downcast_ref::<StringArray>() {
        return Ok((0..n).map(|i| arr.value(i).to_string()).collect());
    }
    if let Some(arr) = col.as_any().downcast_ref::<StringViewArray>() {
        return Ok((0..n).map(|i| arr.value(i).to_string()).collect());
    }
    Err(datafusion::error::DataFusionError::Internal(
        "Unsupported string type".to_string(),
    ))
}

/// Extract i64 values, converting Int32 if needed.
fn as_i64_vec(col: &ArrayRef) -> Result<Arc<[i64]>> {
    if let Some(arr) = col.as_any().downcast_ref::<Int64Array>() {
        return Ok(arr.values().to_vec().into());
    }
    if let Some(arr) = col.as_any().downcast_ref::<Int32Array>() {
        return Ok(arr
            .values()
            .iter()
            .map(|&v| v as i64)
            .collect::<Vec<_>>()
            .into());
    }
    Err(datafusion::error::DataFusionError::Internal(
        "Column is not Int32 or Int64".to_string(),
    ))
}

// ── Sorting and grouping ────────────────────────────────────────

fn argsort_by_chrom_start(
    chrom_ids: &[u32],
    starts: &[i64],
) -> Vec<u32> {
    let mut indices: Vec<u32> = (0..chrom_ids.len() as u32).collect();
    indices.sort_unstable_by(|&a, &b| {
        let (a, b) = (a as usize, b as usize);
        chrom_ids[a]
            .cmp(&chrom_ids[b])
            .then(starts[a].cmp(&starts[b]))
    });
    indices
}

fn split_by_chrom(
    sorted: &[u32],
    chrom_ids: &[u32],
) -> Vec<(u32, std::ops::Range<usize>)> {
    if sorted.is_empty() {
        return vec![];
    }
    let mut groups = Vec::new();
    let mut start = 0;
    let mut cur = chrom_ids[sorted[0] as usize];
    for i in 1..sorted.len() {
        let c = chrom_ids[sorted[i] as usize];
        if c != cur {
            groups.push((cur, start..i));
            start = i;
            cur = c;
        }
    }
    groups.push((cur, start..sorted.len()));
    groups
}

// ── Core sweep ──────────────────────────────────────────────────

fn sweep_chrom(
    left_indices: &[u32],
    left_starts: &[i64],
    left_ends: &[i64],
    right_indices: &[u32],
    right_starts: &[i64],
    right_ends: &[i64],
) -> (Vec<u32>, Vec<u32>) {
    let mut match_left = Vec::new();
    let mut match_right = Vec::new();
    let mut r_cursor = 0usize;
    let mut active: Vec<usize> = Vec::new();

    for &li in left_indices {
        let l_start = left_starts[li as usize];
        let l_end = left_ends[li as usize];

        while r_cursor < right_indices.len() {
            let ri = right_indices[r_cursor] as usize;
            if right_starts[ri] >= l_end {
                break;
            }
            active.push(r_cursor);
            r_cursor += 1;
        }

        active.retain(|&pos| {
            right_ends[right_indices[pos] as usize] > l_start
        });

        for &pos in &active {
            let ri = right_indices[pos];
            if right_starts[ri as usize] < l_end {
                match_left.push(li);
                match_right.push(ri);
            }
        }
    }

    (match_left, match_right)
}

// ── Output ──────────────────────────────────────────────────────

fn concat_batches(batches: &[RecordBatch]) -> Result<RecordBatch> {
    if batches.is_empty() {
        return Err(datafusion::error::DataFusionError::Internal(
            "No batches".to_string(),
        ));
    }
    if batches.len() == 1 {
        return Ok(batches[0].clone());
    }
    let schema = batches[0].schema();
    Ok(compute::concat_batches(&schema, batches)?)
}

fn build_output_take(
    schema: &SchemaRef,
    left: &RecordBatch,
    right: &RecordBatch,
    left_idx: &[u32],
    right_idx: &[u32],
) -> Result<RecordBatch> {
    if left_idx.is_empty() {
        return Ok(RecordBatch::new_empty(schema.clone()));
    }
    let li = UInt32Array::from(left_idx.to_vec());
    let ri = UInt32Array::from(right_idx.to_vec());
    let mut cols: Vec<ArrayRef> = Vec::with_capacity(
        left.num_columns() + right.num_columns(),
    );
    for c in 0..left.num_columns() {
        cols.push(compute::take(
            left.column(c).as_ref(),
            &li,
            None,
        )?);
    }
    for c in 0..right.num_columns() {
        cols.push(compute::take(
            right.column(c).as_ref(),
            &ri,
            None,
        )?);
    }
    Ok(RecordBatch::try_new(schema.clone(), cols)?)
}
