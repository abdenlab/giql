use std::any::Any;
use std::fmt;
use std::ops::Range;
use std::pin::Pin;
use std::sync::Arc;
use std::future::Future;
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
    EquivalenceProperties, LexRequirement, OrderingRequirements,
    Partitioning, PhysicalSortRequirement,
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
    Left,
    Right,
}

/// Streaming sweep-line interval join.
///
/// Materializes the build side (smaller), then streams the probe side
/// batch-by-batch. Each probe batch is swept against the sorted build
/// side per-chromosome, producing output incrementally.
///
/// Requires both inputs sorted by (chrom, start).
#[derive(Debug)]
pub struct SweepLineJoinExec {
    left: Arc<dyn ExecutionPlan>,
    right: Arc<dyn ExecutionPlan>,
    left_cols: IntervalColumns,
    right_cols: IntervalColumns,
    build_side: BuildSide,
    schema: SchemaRef,
    properties: Arc<PlanProperties>,
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
        let properties = Arc::new(PlanProperties::new(
            EquivalenceProperties::new(schema.clone()),
            Partitioning::UnknownPartitioning(1),
            EmissionType::Incremental,
            Boundedness::Bounded,
        ));
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

    fn sort_requirement(cols: &IntervalColumns) -> OrderingRequirements {
        let lex = LexRequirement::new(vec![
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
        .expect("sort requirement should be non-empty");
        OrderingRequirements::new(lex)
    }
}

impl DisplayAs for SweepLineJoinExec {
    fn fmt_as(
        &self,
        _t: DisplayFormatType,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(
            f,
            "SweepLineJoinExec: build={:?}",
            self.build_side
        )
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
    fn properties(&self) -> &Arc<PlanProperties> {
        &self.properties
    }
    fn children(&self) -> Vec<&Arc<dyn ExecutionPlan>> {
        vec![&self.left, &self.right]
    }
    fn required_input_ordering(
        &self,
    ) -> Vec<Option<OrderingRequirements>> {
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
        let (build_plan, probe_plan, build_cols, probe_cols) =
            match self.build_side {
                BuildSide::Right => (
                    self.right.clone(),
                    self.left.clone(),
                    self.right_cols.clone(),
                    self.left_cols.clone(),
                ),
                BuildSide::Left => (
                    self.left.clone(),
                    self.right.clone(),
                    self.left_cols.clone(),
                    self.right_cols.clone(),
                ),
            };

        // Spawn build-side collection as a background task.
        let build_ctx = context.clone();
        let build_fut = tokio::spawn(async move {
            let batches = datafusion::physical_plan::collect(
                build_plan, build_ctx,
            )
            .await?;
            BuildSideData::from_batches(&batches, &build_cols)
        });

        // Open probe-side stream. If the probe plan has multiple
        // partitions (e.g. from RepartitionExec), coalesce them
        // into a single stream so we see all rows.
        use datafusion::physical_plan::coalesce_partitions::CoalescePartitionsExec;
        let probe_plan = if probe_plan
            .properties()
            .partitioning
            .partition_count()
            > 1
        {
            Arc::new(CoalescePartitionsExec::new(probe_plan))
                as Arc<dyn ExecutionPlan>
        } else {
            probe_plan
        };
        let probe_stream = probe_plan.execute(0, context)?;

        Ok(Box::pin(SweepLineStream {
            state: SweepLineState::WaitBuildSide,
            probe_stream,
            build_fut: Some(build_fut),
            build_data: None,
            probe_cols,
            build_side: self.build_side,
            schema: self.schema.clone(),
            metrics: BaselineMetrics::new(&self.metrics, 0),
        }))
    }

    fn metrics(&self) -> Option<MetricsSet> {
        Some(self.metrics.clone_inner())
    }
}

// ── Build side data ─────────────────────────────────────────────

/// Materialized, sorted, and indexed build-side data.
struct BuildSideData {
    batch: RecordBatch,
    starts: Arc<[i64]>,
    ends: Arc<[i64]>,
    /// Sorted indices into the batch, ordered by (chrom, start).
    sorted_order: Vec<u32>,
    /// Per-chromosome ranges into sorted_order, keyed by chrom name.
    chrom_groups: Vec<(String, Range<usize>)>,
}

impl BuildSideData {
    fn from_batches(
        batches: &[RecordBatch],
        cols: &IntervalColumns,
    ) -> Result<Self> {
        let batch = concat_batches(batches)?;
        let chrom_strings =
            extract_all_strings(batch.column(cols.chrom_idx))?;
        let chrom_ids = strings_to_ids(&chrom_strings);
        let starts = as_i64_vec(batch.column(cols.start_idx))?;
        let ends = as_i64_vec(batch.column(cols.end_idx))?;
        let sorted_order = argsort_by_chrom_start(&chrom_ids, &starts);
        // Build groups keyed by chrom name (not ID) so probe batches
        // can match regardless of their own ID assignment.
        let id_groups = split_by_chrom(&sorted_order, &chrom_ids);
        let id_to_name = ids_to_names(&chrom_strings, &chrom_ids);
        let chrom_groups: Vec<(String, Range<usize>)> = id_groups
            .into_iter()
            .map(|(id, range)| (id_to_name[&id].clone(), range))
            .collect();
        Ok(Self { batch, starts, ends, sorted_order, chrom_groups })
    }

    /// Find the sorted index range for a given chromosome name.
    fn chrom_range(&self, chrom: &str) -> Option<Range<usize>> {
        self.chrom_groups
            .iter()
            .find(|(c, _)| c == chrom)
            .map(|(_, r)| r.clone())
    }
}

// ── Stream state machine ────────────────────────────────────────

enum SweepLineState {
    WaitBuildSide,
    FetchProbeBatch,
    Completed,
}

struct SweepLineStream {
    state: SweepLineState,
    probe_stream: SendableRecordBatchStream,
    build_fut: Option<tokio::task::JoinHandle<Result<BuildSideData>>>,
    build_data: Option<Arc<BuildSideData>>,
    probe_cols: IntervalColumns,
    build_side: BuildSide,
    schema: SchemaRef,
    metrics: BaselineMetrics,
}

impl RecordBatchStream for SweepLineStream {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }
}

impl futures::Stream for SweepLineStream {
    type Item = Result<RecordBatch>;

    fn poll_next(
        mut self: Pin<&mut Self>,
        cx: &mut Context<'_>,
    ) -> Poll<Option<Self::Item>> {
        loop {
            match self.state {
                SweepLineState::WaitBuildSide => {
                    // Poll the build-side future.
                    let fut = self.build_fut.as_mut().unwrap();
                    // SAFETY: we only poll this once, and it's
                    // behind an Option that we take after Ready.
                    let fut = unsafe { Pin::new_unchecked(fut) };
                    match fut.poll(cx) {
                        Poll::Pending => return Poll::Pending,
                        Poll::Ready(join_result) => {
                            self.build_fut = None;
                            let build_data = join_result
                                .map_err(|e| {
                                    datafusion::error::DataFusionError::External(
                                        Box::new(e),
                                    )
                                })?
                                ?;
                            self.build_data =
                                Some(Arc::new(build_data));
                            self.state =
                                SweepLineState::FetchProbeBatch;
                        }
                    }
                }

                SweepLineState::FetchProbeBatch => {
                    // Poll the probe stream for the next batch.
                    match self
                        .probe_stream
                        .as_mut()
                        .poll_next(cx)
                    {
                        Poll::Pending => return Poll::Pending,
                        Poll::Ready(None) => {
                            self.state = SweepLineState::Completed;
                            return Poll::Ready(None);
                        }
                        Poll::Ready(Some(Err(e))) => {
                            return Poll::Ready(Some(Err(e)));
                        }
                        Poll::Ready(Some(Ok(probe_batch))) => {
                            if probe_batch.num_rows() == 0 {
                                continue; // skip empty batches
                            }
                            let build = self
                                .build_data
                                .as_ref()
                                .unwrap();
                            let result = process_probe_batch(
                                &self.schema,
                                build,
                                &probe_batch,
                                &self.probe_cols,
                                self.build_side,
                            );
                            if let Ok(ref batch) = result {
                                self.metrics
                                    .record_output(batch.num_rows());
                                if batch.num_rows() == 0 {
                                    continue;
                                }
                            }
                            return Poll::Ready(Some(result));
                        }
                    }
                }

                SweepLineState::Completed => {
                    return Poll::Ready(None);
                }
            }
        }
    }
}

// ── Per-batch processing ────────────────────────────────────────

/// Sweep a single probe batch against the build side, producing
/// matched output rows.
fn process_probe_batch(
    schema: &SchemaRef,
    build: &BuildSideData,
    probe_batch: &RecordBatch,
    probe_cols: &IntervalColumns,
    build_side: BuildSide,
) -> Result<RecordBatch> {
    let probe_chrom_strings =
        extract_all_strings(probe_batch.column(probe_cols.chrom_idx))?;
    let probe_chrom_ids = strings_to_ids(&probe_chrom_strings);
    let probe_starts =
        as_i64_vec(probe_batch.column(probe_cols.start_idx))?;
    let probe_ends =
        as_i64_vec(probe_batch.column(probe_cols.end_idx))?;

    let mut match_build: Vec<u32> = Vec::new();
    let mut match_probe: Vec<u32> = Vec::new();

    // Group probe rows by chromosome, match by name against build
    let probe_order =
        argsort_by_chrom_start(&probe_chrom_ids, &probe_starts);
    let probe_id_groups =
        split_by_chrom(&probe_order, &probe_chrom_ids);
    let id_to_name = ids_to_names(&probe_chrom_strings, &probe_chrom_ids);

    for (probe_cid, p_range) in &probe_id_groups {
        let chrom_name = &id_to_name[probe_cid];
        let b_range = match build.chrom_range(chrom_name) {
            Some(r) => r,
            None => continue,
        };

        let p_idx = &probe_order[p_range.clone()];
        let b_idx = &build.sorted_order[b_range];

        let (bl, pl) = sweep_chrom(
            p_idx,
            &probe_starts,
            &probe_ends,
            b_idx,
            &build.starts,
            &build.ends,
        );

        match_probe.extend_from_slice(&pl);
        match_build.extend_from_slice(&bl);
    }

    if match_probe.is_empty() {
        return Ok(RecordBatch::new_empty(schema.clone()));
    }

    // Build output: column order depends on which side is build/probe.
    // Output must be (left_cols..., right_cols...).
    let (left_batch, right_batch, left_idx, right_idx) =
        match build_side {
            BuildSide::Left => (
                &build.batch,
                probe_batch,
                &match_build,
                &match_probe,
            ),
            BuildSide::Right => (
                probe_batch,
                &build.batch,
                &match_probe,
                &match_build,
            ),
        };

    build_output_take(schema, left_batch, right_batch, left_idx, right_idx)
}

// ── Core sweep ──────────────────────────────────────────────────

/// Sweep probe intervals against build intervals for one chromosome.
/// Returns (build_matches, probe_matches).
fn sweep_chrom(
    probe_indices: &[u32],
    probe_starts: &[i64],
    probe_ends: &[i64],
    build_indices: &[u32],
    build_starts: &[i64],
    build_ends: &[i64],
) -> (Vec<u32>, Vec<u32>) {
    let mut match_build = Vec::new();
    let mut match_probe = Vec::new();
    let mut b_cursor = 0usize;
    let mut active: Vec<usize> = Vec::new();

    for &pi in probe_indices {
        let p_start = probe_starts[pi as usize];
        let p_end = probe_ends[pi as usize];

        // Add build intervals with start < probe_end
        while b_cursor < build_indices.len() {
            let bi = build_indices[b_cursor] as usize;
            if build_starts[bi] >= p_end {
                break;
            }
            active.push(b_cursor);
            b_cursor += 1;
        }

        // Remove expired (end <= probe_start)
        active.retain(|&pos| {
            build_ends[build_indices[pos] as usize] > p_start
        });

        // Emit overlapping pairs
        for &pos in &active {
            let bi = build_indices[pos];
            if build_starts[bi as usize] < p_end {
                match_build.push(bi);
                match_probe.push(pi);
            }
        }
    }

    (match_build, match_probe)
}

// ── Helpers ─────────────────────────────────────────────────────

/// Assign dense u32 IDs to strings (sorted order).
fn strings_to_ids(strings: &[String]) -> Vec<u32> {
    let mut unique: Vec<String> = Vec::new();
    for s in strings {
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
    strings.iter().map(|s| map[s.as_str()]).collect()
}

/// Build reverse mapping from u32 ID → chromosome name.
fn ids_to_names(
    strings: &[String],
    ids: &[u32],
) -> std::collections::HashMap<u32, String> {
    let mut map = std::collections::HashMap::new();
    for (s, &id) in strings.iter().zip(ids.iter()) {
        map.entry(id).or_insert_with(|| s.clone());
    }
    map
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

fn as_i64_vec(col: &ArrayRef) -> Result<Arc<[i64]>> {
    if let Some(arr) = col.as_any().downcast_ref::<Int64Array>() {
        return Ok(arr.values().to_vec().into());
    }
    if let Some(arr) = col.as_any().downcast_ref::<Int32Array>() {
        return Ok(arr.values().iter().map(|&v| v as i64).collect::<Vec<_>>().into());
    }
    Err(datafusion::error::DataFusionError::Internal(
        "Column is not Int32 or Int64".to_string(),
    ))
}

fn argsort_by_chrom_start(chrom_ids: &[u32], starts: &[i64]) -> Vec<u32> {
    let mut idx: Vec<u32> = (0..chrom_ids.len() as u32).collect();
    idx.sort_unstable_by(|&a, &b| {
        let (a, b) = (a as usize, b as usize);
        chrom_ids[a].cmp(&chrom_ids[b]).then(starts[a].cmp(&starts[b]))
    });
    idx
}

fn split_by_chrom(
    sorted: &[u32],
    chrom_ids: &[u32],
) -> Vec<(u32, Range<usize>)> {
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

fn concat_batches(batches: &[RecordBatch]) -> Result<RecordBatch> {
    if batches.is_empty() {
        return Err(datafusion::error::DataFusionError::Internal(
            "No batches".to_string(),
        ));
    }
    if batches.len() == 1 {
        return Ok(batches[0].clone());
    }
    Ok(compute::concat_batches(&batches[0].schema(), batches)?)
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
    let mut cols: Vec<ArrayRef> =
        Vec::with_capacity(left.num_columns() + right.num_columns());
    for c in 0..left.num_columns() {
        cols.push(compute::take(left.column(c).as_ref(), &li, None)?);
    }
    for c in 0..right.num_columns() {
        cols.push(compute::take(right.column(c).as_ref(), &ri, None)?);
    }
    Ok(RecordBatch::try_new(schema.clone(), cols)?)
}
