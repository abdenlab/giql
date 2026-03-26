use std::any::Any;
use std::fmt;
use std::sync::Arc;

use arrow::array::{
    Array, ArrayRef, Int32Array, Int64Array, RecordBatch, StringArray,
    StringViewArray, UInt32Array,
};
use arrow::compute;
use arrow::datatypes::SchemaRef;
use datafusion::common::Result;
use datafusion::execution::SendableRecordBatchStream;
use datafusion::physical_expr::{EquivalenceProperties, Partitioning};
use datafusion::physical_plan::execution_plan::{
    Boundedness, EmissionType,
};
use datafusion::physical_plan::stream::RecordBatchStreamAdapter;
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, PlanProperties,
};

use crate::pattern::IntervalColumns;

/// Sweep-line interval join parallelized by chromosome.
///
/// 1. Collect and concat both inputs into contiguous Arrow arrays
/// 2. Assign integer chromosome IDs (avoids String comparisons)
/// 3. Sort by (chrom_id, start) and split at chromosome boundaries
/// 4. Sweep each chromosome in parallel via tokio::spawn
/// 5. Build output with vectorized compute::take
#[derive(Debug)]
pub struct SweepLineJoinExec {
    left: Arc<dyn ExecutionPlan>,
    right: Arc<dyn ExecutionPlan>,
    left_cols: IntervalColumns,
    right_cols: IntervalColumns,
    skip_sort: bool,
    schema: SchemaRef,
    properties: PlanProperties,
}

impl SweepLineJoinExec {
    pub fn new(
        left: Arc<dyn ExecutionPlan>,
        right: Arc<dyn ExecutionPlan>,
        left_cols: IntervalColumns,
        right_cols: IntervalColumns,
        schema: SchemaRef,
        skip_sort: bool,
    ) -> Self {
        let properties = PlanProperties::new(
            EquivalenceProperties::new(schema.clone()),
            Partitioning::UnknownPartitioning(1),
            EmissionType::Final,
            Boundedness::Bounded,
        );

        Self {
            left,
            right,
            left_cols,
            right_cols,
            skip_sort,
            schema,
            properties,
        }
    }
}

impl DisplayAs for SweepLineJoinExec {
    fn fmt_as(
        &self,
        _t: DisplayFormatType,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(f, "SweepLineJoinExec: skip_sort={}", self.skip_sort)
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
            self.skip_sort,
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
        let schema = self.schema.clone();
        let ctx = context;

        let stream = futures::stream::once(async move {
            let left_batches = datafusion::physical_plan::collect(
                left, ctx.clone(),
            )
            .await?;
            let right_batches = datafusion::physical_plan::collect(
                right, ctx,
            )
            .await?;

            let left_concat = concat_batches(&left_batches)?;
            let right_concat = concat_batches(&right_batches)?;

            if left_concat.num_rows() == 0
                || right_concat.num_rows() == 0
            {
                return Ok(RecordBatch::new_empty(schema));
            }

            // Extract typed columns — no String allocation
            let l_chrom_ids = chrom_to_ids(
                left_concat.column(left_cols.chrom_idx),
            )?;
            let l_starts = as_i64_slice(
                left_concat.column(left_cols.start_idx),
            )?;
            let l_ends = as_i64_slice(
                left_concat.column(left_cols.end_idx),
            )?;

            let r_chrom_ids = chrom_to_ids(
                right_concat.column(right_cols.chrom_idx),
            )?;
            let r_starts = as_i64_slice(
                right_concat.column(right_cols.start_idx),
            )?;
            let r_ends = as_i64_slice(
                right_concat.column(right_cols.end_idx),
            )?;

            // Sort indices by (chrom_id, start)
            let l_order = argsort_by_chrom_start(&l_chrom_ids, &l_starts);
            let r_order = argsort_by_chrom_start(&r_chrom_ids, &r_starts);

            // Split at chromosome boundaries
            let l_groups = split_by_chrom(&l_order, &l_chrom_ids);
            let r_groups = split_by_chrom(&r_order, &r_chrom_ids);

            // Parallel sweep per chromosome
            let mut handles = Vec::with_capacity(l_groups.len());
            for (chrom_id, l_range) in &l_groups {
                let r_range = match r_groups
                    .binary_search_by_key(chrom_id, |(c, _)| *c)
                {
                    Ok(pos) => r_groups[pos].1.clone(),
                    Err(_) => continue,
                };

                let l_idx: Arc<[u32]> =
                    l_order[l_range.clone()].into();
                let r_idx: Arc<[u32]> =
                    r_order[r_range].into();
                let ls = l_starts.clone();
                let le = l_ends.clone();
                let rs = r_starts.clone();
                let re = r_ends.clone();

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

            // Vectorized output
            build_output_take(
                &schema,
                &left_concat,
                &right_concat,
                &all_left,
                &all_right,
            )
        });

        Ok(Box::pin(RecordBatchStreamAdapter::new(
            self.schema.clone(),
            stream,
        )))
    }
}

// ── Column extraction (zero-copy where possible) ────────────────

/// Map chromosome strings to dense integer IDs. Uses a sorted unique
/// list so IDs are consistent across left/right.
fn chrom_to_ids(col: &ArrayRef) -> Result<Vec<u32>> {
    let n = col.len();
    // Build string→id map from unique values
    let mut unique: Vec<String> = Vec::new();
    if let Some(arr) = col.as_any().downcast_ref::<StringArray>() {
        for i in 0..n {
            let s = arr.value(i);
            if !unique.contains(&s.to_string()) {
                unique.push(s.to_string());
            }
        }
        unique.sort();
        let map: std::collections::HashMap<&str, u32> = unique
            .iter()
            .enumerate()
            .map(|(i, s)| (s.as_str(), i as u32))
            .collect();
        return Ok((0..n).map(|i| map[arr.value(i)]).collect());
    }
    if let Some(arr) = col.as_any().downcast_ref::<StringViewArray>()
    {
        for i in 0..n {
            let s = arr.value(i);
            if !unique.contains(&s.to_string()) {
                unique.push(s.to_string());
            }
        }
        unique.sort();
        let map: std::collections::HashMap<&str, u32> = unique
            .iter()
            .enumerate()
            .map(|(i, s)| (s.as_str(), i as u32))
            .collect();
        return Ok((0..n).map(|i| map[arr.value(i)]).collect());
    }
    Err(datafusion::error::DataFusionError::Internal(
        "Unsupported string type".to_string(),
    ))
}

/// Get a reference to the i64 values, converting Int32 if needed.
fn as_i64_slice(col: &ArrayRef) -> Result<Arc<[i64]>> {
    if let Some(arr) = col.as_any().downcast_ref::<Int64Array>() {
        return Ok(arr.values().to_vec().into());
    }
    if let Some(arr) = col.as_any().downcast_ref::<Int32Array>() {
        return Ok(
            arr.values().iter().map(|&v| v as i64).collect::<Vec<_>>().into()
        );
    }
    Err(datafusion::error::DataFusionError::Internal(
        "Column is not Int32 or Int64".to_string(),
    ))
}

// ── Sorting and grouping ────────────────────────────────────────

/// Sort indices by (chrom_id, start).
fn argsort_by_chrom_start(
    chrom_ids: &[u32],
    starts: &[i64],
) -> Vec<u32> {
    let n = chrom_ids.len();
    let mut indices: Vec<u32> = (0..n as u32).collect();
    indices.sort_unstable_by(|&a, &b| {
        let a = a as usize;
        let b = b as usize;
        chrom_ids[a]
            .cmp(&chrom_ids[b])
            .then(starts[a].cmp(&starts[b]))
    });
    indices
}

/// Split sorted indices into contiguous chromosome ranges.
/// Returns Vec<(chrom_id, Range<usize>)> sorted by chrom_id.
fn split_by_chrom(
    sorted_indices: &[u32],
    chrom_ids: &[u32],
) -> Vec<(u32, std::ops::Range<usize>)> {
    if sorted_indices.is_empty() {
        return vec![];
    }

    let mut groups = Vec::new();
    let mut start = 0;
    let mut cur_chrom = chrom_ids[sorted_indices[0] as usize];

    for i in 1..sorted_indices.len() {
        let c = chrom_ids[sorted_indices[i] as usize];
        if c != cur_chrom {
            groups.push((cur_chrom, start..i));
            start = i;
            cur_chrom = c;
        }
    }
    groups.push((cur_chrom, start..sorted_indices.len()));
    groups
}

// ── Core sweep ──────────────────────────────────────────────────

/// Sweep-line for one chromosome. Inputs are sorted index slices.
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

        // Add right intervals with start < l_end
        while r_cursor < right_indices.len() {
            let ri = right_indices[r_cursor] as usize;
            if right_starts[ri] >= l_end {
                break;
            }
            active.push(r_cursor);
            r_cursor += 1;
        }

        // Remove expired (end <= l_start)
        active.retain(|&pos| {
            right_ends[right_indices[pos] as usize] > l_start
        });

        // Emit overlapping pairs
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

// ── Output construction ─────────────────────────────────────────

fn concat_batches(batches: &[RecordBatch]) -> Result<RecordBatch> {
    if batches.is_empty() {
        return Err(datafusion::error::DataFusionError::Internal(
            "No batches to concatenate".to_string(),
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

    let mut columns: Vec<ArrayRef> = Vec::with_capacity(
        left.num_columns() + right.num_columns(),
    );

    for c in 0..left.num_columns() {
        columns
            .push(compute::take(left.column(c).as_ref(), &li, None)?);
    }
    for c in 0..right.num_columns() {
        columns.push(
            compute::take(right.column(c).as_ref(), &ri, None)?,
        );
    }

    Ok(RecordBatch::try_new(schema.clone(), columns)?)
}
