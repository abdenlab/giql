use std::any::Any;
use std::fmt;
use std::sync::Arc;

use arrow::array::{
    Array, ArrayRef, Int32Array, Int64Array, RecordBatch, StringArray,
    StringViewArray,
};
use arrow::datatypes::SchemaRef;
use datafusion::common::Result;
use datafusion::execution::SendableRecordBatchStream;
use datafusion::physical_expr::{EquivalenceProperties, Partitioning};
use datafusion::physical_plan::stream::RecordBatchStreamAdapter;
use datafusion::physical_plan::execution_plan::{
    Boundedness, EmissionType,
};
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, PlanProperties,
};

use crate::pattern::IntervalColumns;

/// Custom execution plan implementing the sweep-line interval join.
///
/// Both inputs are sorted by `(chrom, start)`, then swept left to
/// right. For each left interval, all right intervals whose start is
/// less than the left's end are candidates; those whose end is greater
/// than the left's start are matches.
///
/// Complexity: O((n+m) log(n+m) + k) where k is the output size.
/// If `skip_sort` is true, the sort is assumed already done and the
/// complexity is O(n+m+k).
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
        write!(
            f,
            "SweepLineJoinExec: skip_sort={}",
            self.skip_sort
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
        // Collect ALL partitions from both children. DataFusion may
        // split inputs across multiple partitions (default =
        // num_cpus). We must read every partition to get all rows.
        let left = self.left.clone();
        let right = self.right.clone();
        let left_cols = self.left_cols.clone();
        let right_cols = self.right_cols.clone();
        let schema = self.schema.clone();
        let ctx = context;

        let stream = futures::stream::once(async move {
            let left_batches =
                collect_all_partitions(&left, &ctx).await?;
            let right_batches =
                collect_all_partitions(&right, &ctx).await?;

            let left_intervals =
                extract_intervals(&left_batches, &left_cols)?;
            let right_intervals =
                extract_intervals(&right_batches, &right_cols)?;

            let matches = sweep_line_join(
                &left_intervals,
                &right_intervals,
            );

            build_output_batch(
                &schema,
                &left_batches,
                &right_batches,
                &matches,
            )
        });

        Ok(Box::pin(RecordBatchStreamAdapter::new(
            self.schema.clone(),
            stream,
        )))
    }
}

/// A flattened interval with a pointer back to its batch and row.
#[derive(Debug, Clone)]
struct FlatInterval {
    chrom: String,
    start: i64,
    end: i64,
    batch_idx: usize,
    row_idx: usize,
}

/// Extract all intervals from record batches into a flat sorted vec.
fn extract_intervals(
    batches: &[RecordBatch],
    cols: &IntervalColumns,
) -> Result<Vec<FlatInterval>> {
    let mut intervals = Vec::new();

    for (batch_idx, batch) in batches.iter().enumerate() {
        let chrom_col = batch.column(cols.chrom_idx);
        let start_col = batch.column(cols.start_idx);
        let end_col = batch.column(cols.end_idx);

        for row_idx in 0..batch.num_rows() {
            if chrom_col.is_null(row_idx)
                || start_col.is_null(row_idx)
                || end_col.is_null(row_idx)
            {
                continue;
            }
            let chrom = get_string_value(chrom_col.as_ref(), row_idx)
                .ok_or_else(|| {
                    datafusion::error::DataFusionError::Internal(
                        "Chrom column has unsupported string type"
                            .to_string(),
                    )
                })?;
            let start = get_i64_value(start_col.as_ref(), row_idx)
                .ok_or_else(|| {
                    datafusion::error::DataFusionError::Internal(
                        "Start column is not Int32 or Int64".to_string(),
                    )
                })?;
            let end = get_i64_value(end_col.as_ref(), row_idx)
                .ok_or_else(|| {
                    datafusion::error::DataFusionError::Internal(
                        "End column is not Int32 or Int64".to_string(),
                    )
                })?;
            intervals.push(FlatInterval {
                chrom,
                start,
                end,
                batch_idx,
                row_idx,
            });
        }
    }

    // Sort by (chrom, start)
    intervals.sort_by(|a, b| {
        a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start))
    });

    Ok(intervals)
}

/// Core sweep-line algorithm.
///
/// Both inputs must be sorted by (chrom, start). For each chromosome,
/// maintains an active set of right intervals and sweeps left to right.
fn sweep_line_join(
    left: &[FlatInterval],
    right: &[FlatInterval],
) -> Vec<(usize, usize, usize, usize)> {
    // (left_batch, left_row, right_batch, right_row)
    let mut matches = Vec::new();

    let mut right_idx = 0;
    let mut active: Vec<usize> = Vec::new(); // indices into right

    for l in left {

        // Advance right_idx to add right intervals that START before
        // ANY future left interval could end. Since left is sorted by
        // start, and a future left interval's end could be arbitrarily
        // large, we add all right intervals with start < l.end.
        // This is conservative — we check exact overlap below.
        while right_idx < right.len() {
            let r = &right[right_idx];
            if r.chrom < l.chrom {
                right_idx += 1;
                continue;
            }
            if r.chrom > l.chrom {
                break;
            }
            if r.start >= l.end {
                break;
            }
            active.push(right_idx);
            right_idx += 1;
        }

        // Remove only truly expired intervals: those whose END is at
        // or before the current left START. Since left is sorted by
        // start, any interval with end <= l.start can never overlap
        // any future left interval either.
        //
        // We do NOT filter on r.start < l.end here because a wide
        // right interval (r.start before a previous narrow left's
        // end) may still be needed by a LATER wider left interval.
        active.retain(|&ri| {
            let r = &right[ri];
            r.chrom == l.chrom && r.end > l.start
        });

        // Emit matches: check the full overlap condition inline.
        for &ri in &active {
            let r = &right[ri];
            if r.start < l.end {
                matches.push((
                    l.batch_idx,
                    l.row_idx,
                    r.batch_idx,
                    r.row_idx,
                ));
            }
        }
    }

    matches
}

/// Build the output RecordBatch from matched pairs.
fn build_output_batch(
    schema: &SchemaRef,
    left_batches: &[RecordBatch],
    right_batches: &[RecordBatch],
    matches: &[(usize, usize, usize, usize)],
) -> Result<RecordBatch> {
    if matches.is_empty() {
        return Ok(RecordBatch::new_empty(schema.clone()));
    }

    let left_schema = left_batches[0].schema();
    let right_schema = right_batches[0].schema();
    let num_left_cols = left_schema.fields().len();
    let num_right_cols = right_schema.fields().len();

    let mut columns: Vec<ArrayRef> =
        Vec::with_capacity(num_left_cols + num_right_cols);

    for col_idx in 0..num_left_cols {
        let values: Vec<ArrayRef> = matches
            .iter()
            .map(|&(lb, lr, _, _)| {
                left_batches[lb].column(col_idx).slice(lr, 1)
            })
            .collect();

        let refs: Vec<&dyn Array> =
            values.iter().map(|a| a.as_ref()).collect();
        columns.push(arrow::compute::concat(&refs)?);
    }

    for col_idx in 0..num_right_cols {
        let values: Vec<ArrayRef> = matches
            .iter()
            .map(|&(_, _, rb, rr)| {
                right_batches[rb].column(col_idx).slice(rr, 1)
            })
            .collect();

        let refs: Vec<&dyn Array> =
            values.iter().map(|a| a.as_ref()).collect();
        columns.push(arrow::compute::concat(&refs)?);
    }

    Ok(RecordBatch::try_new(schema.clone(), columns)?)
}

/// Collect all record batches from all partitions of an execution
/// plan. Uses DataFusion's `collect` which spawns partition tasks
/// concurrently — required because RepartitionExec uses shared
/// channels that break under sequential execution.
async fn collect_all_partitions(
    plan: &Arc<dyn ExecutionPlan>,
    context: &Arc<datafusion::execution::TaskContext>,
) -> Result<Vec<RecordBatch>> {
    datafusion::physical_plan::collect(
        plan.clone(),
        context.clone(),
    )
    .await
}

/// Extract an i64 value from an array that may be Int32Array or
/// Int64Array.
fn get_i64_value(array: &dyn Array, idx: usize) -> Option<i64> {
    array
        .as_any()
        .downcast_ref::<Int64Array>()
        .map(|arr| arr.value(idx))
        .or_else(|| {
            array
                .as_any()
                .downcast_ref::<Int32Array>()
                .map(|arr| arr.value(idx) as i64)
        })
}

/// Extract a string value from an array that may be StringArray or
/// StringViewArray (DataFusion v47+ uses StringViewArray by default).
fn get_string_value(
    array: &dyn Array,
    idx: usize,
) -> Option<String> {
    array
        .as_any()
        .downcast_ref::<StringArray>()
        .map(|arr| arr.value(idx).to_string())
        .or_else(|| {
            array
                .as_any()
                .downcast_ref::<StringViewArray>()
                .map(|arr| arr.value(idx).to_string())
        })
}
