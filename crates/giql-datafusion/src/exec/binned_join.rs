use std::any::Any;
use std::collections::HashMap;
use std::fmt;
use std::sync::Arc;

use arrow::array::{Array, ArrayRef, Int64Array, RecordBatch, StringArray};
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

/// Custom execution plan implementing the binned equi-join for
/// interval overlap.
///
/// Each interval is expanded into genome bins of fixed size. A hash
/// map is built from the right side keyed by `(chrom, bin_id)`. The
/// left side probes the map, and a post-filter removes false positives.
/// Output is deduplicated to avoid emitting duplicate pairs when an
/// interval spans multiple bins.
///
/// Complexity: O((n+m) * avg_replication + k) where avg_replication
/// is mean_width / bin_size + 1.
#[derive(Debug)]
pub struct BinnedJoinExec {
    left: Arc<dyn ExecutionPlan>,
    right: Arc<dyn ExecutionPlan>,
    left_cols: IntervalColumns,
    right_cols: IntervalColumns,
    bin_size: usize,
    schema: SchemaRef,
    properties: PlanProperties,
}

impl BinnedJoinExec {
    pub fn new(
        left: Arc<dyn ExecutionPlan>,
        right: Arc<dyn ExecutionPlan>,
        left_cols: IntervalColumns,
        right_cols: IntervalColumns,
        schema: SchemaRef,
        bin_size: usize,
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
            bin_size,
            schema,
            properties,
        }
    }
}

impl DisplayAs for BinnedJoinExec {
    fn fmt_as(
        &self,
        _t: DisplayFormatType,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(f, "BinnedJoinExec: bin_size={}", self.bin_size)
    }
}

impl ExecutionPlan for BinnedJoinExec {
    fn name(&self) -> &str {
        "BinnedJoinExec"
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
        Ok(Arc::new(BinnedJoinExec::new(
            children[0].clone(),
            children[1].clone(),
            self.left_cols.clone(),
            self.right_cols.clone(),
            self.schema.clone(),
            self.bin_size,
        )))
    }

    fn execute(
        &self,
        partition: usize,
        context: Arc<datafusion::execution::TaskContext>,
    ) -> Result<SendableRecordBatchStream> {
        let left_stream =
            self.left.execute(partition, context.clone())?;
        let right_stream = self.right.execute(partition, context)?;

        let left_cols = self.left_cols.clone();
        let right_cols = self.right_cols.clone();
        let schema = self.schema.clone();
        let bin_size = self.bin_size;

        let stream = futures::stream::once(async move {
            let left_batches = collect_batches(left_stream).await?;
            let right_batches =
                collect_batches(right_stream).await?;

            binned_join(
                &schema,
                &left_batches,
                &right_batches,
                &left_cols,
                &right_cols,
                bin_size,
            )
        });

        Ok(Box::pin(RecordBatchStreamAdapter::new(
            self.schema.clone(),
            stream,
        )))
    }
}

/// A reference to a specific row in a batch.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
struct RowRef {
    batch_idx: usize,
    row_idx: usize,
}

/// Interval data extracted from a row.
struct IntervalRow {
    chrom: String,
    start: i64,
    end: i64,
    row_ref: RowRef,
}

/// Core binned join algorithm.
fn binned_join(
    schema: &SchemaRef,
    left_batches: &[RecordBatch],
    right_batches: &[RecordBatch],
    left_cols: &IntervalColumns,
    right_cols: &IntervalColumns,
    bin_size: usize,
) -> Result<RecordBatch> {
    let right_intervals =
        extract_interval_rows(right_batches, right_cols)?;
    let mut right_map: HashMap<(String, i64), Vec<usize>> =
        HashMap::new();

    for (idx, interval) in right_intervals.iter().enumerate() {
        let start_bin = interval.start / bin_size as i64;
        let end_bin = (interval.end - 1) / bin_size as i64;
        for bin in start_bin..=end_bin {
            right_map
                .entry((interval.chrom.clone(), bin))
                .or_default()
                .push(idx);
        }
    }

    let left_intervals =
        extract_interval_rows(left_batches, left_cols)?;

    let mut seen = std::collections::HashSet::new();
    let mut matches: Vec<(RowRef, RowRef)> = Vec::new();

    for (li, l) in left_intervals.iter().enumerate() {
        let start_bin = l.start / bin_size as i64;
        let end_bin = (l.end - 1) / bin_size as i64;

        for bin in start_bin..=end_bin {
            let key = (l.chrom.clone(), bin);
            if let Some(right_indices) = right_map.get(&key) {
                for &ri in right_indices {
                    if seen.contains(&(li, ri)) {
                        continue;
                    }

                    let r = &right_intervals[ri];
                    if l.start < r.end && l.end > r.start {
                        seen.insert((li, ri));
                        matches.push((
                            l.row_ref.clone(),
                            r.row_ref.clone(),
                        ));
                    }
                }
            }
        }
    }

    build_output(schema, left_batches, right_batches, &matches)
}

/// Extract interval rows from batches.
fn extract_interval_rows(
    batches: &[RecordBatch],
    cols: &IntervalColumns,
) -> Result<Vec<IntervalRow>> {
    let mut rows = Vec::new();

    for (batch_idx, batch) in batches.iter().enumerate() {
        let chroms = batch
            .column(cols.chrom_idx)
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| {
                datafusion::error::DataFusionError::Internal(
                    "Chrom column is not StringArray".to_string(),
                )
            })?;

        let starts = batch
            .column(cols.start_idx)
            .as_any()
            .downcast_ref::<Int64Array>()
            .ok_or_else(|| {
                datafusion::error::DataFusionError::Internal(
                    "Start column is not Int64Array".to_string(),
                )
            })?;

        let ends = batch
            .column(cols.end_idx)
            .as_any()
            .downcast_ref::<Int64Array>()
            .ok_or_else(|| {
                datafusion::error::DataFusionError::Internal(
                    "End column is not Int64Array".to_string(),
                )
            })?;

        for row_idx in 0..batch.num_rows() {
            if chroms.is_null(row_idx)
                || starts.is_null(row_idx)
                || ends.is_null(row_idx)
            {
                continue;
            }
            rows.push(IntervalRow {
                chrom: chroms.value(row_idx).to_string(),
                start: starts.value(row_idx),
                end: ends.value(row_idx),
                row_ref: RowRef {
                    batch_idx,
                    row_idx,
                },
            });
        }
    }

    Ok(rows)
}

/// Build output RecordBatch from matched row pairs.
fn build_output(
    schema: &SchemaRef,
    left_batches: &[RecordBatch],
    right_batches: &[RecordBatch],
    matches: &[(RowRef, RowRef)],
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
            .map(|(lr, _)| {
                left_batches[lr.batch_idx]
                    .column(col_idx)
                    .slice(lr.row_idx, 1)
            })
            .collect();

        let refs: Vec<&dyn Array> =
            values.iter().map(|a| a.as_ref()).collect();
        columns.push(arrow::compute::concat(&refs)?);
    }

    for col_idx in 0..num_right_cols {
        let values: Vec<ArrayRef> = matches
            .iter()
            .map(|(_, rr)| {
                right_batches[rr.batch_idx]
                    .column(col_idx)
                    .slice(rr.row_idx, 1)
            })
            .collect();

        let refs: Vec<&dyn Array> =
            values.iter().map(|a| a.as_ref()).collect();
        columns.push(arrow::compute::concat(&refs)?);
    }

    Ok(RecordBatch::try_new(schema.clone(), columns)?)
}

/// Collect all batches from a stream.
async fn collect_batches(
    stream: SendableRecordBatchStream,
) -> Result<Vec<RecordBatch>> {
    use futures::StreamExt;

    let mut batches = Vec::new();
    let mut stream = stream;

    while let Some(batch) = stream.next().await {
        batches.push(batch?);
    }

    Ok(batches)
}
