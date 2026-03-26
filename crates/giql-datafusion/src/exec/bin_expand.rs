use std::any::Any;
use std::fmt;
use std::sync::Arc;

use arrow::array::{
    Array, ArrayRef, Int32Array, Int64Array, RecordBatch,
    UInt64Array,
};
use arrow::datatypes::{DataType, Field, SchemaRef};
use datafusion::common::Result;
use datafusion::execution::SendableRecordBatchStream;
use datafusion::physical_expr::EquivalenceProperties;
use datafusion::physical_plan::execution_plan::{
    Boundedness, EmissionType,
};
use datafusion::physical_plan::stream::RecordBatchStreamAdapter;
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, Partitioning,
    PlanProperties,
};

/// Stateless per-partition exec that expands each interval into rows
/// for each genome bin it touches.
///
/// Appends two columns:
/// - `__giql_bin`: the bin ID for this expanded row
/// - `__giql_first_bin`: `start / bin_size` — the interval's first bin
///
/// The `__giql_first_bin` column enables the canonical-bin
/// deduplication trick: after the hash join, a filter keeps only
/// matches where `__giql_bin == max(left.__giql_first_bin,
/// right.__giql_first_bin)`, ensuring each pair is emitted exactly
/// once.
#[derive(Debug)]
pub struct BinExpandExec {
    input: Arc<dyn ExecutionPlan>,
    start_col_idx: usize,
    end_col_idx: usize,
    bin_size: usize,
    schema: SchemaRef,
    properties: Arc<PlanProperties>,
}

impl BinExpandExec {
    pub fn new(
        input: Arc<dyn ExecutionPlan>,
        start_col_idx: usize,
        end_col_idx: usize,
        bin_size: usize,
    ) -> Self {
        let input_schema = input.schema();
        let mut fields: Vec<Arc<Field>> =
            input_schema.fields().iter().cloned().collect();
        fields.push(Arc::new(Field::new(
            "__giql_bin",
            DataType::Int64,
            false,
        )));
        fields.push(Arc::new(Field::new(
            "__giql_first_bin",
            DataType::Int64,
            false,
        )));
        let schema =
            Arc::new(arrow::datatypes::Schema::new(fields));

        let properties = Arc::new(PlanProperties::new(
            EquivalenceProperties::new(schema.clone()),
            input.properties().partitioning.clone(),
            EmissionType::Incremental,
            Boundedness::Bounded,
        ));

        Self {
            input,
            start_col_idx,
            end_col_idx,
            bin_size,
            schema,
            properties,
        }
    }

    /// Number of columns added by this exec (bin + first_bin).
    pub const EXTRA_COLS: usize = 2;
}

impl DisplayAs for BinExpandExec {
    fn fmt_as(
        &self,
        _t: DisplayFormatType,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(
            f,
            "BinExpandExec: start_col={}, end_col={}, bin_size={}",
            self.start_col_idx, self.end_col_idx, self.bin_size
        )
    }
}

impl ExecutionPlan for BinExpandExec {
    fn name(&self) -> &str {
        "BinExpandExec"
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
        vec![&self.input]
    }

    fn with_new_children(
        self: Arc<Self>,
        children: Vec<Arc<dyn ExecutionPlan>>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        Ok(Arc::new(BinExpandExec::new(
            children[0].clone(),
            self.start_col_idx,
            self.end_col_idx,
            self.bin_size,
        )))
    }

    fn execute(
        &self,
        partition: usize,
        context: Arc<datafusion::execution::TaskContext>,
    ) -> Result<SendableRecordBatchStream> {
        let input_stream =
            self.input.execute(partition, context)?;
        let start_idx = self.start_col_idx;
        let end_idx = self.end_col_idx;
        let bin_size = self.bin_size as i64;
        let schema = self.schema.clone();

        let output_stream = futures::stream::unfold(
            input_stream,
            move |mut stream| {
                let schema = schema.clone();
                async move {
                    use futures::StreamExt;
                    match stream.next().await {
                        Some(Ok(batch)) => {
                            let result = expand_batch(
                                &batch, start_idx, end_idx,
                                bin_size, &schema,
                            );
                            Some((result, stream))
                        }
                        Some(Err(e)) => Some((Err(e), stream)),
                        None => None,
                    }
                }
            },
        );

        Ok(Box::pin(RecordBatchStreamAdapter::new(
            self.schema.clone(),
            output_stream,
        )))
    }
}

/// Expand a single batch: replicate each row for each bin it touches,
/// adding `__giql_bin` and `__giql_first_bin` columns.
fn expand_batch(
    batch: &RecordBatch,
    start_idx: usize,
    end_idx: usize,
    bin_size: i64,
    output_schema: &SchemaRef,
) -> Result<RecordBatch> {
    let num_rows = batch.num_rows();
    if num_rows == 0 {
        return Ok(RecordBatch::new_empty(output_schema.clone()));
    }

    let start_col = batch.column(start_idx);
    let end_col = batch.column(end_idx);

    let mut bin_ranges: Vec<(i64, i64)> = Vec::with_capacity(num_rows);
    let mut total_expanded = 0usize;

    for row in 0..num_rows {
        let start = get_i64(start_col.as_ref(), row).unwrap_or(0);
        let end = get_i64(end_col.as_ref(), row).unwrap_or(0);
        let first_bin = start / bin_size;
        let last_bin = (end - 1).max(0) / bin_size;
        let n_bins = (last_bin - first_bin + 1).max(0) as usize;
        bin_ranges.push((first_bin, last_bin));
        total_expanded += n_bins;
    }

    let mut row_indices: Vec<u64> =
        Vec::with_capacity(total_expanded);
    let mut bin_ids: Vec<i64> = Vec::with_capacity(total_expanded);
    let mut first_bins: Vec<i64> =
        Vec::with_capacity(total_expanded);

    for (row, &(first_bin, last_bin)) in
        bin_ranges.iter().enumerate()
    {
        for bin in first_bin..=last_bin {
            row_indices.push(row as u64);
            bin_ids.push(bin);
            first_bins.push(first_bin);
        }
    }

    let row_indices_arr = UInt64Array::from(row_indices);
    let mut columns: Vec<ArrayRef> =
        Vec::with_capacity(batch.num_columns() + BinExpandExec::EXTRA_COLS);

    for col_idx in 0..batch.num_columns() {
        let col = batch.column(col_idx);
        let taken =
            arrow::compute::take(col.as_ref(), &row_indices_arr, None)?;
        columns.push(taken);
    }

    columns.push(Arc::new(Int64Array::from(bin_ids)));
    columns.push(Arc::new(Int64Array::from(first_bins)));

    Ok(RecordBatch::try_new(output_schema.clone(), columns)?)
}

fn get_i64(array: &dyn Array, idx: usize) -> Option<i64> {
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
