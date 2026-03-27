use std::any::Any;
use std::fmt;
use std::sync::Arc;

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

/// Binned interval join that delegates to DataFusion's SQL engine.
///
/// Collects both inputs, registers them as temporary tables, and
/// executes the binned equi-join as a SQL query through a fresh
/// SessionContext. This lets DataFusion's built-in UNNEST and
/// DISTINCT optimizations handle the bin expansion and dedup,
/// matching the performance of hand-written SQL.
#[derive(Debug)]
pub struct BinnedSqlExec {
    left: Arc<dyn ExecutionPlan>,
    right: Arc<dyn ExecutionPlan>,
    left_cols: IntervalColumns,
    right_cols: IntervalColumns,
    bin_size: usize,
    schema: SchemaRef,
    properties: Arc<PlanProperties>,
}

impl BinnedSqlExec {
    pub fn new(
        left: Arc<dyn ExecutionPlan>,
        right: Arc<dyn ExecutionPlan>,
        left_cols: IntervalColumns,
        right_cols: IntervalColumns,
        schema: SchemaRef,
        bin_size: usize,
    ) -> Self {
        let properties = Arc::new(PlanProperties::new(
            EquivalenceProperties::new(schema.clone()),
            Partitioning::UnknownPartitioning(1),
            EmissionType::Final,
            Boundedness::Bounded,
        ));

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

    fn build_sql(&self) -> String {
        let b = self.bin_size;
        let lc = &self.left_cols;
        let rc = &self.right_cols;

        format!(
            "WITH __giql_left AS (\
             SELECT *, UNNEST(range(\
               CAST(\"{ls}\" / {b} AS BIGINT), \
               CAST((\"{le}\" - 1) / {b} + 1 AS BIGINT)\
             )) AS __giql_bin FROM __giql_a), \
             __giql_right AS (\
             SELECT *, UNNEST(range(\
               CAST(\"{rs}\" / {b} AS BIGINT), \
               CAST((\"{re}\" - 1) / {b} + 1 AS BIGINT)\
             )) AS __giql_bin FROM __giql_b) \
             SELECT DISTINCT \
               l.\"{lch}\", l.\"{ls}\", l.\"{le}\", \
               r.\"{rch}\" AS \"{rch}\", \
               r.\"{rs}\" AS \"{rs}\", \
               r.\"{re}\" AS \"{re}\" \
             FROM __giql_left AS l \
             JOIN __giql_right AS r \
             ON l.\"{lch}\" = r.\"{rch}\" \
             AND l.__giql_bin = r.__giql_bin \
             WHERE l.\"{ls}\" < r.\"{re}\" \
             AND l.\"{le}\" > r.\"{rs}\"",
            b = b,
            lch = lc.chrom_col,
            ls = lc.start_col,
            le = lc.end_col,
            rch = rc.chrom_col,
            rs = rc.start_col,
            re = rc.end_col,
        )
    }
}

impl DisplayAs for BinnedSqlExec {
    fn fmt_as(
        &self,
        _t: DisplayFormatType,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(f, "BinnedSqlExec: bin_size={}", self.bin_size)
    }
}

impl ExecutionPlan for BinnedSqlExec {
    fn name(&self) -> &str {
        "BinnedSqlExec"
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
        Ok(Arc::new(BinnedSqlExec::new(
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
        _partition: usize,
        context: Arc<datafusion::execution::TaskContext>,
    ) -> Result<SendableRecordBatchStream> {
        let left_plan = self.left.clone();
        let right_plan = self.right.clone();
        let sql = self.build_sql();
        let schema = self.schema.clone();

        let stream = futures::stream::once(async move {
            use datafusion::prelude::SessionContext;

            // Collect both sides
            let left_batches = datafusion::physical_plan::collect(
                left_plan,
                context.clone(),
            )
            .await?;
            let right_batches = datafusion::physical_plan::collect(
                right_plan, context,
            )
            .await?;

            // Create a fresh context and register as memory tables
            let ctx = SessionContext::new();
            let left_table =
                datafusion::datasource::MemTable::try_new(
                    left_batches[0].schema(),
                    vec![left_batches],
                )?;
            let right_table =
                datafusion::datasource::MemTable::try_new(
                    right_batches[0].schema(),
                    vec![right_batches],
                )?;
            ctx.register_table("__giql_a", Arc::new(left_table))?;
            ctx.register_table("__giql_b", Arc::new(right_table))?;

            // Execute the binned SQL
            let df = ctx.sql(&sql).await?;
            let batches = df.collect().await?;

            if batches.is_empty() {
                return Ok(
                    arrow::record_batch::RecordBatch::new_empty(schema),
                );
            }

            // Concat all result batches
            Ok(arrow::compute::concat_batches(
                &batches[0].schema(),
                &batches,
            )?)
        });

        Ok(Box::pin(RecordBatchStreamAdapter::new(
            self.schema.clone(),
            stream,
        )))
    }
}
