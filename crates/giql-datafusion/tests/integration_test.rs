//! Integration tests for the INTERSECTS join optimizer.
//!
//! These tests exercise the full pipeline: create Parquet files →
//! register with DataFusion → add optimizer rule → execute
//! INTERSECTS join SQL → verify results.

use std::path::Path;
use std::sync::Arc;

use arrow::array::{Int64Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use datafusion::execution::SessionStateBuilder;
use datafusion::prelude::*;
use parquet::arrow::ArrowWriter;
use tempfile::TempDir;

use giql_datafusion::{IntersectsOptimizerConfig, register_optimizer};

/// Write a Parquet file with the given genomic intervals.
fn write_intervals_parquet(
    dir: &Path,
    filename: &str,
    chroms: &[&str],
    starts: &[i64],
    ends: &[i64],
) -> std::path::PathBuf {
    let schema = Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::Int64, false),
        Field::new("end", DataType::Int64, false),
    ]));

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(StringArray::from(chroms.to_vec())),
            Arc::new(Int64Array::from(starts.to_vec())),
            Arc::new(Int64Array::from(ends.to_vec())),
        ],
    )
    .unwrap();

    let path = dir.join(filename);
    let file = std::fs::File::create(&path).unwrap();
    let mut writer =
        ArrowWriter::try_new(file, schema, None).unwrap();
    writer.write(&batch).unwrap();
    writer.close().unwrap();
    path
}

/// Create a SessionContext with the INTERSECTS optimizer registered.
fn make_ctx_with_optimizer() -> SessionContext {
    let config = IntersectsOptimizerConfig::default();
    let state = SessionStateBuilder::new()
        .with_default_features()
        .build();
    let state = register_optimizer(state, config);
    SessionContext::from(state)
}

/// The standard INTERSECTS join SQL using the chrom/start/end
/// predicate pattern that the optimizer detects.
const INTERSECTS_SQL: &str = "\
    SELECT a.chrom, a.start, a.\"end\", \
           b.chrom AS chrom_b, b.start AS start_b, b.\"end\" AS end_b \
    FROM a JOIN b \
    ON a.chrom = b.chrom \
    AND a.start < b.\"end\" \
    AND a.\"end\" > b.start";

// ── Correctness tests ──────────────────────────────────────────

#[tokio::test]
async fn test_overlapping_intervals_returns_pairs() {
    let dir = TempDir::new().unwrap();
    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &["chr1", "chr1", "chr1"],
        &[100, 300, 600],
        &[250, 500, 800],
    );
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &["chr1", "chr1"],
        &[200, 700],
        &[400, 900],
    );

    let ctx = make_ctx_with_optimizer();
    ctx.register_parquet("a", path_a.to_str().unwrap(), Default::default())
        .await
        .unwrap();
    ctx.register_parquet("b", path_b.to_str().unwrap(), Default::default())
        .await
        .unwrap();

    let result = ctx.sql(INTERSECTS_SQL).await.unwrap();
    let batches = result.collect().await.unwrap();
    let total_rows: usize =
        batches.iter().map(|b| b.num_rows()).sum();

    // Expected overlaps:
    // a[100,250) overlaps b[200,400) ✓
    // a[300,500) overlaps b[200,400) ✓
    // a[600,800) overlaps b[700,900) ✓
    assert_eq!(total_rows, 3);
}

#[tokio::test]
async fn test_no_overlapping_intervals_returns_empty() {
    let dir = TempDir::new().unwrap();
    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &["chr1", "chr1"],
        &[100, 300],
        &[200, 400],
    );
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &["chr1", "chr1"],
        &[500, 700],
        &[600, 800],
    );

    let ctx = make_ctx_with_optimizer();
    ctx.register_parquet("a", path_a.to_str().unwrap(), Default::default())
        .await
        .unwrap();
    ctx.register_parquet("b", path_b.to_str().unwrap(), Default::default())
        .await
        .unwrap();

    let result = ctx.sql(INTERSECTS_SQL).await.unwrap();
    let batches = result.collect().await.unwrap();
    let total_rows: usize =
        batches.iter().map(|b| b.num_rows()).sum();

    assert_eq!(total_rows, 0);
}

#[tokio::test]
async fn test_different_chromosomes_returns_empty() {
    let dir = TempDir::new().unwrap();
    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &["chr1", "chr1"],
        &[100, 300],
        &[500, 600],
    );
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &["chr2", "chr2"],
        &[100, 300],
        &[500, 600],
    );

    let ctx = make_ctx_with_optimizer();
    ctx.register_parquet("a", path_a.to_str().unwrap(), Default::default())
        .await
        .unwrap();
    ctx.register_parquet("b", path_b.to_str().unwrap(), Default::default())
        .await
        .unwrap();

    let result = ctx.sql(INTERSECTS_SQL).await.unwrap();
    let batches = result.collect().await.unwrap();
    let total_rows: usize =
        batches.iter().map(|b| b.num_rows()).sum();

    assert_eq!(total_rows, 0);
}

#[tokio::test]
async fn test_adjacent_intervals_no_overlap() {
    // Half-open interval semantics: [100,200) and [200,300) do NOT
    // overlap because 200 is not < 200.
    let dir = TempDir::new().unwrap();
    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &["chr1"],
        &[100],
        &[200],
    );
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &["chr1"],
        &[200],
        &[300],
    );

    let ctx = make_ctx_with_optimizer();
    ctx.register_parquet("a", path_a.to_str().unwrap(), Default::default())
        .await
        .unwrap();
    ctx.register_parquet("b", path_b.to_str().unwrap(), Default::default())
        .await
        .unwrap();

    let result = ctx.sql(INTERSECTS_SQL).await.unwrap();
    let batches = result.collect().await.unwrap();
    let total_rows: usize =
        batches.iter().map(|b| b.num_rows()).sum();

    assert_eq!(total_rows, 0);
}

#[tokio::test]
async fn test_containment_counts_as_overlap() {
    // [100,500) fully contains [200,300)
    let dir = TempDir::new().unwrap();
    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &["chr1"],
        &[100],
        &[500],
    );
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &["chr1"],
        &[200],
        &[300],
    );

    let ctx = make_ctx_with_optimizer();
    ctx.register_parquet("a", path_a.to_str().unwrap(), Default::default())
        .await
        .unwrap();
    ctx.register_parquet("b", path_b.to_str().unwrap(), Default::default())
        .await
        .unwrap();

    let result = ctx.sql(INTERSECTS_SQL).await.unwrap();
    let batches = result.collect().await.unwrap();
    let total_rows: usize =
        batches.iter().map(|b| b.num_rows()).sum();

    assert_eq!(total_rows, 1);
}

// ── Optimizer registration test ────────────────────────────────

#[tokio::test]
async fn test_register_optimizer_present_in_state() {
    let state = SessionStateBuilder::new()
        .with_default_features()
        .build();
    let n_before = state.physical_optimizers().len();

    let config = IntersectsOptimizerConfig::default();
    let state = register_optimizer(state, config);

    assert_eq!(state.physical_optimizers().len(), n_before + 1);
    assert_eq!(
        state.physical_optimizers().last().unwrap().name(),
        "intersects_optimizer"
    );
}
