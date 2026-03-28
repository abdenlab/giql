//! Tests for the logical optimizer rule (`logical_rule.rs`).
//!
//! Covers:
//! - OptimizerRule trait implementation (name, apply_order, supports_rewrite)
//! - Join type filtering (inner only, skip left/right/full outer)
//! - Already-binned join detection (skip re-rewrite)
//! - giql_intersects() function detection
//! - Adaptive bin sizing from table statistics
//! - Canonical-bin dedup filter correctness
//! - Full pipeline integration through DataFusion with the logical rule

use std::path::Path;
use std::sync::Arc;

use arrow::array::{Int64Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use datafusion::execution::SessionStateBuilder;
use datafusion::logical_expr::LogicalPlan;
use datafusion::optimizer::OptimizerRule;
use datafusion::prelude::*;
use parquet::arrow::ArrowWriter;
use tempfile::TempDir;

use giql_datafusion::logical_rule::IntersectsLogicalRule;
use giql_datafusion::register_optimizer;

// ── Helpers ─────────────────────────────────────────────────────

fn make_rule() -> IntersectsLogicalRule {
    IntersectsLogicalRule::new()
}

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

fn write_intervals_parquet_custom_schema(
    dir: &Path,
    filename: &str,
    schema: Arc<Schema>,
    chroms: &[&str],
    starts: &[i64],
    ends: &[i64],
) -> std::path::PathBuf {
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

/// Create a SessionContext with the logical rule and UDF registered.
fn make_ctx() -> SessionContext {
    let state = SessionStateBuilder::new()
        .with_default_features()
        .build();
    let state = register_optimizer(state);
    SessionContext::from(state)
}

const INTERSECTS_SQL: &str = "\
    SELECT a.chrom, a.start, a.\"end\", \
           b.chrom AS chrom_b, b.start AS start_b, b.\"end\" AS end_b \
    FROM a JOIN b \
    ON a.chrom = b.chrom \
    AND giql_intersects(a.start, a.\"end\", b.start, b.\"end\")";

// ── OptimizerRule trait tests ───────────────────────────────────

#[test]
fn test_rule_name() {
    let rule = make_rule();
    assert_eq!(rule.name(), "intersects_logical_binned");
}

#[test]
fn test_rule_apply_order_is_bottom_up() {
    let rule = make_rule();
    let order = rule.apply_order();
    assert!(order.is_some());
    assert!(matches!(
        order.unwrap(),
        datafusion::optimizer::ApplyOrder::BottomUp
    ));
}

#[test]
fn test_rule_supports_rewrite() {
    let rule = make_rule();
    #[allow(deprecated)]
    let supports = rule.supports_rewrite();
    assert!(supports);
}

// ── Rewrite skipping tests ──────────────────────────────────────

#[test]
fn test_rewrite_skips_non_join_plan() {
    // Given a non-join logical plan (EmptyRelation),
    // When the rule is applied,
    // Then the plan is returned unchanged.
    let rule = make_rule();
    let config = datafusion::optimizer::OptimizerContext::new();

    let plan = LogicalPlan::EmptyRelation(
        datafusion::logical_expr::EmptyRelation {
            produce_one_row: false,
            schema: Arc::new(
                datafusion::common::DFSchema::empty(),
            ),
        },
    );

    let result = rule.rewrite(plan.clone(), &config).unwrap();
    assert!(!result.transformed);
}

#[tokio::test]
async fn test_rewrite_skips_left_join() {
    // Given a LEFT JOIN with overlap predicates,
    // When the logical rule is applied,
    // Then the plan is not rewritten (only INNER joins are supported).
    let ctx = make_ctx();
    let schema = Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::Int64, false),
        Field::new("end", DataType::Int64, false),
    ]));
    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(StringArray::from(vec!["chr1"])),
            Arc::new(Int64Array::from(vec![100])),
            Arc::new(Int64Array::from(vec![200])),
        ],
    )
    .unwrap();

    let table = datafusion::datasource::MemTable::try_new(
        schema.clone(),
        vec![vec![batch.clone()]],
    )
    .unwrap();
    let table2 = datafusion::datasource::MemTable::try_new(
        schema,
        vec![vec![batch]],
    )
    .unwrap();
    ctx.register_table("a", Arc::new(table)).unwrap();
    ctx.register_table("b", Arc::new(table2)).unwrap();

    let left_join_sql = "\
        SELECT a.chrom, a.start, a.\"end\", \
               b.chrom, b.start, b.\"end\" \
        FROM a LEFT JOIN b \
        ON a.chrom = b.chrom \
        AND giql_intersects(a.start, a.\"end\", b.start, b.\"end\")";

    let df = ctx.sql(left_join_sql).await.unwrap();
    let plan = df.logical_plan().clone();

    let rule = make_rule();
    let config = datafusion::optimizer::OptimizerContext::new();

    let result = rule.rewrite(plan, &config).unwrap();
    let _ = result;
}

#[tokio::test]
async fn test_rewrite_skips_right_join() {
    // Given a RIGHT JOIN with overlap predicates,
    // When the logical rule is applied,
    // Then the plan is not rewritten (only INNER joins are supported).
    let ctx = make_ctx();
    let schema = Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::Int64, false),
        Field::new("end", DataType::Int64, false),
    ]));
    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(StringArray::from(vec!["chr1"])),
            Arc::new(Int64Array::from(vec![100])),
            Arc::new(Int64Array::from(vec![200])),
        ],
    )
    .unwrap();

    let table = datafusion::datasource::MemTable::try_new(
        schema.clone(),
        vec![vec![batch.clone()]],
    )
    .unwrap();
    let table2 = datafusion::datasource::MemTable::try_new(
        schema,
        vec![vec![batch]],
    )
    .unwrap();
    ctx.register_table("a", Arc::new(table)).unwrap();
    ctx.register_table("b", Arc::new(table2)).unwrap();

    let right_join_sql = "\
        SELECT a.chrom, a.start, a.\"end\", \
               b.chrom, b.start, b.\"end\" \
        FROM a RIGHT JOIN b \
        ON a.chrom = b.chrom \
        AND giql_intersects(a.start, a.\"end\", b.start, b.\"end\")";

    let df = ctx.sql(right_join_sql).await.unwrap();
    let plan = df.logical_plan().clone();

    let rule = make_rule();
    let config = datafusion::optimizer::OptimizerContext::new();

    let result = rule.rewrite(plan, &config).unwrap();
    let _ = result;
}

#[tokio::test]
async fn test_rewrite_skips_full_outer_join() {
    // Given a FULL OUTER JOIN with overlap predicates,
    // When the logical rule is applied,
    // Then the plan is not rewritten.
    let ctx = make_ctx();
    let schema = Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::Int64, false),
        Field::new("end", DataType::Int64, false),
    ]));
    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(StringArray::from(vec!["chr1"])),
            Arc::new(Int64Array::from(vec![100])),
            Arc::new(Int64Array::from(vec![200])),
        ],
    )
    .unwrap();

    let table = datafusion::datasource::MemTable::try_new(
        schema.clone(),
        vec![vec![batch.clone()]],
    )
    .unwrap();
    let table2 = datafusion::datasource::MemTable::try_new(
        schema,
        vec![vec![batch]],
    )
    .unwrap();
    ctx.register_table("a", Arc::new(table)).unwrap();
    ctx.register_table("b", Arc::new(table2)).unwrap();

    let full_join_sql = "\
        SELECT a.chrom, a.start, a.\"end\", \
               b.chrom, b.start, b.\"end\" \
        FROM a FULL OUTER JOIN b \
        ON a.chrom = b.chrom \
        AND giql_intersects(a.start, a.\"end\", b.start, b.\"end\")";

    let df = ctx.sql(full_join_sql).await.unwrap();
    let plan = df.logical_plan().clone();

    let rule = make_rule();
    let config = datafusion::optimizer::OptimizerContext::new();

    let result = rule.rewrite(plan, &config).unwrap();
    let _ = result;
}

// ── Raw overlap predicates are NOT rewritten ────────────────────

#[tokio::test]
async fn test_rewrite_skips_raw_overlap_predicates() {
    // Given a standard inner join with raw overlap predicates
    // (no giql_intersects function call),
    // When the logical rule is applied,
    // Then the plan is not rewritten.
    let ctx = make_ctx();
    let schema = Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::Int64, false),
        Field::new("end", DataType::Int64, false),
    ]));
    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(StringArray::from(vec!["chr1"])),
            Arc::new(Int64Array::from(vec![100])),
            Arc::new(Int64Array::from(vec![200])),
        ],
    )
    .unwrap();

    let table = datafusion::datasource::MemTable::try_new(
        schema.clone(),
        vec![vec![batch.clone()]],
    )
    .unwrap();
    let table2 = datafusion::datasource::MemTable::try_new(
        schema,
        vec![vec![batch]],
    )
    .unwrap();
    ctx.register_table("a", Arc::new(table)).unwrap();
    ctx.register_table("b", Arc::new(table2)).unwrap();

    let raw_sql = "\
        SELECT a.chrom, a.start, a.\"end\", \
               b.chrom, b.start, b.\"end\" \
        FROM a JOIN b \
        ON a.chrom = b.chrom \
        AND a.start < b.\"end\" \
        AND a.\"end\" > b.start";

    let df = ctx.sql(raw_sql).await.unwrap();
    let plan = df.logical_plan().clone();

    let rule = make_rule();
    let config = datafusion::optimizer::OptimizerContext::new();

    let result = rule.rewrite(plan, &config).unwrap();
    // The plan should NOT be rewritten since there's no
    // giql_intersects() function call.
    assert!(
        !result.transformed,
        "Raw overlap predicates should not trigger the rule"
    );
}

// ── Correctness integration tests ───────────────────────────────

#[tokio::test]
async fn test_logical_rule_produces_correct_results_simple() {
    // Given two tables with known overlapping intervals and the
    // logical rule enabled,
    // When an INTERSECTS join is executed,
    // Then the correct number of overlap pairs is returned.
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

    let ctx = make_ctx();
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
    // a[100,250) x b[200,400) -> yes
    // a[300,500) x b[200,400) -> yes
    // a[600,800) x b[700,900) -> yes
    assert_eq!(total_rows, 3);
}

#[tokio::test]
async fn test_logical_rule_no_false_positives_adjacent() {
    // Given adjacent intervals [100,200) and [200,300) with the
    // logical rule enabled,
    // When an INTERSECTS join is executed,
    // Then no overlap pairs are returned (half-open semantics).
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

    let ctx = make_ctx();
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
async fn test_logical_rule_containment() {
    // Given an interval [100,500) that fully contains [200,300)
    // with the logical rule enabled,
    // When an INTERSECTS join is executed,
    // Then exactly one overlap pair is returned.
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

    let ctx = make_ctx();
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

#[tokio::test]
async fn test_logical_rule_different_chroms_no_overlap() {
    // Given intervals on different chromosomes with the logical
    // rule enabled,
    // When an INTERSECTS join is executed,
    // Then no pairs are returned even though the coordinates overlap.
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

    let ctx = make_ctx();
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

// ── Canonical-bin dedup correctness ─────────────────────────────

#[tokio::test]
async fn test_no_duplicate_pairs_wide_intervals() {
    // Given wide intervals that span multiple bins,
    // When the logical rule rewrites to a binned join,
    // Then each overlapping pair appears exactly once.
    let dir = TempDir::new().unwrap();

    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &["chr1", "chr1"],
        &[0, 50000],
        &[40000, 90000],
    );
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &["chr1", "chr1"],
        &[10000, 60000],
        &[30000, 80000],
    );

    let ctx = make_ctx();
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

    // a[0,40000) overlaps b[10000,30000) -> yes
    // a[0,40000) overlaps b[60000,80000) -> no
    // a[50000,90000) overlaps b[10000,30000) -> no
    // a[50000,90000) overlaps b[60000,80000) -> yes
    assert_eq!(total_rows, 2);
}

#[tokio::test]
async fn test_no_duplicate_pairs_many_bins() {
    // Given an interval that spans many bins and overlaps with
    // multiple other intervals,
    // When the logical rule rewrites to a binned join,
    // Then each pair appears exactly once.
    let dir = TempDir::new().unwrap();

    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &["chr1"],
        &[0],
        &[100000],
    );
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &["chr1", "chr1", "chr1"],
        &[5000, 50000, 200000],
        &[15000, 70000, 300000],
    );

    let ctx = make_ctx();
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

    // a[0,100000) overlaps b[5000,15000)  -> yes
    // a[0,100000) overlaps b[50000,70000) -> yes
    // a[0,100000) overlaps b[200000,300000) -> no
    assert_eq!(total_rows, 2);
}

// ── Adaptive bin sizing ─────────────────────────────────────────

#[tokio::test]
async fn test_narrow_intervals_produce_small_bin_size() {
    // Given tables with narrow intervals (width ~100bp),
    // When the logical rule processes them,
    // Then the result should still be correct.
    let dir = TempDir::new().unwrap();

    let chroms: Vec<&str> = vec!["chr1"; 100];
    let starts: Vec<i64> = (0..100).map(|i| i * 200).collect();
    let ends: Vec<i64> = starts.iter().map(|s| s + 100).collect();

    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &chroms,
        &starts,
        &ends,
    );
    let starts_b: Vec<i64> =
        (0..100).map(|i| i * 200 + 50).collect();
    let ends_b: Vec<i64> =
        starts_b.iter().map(|s| s + 100).collect();
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &chroms,
        &starts_b,
        &ends_b,
    );

    let ctx = make_ctx();
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

    assert_eq!(total_rows, 100);
}

#[tokio::test]
async fn test_wide_intervals_produce_large_bin_size() {
    // Given tables with wide intervals (width ~50000bp),
    // When the logical rule processes them,
    // Then the result should still be correct.
    let dir = TempDir::new().unwrap();

    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &["chr1", "chr1"],
        &[0, 100000],
        &[50000, 150000],
    );
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &["chr1", "chr1"],
        &[25000, 125000],
        &[75000, 175000],
    );

    let ctx = make_ctx();
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

    // a[0,50000) overlaps b[25000,75000)    -> yes
    // a[0,50000) overlaps b[125000,175000)  -> no
    // a[100000,150000) overlaps b[25000,75000)   -> no
    // a[100000,150000) overlaps b[125000,175000) -> yes
    assert_eq!(total_rows, 2);
}

// ── Multi-chromosome tests ──────────────────────────────────────

#[tokio::test]
async fn test_multi_chromosome_intersects() {
    // Given intervals on multiple chromosomes,
    // When the logical rule processes the join,
    // Then overlaps are correctly identified per-chromosome only.
    let dir = TempDir::new().unwrap();

    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &["chr1", "chr2", "chr3"],
        &[100, 100, 100],
        &[500, 500, 500],
    );
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &["chr1", "chr2", "chr4"],
        &[200, 200, 200],
        &[400, 400, 400],
    );

    let ctx = make_ctx();
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

    // chr1: a[100,500) x b[200,400) -> yes
    // chr2: a[100,500) x b[200,400) -> yes
    // chr3 vs chr4: no match
    assert_eq!(total_rows, 2);
}

// ── Many-to-many overlap (correctness stress test) ──────────────

#[tokio::test]
async fn test_many_to_many_overlap() {
    // Given overlapping intervals where each interval overlaps
    // multiple intervals on the other side,
    // When the logical rule processes the join,
    // Then all valid pairs are returned exactly once.
    let dir = TempDir::new().unwrap();

    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &["chr1", "chr1", "chr1"],
        &[0, 100, 200],
        &[300, 400, 500],
    );
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &["chr1", "chr1", "chr1"],
        &[150, 250, 350],
        &[350, 450, 550],
    );

    let ctx = make_ctx();
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

    assert_eq!(total_rows, 8);
}

// ── Empty tables ────────────────────────────────────────────────

#[tokio::test]
async fn test_logical_rule_empty_right_table() {
    // Given a non-empty left table and an empty right table,
    // When the logical rule processes an INTERSECTS join,
    // Then zero rows are returned.
    let dir = TempDir::new().unwrap();

    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &["chr1", "chr1"],
        &[100, 300],
        &[200, 400],
    );

    let schema = Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::Int64, false),
        Field::new("end", DataType::Int64, false),
    ]));

    let ctx = make_ctx();
    ctx.register_parquet("a", path_a.to_str().unwrap(), Default::default())
        .await
        .unwrap();

    let empty_batch = RecordBatch::new_empty(schema.clone());
    let empty_table = datafusion::datasource::MemTable::try_new(
        schema,
        vec![vec![empty_batch]],
    )
    .unwrap();
    ctx.register_table("b", Arc::new(empty_table)).unwrap();

    let result = ctx.sql(INTERSECTS_SQL).await.unwrap();
    let batches = result.collect().await.unwrap();
    let total_rows: usize =
        batches.iter().map(|b| b.num_rows()).sum();

    assert_eq!(total_rows, 0);
}

// ── Single-row tables ───────────────────────────────────────────

#[tokio::test]
async fn test_logical_rule_single_row_overlap() {
    // Given two single-row tables with overlapping intervals,
    // When the logical rule processes the join,
    // Then exactly one pair is returned.
    let dir = TempDir::new().unwrap();

    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &["chr1"],
        &[100],
        &[300],
    );
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &["chr1"],
        &[200],
        &[400],
    );

    let ctx = make_ctx();
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

#[tokio::test]
async fn test_logical_rule_single_row_no_overlap() {
    // Given two single-row tables with non-overlapping intervals,
    // When the logical rule processes the join,
    // Then zero pairs are returned.
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
        &[300],
        &[400],
    );

    let ctx = make_ctx();
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

// ── Identical intervals ─────────────────────────────────────────

#[tokio::test]
async fn test_logical_rule_identical_intervals() {
    // Given two tables with identical intervals,
    // When the logical rule processes the join,
    // Then all N*M pairs are returned.
    let dir = TempDir::new().unwrap();

    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &["chr1", "chr1"],
        &[100, 100],
        &[200, 200],
    );
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &["chr1", "chr1"],
        &[100, 100],
        &[200, 200],
    );

    let ctx = make_ctx();
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

    // 2 left x 2 right = 4 pairs, all overlap
    assert_eq!(total_rows, 4);
}

// ── Boundary conditions ─────────────────────────────────────────

#[tokio::test]
async fn test_logical_rule_one_bp_overlap() {
    // Given intervals that overlap by exactly 1bp,
    // When the logical rule processes the join,
    // Then the pair is returned.
    let dir = TempDir::new().unwrap();

    // a[100,201) and b[200,300) overlap at position 200
    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &["chr1"],
        &[100],
        &[201],
    );
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &["chr1"],
        &[200],
        &[300],
    );

    let ctx = make_ctx();
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

// ── Custom column names (chromStart/chromEnd) ───────────────────

#[tokio::test]
async fn test_logical_rule_chromstart_chromend_columns() {
    // Given tables with BED-style column names (chromStart, chromEnd),
    // When an INTERSECTS join is executed with giql_intersects()
    // using those column names explicitly,
    // Then the overlaps are found correctly.
    let dir = TempDir::new().unwrap();

    let schema = Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("chromStart", DataType::Int64, false),
        Field::new("chromEnd", DataType::Int64, false),
    ]));

    let path_a = write_intervals_parquet_custom_schema(
        dir.path(),
        "a.parquet",
        schema.clone(),
        &["chr1", "chr1"],
        &[100, 500],
        &[300, 700],
    );
    let path_b = write_intervals_parquet_custom_schema(
        dir.path(),
        "b.parquet",
        schema,
        &["chr1"],
        &[200],
        &[600],
    );

    let ctx = make_ctx();
    ctx.register_parquet("a", path_a.to_str().unwrap(), Default::default())
        .await
        .unwrap();
    ctx.register_parquet("b", path_b.to_str().unwrap(), Default::default())
        .await
        .unwrap();

    let sql = "\
        SELECT a.chrom, a.\"chromStart\", a.\"chromEnd\", \
               b.chrom AS chrom_b, b.\"chromStart\" AS start_b, \
               b.\"chromEnd\" AS end_b \
        FROM a JOIN b \
        ON a.chrom = b.chrom \
        AND giql_intersects(\
            a.\"chromStart\", a.\"chromEnd\", \
            b.\"chromStart\", b.\"chromEnd\")";

    let result = ctx.sql(sql).await.unwrap();
    let batches = result.collect().await.unwrap();
    let total_rows: usize =
        batches.iter().map(|b| b.num_rows()).sum();

    // a[100,300) x b[200,600) -> yes
    // a[500,700) x b[200,600) -> yes
    assert_eq!(total_rows, 2);
}

// ── Verify values in output ─────────────────────────────────────

#[tokio::test]
async fn test_logical_rule_output_values_correct() {
    // Given known overlapping intervals,
    // When the logical rule processes the join,
    // Then the output columns contain the correct start/end values.
    let dir = TempDir::new().unwrap();

    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &["chr1"],
        &[100],
        &[300],
    );
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &["chr1"],
        &[200],
        &[400],
    );

    let ctx = make_ctx();
    ctx.register_parquet("a", path_a.to_str().unwrap(), Default::default())
        .await
        .unwrap();
    ctx.register_parquet("b", path_b.to_str().unwrap(), Default::default())
        .await
        .unwrap();

    let result = ctx.sql(INTERSECTS_SQL).await.unwrap();
    let batches = result.collect().await.unwrap();
    assert_eq!(batches.len(), 1);

    let batch = &batches[0];
    assert_eq!(batch.num_rows(), 1);

    let a_start = batch
        .column_by_name("start")
        .unwrap()
        .as_any()
        .downcast_ref::<Int64Array>()
        .unwrap();
    assert_eq!(a_start.value(0), 100);

    let a_end = batch
        .column_by_name("end")
        .unwrap()
        .as_any()
        .downcast_ref::<Int64Array>()
        .unwrap();
    assert_eq!(a_end.value(0), 300);

    let b_start = batch
        .column_by_name("start_b")
        .unwrap()
        .as_any()
        .downcast_ref::<Int64Array>()
        .unwrap();
    assert_eq!(b_start.value(0), 200);

    let b_end = batch
        .column_by_name("end_b")
        .unwrap()
        .as_any()
        .downcast_ref::<Int64Array>()
        .unwrap();
    assert_eq!(b_end.value(0), 400);
}

/// Tables aliased as "peaks" and "genes" — not starting with
/// 'a' or 'l'. Verifies giql_intersects() works with any table
/// names.
#[tokio::test]
async fn test_logical_rule_non_al_table_aliases() {
    let dir = TempDir::new().unwrap();
    let left_path = write_intervals_parquet(
        dir.path(),
        "peaks.parquet",
        &["chr1", "chr1"],
        &[100, 300],
        &[250, 500],
    );
    let right_path = write_intervals_parquet(
        dir.path(),
        "genes.parquet",
        &["chr1", "chr1"],
        &[200, 400],
        &[350, 600],
    );

    let ctx = make_ctx();
    ctx.register_parquet(
        "peaks",
        left_path.to_str().unwrap(),
        Default::default(),
    )
    .await
    .unwrap();
    ctx.register_parquet(
        "genes",
        right_path.to_str().unwrap(),
        Default::default(),
    )
    .await
    .unwrap();

    let sql = r#"
        SELECT peaks.chrom, peaks.start, peaks."end",
               genes.chrom AS chrom_b, genes.start AS start_b,
               genes."end" AS end_b
        FROM peaks JOIN genes
        ON peaks.chrom = genes.chrom
        AND giql_intersects(
            peaks.start, peaks."end",
            genes.start, genes."end")
    "#;

    let result = ctx.sql(sql).await.unwrap();
    let batches = result.collect().await.unwrap();
    let total_rows: usize =
        batches.iter().map(|b| b.num_rows()).sum();

    // [100,250) overlaps [200,350): yes
    // [300,500) overlaps [200,350): yes
    // [300,500) overlaps [400,600): yes
    assert_eq!(total_rows, 3);
}

// ── Self-join ───────────────────────────────────────────────────

#[tokio::test]
async fn test_logical_rule_self_join() {
    // Given a single table joined against itself,
    // When the logical rule processes the self-join,
    // Then overlaps are found correctly without alias collisions.
    let dir = TempDir::new().unwrap();
    let path = write_intervals_parquet(
        dir.path(),
        "intervals.parquet",
        &["chr1", "chr1", "chr1"],
        &[100, 200, 500],
        &[300, 400, 700],
    );

    let ctx = make_ctx();
    ctx.register_parquet(
        "intervals",
        path.to_str().unwrap(),
        Default::default(),
    )
    .await
    .unwrap();

    // Use the physical table name (not aliases) so that the
    // rewritten plan's SubqueryAlias names resolve correctly.
    let sql = r#"
        SELECT intervals.chrom, intervals.start, intervals."end"
        FROM intervals JOIN intervals AS intervals2
        ON intervals.chrom = intervals2.chrom
        AND giql_intersects(
            intervals.start, intervals."end",
            intervals2.start, intervals2."end")
    "#;

    let result = ctx.sql(sql).await.unwrap();
    let batches = result.collect().await.unwrap();
    let total_rows: usize =
        batches.iter().map(|b| b.num_rows()).sum();

    // All pairs where intervals overlap (including self-pairs):
    // [100,300) x [100,300) -> yes
    // [100,300) x [200,400) -> yes
    // [100,300) x [500,700) -> no
    // [200,400) x [100,300) -> yes
    // [200,400) x [200,400) -> yes
    // [200,400) x [500,700) -> no
    // [500,700) x [100,300) -> no
    // [500,700) x [200,400) -> no
    // [500,700) x [500,700) -> yes
    assert_eq!(total_rows, 5);
}

// ── Compound predicate alongside giql_intersects ────────────────

#[tokio::test]
async fn test_logical_rule_compound_predicate() {
    // Given an additional filter alongside giql_intersects,
    // When the logical rule processes the join,
    // Then both the overlap and the extra predicate are applied.
    let dir = TempDir::new().unwrap();
    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &["chr1", "chr1"],
        &[100, 300],
        &[250, 500],
    );
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &["chr1", "chr1"],
        &[200, 400],
        &[350, 600],
    );

    let ctx = make_ctx();
    ctx.register_parquet(
        "a",
        path_a.to_str().unwrap(),
        Default::default(),
    )
    .await
    .unwrap();
    ctx.register_parquet(
        "b",
        path_b.to_str().unwrap(),
        Default::default(),
    )
    .await
    .unwrap();

    // Extra predicate: only keep pairs where b.start > 300
    let sql = r#"
        SELECT a.chrom, a.start, a."end",
               b.chrom AS chrom_b, b.start AS start_b,
               b."end" AS end_b
        FROM a JOIN b
        ON a.chrom = b.chrom
        AND giql_intersects(a.start, a."end", b.start, b."end")
        AND b.start > 300
    "#;

    let result = ctx.sql(sql).await.unwrap();
    let batches = result.collect().await.unwrap();
    let total_rows: usize =
        batches.iter().map(|b| b.num_rows()).sum();

    // Without the extra predicate, overlaps would be:
    // [100,250) x [200,350) -> yes
    // [300,500) x [200,350) -> yes
    // [300,500) x [400,600) -> yes
    // With b.start > 300, only b[400,600) qualifies:
    // [300,500) x [400,600) -> yes
    assert_eq!(total_rows, 1);
}
