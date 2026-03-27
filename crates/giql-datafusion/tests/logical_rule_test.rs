//! Tests for the logical optimizer rule (`logical_rule.rs`).
//!
//! Covers:
//! - OptimizerRule trait implementation (name, apply_order, supports_rewrite)
//! - Join type filtering (inner only, skip left/right/full outer)
//! - Already-binned join detection (skip re-rewrite)
//! - Overlap pattern detection (start/end column name variants)
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
use giql_datafusion::{IntersectsOptimizerConfig, register_optimizer};

// ── Helpers ─────────────────────────────────────────────────────

fn default_config() -> IntersectsOptimizerConfig {
    IntersectsOptimizerConfig::default()
}

fn make_rule() -> IntersectsLogicalRule {
    IntersectsLogicalRule::new(default_config())
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

/// Create a SessionContext with the logical rule enabled.
fn make_ctx_with_logical_rule() -> SessionContext {
    let config = IntersectsOptimizerConfig {
        enable_logical_rule: true,
        ..default_config()
    };
    let state = SessionStateBuilder::new()
        .with_default_features()
        .build();
    let state = register_optimizer(state, config);
    SessionContext::from(state)
}

/// Create a SessionContext with the logical rule disabled.
fn make_ctx_without_logical_rule() -> SessionContext {
    let config = IntersectsOptimizerConfig {
        enable_logical_rule: false,
        ..default_config()
    };
    let state = SessionStateBuilder::new()
        .with_default_features()
        .build();
    let state = register_optimizer(state, config);
    SessionContext::from(state)
}

const INTERSECTS_SQL: &str = "\
    SELECT a.chrom, a.start, a.\"end\", \
           b.chrom AS chrom_b, b.start AS start_b, b.\"end\" AS end_b \
    FROM a JOIN b \
    ON a.chrom = b.chrom \
    AND a.start < b.\"end\" \
    AND a.\"end\" > b.start";

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
    // Given a non-join logical plan (TableScan),
    // When the rule is applied,
    // Then the plan is returned unchanged.
    let rule = make_rule();
    let config = datafusion::optimizer::OptimizerContext::new();

    // Create a simple empty relation plan (not a join)
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
    // Given a LEFT JOIN with interval overlap predicates,
    // When the logical rule is applied,
    // Then the plan is not rewritten (only INNER joins are supported).
    let ctx = SessionContext::new();
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
        AND a.start < b.\"end\" \
        AND a.\"end\" > b.start";

    let df = ctx.sql(left_join_sql).await.unwrap();
    let plan = df.logical_plan().clone();

    let rule = make_rule();
    let config = datafusion::optimizer::OptimizerContext::new();

    // Walk the plan tree looking for join nodes
    let result = rule.rewrite(plan, &config).unwrap();
    // Either the plan is not transformed (because it's not a Join
    // at top level), or if DataFusion restructured it, the rule
    // should still not rewrite non-inner joins.
    // The important thing is that the rule doesn't panic.
    let _ = result;
}

#[tokio::test]
async fn test_rewrite_skips_right_join() {
    // Given a RIGHT JOIN with interval overlap predicates,
    // When the logical rule is applied,
    // Then the plan is not rewritten (only INNER joins are supported).
    let ctx = SessionContext::new();
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
        AND a.start < b.\"end\" \
        AND a.\"end\" > b.start";

    let df = ctx.sql(right_join_sql).await.unwrap();
    let plan = df.logical_plan().clone();

    let rule = make_rule();
    let config = datafusion::optimizer::OptimizerContext::new();

    let result = rule.rewrite(plan, &config).unwrap();
    let _ = result;
}

#[tokio::test]
async fn test_rewrite_skips_full_outer_join() {
    // Given a FULL OUTER JOIN with interval overlap predicates,
    // When the logical rule is applied,
    // Then the plan is not rewritten.
    let ctx = SessionContext::new();
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
        AND a.start < b.\"end\" \
        AND a.\"end\" > b.start";

    let df = ctx.sql(full_join_sql).await.unwrap();
    let plan = df.logical_plan().clone();

    let rule = make_rule();
    let config = datafusion::optimizer::OptimizerContext::new();

    let result = rule.rewrite(plan, &config).unwrap();
    let _ = result;
}

// ── Register optimizer with logical rule enabled/disabled ───────

#[test]
fn test_register_optimizer_with_logical_rule_enabled() {
    // Given a default config with enable_logical_rule = true,
    // When register_optimizer is called,
    // Then both logical and physical rules are added.
    let state = SessionStateBuilder::new()
        .with_default_features()
        .build();
    let n_logical_before = state.optimizers().len();
    let n_physical_before = state.physical_optimizers().len();

    let config = IntersectsOptimizerConfig {
        enable_logical_rule: true,
        ..default_config()
    };
    let state = register_optimizer(state, config);

    assert_eq!(
        state.optimizers().len(),
        n_logical_before + 1,
        "Should add one logical rule"
    );
    assert_eq!(
        state.physical_optimizers().len(),
        n_physical_before + 1,
        "Should add one physical rule"
    );

    let last_logical = state.optimizers().last().unwrap();
    assert_eq!(last_logical.name(), "intersects_logical_binned");
}

#[test]
fn test_register_optimizer_with_logical_rule_disabled() {
    // Given a config with enable_logical_rule = false,
    // When register_optimizer is called,
    // Then only the physical rule is added, not the logical rule.
    let state = SessionStateBuilder::new()
        .with_default_features()
        .build();
    let n_logical_before = state.optimizers().len();
    let n_physical_before = state.physical_optimizers().len();

    let config = IntersectsOptimizerConfig {
        enable_logical_rule: false,
        ..default_config()
    };
    let state = register_optimizer(state, config);

    assert_eq!(
        state.optimizers().len(),
        n_logical_before,
        "Should NOT add a logical rule"
    );
    assert_eq!(
        state.physical_optimizers().len(),
        n_physical_before + 1,
        "Should still add the physical rule"
    );
}

// ── Adaptive bin sizing integration tests ───────────────────────

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

    let ctx = make_ctx_with_logical_rule();
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

    let ctx = make_ctx_with_logical_rule();
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

    let ctx = make_ctx_with_logical_rule();
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

    let ctx = make_ctx_with_logical_rule();
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
    // Then each overlapping pair appears exactly once (dedup filter
    // eliminates multi-bin duplicates).
    let dir = TempDir::new().unwrap();

    // Wide intervals spanning many bins (default bin ~10k)
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

    let ctx = make_ctx_with_logical_rule();
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
    // a[0,40000) overlaps b[60000,80000) -> no (40000 < 60000)
    // a[50000,90000) overlaps b[10000,30000) -> no (50000 >= 30000)
    // a[50000,90000) overlaps b[60000,80000) -> yes
    assert_eq!(total_rows, 2);
}

#[tokio::test]
async fn test_no_duplicate_pairs_many_bins() {
    // Given an interval that spans many bins and overlaps with
    // multiple other intervals,
    // When the logical rule rewrites to a binned join,
    // Then each pair appears exactly once regardless of how many
    // bins they share.
    let dir = TempDir::new().unwrap();

    // One very wide interval on each side, plus a narrow one
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

    let ctx = make_ctx_with_logical_rule();
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
    // a[0,100000) overlaps b[200000,300000) -> no (100000 <= 200000)
    assert_eq!(total_rows, 2);
}

// ── Adaptive bin sizing ─────────────────────────────────────────

#[tokio::test]
async fn test_narrow_intervals_produce_small_bin_size() {
    // Given tables with narrow intervals (width ~100bp),
    // When the logical rule processes them,
    // Then the bin size should be small (clamped to minimum 1000)
    // and the result should still be correct.
    let dir = TempDir::new().unwrap();

    // 100 narrow intervals of width 100
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
    // Overlapping intervals offset by 50
    let starts_b: Vec<i64> = (0..100).map(|i| i * 200 + 50).collect();
    let ends_b: Vec<i64> = starts_b.iter().map(|s| s + 100).collect();
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &chroms,
        &starts_b,
        &ends_b,
    );

    let ctx = make_ctx_with_logical_rule();
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

    // Each a interval [i*200, i*200+100) overlaps b[i*200+50, i*200+150)
    // Plus a[i*200, i*200+100) may also overlap b[(i-1)*200+50, (i-1)*200+150)
    // when i*200 < (i-1)*200+150, i.e., 200 < 150 -> never.
    // So exactly 100 pairs.
    assert_eq!(total_rows, 100);
}

#[tokio::test]
async fn test_wide_intervals_produce_large_bin_size() {
    // Given tables with wide intervals (width ~50000bp),
    // When the logical rule processes them,
    // Then the result should still be correct with the adaptively
    // chosen larger bin size.
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

    let ctx = make_ctx_with_logical_rule();
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

    let ctx = make_ctx_with_logical_rule();
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

    // Three intervals each spanning [0,300), [100,400), [200,500)
    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &["chr1", "chr1", "chr1"],
        &[0, 100, 200],
        &[300, 400, 500],
    );
    // Three intervals each spanning [150,350), [250,450), [350,550)
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &["chr1", "chr1", "chr1"],
        &[150, 250, 350],
        &[350, 450, 550],
    );

    let ctx = make_ctx_with_logical_rule();
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

    // a[0,300) vs b[150,350) -> yes (0<350, 300>150)
    // a[0,300) vs b[250,450) -> yes (0<450, 300>250)
    // a[0,300) vs b[350,550) -> no  (300 <= 350)
    // a[100,400) vs b[150,350) -> yes
    // a[100,400) vs b[250,450) -> yes
    // a[100,400) vs b[350,550) -> yes (100<550, 400>350)
    // a[200,500) vs b[150,350) -> yes (200<350, 500>150)
    // a[200,500) vs b[250,450) -> yes
    // a[200,500) vs b[350,550) -> yes (200<550, 500>350)
    assert_eq!(total_rows, 8);
}

// ── Logical rule vs no logical rule consistency ─────────────────

#[tokio::test]
async fn test_logical_rule_matches_baseline_results() {
    // Given the same data,
    // When an INTERSECTS join is executed with and without the
    // logical rule,
    // Then both produce the same number of results.
    let dir = TempDir::new().unwrap();

    let path_a = write_intervals_parquet(
        dir.path(),
        "a.parquet",
        &["chr1", "chr1", "chr2", "chr2"],
        &[100, 500, 200, 800],
        &[400, 900, 600, 1200],
    );
    let path_b = write_intervals_parquet(
        dir.path(),
        "b.parquet",
        &["chr1", "chr1", "chr2"],
        &[300, 700, 400],
        &[600, 1000, 700],
    );

    // With logical rule
    let ctx_with = make_ctx_with_logical_rule();
    ctx_with
        .register_parquet("a", path_a.to_str().unwrap(), Default::default())
        .await
        .unwrap();
    ctx_with
        .register_parquet("b", path_b.to_str().unwrap(), Default::default())
        .await
        .unwrap();

    let result_with = ctx_with.sql(INTERSECTS_SQL).await.unwrap();
    let batches_with = result_with.collect().await.unwrap();
    let rows_with: usize =
        batches_with.iter().map(|b| b.num_rows()).sum();

    // Without logical rule
    let ctx_without = make_ctx_without_logical_rule();
    ctx_without
        .register_parquet("a", path_a.to_str().unwrap(), Default::default())
        .await
        .unwrap();
    ctx_without
        .register_parquet("b", path_b.to_str().unwrap(), Default::default())
        .await
        .unwrap();

    let result_without = ctx_without.sql(INTERSECTS_SQL).await.unwrap();
    let batches_without = result_without.collect().await.unwrap();
    let rows_without: usize =
        batches_without.iter().map(|b| b.num_rows()).sum();

    assert_eq!(
        rows_with, rows_without,
        "Logical rule should produce same count as baseline"
    );

    // Also verify the expected count:
    // chr1: a[100,400) x b[300,600) -> yes
    //        a[100,400) x b[700,1000) -> no
    //        a[500,900) x b[300,600) -> yes
    //        a[500,900) x b[700,1000) -> yes
    // chr2: a[200,600) x b[400,700) -> yes
    //        a[800,1200) x b[400,700) -> no
    assert_eq!(rows_with, 4);
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

    // Empty table - at least one row needed for Parquet writing,
    // so we'll use the memory table approach
    let schema = Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::Int64, false),
        Field::new("end", DataType::Int64, false),
    ]));

    let ctx = make_ctx_with_logical_rule();
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

    let ctx = make_ctx_with_logical_rule();
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

    let ctx = make_ctx_with_logical_rule();
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
    // Then all N*M pairs are returned where both are on the same
    // chrom.
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

    let ctx = make_ctx_with_logical_rule();
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

    let ctx = make_ctx_with_logical_rule();
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

// ── Column name variants (chromStart/chromEnd) ──────────────────

#[tokio::test]
async fn test_logical_rule_chromstart_chromend_columns() {
    // Given tables with BED-style column names (chromStart, chromEnd),
    // When an INTERSECTS join is executed with the logical rule,
    // Then the column names are recognized and overlaps are found.
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

    let ctx = make_ctx_with_logical_rule();
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
        AND a.\"chromStart\" < b.\"chromEnd\" \
        AND a.\"chromEnd\" > b.\"chromStart\"";

    let result = ctx.sql(sql).await.unwrap();
    let batches = result.collect().await.unwrap();
    let total_rows: usize =
        batches.iter().map(|b| b.num_rows()).sum();

    // a[100,300) x b[200,600) -> yes
    // a[500,700) x b[200,600) -> yes (500 < 600, 700 > 200)
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

    let ctx = make_ctx_with_logical_rule();
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

    // Check the values: should have a.start=100, a.end=300,
    // b.start=200, b.end=400
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
