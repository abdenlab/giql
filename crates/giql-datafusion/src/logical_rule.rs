use std::sync::Arc;

use datafusion::common::tree_node::Transformed;
use datafusion::common::{Column, Result, ScalarValue};
use datafusion::datasource::listing::ListingTable;
use datafusion::datasource::source_as_provider;
use datafusion::logical_expr::expr::ScalarFunction;
use datafusion::logical_expr::{
    BinaryExpr, Expr, Join, JoinType, LogicalPlan, LogicalPlanBuilder,
    Operator,
};
use datafusion::optimizer::{OptimizerConfig, OptimizerRule};
use datafusion::prelude::*;

/// Logical optimizer rule that rewrites interval overlap joins into
/// binned equi-joins using UNNEST.
///
/// Detects `giql_intersects(start_a, end_a, start_b, end_b)` function
/// calls in join filters (emitted by the GIQL transpiler's
/// `"datafusion"` dialect) and rewrites them to:
///
///   `SELECT ... FROM Unnest(a + bins) JOIN Unnest(b + bins)
///    ON chrom = chrom AND bin = bin WHERE start < end AND end > start`
///
/// DataFusion handles UNNEST, hash join, and dedup natively with
/// full parallelism.
/// Configuration for the INTERSECTS logical optimizer rule.
#[derive(Debug, Clone, Default)]
pub struct IntersectsConfig {
    /// Force the binned equi-join strategy instead of the default
    /// COI tree join. The bin size is chosen adaptively from
    /// Parquet sampling or column-level statistics. This is an
    /// escape hatch for benchmarking; COI tree is faster in all
    /// tested distributions.
    pub force_binned: bool,
}

#[derive(Debug)]
pub struct IntersectsLogicalRule {
    config: IntersectsConfig,
}

impl IntersectsLogicalRule {
    pub fn new() -> Self {
        Self {
            config: IntersectsConfig::default(),
        }
    }

    pub fn with_config(config: IntersectsConfig) -> Self {
        Self { config }
    }
}

impl Default for IntersectsLogicalRule {
    fn default() -> Self {
        Self::new()
    }
}

impl OptimizerRule for IntersectsLogicalRule {
    fn name(&self) -> &str {
        "intersects_logical_binned"
    }

    fn apply_order(&self) -> Option<datafusion::optimizer::ApplyOrder> {
        Some(datafusion::optimizer::ApplyOrder::BottomUp)
    }

    fn supports_rewrite(&self) -> bool {
        true
    }

    fn rewrite(
        &self,
        plan: LogicalPlan,
        _config: &dyn OptimizerConfig,
    ) -> Result<Transformed<LogicalPlan>> {
        let LogicalPlan::Join(ref join) = plan else {
            return Ok(Transformed::no(plan));
        };

        if join.join_type != JoinType::Inner {
            return Ok(Transformed::no(plan));
        }

        // Skip if already rewritten (has __giql_bins equi-keys)
        let already_binned = join.on.iter().any(|(l, _)| {
            if let Expr::Column(c) = l {
                c.name.starts_with("__giql_bins")
            } else {
                false
            }
        });
        if already_binned {
            return Ok(Transformed::no(plan));
        }

        // Detect giql_intersects() function call in the filter
        let overlap = match &join.filter {
            Some(filter) => detect_giql_intersects(filter),
            None => None,
        };

        let Some((start_a, end_a, start_b, end_b)) = overlap else {
            return Ok(Transformed::no(plan));
        };

        if self.config.force_binned {
            // Binned equi-join path (escape hatch for benchmarking).
            let left_stats = get_table_stats(
                &join.left, &start_a.name, &end_a.name,
            );
            let right_stats = get_table_stats(
                &join.right, &start_b.name, &end_b.name,
            );
            let bin_size =
                choose_bin_size(&left_stats, &right_stats);

            log::debug!(
                "INTERSECTS logical rule: rewriting to \
                 binned join (force_binned), bin_size={bin_size}"
            );

            let rewritten_filter =
                join.filter.as_ref().map(|f| {
                    replace_giql_intersects(
                        f, &start_a, &end_a, &start_b, &end_b,
                    )
                });

            let rewritten = rewrite_to_binned(
                join,
                bin_size,
                &start_a,
                &end_a,
                &start_b,
                &end_b,
                rewritten_filter.as_ref(),
            )?;
            return Ok(Transformed::yes(rewritten));
        }

        // Default: COI tree join — faster across all tested
        // distributions, no bin replication overhead.
        use crate::coitree::COITreeJoinNode;
        use datafusion::logical_expr::Extension;

        log::debug!("INTERSECTS logical rule: using COI tree join");

        let on: Vec<(Column, Column)> = join
            .on
            .iter()
            .map(|(l, r)| {
                (
                    extract_column(l).unwrap_or_else(|| {
                        Column::new(None::<&str>, "chrom")
                    }),
                    extract_column(r).unwrap_or_else(|| {
                        Column::new(None::<&str>, "chrom")
                    }),
                )
            })
            .collect();

        let node = COITreeJoinNode {
            left: Arc::new((*join.left).clone()),
            right: Arc::new((*join.right).clone()),
            on,
            start_a,
            end_a,
            start_b,
            end_b,
            schema: join.schema.clone(),
        };

        Ok(Transformed::yes(LogicalPlan::Extension(Extension {
            node: Arc::new(node),
        })))
    }
}

// ── Pattern detection ───────────────────────────────────────────

/// Detect `giql_intersects(start_a, end_a, start_b, end_b)` in a
/// filter expression. Searches through AND-combined predicates.
///
/// Returns `(start_a, end_a, start_b, end_b)` column references.
fn detect_giql_intersects(
    expr: &Expr,
) -> Option<(Column, Column, Column, Column)> {
    match expr {
        Expr::ScalarFunction(func)
            if func.name() == "giql_intersects"
                && func.args.len() == 4 =>
        {
            let start_a = extract_column(&func.args[0])?;
            let end_a = extract_column(&func.args[1])?;
            let start_b = extract_column(&func.args[2])?;
            let end_b = extract_column(&func.args[3])?;
            Some((start_a, end_a, start_b, end_b))
        }
        Expr::BinaryExpr(BinaryExpr {
            left,
            op: Operator::And,
            right,
        }) => detect_giql_intersects(left)
            .or_else(|| detect_giql_intersects(right)),
        _ => None,
    }
}

/// Extract a Column from an Expr, handling TryCast/Cast wrappers
/// that DataFusion may insert during type coercion.
fn extract_column(expr: &Expr) -> Option<Column> {
    match expr {
        Expr::Column(c) => Some(c.clone()),
        Expr::Cast(cast) => extract_column(&cast.expr),
        Expr::TryCast(tc) => extract_column(&tc.expr),
        _ => None,
    }
}

/// Replace `giql_intersects(start_a, end_a, start_b, end_b)` in an
/// expression tree with `start_a < end_b AND end_a > start_b`.
fn replace_giql_intersects(
    expr: &Expr,
    start_a: &Column,
    end_a: &Column,
    start_b: &Column,
    end_b: &Column,
) -> Expr {
    match expr {
        Expr::ScalarFunction(func)
            if func.name() == "giql_intersects" =>
        {
            // start_a < end_b AND end_a > start_b
            Expr::BinaryExpr(BinaryExpr {
                left: Box::new(Expr::BinaryExpr(BinaryExpr {
                    left: Box::new(Expr::Column(start_a.clone())),
                    op: Operator::Lt,
                    right: Box::new(Expr::Column(end_b.clone())),
                })),
                op: Operator::And,
                right: Box::new(Expr::BinaryExpr(BinaryExpr {
                    left: Box::new(Expr::Column(end_a.clone())),
                    op: Operator::Gt,
                    right: Box::new(Expr::Column(start_b.clone())),
                })),
            })
        }
        Expr::BinaryExpr(BinaryExpr { left, op, right }) => {
            Expr::BinaryExpr(BinaryExpr {
                left: Box::new(replace_giql_intersects(
                    left, start_a, end_a, start_b, end_b,
                )),
                op: *op,
                right: Box::new(replace_giql_intersects(
                    right, start_a, end_a, start_b, end_b,
                )),
            })
        }
        other => other.clone(),
    }
}

// ── Stats collection ────────────────────────────────────────────

struct LogicalStats {
    /// Reserved for future use (e.g., skipping binning on tiny
    /// tables where a nested-loop join would be cheaper).
    #[allow(dead_code)]
    row_count: Option<usize>,
    start_min: Option<i64>,
    start_max: Option<i64>,
    end_min: Option<i64>,
    end_max: Option<i64>,
    /// Sampled width statistics from Parquet row groups.
    sampled: Option<SampledWidthStats>,
}

/// Width statistics computed from sampled Parquet row groups.
struct SampledWidthStats {
    /// Smallest bin size with mean replication <= 2.0.
    optimal_bin: i64,
    /// Median interval width.
    median: i64,
}

fn get_table_stats(
    plan: &LogicalPlan,
    start_col_name: &str,
    end_col_name: &str,
) -> Option<LogicalStats> {
    match plan {
        LogicalPlan::TableScan(ts) => {
            let provider = source_as_provider(&ts.source).ok()?;

            let mut stats = provider
                .statistics()
                .and_then(|s| {
                    stats_to_logical(
                        &s,
                        &ts.source.schema(),
                        start_col_name,
                        end_col_name,
                    )
                })
                .unwrap_or(LogicalStats {
                    row_count: None,
                    start_min: None,
                    start_max: None,
                    end_min: None,
                    end_max: None,
                    sampled: None,
                });

            // Try lightweight Parquet sampling for accurate width
            stats.sampled = try_sample_from_listing(
                provider.as_ref(),
                start_col_name,
                end_col_name,
            );

            Some(stats)
        }
        _ => plan.inputs().first().and_then(|child| {
            get_table_stats(child, start_col_name, end_col_name)
        }),
    }
}

fn stats_to_logical(
    stats: &datafusion::common::Statistics,
    schema: &arrow::datatypes::SchemaRef,
    start_col_name: &str,
    end_col_name: &str,
) -> Option<LogicalStats> {
    let row_count = match stats.num_rows {
        datafusion::common::stats::Precision::Exact(n) => Some(n),
        datafusion::common::stats::Precision::Inexact(n) => Some(n),
        _ => None,
    };
    let start_idx = schema
        .fields()
        .iter()
        .position(|f| f.name() == start_col_name)?;
    let end_idx = schema
        .fields()
        .iter()
        .position(|f| f.name() == end_col_name)?;
    let col_stats = &stats.column_statistics;
    let start_stats = col_stats.get(start_idx)?;
    let end_stats = col_stats.get(end_idx)?;
    Some(LogicalStats {
        row_count,
        start_min: scalar_to_i64(&start_stats.min_value),
        start_max: scalar_to_i64(&start_stats.max_value),
        end_min: scalar_to_i64(&end_stats.min_value),
        end_max: scalar_to_i64(&end_stats.max_value),
        sampled: None, // filled by try_sample_from_listing
    })
}

fn scalar_to_i64(
    precision: &datafusion::common::stats::Precision<ScalarValue>,
) -> Option<i64> {
    match precision {
        datafusion::common::stats::Precision::Exact(v)
        | datafusion::common::stats::Precision::Inexact(v) => match v {
            ScalarValue::Int32(Some(n)) => Some(*n as i64),
            ScalarValue::Int64(Some(n)) => Some(*n),
            _ => None,
        },
        _ => None,
    }
}

// ── Lightweight Parquet sampling ─────────────────────────────────

/// Try to sample interval widths from a Parquet-backed ListingTable.
///
/// Returns `None` silently if the provider is not a ListingTable,
/// the file is not Parquet, or any I/O error occurs.
fn try_sample_from_listing(
    provider: &dyn datafusion::catalog::TableProvider,
    start_col: &str,
    end_col: &str,
) -> Option<SampledWidthStats> {
    let listing = provider.as_any().downcast_ref::<ListingTable>()?;
    let path_str = listing.table_paths().first()?.as_str();

    // ListingTableUrl stores file:// URLs. For remote sources
    // (s3://, gs://, etc.) the else branch produces a path that
    // won't exist on disk — File::open fails and we fall back to
    // column-level stats gracefully.
    let fs_path = if let Some(p) = path_str.strip_prefix("file://") {
        std::path::PathBuf::from(p)
    } else {
        std::path::PathBuf::from(format!("/{path_str}"))
    };

    match sample_width_stats(&fs_path, start_col, end_col) {
        Some(stats) => {
            log::debug!(
                "INTERSECTS logical: sampled optimal_bin={}, \
                 median={} from {path_str}",
                stats.optimal_bin,
                stats.median
            );
            Some(stats)
        }
        None => {
            log::debug!(
                "INTERSECTS logical: Parquet sampling failed \
                 for {path_str}"
            );
            None
        }
    }
}

/// Read start/end columns from 1–3 representative Parquet row groups
/// and choose the optimal bin size by minimizing replication cost.
///
/// Binary searches for the smallest bin size B such that the mean
/// replication factor `mean(ceil(w_i / B))` is at most
/// `TARGET_MEAN_REPLICATION`. This naturally handles all width
/// distributions: the few wide outliers pull B upward in proportion
/// to their actual cost, without over-sizing bins for the majority.
fn sample_width_stats(
    path: &std::path::Path,
    start_col: &str,
    end_col: &str,
) -> Option<SampledWidthStats> {
    use arrow::array::{Array, Int64Array};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
    use parquet::arrow::ProjectionMask;

    let file = std::fs::File::open(path).ok()?;
    let builder =
        ParquetRecordBatchReaderBuilder::try_new(file).ok()?;

    let parquet_schema = builder.parquet_schema().clone();
    let num_row_groups = builder.metadata().num_row_groups();
    if num_row_groups == 0 {
        return None;
    }

    // Find column indices in the Parquet schema
    let start_idx = parquet_schema
        .columns()
        .iter()
        .position(|c| c.name() == start_col)?;
    let end_idx = parquet_schema
        .columns()
        .iter()
        .position(|c| c.name() == end_col)?;

    // Select representative row groups: first, middle, last
    let mut rg_indices = vec![0];
    if num_row_groups > 2 {
        rg_indices.push(num_row_groups / 2);
    }
    if num_row_groups > 1 {
        rg_indices.push(num_row_groups - 1);
    }

    let mask = ProjectionMask::leaves(
        &parquet_schema,
        [start_idx, end_idx],
    );

    // Cap batch size to bound memory for very large row groups.
    let reader = builder
        .with_projection(mask)
        .with_row_groups(rg_indices)
        .with_batch_size(100_000)
        .build()
        .ok()?;

    let mut widths: Vec<i64> = Vec::new();
    const MAX_SAMPLES: usize = 300_000; // ~3 row groups × 100K

    for batch in reader {
        let batch = batch.ok()?;
        // ProjectionMask preserves original column order, so
        // column 0 is start and column 1 is end (assuming
        // start_idx < end_idx, which holds for all standard
        // genomic schemas).
        let starts = batch
            .column(0)
            .as_any()
            .downcast_ref::<Int64Array>()?;
        let ends = batch
            .column(1)
            .as_any()
            .downcast_ref::<Int64Array>()?;

        for i in 0..batch.num_rows() {
            if !starts.is_null(i) && !ends.is_null(i) {
                let w = ends.value(i) - starts.value(i);
                // Skip malformed intervals where end < start
                if w > 0 {
                    widths.push(w);
                }
            }
        }

        if widths.len() >= MAX_SAMPLES {
            break;
        }
    }

    if widths.is_empty() {
        return None;
    }

    widths.sort_unstable();
    let median = widths[widths.len() / 2];
    let optimal_bin = find_optimal_bin_size(&widths);
    Some(SampledWidthStats {
        optimal_bin,
        median,
    })
}

/// Find the smallest bin size B such that the mean replication
/// factor across all sampled widths is at most TARGET.
///
/// Binary searches over [1, max_width]. For each candidate B,
/// `mean_replication(B) = sum(ceil(w_i / B)) / N`. Since widths
/// are sorted, all w_i <= B contribute ceil=1 — a binary search
/// finds the cutoff index, making each evaluation O(N_above_B)
/// rather than O(N).
fn find_optimal_bin_size(sorted_widths: &[i64]) -> i64 {
    /// Target mean replication factor. Each interval is copied
    /// into ceil(width / bin_size) bins on average. 2.0 means
    /// the average interval spans at most 2 bins — good
    /// selectivity with bounded replication.
    const TARGET_MEAN_REPL: f64 = 2.0;

    let max_width = *sorted_widths.last().unwrap();

    // Binary search: lo always has mean_repl > target,
    // hi always has mean_repl <= target.
    let mut lo: i64 = 1;
    let mut hi: i64 = max_width;

    // At hi = max_width, every interval fits in 1-2 bins, so
    // mean_repl <= 2.0. At lo = 1, mean_repl = mean(widths).
    // If even max_width doesn't meet the target (shouldn't happen
    // since ceil(w/w) <= 2 for all w), return max_width.
    if mean_replication(sorted_widths, hi) > TARGET_MEAN_REPL {
        return hi;
    }

    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if mid == lo {
            break;
        }
        if mean_replication(sorted_widths, mid) <= TARGET_MEAN_REPL {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    hi
}

/// Compute mean(ceil(w_i / B)) for a sorted widths array.
///
/// All widths <= B have ceil = 1. Binary search finds the cutoff,
/// then only iterate the tail above B.
fn mean_replication(sorted_widths: &[i64], bin_size: i64) -> f64 {
    let n = sorted_widths.len();
    // Find first index where width > bin_size
    let cutoff =
        sorted_widths.partition_point(|&w| w <= bin_size);
    // All [0..cutoff) contribute 1 each
    let mut total: i64 = cutoff as i64;
    // [cutoff..n) contribute ceil(w / bin_size) each
    for &w in &sorted_widths[cutoff..] {
        total += (w + bin_size - 1) / bin_size; // ceil division
    }
    total as f64 / n as f64
}

// ── Strategy decision ───────────────────────────────────────────

/// Default bin size when stats are unavailable.
const DEFAULT_BIN_SIZE: usize = 10_000;

/// Choose a bin size for the forced-binning path.
///
/// Uses sampled p95 width from Parquet when available, falling
/// back to a column-level min/max heuristic.
fn choose_bin_size(
    left: &Option<LogicalStats>,
    right: &Option<LogicalStats>,
) -> usize {
    // Tier 1: cost-optimal bin size from Parquet sampling.
    let sampled: Option<&SampledWidthStats> = [left, right]
        .iter()
        .filter_map(|s| s.as_ref()?.sampled.as_ref())
        .next();

    if let Some(stats) = sampled {
        let bin_size = (stats.optimal_bin.max(1) as usize)
            .clamp(1_000, 1_000_000);
        log::debug!(
            "INTERSECTS logical: bin_size={bin_size} \
             (from sampled optimal_bin={})",
            stats.optimal_bin
        );
        return bin_size;
    }

    // Tier 2: column-level min/max heuristic.
    let width_from_stats = |s: &LogicalStats| -> Option<i64> {
        let min_start = s.start_min?;
        let max_start = s.start_max?;
        let min_end = s.end_min?;
        let max_end = s.end_max?;
        let w1 = min_end - min_start;
        let w2 = max_end - max_start;
        Some(w1.max(w2).max(1))
    };

    let l_width = left.as_ref().and_then(width_from_stats);
    let r_width = right.as_ref().and_then(width_from_stats);

    let bin_size = match (l_width, r_width) {
        (Some(l), Some(r)) => {
            (l.max(r).max(1) as usize).clamp(1_000, 1_000_000)
        }
        (Some(w), None) | (None, Some(w)) => {
            (w.max(1) as usize).clamp(1_000, 1_000_000)
        }
        (None, None) => DEFAULT_BIN_SIZE,
    };

    log::debug!(
        "INTERSECTS logical: bin_size={bin_size} \
         (column-level heuristic)"
    );
    bin_size
}

// ── Plan rewrite ────────────────────────────────────────────────

/// Extract the table qualifier from a plan's schema.
///
/// Uses the qualifier of the first column in the plan's output
/// schema, which reflects SQL aliases (e.g., `intervals2` in
/// `FROM intervals JOIN intervals AS intervals2`). This is more
/// robust than walking to the TableScan, which would return the
/// physical table name and miss SQL aliases.
fn get_plan_qualifier(plan: &LogicalPlan) -> Option<String> {
    plan.schema()
        .columns()
        .first()
        .and_then(|c| c.relation.as_ref())
        .map(|r| r.table().to_string())
}

fn rewrite_to_binned(
    join: &Join,
    bin_size: usize,
    start_a: &Column,
    end_a: &Column,
    start_b: &Column,
    end_b: &Column,
    rewritten_filter: Option<&Expr>,
) -> Result<LogicalPlan> {
    let bs = bin_size as i64;

    // Get table qualifiers for aliasing after UNNEST. Uses the
    // schema qualifier (which reflects SQL aliases) so column
    // references in the filter resolve correctly after the rewrite.
    let left_alias = get_plan_qualifier(&join.left)
        .unwrap_or_else(|| "l".to_string());
    let right_alias = get_plan_qualifier(&join.right)
        .unwrap_or_else(|| "r".to_string());

    let left_expanded = expand_with_bins(
        (*join.left).clone(),
        "__giql_bins_l",
        bs,
        &left_alias,
        start_a,
        end_a,
    )?;
    let right_expanded = expand_with_bins(
        (*join.right).clone(),
        "__giql_bins_r",
        bs,
        &right_alias,
        start_b,
        end_b,
    )?;

    // Equi-keys: original keys re-qualified with the aliases +
    // bin columns
    let mut left_keys: Vec<Expr> = join
        .on
        .iter()
        .map(|(l, _)| {
            if let Expr::Column(c) = l {
                Expr::Column(Column::new(
                    Some(left_alias.clone()),
                    &c.name,
                ))
            } else {
                l.clone()
            }
        })
        .collect();
    let mut right_keys: Vec<Expr> = join
        .on
        .iter()
        .map(|(_, r)| {
            if let Expr::Column(c) = r {
                Expr::Column(Column::new(
                    Some(right_alias.clone()),
                    &c.name,
                ))
            } else {
                r.clone()
            }
        })
        .collect();
    left_keys
        .push(col(format!("{left_alias}.__giql_bins_l")));
    right_keys
        .push(col(format!("{right_alias}.__giql_bins_r")));

    // Build the join with the rewritten filter (giql_intersects
    // replaced with real overlap predicates) and extra bin
    // equi-keys.
    let joined = LogicalPlanBuilder::from(left_expanded)
        .join_with_expr_keys(
            right_expanded,
            JoinType::Inner,
            (left_keys, right_keys),
            rewritten_filter.cloned(),
        )?
        .build()?;

    // Add canonical-bin filter to eliminate duplicates from
    // multi-bin matches. For each pair, only emit from the bin
    // that equals the GREATER of the two intervals' first bins.
    // This makes DISTINCT unnecessary.
    let left_first_bin = cast(
        Expr::Column(Column::new(
            Some(left_alias.clone()),
            &start_a.name,
        )),
        arrow::datatypes::DataType::Int64,
    ) / lit(bs);
    let right_first_bin = cast(
        Expr::Column(Column::new(
            Some(right_alias.clone()),
            &start_b.name,
        )),
        arrow::datatypes::DataType::Int64,
    ) / lit(bs);

    // canonical_bin = CASE WHEN l_fb >= r_fb THEN l_fb ELSE r_fb
    let canonical_bin =
        Expr::Case(datafusion::logical_expr::expr::Case {
            expr: None,
            when_then_expr: vec![(
                Box::new(
                    left_first_bin.clone().gt_eq(right_first_bin.clone()),
                ),
                Box::new(left_first_bin),
            )],
            else_expr: Some(Box::new(right_first_bin)),
        });

    let bins_col =
        col(format!("{left_alias}.__giql_bins_l"));
    let dedup_filter = bins_col.eq(canonical_bin);

    let filtered = LogicalPlanBuilder::from(joined)
        .filter(dedup_filter)?
        .build()?;

    // Project to strip bin columns. The canonical-bin filter above
    // already ensures each pair appears exactly once, so DISTINCT
    // is unnecessary.
    let output_exprs: Vec<Expr> = filtered
        .schema()
        .columns()
        .into_iter()
        .filter(|c| {
            c.name != "__giql_bins_l" && c.name != "__giql_bins_r"
        })
        .map(|c| Expr::Column(c))
        .collect();

    LogicalPlanBuilder::from(filtered)
        .project(output_exprs)?
        .build()
}

/// Add a `range(start/B, (end-1)/B + 1)` column, unnest it, and
/// wrap in a SubqueryAlias to preserve the table qualifier for
/// the join filter.
fn expand_with_bins(
    input: LogicalPlan,
    bin_col_name: &str,
    bin_size: i64,
    table_alias: &str,
    start_col: &Column,
    end_col: &Column,
) -> Result<LogicalPlan> {
    let schema = input.schema().clone();

    // Cast start/end to Int64 first, then compute bin boundaries:
    //   range(CAST(start AS BIGINT) / B, (CAST(end AS BIGINT) - 1)
    //   / B + 1)
    let start_i64 = cast(
        Expr::Column(start_col.clone()),
        arrow::datatypes::DataType::Int64,
    );
    let end_i64 = cast(
        Expr::Column(end_col.clone()),
        arrow::datatypes::DataType::Int64,
    );
    let start_bin = start_i64 / lit(bin_size);
    let end_bin =
        (end_i64 - lit(1i64)) / lit(bin_size) + lit(1i64);

    // Build: SELECT *, range(start_bin, end_bin) AS __giql_bins
    // FROM input. Then UNNEST(__giql_bins)
    let range_expr =
        Expr::ScalarFunction(ScalarFunction::new_udf(
            datafusion::functions_nested::range::range_udf(),
            vec![start_bin, end_bin],
        ))
        .alias(bin_col_name);

    let mut proj_exprs: Vec<Expr> = schema
        .columns()
        .into_iter()
        .map(|c| Expr::Column(c))
        .collect();
    proj_exprs.push(range_expr);

    let with_bins = LogicalPlanBuilder::from(input)
        .project(proj_exprs)?
        .build()?;

    // Unnest the bin list column, then re-apply the table alias
    // so that qualified column references (e.g., a.start) in the
    // join filter resolve correctly against this side.
    LogicalPlanBuilder::from(with_bins)
        .unnest_column(bin_col_name)?
        .alias(table_alias)?
        .build()
}
