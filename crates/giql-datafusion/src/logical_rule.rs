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
#[derive(Debug)]
pub struct IntersectsLogicalRule;

impl IntersectsLogicalRule {
    pub fn new() -> Self {
        Self
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

        // Get stats from TableScan children for adaptive bin sizing.
        let left_stats =
            get_table_stats(&join.left, &start_a.name, &end_a.name);
        let right_stats =
            get_table_stats(&join.right, &start_b.name, &end_b.name);

        let bin_size = choose_bin_size(&left_stats, &right_stats);

        log::debug!(
            "INTERSECTS logical rule: rewriting to binned join, \
             bin_size={bin_size}"
        );

        // Replace giql_intersects() with real overlap predicates
        // before building the binned join, since the placeholder
        // UDF cannot execute.
        let rewritten_filter = join.filter.as_ref().map(|f| {
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
        Ok(Transformed::yes(rewritten))
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
    #[allow(dead_code)]
    row_count: Option<usize>,
    start_min: Option<i64>,
    start_max: Option<i64>,
    end_min: Option<i64>,
    end_max: Option<i64>,
    /// p95 interval width from lightweight Parquet sampling.
    sampled_width_p95: Option<i64>,
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
                    sampled_width_p95: None,
                });

            // Try lightweight Parquet sampling for accurate width
            stats.sampled_width_p95 = try_sample_from_listing(
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
        sampled_width_p95: None, // filled by try_sample_from_listing
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
) -> Option<i64> {
    let listing = provider.as_any().downcast_ref::<ListingTable>()?;
    let path_str = listing.table_paths().first()?.as_str();

    // ListingTableUrl stores file:// URLs
    let fs_path = if let Some(p) = path_str.strip_prefix("file://") {
        std::path::PathBuf::from(p)
    } else {
        std::path::PathBuf::from(format!("/{path_str}"))
    };

    match sample_width_p95(&fs_path, start_col, end_col) {
        Some(p95) => {
            log::debug!(
                "INTERSECTS logical: sampled p95 width={p95} \
                 from {path_str}"
            );
            Some(p95)
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
/// and return the p95 interval width.
fn sample_width_p95(
    path: &std::path::Path,
    start_col: &str,
    end_col: &str,
) -> Option<i64> {
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

    let reader = builder
        .with_projection(mask)
        .with_row_groups(rg_indices)
        .build()
        .ok()?;

    let mut widths: Vec<i64> = Vec::new();

    for batch in reader {
        let batch = batch.ok()?;
        // After projection, columns are at indices 0 and 1
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
                widths.push(ends.value(i) - starts.value(i));
            }
        }
    }

    if widths.is_empty() {
        return None;
    }

    widths.sort_unstable();
    let p95_idx = (widths.len() * 95 / 100).min(widths.len() - 1);
    Some(widths[p95_idx])
}

// ── Strategy decision ───────────────────────────────────────────

/// Default bin size when stats are unavailable.
const DEFAULT_BIN_SIZE: usize = 10_000;

/// Choose a bin size from table statistics.
///
/// Estimates representative interval width from column-level min/max
/// stats using two independent signals:
///
/// - `min(end) - min(start)`: width estimate from the leftmost
///   interval (the row with the smallest start likely ends near
///   `min(end)`)
/// - `max(end) - max(start)`: width estimate from the rightmost
///   interval (the row with the largest start likely ends near
///   `max(end)`)
///
/// Takes the max of both to be conservative — a larger bin size
/// means fewer bins per interval, avoiding replication blowup at the
/// cost of more false-positive bin matches (filtered by the overlap
/// predicate).
fn choose_bin_size(
    left: &Option<LogicalStats>,
    right: &Option<LogicalStats>,
) -> usize {
    // Tier 1: Use sampled p95 width if available from either side.
    // This reads actual interval widths from Parquet row groups and
    // is robust against all endpoint distributions.
    let sampled = [left, right]
        .iter()
        .filter_map(|s| s.as_ref()?.sampled_width_p95)
        .max();

    if let Some(p95) = sampled {
        let bin_size =
            (p95.max(1) as usize).clamp(1_000, 1_000_000);
        log::debug!(
            "INTERSECTS logical: bin_size={bin_size} \
             (from sampled p95 width={p95})"
        );
        return bin_size;
    }

    // Tier 2: Fall back to column-level min/max heuristic.
    let width_from_stats = |s: &LogicalStats| -> Option<i64> {
        let min_start = s.start_min?;
        let max_start = s.start_max?;
        let min_end = s.end_min?;
        let max_end = s.end_max?;
        // Two independent width estimates; take the max.
        let w1 = min_end - min_start;
        let w2 = max_end - max_start;
        Some(w1.max(w2).max(1))
    };

    let l_width = left.as_ref().and_then(width_from_stats);
    let r_width = right.as_ref().and_then(width_from_stats);

    match (l_width, r_width) {
        (Some(l), Some(r)) => {
            let w = l.max(r).max(1) as usize;
            let bin_size = w.clamp(1_000, 1_000_000);
            log::debug!(
                "INTERSECTS logical: adaptive bin_size={bin_size} \
                 (from widths l={l}, r={r})"
            );
            bin_size
        }
        (Some(w), None) | (None, Some(w)) => {
            let bin_size = (w.max(1) as usize).clamp(1_000, 1_000_000);
            log::debug!(
                "INTERSECTS logical: adaptive bin_size={bin_size} \
                 (partial stats, width={w})"
            );
            bin_size
        }
        (None, None) => {
            log::debug!(
                "INTERSECTS logical: using default \
                 bin_size={DEFAULT_BIN_SIZE}"
            );
            DEFAULT_BIN_SIZE
        }
    }
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
