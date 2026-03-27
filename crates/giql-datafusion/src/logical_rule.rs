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

use crate::IntersectsOptimizerConfig;

/// Logical optimizer rule that rewrites interval overlap joins into
/// binned equi-joins using UNNEST.
///
/// Detects:
///   `JOIN ON a.chrom = b.chrom WHERE a.start < b.end AND a.end > b.start`
///
/// Rewrites to:
///   `SELECT DISTINCT ... FROM Unnest(a + bins) JOIN Unnest(b + bins)
///    ON chrom = chrom AND bin = bin WHERE start < end AND end > start`
///
/// DataFusion handles UNNEST, hash join, and DISTINCT natively with
/// full parallelism. The physical optimizer rule handles sweep-line
/// for heavy-tailed distributions; this rule handles the binned case.
#[derive(Debug)]
pub struct IntersectsLogicalRule {
    config: IntersectsOptimizerConfig,
}

impl IntersectsLogicalRule {
    pub fn new(config: IntersectsOptimizerConfig) -> Self {
        Self { config }
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

        // Detect interval overlap pattern in the filter
        let left_schema = join.left.schema();
        let overlap = match &join.filter {
            Some(filter) => {
                detect_overlap_columns(filter, &left_schema)
            }
            None => None,
        };

        let Some((start_a, end_a, start_b, end_b)) = overlap else {
            return Ok(Transformed::no(plan));
        };

        // Get stats from TableScan children to decide strategy.
        // If stats aren't available (common for ListingTable without
        // collect_statistics), default to binned with a reasonable
        // bin size. The physical optimizer rule will still catch
        // heavy-tailed distributions via Parquet metadata sampling
        // and replace with sweep-line if needed.
        let left_stats = get_table_stats(&join.left);
        let right_stats = get_table_stats(&join.right);

        let bin_size = choose_bin_size(&left_stats, &right_stats);

        log::debug!(
            "INTERSECTS logical rule: rewriting to binned join, \
             bin_size={bin_size}"
        );

        let rewritten =
            rewrite_to_binned(join, bin_size, &start_a, &start_b)?;
        Ok(Transformed::yes(rewritten))
    }
}

// ── Pattern detection ───────────────────────────────────────────

/// Detect interval overlap predicates in a filter expression.
///
/// Looks for `col_a < col_b AND col_c > col_d` where:
/// - One comparison has a "start" col on one side and "end" on the other
/// - The other has "end" on one side and "start" on the other
///
/// Returns `(start_a, end_a, start_b, end_b)` column names.
/// Detect interval overlap predicates in a filter expression.
///
/// Checks the join's left child schema to determine which columns
/// belong to which side — no heuristics based on table name.
fn detect_overlap_columns(
    expr: &Expr,
    left_schema: &datafusion::common::DFSchemaRef,
) -> Option<(Column, Column, Column, Column)> {
    let Expr::BinaryExpr(BinaryExpr {
        left,
        op: Operator::And,
        right,
    }) = expr
    else {
        return None;
    };

    try_extract_overlap(left, right, left_schema)
        .or_else(|| try_extract_overlap(right, left, left_schema))
}

fn try_extract_overlap(
    pred_a: &Expr,
    pred_b: &Expr,
    left_schema: &datafusion::common::DFSchemaRef,
) -> Option<(Column, Column, Column, Column)> {
    let (lt_left, lt_right) = extract_comparison(pred_a, Operator::Lt)?;
    let (gt_left, gt_right) = extract_comparison(pred_b, Operator::Gt)?;

    let all = [&lt_left, &lt_right, &gt_left, &gt_right];

    let is_left = |c: &Column| column_in_schema(c, left_schema);

    let left_start = all.iter().find(|c| is_start(&c.name) && is_left(c));
    let left_end = all.iter().find(|c| is_end(&c.name) && is_left(c));
    let right_start = all.iter().find(|c| is_start(&c.name) && !is_left(c));
    let right_end = all.iter().find(|c| is_end(&c.name) && !is_left(c));

    Some((
        (*left_start?).clone(),
        (*left_end?).clone(),
        (*right_start?).clone(),
        (*right_end?).clone(),
    ))
}

fn extract_comparison(
    expr: &Expr,
    expected_op: Operator,
) -> Option<(Column, Column)> {
    let Expr::BinaryExpr(BinaryExpr { left, op, right }) = expr else {
        return None;
    };
    if *op != expected_op {
        return None;
    }
    let left_col = extract_column(left)?;
    let right_col = extract_column(right)?;
    Some((left_col, right_col))
}

/// Extract a Column from an Expr, handling TryCast/Cast wrappers
/// that DataFusion may insert.
fn extract_column(expr: &Expr) -> Option<Column> {
    match expr {
        Expr::Column(c) => Some(c.clone()),
        Expr::Cast(cast) => extract_column(&cast.expr),
        Expr::TryCast(tc) => extract_column(&tc.expr),
        _ => None,
    }
}

fn is_start(name: &str) -> bool {
    let n = name.to_lowercase();
    n == "start" || n == "chromstart" || n == "pos_start" || n == "begin"
}

fn is_end(name: &str) -> bool {
    let n = name.to_lowercase();
    n == "end" || n == "chromend" || n == "pos_end" || n == "stop"
}

/// Check whether a column belongs to a schema by matching its
/// qualified name against the schema's columns.
fn column_in_schema(
    col: &Column,
    schema: &datafusion::common::DFSchemaRef,
) -> bool {
    schema.has_column(col)
}

// ── Stats collection ────────────────────────────────────────────

struct LogicalStats {
    row_count: Option<usize>,
    start_min: Option<i64>,
    start_max: Option<i64>,
    end_min: Option<i64>,
    end_max: Option<i64>,
}

fn get_table_stats(plan: &LogicalPlan) -> Option<LogicalStats> {
    match plan {
        LogicalPlan::TableScan(ts) => {
            let provider = match source_as_provider(&ts.source) {
                Ok(p) => p,
                Err(e) => {
                    log::debug!(
                        "  get_table_stats: source_as_provider failed: {e}"
                    );
                    return None;
                }
            };

            // Try TableProvider::statistics() first
            if let Some(stats) = provider.statistics() {
                return stats_to_logical(&stats, &ts.source.schema());
            }

            // Fall back to reading Parquet metadata directly
            let listing = match provider
                .as_any()
                .downcast_ref::<ListingTable>()
            {
                Some(lt) => lt,
                None => {
                    log::debug!(
                        "  get_table_stats: not a ListingTable: {}",
                        std::any::type_name_of_val(provider.as_ref()),
                    );
                    return None;
                }
            };
            let paths = listing.table_paths();
            let path = paths.first()?;
            let path_str = path.as_str();

            // ListingTableUrl stores file:// URLs; extract the
            // filesystem path
            let fs_path = if let Some(p) = path_str.strip_prefix("file://") {
                std::path::PathBuf::from(p)
            } else {
                std::path::PathBuf::from(format!("/{path_str}"))
            };

            let schema = ts.source.schema();
            let start_col = schema
                .fields()
                .iter()
                .find(|f| is_start(f.name()))?
                .name()
                .as_str();
            let end_col = schema
                .fields()
                .iter()
                .find(|f| is_end(f.name()))?
                .name()
                .as_str();

            // Read Parquet metadata (file footer only — fast)
            let meta = crate::stats::metadata::collect_metadata(
                &fs_path, start_col, end_col,
            )
            .ok()?;

            // Aggregate per-row-group bounds
            let start_min = meta.row_group_bounds.iter()
                .map(|rg| rg.min_start)
                .min();
            let start_max = meta.row_group_bounds.iter()
                .map(|rg| rg.max_start)
                .max();
            let end_min = meta.row_group_bounds.iter()
                .map(|rg| rg.min_end)
                .min();
            let end_max = meta.row_group_bounds.iter()
                .map(|rg| rg.max_end)
                .max();

            Some(LogicalStats {
                row_count: Some(meta.total_rows),
                start_min,
                start_max,
                end_min,
                end_max,
            })
        }
        _ => plan
            .inputs()
            .first()
            .and_then(|child| get_table_stats(child)),
    }
}

fn stats_to_logical(
    stats: &datafusion::common::Statistics,
    schema: &arrow::datatypes::SchemaRef,
) -> Option<LogicalStats> {
    let row_count = match stats.num_rows {
        datafusion::common::stats::Precision::Exact(n) => Some(n),
        datafusion::common::stats::Precision::Inexact(n) => Some(n),
        _ => None,
    };
    let start_idx = schema
        .fields()
        .iter()
        .position(|f| is_start(f.name()))?;
    let end_idx = schema
        .fields()
        .iter()
        .position(|f| is_end(f.name()))?;
    let col_stats = &stats.column_statistics;
    let start_stats = col_stats.get(start_idx)?;
    let end_stats = col_stats.get(end_idx)?;
    Some(LogicalStats {
        row_count,
        start_min: scalar_to_i64(&start_stats.min_value),
        start_max: scalar_to_i64(&start_stats.max_value),
        end_min: scalar_to_i64(&end_stats.min_value),
        end_max: scalar_to_i64(&end_stats.max_value),
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

// ── Strategy decision ───────────────────────────────────────────

/// Default bin size when stats are unavailable.
const DEFAULT_BIN_SIZE: usize = 10_000;

/// Choose a bin size from Parquet metadata.
///
/// The width signal `max(end) - max(start)` approximates the width
/// of the widest intervals. We use this as the bin size so that most
/// intervals fit in a single bin (replication factor ≈ 1).
///
/// Always returns `Some(bin_size)` — binned joins are correct for
/// all distributions, and DataFusion's native UNNEST + hash join is
/// fast enough that adaptive bin sizing is the only knob needed.
fn choose_bin_size(
    left: &Option<LogicalStats>,
    right: &Option<LogicalStats>,
) -> usize {
    let width_from_stats = |s: &LogicalStats| -> Option<i64> {
        let max_start = s.start_max?;
        let max_end = s.end_max?;
        Some((max_end - max_start).max(1))
    };

    let l_width = left.as_ref().and_then(width_from_stats);
    let r_width = right.as_ref().and_then(width_from_stats);

    match (l_width, r_width) {
        (Some(l), Some(r)) => {
            let w = l.max(r) as usize;
            let bin_size = w.clamp(1_000, 1_000_000);
            log::debug!(
                "INTERSECTS logical: adaptive bin_size={bin_size} \
                 (from widths l={l}, r={r})"
            );
            bin_size
        }
        (Some(w), None) | (None, Some(w)) => {
            let bin_size = (w as usize).clamp(1_000, 1_000_000);
            log::debug!(
                "INTERSECTS logical: adaptive bin_size={bin_size} \
                 (partial stats, width={w})"
            );
            bin_size
        }
        (None, None) => {
            log::debug!(
                "INTERSECTS logical: using default bin_size={DEFAULT_BIN_SIZE}"
            );
            DEFAULT_BIN_SIZE
        }
    }
}

// ── Plan rewrite ────────────────────────────────────────────────

/// Extract table name from a logical plan (walks to TableScan).
fn get_table_name(plan: &LogicalPlan) -> Option<String> {
    match plan {
        LogicalPlan::TableScan(ts) => {
            Some(ts.table_name.table().to_string())
        }
        _ => plan
            .inputs()
            .first()
            .and_then(|child| get_table_name(child)),
    }
}

fn rewrite_to_binned(
    join: &Join,
    bin_size: usize,
    start_a: &Column,
    start_b: &Column,
) -> Result<LogicalPlan> {
    let bs = bin_size as i64;

    // Get table names for aliasing after UNNEST
    let left_alias = get_table_name(&join.left)
        .unwrap_or_else(|| "l".to_string());
    let right_alias = get_table_name(&join.right)
        .unwrap_or_else(|| "r".to_string());

    let left_expanded = expand_with_bins(
        (*join.left).clone(),
        "__giql_bins_l",
        bs,
        &left_alias,
    )?;
    let right_expanded = expand_with_bins(
        (*join.right).clone(),
        "__giql_bins_r",
        bs,
        &right_alias,
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
    left_keys.push(col(format!(
        "{left_alias}.__giql_bins_l"
    )));
    right_keys.push(col(format!(
        "{right_alias}.__giql_bins_r"
    )));

    // Build the join with the original filter and extra bin equi-keys.
    // Wrap in a subquery alias to isolate the schema, then project
    // away the bin columns and add DISTINCT.
    let joined = LogicalPlanBuilder::from(left_expanded)
        .join_with_expr_keys(
            right_expanded,
            JoinType::Inner,
            (left_keys, right_keys),
            join.filter.clone(),
        )?
        .build()?;

    // Add canonical-bin filter to eliminate duplicates from
    // multi-bin matches. For each pair, only emit from the bin
    // that equals the GREATER of the two intervals' first bins.
    // This makes DISTINCT unnecessary.
    //
    // We use GREATEST(left_first_bin, right_first_bin) but since
    // DataFusion doesn't have GREATEST, we use:
    //   CASE WHEN left_first_bin >= right_first_bin
    //        THEN left_first_bin
    //        ELSE right_first_bin END
    //
    // Actually, we use the simpler approach: just keep the row
    // where __giql_bins equals the max of (start_a/B, start_b/B).
    // Since we already have the start columns, we can compute this.
    let left_first_bin = cast(
        Expr::Column(start_a.clone()),
        arrow::datatypes::DataType::Int64,
    ) / lit(bs);
    let right_first_bin = cast(
        Expr::Column(start_b.clone()),
        arrow::datatypes::DataType::Int64,
    ) / lit(bs);

    // canonical_bin = CASE WHEN l_fb >= r_fb THEN l_fb ELSE r_fb END
    let canonical_bin = Expr::Case(datafusion::logical_expr::expr::Case {
        expr: None,
        when_then_expr: vec![(
            Box::new(left_first_bin.clone().gt_eq(right_first_bin.clone())),
            Box::new(left_first_bin),
        )],
        else_expr: Some(Box::new(right_first_bin)),
    });

    // We need the bins column from the left side. After the join
    // with aliases, the left bin column is qualified.
    let bins_col = col(format!("{left_alias}.__giql_bins_l"));
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
        .filter(|c| !c.name.starts_with("__giql_bins"))
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
) -> Result<LogicalPlan> {
    let schema = input.schema().clone();

    // Find start and end columns
    let start_col = schema
        .columns()
        .into_iter()
        .find(|c| is_start(&c.name))
        .ok_or_else(|| {
            datafusion::error::DataFusionError::Plan(
                "No start column found".to_string(),
            )
        })?;
    let end_col = schema
        .columns()
        .into_iter()
        .find(|c| is_end(&c.name))
        .ok_or_else(|| {
            datafusion::error::DataFusionError::Plan(
                "No end column found".to_string(),
            )
        })?;

    // Cast start/end to Int64 first, then compute bin boundaries:
    //   range(CAST(start AS BIGINT) / B, (CAST(end AS BIGINT) - 1) / B + 1)
    let start_i64 = cast(
        Expr::Column(start_col),
        arrow::datatypes::DataType::Int64,
    );
    let end_i64 = cast(
        Expr::Column(end_col),
        arrow::datatypes::DataType::Int64,
    );
    let start_bin = start_i64 / lit(bin_size);
    let end_bin = (end_i64 - lit(1i64)) / lit(bin_size) + lit(1i64);

    // Build: SELECT *, range(start_bin, end_bin) AS __giql_bins FROM input
    // Then UNNEST(__giql_bins)
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
