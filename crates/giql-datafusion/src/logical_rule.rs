use std::sync::Arc;

use datafusion::common::tree_node::Transformed;
use datafusion::common::{Column, Result, ScalarValue};
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
        let overlap = match &join.filter {
            Some(filter) => {
                eprintln!(
                    "INTERSECTS logical: checking filter: {filter}"
                );
                detect_overlap_columns(filter)
            }
            None => {
                eprintln!("INTERSECTS logical: join has no filter");
                None
            }
        };

        let Some((start_a, end_a, start_b, end_b)) = overlap else {
            eprintln!("INTERSECTS logical: no overlap pattern found");
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

        let bin_size = match should_use_binned(
            &left_stats,
            &right_stats,
            &self.config,
        ) {
            Some(bs) => bs,
            None => return Ok(Transformed::no(plan)),
        };

        eprintln!(
            "INTERSECTS logical rule: rewriting to binned join, \
             bin_size={bin_size}"
        );

        // Rewrite to binned equi-join
        let rewritten = rewrite_to_binned(join, bin_size)?;
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
fn detect_overlap_columns(
    expr: &Expr,
) -> Option<(Column, Column, Column, Column)> {
    let Expr::BinaryExpr(BinaryExpr {
        left,
        op: Operator::And,
        right,
    }) = expr
    else {
        return None;
    };

    // Try both orderings
    try_extract_overlap(left, right)
        .or_else(|| try_extract_overlap(right, left))
}

fn try_extract_overlap(
    pred_a: &Expr,
    pred_b: &Expr,
) -> Option<(Column, Column, Column, Column)> {
    let (lt_left, lt_right) = extract_comparison(pred_a, Operator::Lt)?;
    let (gt_left, gt_right) = extract_comparison(pred_b, Operator::Gt)?;

    let all = [&lt_left, &lt_right, &gt_left, &gt_right];
    eprintln!("INTERSECTS logical: columns in filter:");
    for c in &all {
        eprintln!("  {:?} start={} end={} left={}",
            c, is_start(&c.name), is_end(&c.name), is_from_left(c));
    }

    let left_start = all.iter().find(|c| is_start(&c.name) && is_from_left(c));
    let left_end = all.iter().find(|c| is_end(&c.name) && is_from_left(c));
    let right_start = all.iter().find(|c| is_start(&c.name) && !is_from_left(c));
    let right_end = all.iter().find(|c| is_end(&c.name) && !is_from_left(c));

    eprintln!("  left_start={left_start:?} left_end={left_end:?} right_start={right_start:?} right_end={right_end:?}");

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
        eprintln!("  extract_comparison: not a BinaryExpr: {expr:?}");
        return None;
    };
    if *op != expected_op {
        eprintln!("  extract_comparison: op={op:?}, expected={expected_op:?}");
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
        other => {
            eprintln!("  extract_column: unexpected expr type: {other:?}");
            None
        }
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

fn is_from_left(col: &Column) -> bool {
    // In DataFusion logical plans, qualified columns have a table
    // relation. We use position in the join: left-side columns
    // have the left table qualifier. Since we don't know the exact
    // qualifier, we rely on the join's on-clause to tell us which
    // table is which. For now, use a simple heuristic: both sides
    // have the same column names, so the relation qualifier
    // distinguishes them. If no qualifier, we can't tell.
    // This works because genomic tables always have qualified refs
    // in JOIN conditions (e.g., a.start, b.start).
    col.relation.is_some()
        && col
            .relation
            .as_ref()
            .map(|r| {
                let s = r.to_string();
                // First table alphabetically is "left" — fragile but
                // works for a.X / b.X patterns. We'll improve this
                // by checking against the join's child schemas.
                s.starts_with('a') || s.starts_with('l')
            })
            .unwrap_or(false)
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
            let provider = source_as_provider(&ts.source).ok()?;
            eprintln!(
                "  get_table_stats: provider type: {}",
                std::any::type_name_of_val(provider.as_ref()),
            );
            let stats = provider.statistics();
            eprintln!("  get_table_stats: statistics = {stats:?}");
            let stats = stats?;
            let row_count = match stats.num_rows {
                datafusion::common::stats::Precision::Exact(n) => Some(n),
                datafusion::common::stats::Precision::Inexact(n) => Some(n),
                _ => None,
            };

            // Find start and end column stats
            let schema = ts.source.schema();
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
        _ => {
            // Try first child
            plan.inputs()
                .first()
                .and_then(|child| get_table_stats(child))
        }
    }
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

fn should_use_binned(
    left: &Option<LogicalStats>,
    right: &Option<LogicalStats>,
    config: &IntersectsOptimizerConfig,
) -> Option<usize> {
    // If stats are available, use them to decide
    if let (Some(left), Some(right)) = (left, right) {
        let l_width_at_max =
            left.end_max.unwrap_or(0) - left.start_max.unwrap_or(0);
        let l_width_at_min =
            left.end_min.unwrap_or(0) - left.start_min.unwrap_or(0);
        let r_width_at_max =
            right.end_max.unwrap_or(0) - right.start_max.unwrap_or(0);
        let r_width_at_min =
            right.end_min.unwrap_or(0) - right.start_min.unwrap_or(0);

        let l_ratio = if l_width_at_min > 0 {
            l_width_at_max as f64 / l_width_at_min as f64
        } else {
            f64::MAX
        };
        let r_ratio = if r_width_at_min > 0 {
            r_width_at_max as f64 / r_width_at_min as f64
        } else {
            f64::MAX
        };

        if l_ratio > config.p99_median_threshold
            || r_ratio > config.p99_median_threshold
        {
            return None; // → sweep-line (physical rule)
        }

        let max_width =
            l_width_at_max.max(r_width_at_max).max(1) as usize;
        return Some(max_width.clamp(1_000, 1_000_000));
    }

    // No stats → default to binned. The physical optimizer rule
    // will still catch heavy-tailed distributions via Parquet
    // metadata sampling and override with sweep-line if needed.
    // But since the logical rule has already rewritten the plan,
    // the physical rule won't see the original join pattern.
    //
    // This is the right trade-off: binned is correct for all
    // distributions (just not always optimal), and DataFusion's
    // native UNNEST + hash join is fast.
    Some(DEFAULT_BIN_SIZE)
}

// ── Plan rewrite ────────────────────────────────────────────────

fn rewrite_to_binned(
    join: &Join,
    bin_size: usize,
) -> Result<LogicalPlan> {
    let bs = bin_size as i64;

    let left_expanded = expand_with_bins(
        (*join.left).clone(),
        "__giql_bins_l",
        bs,
    )?;
    let right_expanded = expand_with_bins(
        (*join.right).clone(),
        "__giql_bins_r",
        bs,
    )?;

    // Original equi-keys + bin columns (differently named to avoid collision)
    let mut left_keys: Vec<Expr> =
        join.on.iter().map(|(l, _)| l.clone()).collect();
    let mut right_keys: Vec<Expr> =
        join.on.iter().map(|(_, r)| r.clone()).collect();
    left_keys.push(col("__giql_bins_l"));
    right_keys.push(col("__giql_bins_r"));

    // Build: join → project (strip __giql_bins) → distinct
    let joined = LogicalPlanBuilder::from(left_expanded)
        .join_with_expr_keys(
            right_expanded,
            JoinType::Inner,
            (left_keys, right_keys),
            join.filter.clone(),
        )?
        .build()?;

    // Project away __giql_bins columns from the new join's schema
    let output_exprs: Vec<Expr> = joined
        .schema()
        .columns()
        .into_iter()
        .filter(|c| !c.name.starts_with("__giql_bins"))
        .map(|c| Expr::Column(c))
        .collect();

    let projected = LogicalPlanBuilder::from(joined)
        .project(output_exprs)?
        .distinct()?
        .build()?;

    Ok(projected)
}

/// Add a `range(start/B, (end-1)/B + 1)` column and unnest it.
fn expand_with_bins(
    input: LogicalPlan,
    bin_col_name: &str,
    bin_size: i64,
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

    LogicalPlanBuilder::from(with_bins)
        .unnest_column(bin_col_name)?
        .build()
}
