use std::path::PathBuf;
use std::sync::Arc;

use arrow::datatypes::SchemaRef;
use datafusion::common::{JoinSide, Result};
use datafusion::physical_plan::joins::utils::ColumnIndex;
use datafusion::physical_plan::ExecutionPlan;

/// Column indices for the genomic interval columns on one side of a
/// join.
#[derive(Debug, Clone)]
pub struct IntervalColumns {
    /// Column name for chromosome.
    pub chrom_col: String,
    /// Column name for interval start.
    pub start_col: String,
    /// Column name for interval end.
    pub end_col: String,
    /// Column index for chromosome in the side's schema.
    pub chrom_idx: usize,
    /// Column index for start in the side's schema.
    pub start_idx: usize,
    /// Column index for end in the side's schema.
    pub end_idx: usize,
}

/// A detected interval overlap join pattern in the physical plan.
#[derive(Debug)]
pub struct IntervalJoinMatch {
    /// The left input execution plan.
    pub left: Arc<dyn ExecutionPlan>,
    /// The right input execution plan.
    pub right: Arc<dyn ExecutionPlan>,
    /// Interval column indices on the left side.
    pub left_cols: IntervalColumns,
    /// Interval column indices on the right side.
    pub right_cols: IntervalColumns,
    /// Output schema of the join node being replaced.
    pub output_schema: SchemaRef,
    /// Parquet file paths from the left source, if available.
    pub left_parquet_paths: Vec<PathBuf>,
    /// Parquet file paths from the right source, if available.
    pub right_parquet_paths: Vec<PathBuf>,
}

/// Attempt to detect an interval overlap join in the given execution
/// plan node.
///
/// Looks for join nodes (hash join, nested loop join) with predicates
/// matching the pattern:
///
/// ```text
/// left.chrom = right.chrom
///   AND left.start < right.end
///   AND left.end > right.start
/// ```
pub fn detect_interval_join(
    plan: &Arc<dyn ExecutionPlan>,
) -> Result<Option<IntervalJoinMatch>> {
    let plan_name = plan.name();

    match plan_name {
        "NestedLoopJoinExec" => detect_from_nested_loop_join(plan),
        "HashJoinExec" => detect_from_hash_join(plan),
        _ => Ok(None),
    }
}

/// Detect interval overlap in a NestedLoopJoinExec.
fn detect_from_nested_loop_join(
    _plan: &Arc<dyn ExecutionPlan>,
) -> Result<Option<IntervalJoinMatch>> {
    // NestedLoopJoinExec detection requires walking the full filter
    // expression tree to find all three predicates (chrom equality +
    // range overlap). This is deferred to a future iteration — the
    // HashJoinExec path handles the common case where DataFusion
    // separates the equi-key (chrom) from the range filter.
    Ok(None)
}

/// Detect interval overlap in a HashJoinExec.
fn detect_from_hash_join(
    plan: &Arc<dyn ExecutionPlan>,
) -> Result<Option<IntervalJoinMatch>> {
    use datafusion::physical_plan::joins::HashJoinExec;

    let hj = match plan.as_any().downcast_ref::<HashJoinExec>() {
        Some(hj) => hj,
        None => return Ok(None),
    };

    let filter = match hj.filter() {
        Some(f) => f,
        None => return Ok(None),
    };

    let left = hj.left().clone();
    let right = hj.right().clone();
    let left_schema = left.schema();
    let right_schema = right.schema();

    // Extract the equi-join key — should be a chromosome column
    let on = hj.on();
    if on.is_empty() {
        return Ok(None);
    }

    use datafusion::physical_expr::expressions::Column;

    let (left_chrom_key, right_chrom_key) = &on[0];
    let left_chrom_col =
        match left_chrom_key.as_any().downcast_ref::<Column>() {
            Some(c) => c,
            None => return Ok(None),
        };
    let right_chrom_col =
        match right_chrom_key.as_any().downcast_ref::<Column>() {
            Some(c) => c,
            None => return Ok(None),
        };
    let left_chrom_name = left_chrom_col.name().to_string();
    let right_chrom_name = right_chrom_col.name().to_string();

    let left_chrom_idx = match left_schema
        .fields()
        .iter()
        .position(|f| f.name() == &left_chrom_name)
    {
        Some(idx) => idx,
        None => return Ok(None),
    };
    let right_chrom_idx = match right_schema
        .fields()
        .iter()
        .position(|f| f.name() == &right_chrom_name)
    {
        Some(idx) => idx,
        None => return Ok(None),
    };

    // Extract start/end columns from the filter expression
    let filter_expr = filter.expression();
    let column_indices = filter.column_indices();
    let range_cols = match extract_range_columns_from_filter(
        filter_expr,
        column_indices,
        &left_schema,
        &right_schema,
    ) {
        Some(cols) => cols,
        None => return Ok(None),
    };

    let left_start_idx = match left_schema
        .fields()
        .iter()
        .position(|f| f.name() == &range_cols.0)
    {
        Some(idx) => idx,
        None => return Ok(None),
    };
    let left_end_idx = match left_schema
        .fields()
        .iter()
        .position(|f| f.name() == &range_cols.1)
    {
        Some(idx) => idx,
        None => return Ok(None),
    };
    let right_start_idx = match right_schema
        .fields()
        .iter()
        .position(|f| f.name() == &range_cols.2)
    {
        Some(idx) => idx,
        None => return Ok(None),
    };
    let right_end_idx = match right_schema
        .fields()
        .iter()
        .position(|f| f.name() == &range_cols.3)
    {
        Some(idx) => idx,
        None => return Ok(None),
    };

    let left_cols = IntervalColumns {
        chrom_col: left_chrom_name,
        start_col: range_cols.0,
        end_col: range_cols.1,
        chrom_idx: left_chrom_idx,
        start_idx: left_start_idx,
        end_idx: left_end_idx,
    };

    let right_cols = IntervalColumns {
        chrom_col: right_chrom_name,
        start_col: range_cols.2,
        end_col: range_cols.3,
        chrom_idx: right_chrom_idx,
        start_idx: right_start_idx,
        end_idx: right_end_idx,
    };

    let left_parquet_paths = find_parquet_paths(&left);
    let right_parquet_paths = find_parquet_paths(&right);

    Ok(Some(IntervalJoinMatch {
        left,
        right,
        left_cols,
        right_cols,
        output_schema: plan.schema(),
        left_parquet_paths,
        right_parquet_paths,
    }))
}

/// Extract range column names from a HashJoin filter expression.
///
/// The filter should contain `left.start < right.end AND left.end >
/// right.start`. Returns `Some((left_start, left_end, right_start,
/// right_end))` column names if the pattern matches.
fn extract_range_columns_from_filter(
    expr: &Arc<dyn datafusion::physical_plan::PhysicalExpr>,
    column_indices: &[ColumnIndex],
    left_schema: &SchemaRef,
    right_schema: &SchemaRef,
) -> Option<(String, String, String, String)> {
    use datafusion::logical_expr::Operator;
    use datafusion::physical_expr::expressions::BinaryExpr;

    let binary = expr.as_any().downcast_ref::<BinaryExpr>()?;

    if *binary.op() != Operator::And {
        return None;
    }

    let left_pred = binary.left();
    let right_pred = binary.right();

    // Try both orderings of the two predicates
    

    try_extract_range_pair(
        left_pred,
        right_pred,
        column_indices,
        left_schema,
        right_schema,
    )
    .or_else(|| {
        try_extract_range_pair(
            right_pred,
            left_pred,
            column_indices,
            left_schema,
            right_schema,
        )
    })
}

/// Try to extract (left_start, left_end, right_start, right_end) from
/// a pair of predicates where one is Lt and one is Gt.
fn try_extract_range_pair(
    pred_a: &Arc<dyn datafusion::physical_plan::PhysicalExpr>,
    pred_b: &Arc<dyn datafusion::physical_plan::PhysicalExpr>,
    column_indices: &[ColumnIndex],
    left_schema: &SchemaRef,
    right_schema: &SchemaRef,
) -> Option<(String, String, String, String)> {
    let (lt_left, lt_right) = extract_lt_columns(
        pred_a,
        column_indices,
        left_schema,
        right_schema,
    )?;
    let (gt_left, gt_right) = extract_gt_columns(
        pred_b,
        column_indices,
        left_schema,
        right_schema,
    )?;

    // lt pattern: left.start < right.end
    // gt pattern: left.end > right.start
    Some((lt_left, gt_left, gt_right, lt_right))
}

/// Extract columns from a `<` comparison.
fn extract_lt_columns(
    expr: &Arc<dyn datafusion::physical_plan::PhysicalExpr>,
    column_indices: &[ColumnIndex],
    left_schema: &SchemaRef,
    right_schema: &SchemaRef,
) -> Option<(String, String)> {
    use datafusion::logical_expr::Operator;
    use datafusion::physical_expr::expressions::{BinaryExpr, Column};

    let binary = expr.as_any().downcast_ref::<BinaryExpr>()?;
    if *binary.op() != Operator::Lt {
        return None;
    }

    let left_col = binary.left().as_any().downcast_ref::<Column>()?;
    let right_col =
        binary.right().as_any().downcast_ref::<Column>()?;

    let left_name = resolve_column_name(
        left_col.index(),
        column_indices,
        left_schema,
        right_schema,
        true,
    )?;
    let right_name = resolve_column_name(
        right_col.index(),
        column_indices,
        left_schema,
        right_schema,
        false,
    )?;

    Some((left_name, right_name))
}

/// Extract columns from a `>` comparison.
fn extract_gt_columns(
    expr: &Arc<dyn datafusion::physical_plan::PhysicalExpr>,
    column_indices: &[ColumnIndex],
    left_schema: &SchemaRef,
    right_schema: &SchemaRef,
) -> Option<(String, String)> {
    use datafusion::logical_expr::Operator;
    use datafusion::physical_expr::expressions::{BinaryExpr, Column};

    let binary = expr.as_any().downcast_ref::<BinaryExpr>()?;
    if *binary.op() != Operator::Gt {
        return None;
    }

    let left_col = binary.left().as_any().downcast_ref::<Column>()?;
    let right_col =
        binary.right().as_any().downcast_ref::<Column>()?;

    let left_name = resolve_column_name(
        left_col.index(),
        column_indices,
        left_schema,
        right_schema,
        true,
    )?;
    let right_name = resolve_column_name(
        right_col.index(),
        column_indices,
        left_schema,
        right_schema,
        false,
    )?;

    Some((left_name, right_name))
}

/// Resolve a filter-local column index to a column name.
fn resolve_column_name(
    filter_idx: usize,
    column_indices: &[ColumnIndex],
    left_schema: &SchemaRef,
    right_schema: &SchemaRef,
    expect_left: bool,
) -> Option<String> {
    if filter_idx >= column_indices.len() {
        return None;
    }

    let col_idx = &column_indices[filter_idx];
    let is_left = matches!(col_idx.side, JoinSide::Left);

    if is_left != expect_left {
        return None;
    }

    let schema = if is_left { left_schema } else { right_schema };
    if col_idx.index >= schema.fields().len() {
        return None;
    }
    let field = schema.field(col_idx.index);
    Some(field.name().clone())
}

/// Recursively find Parquet file paths in the plan tree.
fn find_parquet_paths(plan: &Arc<dyn ExecutionPlan>) -> Vec<PathBuf> {
    use datafusion::datasource::source::DataSourceExec;
    use datafusion::datasource::physical_plan::parquet::source::ParquetSource;

    let mut paths = Vec::new();

    if let Some(ds_exec) =
        plan.as_any().downcast_ref::<DataSourceExec>()
    {
        if let Some((file_config, _parquet_source)) =
            ds_exec.downcast_to_file_source::<ParquetSource>()
        {
            for group in &file_config.file_groups {
                for file in group.iter() {
                    paths.push(PathBuf::from(
                        file.object_meta.location.as_ref(),
                    ));
                }
            }
            return paths;
        }
    }

    for child in plan.children() {
        paths.extend(find_parquet_paths(child));
    }

    paths
}
