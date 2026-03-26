use std::fs::File;
use std::path::Path;

use datafusion::common::Result;
use parquet::file::reader::FileReader;
use parquet::file::serialized_reader::SerializedFileReader;
use parquet::file::statistics::Statistics;

use super::types::RowGroupBounds;

/// Statistics extracted from Parquet file footer metadata only (no data
/// pages read). This is effectively free — it reads the file footer
/// which is already cached by the OS after open.
#[derive(Debug, Clone)]
pub struct MetadataStats {
    /// Per-row-group bounds for start/end columns.
    pub row_group_bounds: Vec<RowGroupBounds>,
    /// Total row count across all row groups.
    pub total_rows: usize,
    /// Global minimum start position.
    pub domain_min: i64,
    /// Global maximum end position.
    pub domain_max: i64,
    /// Whether the file declares itself sorted by the start column.
    pub is_sorted_by_start: bool,
    /// Whether page-level column index is present.
    pub has_page_index: bool,
}

/// Collect metadata-tier statistics from a Parquet file footer.
///
/// Reads only the file footer — no data pages are touched. Returns
/// per-row-group column statistics for the start and end columns, plus
/// file-level properties like sort order and page index presence.
pub fn collect_metadata(
    path: &Path,
    start_col: &str,
    end_col: &str,
) -> Result<MetadataStats> {
    let file = File::open(path).map_err(|e| {
        datafusion::error::DataFusionError::External(Box::new(e))
    })?;
    let reader = SerializedFileReader::new(file).map_err(|e| {
        datafusion::error::DataFusionError::External(Box::new(e))
    })?;

    let parquet_metadata = reader.metadata();
    let file_metadata = parquet_metadata.file_metadata();
    let schema = file_metadata.schema_descr();

    // Find column indices for start and end
    let start_idx = find_column_index(schema, start_col)?;
    let end_idx = find_column_index(schema, end_col)?;

    // Check sort order from file metadata
    let is_sorted_by_start = check_sort_order(file_metadata, start_col, schema);

    let num_row_groups = parquet_metadata.num_row_groups();
    let mut row_group_bounds = Vec::with_capacity(num_row_groups);
    let mut total_rows: usize = 0;
    let mut global_min_start = i64::MAX;
    let mut global_max_end = i64::MIN;
    let mut has_page_index = false;

    for rg_idx in 0..num_row_groups {
        let rg_metadata = parquet_metadata.row_group(rg_idx);
        let row_count = rg_metadata.num_rows() as usize;
        total_rows += row_count;

        let start_col_meta = rg_metadata.column(start_idx);
        let end_col_meta = rg_metadata.column(end_idx);

        // Check page index on first row group
        if rg_idx == 0 {
            has_page_index = start_col_meta.offset_index_offset().is_some()
                || start_col_meta.column_index_offset().is_some();
        }

        let (min_start, max_start) = extract_i64_min_max(
            start_col_meta.statistics(),
            start_col,
        )?;
        let (min_end, max_end) = extract_i64_min_max(
            end_col_meta.statistics(),
            end_col,
        )?;

        if min_start < global_min_start {
            global_min_start = min_start;
        }
        if max_end > global_max_end {
            global_max_end = max_end;
        }

        row_group_bounds.push(RowGroupBounds {
            min_start,
            max_start,
            min_end,
            max_end,
            row_count,
        });
    }

    Ok(MetadataStats {
        row_group_bounds,
        total_rows,
        domain_min: global_min_start,
        domain_max: global_max_end,
        is_sorted_by_start,
        has_page_index,
    })
}

/// Find the column index for a named column in the Parquet schema.
fn find_column_index(
    schema: &parquet::schema::types::SchemaDescriptor,
    col_name: &str,
) -> Result<usize> {
    for i in 0..schema.num_columns() {
        if schema.column(i).name() == col_name {
            return Ok(i);
        }
    }
    Err(datafusion::error::DataFusionError::Plan(format!(
        "Column '{col_name}' not found in Parquet schema"
    )))
}

/// Check whether the file declares itself sorted by the start column.
fn check_sort_order(
    file_metadata: &parquet::file::metadata::FileMetaData,
    start_col: &str,
    schema: &parquet::schema::types::SchemaDescriptor,
) -> bool {
    // Check key-value metadata for sorting_columns
    if let Some(kv_metadata) = file_metadata.key_value_metadata() {
        for kv in kv_metadata {
            if kv.key == "sorting_columns" || kv.key == "pandas.sort_columns" {
                if let Some(ref value) = kv.value {
                    if value.contains(start_col) {
                        return true;
                    }
                }
            }
        }
    }

    // Check if column order metadata indicates ascending on start column
    if let Some(sort_order) = file_metadata.column_orders() {
        if let Ok(start_idx) = find_column_index(schema, start_col) {
            if start_idx < sort_order.len() {
                // TypeDefinedOrder means the natural ordering applies,
                // which combined with sorted row groups suggests sorted data.
                // However, this only indicates comparison semantics, not
                // that data is actually sorted. We need sorting_columns
                // metadata for a definitive answer. Return false here.
            }
        }
    }

    false
}

/// Extract min and max i64 values from Parquet column statistics.
fn extract_i64_min_max(
    statistics: Option<&Statistics>,
    col_name: &str,
) -> Result<(i64, i64)> {
    match statistics {
        Some(Statistics::Int32(stats)) => {
            let min_val = stats.min_opt().ok_or_else(|| {
                datafusion::error::DataFusionError::Plan(format!(
                    "Column '{col_name}' Int32 stats missing min"
                ))
            })?;
            let max_val = stats.max_opt().ok_or_else(|| {
                datafusion::error::DataFusionError::Plan(format!(
                    "Column '{col_name}' Int32 stats missing max"
                ))
            })?;
            Ok((*min_val as i64, *max_val as i64))
        }
        Some(Statistics::Int64(stats)) => {
            let min_val = stats.min_opt().ok_or_else(|| {
                datafusion::error::DataFusionError::Plan(format!(
                    "Column '{col_name}' Int64 stats missing min"
                ))
            })?;
            let max_val = stats.max_opt().ok_or_else(|| {
                datafusion::error::DataFusionError::Plan(format!(
                    "Column '{col_name}' Int64 stats missing max"
                ))
            })?;
            Ok((*min_val, *max_val))
        }
        Some(_) => Err(datafusion::error::DataFusionError::Plan(format!(
            "Column '{col_name}' has unsupported statistics type for interval bounds"
        ))),
        None => Err(datafusion::error::DataFusionError::Plan(format!(
            "Column '{col_name}' has no statistics in Parquet metadata"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::array::Int64Array;
    use arrow::datatypes::{DataType, Field, Schema};
    use arrow::record_batch::RecordBatch;
    use parquet::arrow::ArrowWriter;
    use std::sync::Arc;
    use tempfile::NamedTempFile;

    fn write_test_parquet(
        starts: &[i64],
        ends: &[i64],
    ) -> NamedTempFile {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let chroms: Vec<&str> = vec!["chr1"; starts.len()];
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(arrow::array::StringArray::from(chroms)),
                Arc::new(Int64Array::from(starts.to_vec())),
                Arc::new(Int64Array::from(ends.to_vec())),
            ],
        )
        .unwrap();

        let file = NamedTempFile::new().unwrap();
        let mut writer =
            ArrowWriter::try_new(file.reopen().unwrap(), schema, None)
                .unwrap();
        writer.write(&batch).unwrap();
        writer.close().unwrap();
        file
    }

    #[test]
    fn test_collect_metadata_basic() {
        let starts = vec![100, 200, 300, 400, 500];
        let ends = vec![150, 280, 350, 450, 600];
        let file = write_test_parquet(&starts, &ends);

        let stats =
            collect_metadata(file.path(), "start", "end").unwrap();

        assert_eq!(stats.total_rows, 5);
        assert_eq!(stats.domain_min, 100);
        assert_eq!(stats.domain_max, 600);
        assert_eq!(stats.row_group_bounds.len(), 1);

        let rg = &stats.row_group_bounds[0];
        assert_eq!(rg.min_start, 100);
        assert_eq!(rg.max_start, 500);
        assert_eq!(rg.min_end, 150);
        assert_eq!(rg.max_end, 600);
        assert_eq!(rg.row_count, 5);
    }

    #[test]
    fn test_collect_metadata_missing_column() {
        let file = write_test_parquet(&[100], &[200]);
        let result =
            collect_metadata(file.path(), "nonexistent", "end");
        assert!(result.is_err());
    }
}
