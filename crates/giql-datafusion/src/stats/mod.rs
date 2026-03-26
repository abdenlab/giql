pub mod metadata;
pub mod sampler;
pub mod types;

pub use metadata::{collect_metadata, MetadataStats};
pub use sampler::sample_widths;
pub use types::{IntervalStats, RowGroupBounds, WidthStats};

use datafusion::common::Result;
use std::path::Path;

/// Collect full interval statistics from a Parquet file by combining
/// metadata (tier 1, free) with lightweight sampling (tier 2,
/// milliseconds).
///
/// Returns `None` if the file cannot be read or lacks the required
/// columns.
pub fn collect_parquet_stats(
    path: &Path,
    start_col: &str,
    end_col: &str,
    max_sample_row_groups: usize,
) -> Result<IntervalStats> {
    // Tier 1: file footer metadata
    let meta = collect_metadata(path, start_col, end_col)?;

    // Select representative row groups for sampling
    let sample_indices =
        select_representative_row_groups(meta.row_group_bounds.len(), max_sample_row_groups);

    // Tier 2: lightweight sampling
    let width = sample_widths(path, start_col, end_col, &sample_indices)?;

    Ok(IntervalStats {
        row_count: meta.total_rows,
        domain_min: meta.domain_min,
        domain_max: meta.domain_max,
        is_sorted_by_start: meta.is_sorted_by_start,
        row_group_bounds: meta.row_group_bounds,
        width,
    })
}

/// Pick representative row groups for sampling: first, middle, last.
fn select_representative_row_groups(
    num_row_groups: usize,
    max_sample: usize,
) -> Vec<usize> {
    if num_row_groups == 0 {
        return vec![];
    }
    if num_row_groups == 1 || max_sample == 1 {
        return vec![0];
    }

    let last = num_row_groups - 1;
    if num_row_groups == 2 || max_sample == 2 {
        return vec![0, last];
    }

    // First, middle, last
    let mid = num_row_groups / 2;
    let mut indices = vec![0, mid, last];
    indices.truncate(max_sample);
    indices
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_select_representative_single() {
        assert_eq!(select_representative_row_groups(1, 3), vec![0]);
    }

    #[test]
    fn test_select_representative_two() {
        assert_eq!(select_representative_row_groups(2, 3), vec![0, 1]);
    }

    #[test]
    fn test_select_representative_many() {
        assert_eq!(
            select_representative_row_groups(10, 3),
            vec![0, 5, 9]
        );
    }

    #[test]
    fn test_select_representative_max_one() {
        assert_eq!(select_representative_row_groups(10, 1), vec![0]);
    }

    #[test]
    fn test_collect_parquet_stats_uniform() {
        use arrow::array::{Int64Array, StringArray};
        use arrow::datatypes::{DataType, Field, Schema};
        use arrow::record_batch::RecordBatch;
        use parquet::arrow::ArrowWriter;
        use std::sync::Arc;
        use tempfile::NamedTempFile;

        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let starts: Vec<i64> = (0..50).map(|i| i * 200).collect();
        let ends: Vec<i64> = starts.iter().map(|s| s + 100).collect();
        let chroms: Vec<&str> = vec!["chr1"; 50];
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(chroms)),
                Arc::new(Int64Array::from(starts.clone())),
                Arc::new(Int64Array::from(ends.clone())),
            ],
        )
        .unwrap();

        let file = NamedTempFile::new().unwrap();
        let mut writer =
            ArrowWriter::try_new(file.reopen().unwrap(), schema, None)
                .unwrap();
        writer.write(&batch).unwrap();
        writer.close().unwrap();

        let stats =
            collect_parquet_stats(file.path(), "start", "end", 3)
                .unwrap();

        assert_eq!(stats.row_count, 50);
        assert_eq!(stats.domain_min, 0);
        assert_eq!(stats.domain_max, *ends.last().unwrap());
        assert!((stats.width.median - 100.0).abs() < 1e-6);
        assert!(stats.width.cv < 0.01);
    }

    #[test]
    fn test_collect_parquet_stats_nonexistent_file() {
        let result = collect_parquet_stats(
            Path::new("/tmp/nonexistent_file.parquet"),
            "start",
            "end",
            3,
        );
        assert!(result.is_err());
    }
}
