use std::fs::File;
use std::path::Path;

use arrow::array::{Array, Int64Array};
use datafusion::common::Result;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ProjectionMask;

use super::types::WidthStats;

/// Read start and end columns from selected row groups and compute
/// width distribution statistics.
///
/// Only two columns are read from each row group — no other data is
/// touched. For a multi-GB dataset this typically completes in
/// milliseconds.
pub fn sample_widths(
    path: &Path,
    start_col: &str,
    end_col: &str,
    row_group_indices: &[usize],
) -> Result<WidthStats> {
    if row_group_indices.is_empty() {
        return Err(datafusion::error::DataFusionError::Plan(
            "No row groups to sample".to_string(),
        ));
    }

    let file = File::open(path).map_err(|e| {
        datafusion::error::DataFusionError::External(Box::new(e))
    })?;

    let builder =
        ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| {
            datafusion::error::DataFusionError::External(Box::new(e))
        })?;

    // Find column indices in the Arrow schema
    let arrow_schema = builder.schema().clone();
    let start_idx = arrow_schema
        .fields()
        .iter()
        .position(|f| f.name() == start_col)
        .ok_or_else(|| {
            datafusion::error::DataFusionError::Plan(format!(
                "Column '{start_col}' not found in Parquet schema"
            ))
        })?;
    let end_idx = arrow_schema
        .fields()
        .iter()
        .position(|f| f.name() == end_col)
        .ok_or_else(|| {
            datafusion::error::DataFusionError::Plan(format!(
                "Column '{end_col}' not found in Parquet schema"
            ))
        })?;

    // Project only the start and end columns using the Parquet schema
    let parquet_schema = builder.parquet_schema();
    let projection =
        ProjectionMask::leaves(parquet_schema, vec![start_idx, end_idx]);

    let reader = builder
        .with_projection(projection)
        .with_row_groups(row_group_indices.to_vec())
        .build()
        .map_err(|e| {
            datafusion::error::DataFusionError::External(Box::new(e))
        })?;

    let mut widths: Vec<f64> = Vec::new();

    for batch_result in reader {
        let batch = batch_result.map_err(|e| {
            datafusion::error::DataFusionError::External(Box::new(e))
        })?;

        // Columns are projected, so index 0 = start, 1 = end
        let starts = extract_i64_column(&batch, 0, start_col)?;
        let ends = extract_i64_column(&batch, 1, end_col)?;

        for i in 0..batch.num_rows() {
            if !starts.is_null(i) && !ends.is_null(i) {
                let w = (ends.value(i) - starts.value(i)) as f64;
                widths.push(w);
            }
        }
    }

    if widths.is_empty() {
        return Err(datafusion::error::DataFusionError::Plan(
            "No valid intervals found in sampled row groups".to_string(),
        ));
    }

    Ok(compute_width_stats(&mut widths))
}

/// Compute width distribution statistics from a vector of widths.
pub(crate) fn compute_width_stats(widths: &mut [f64]) -> WidthStats {
    widths.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let n = widths.len();
    let median = percentile_sorted(widths, 0.5);
    let p95 = percentile_sorted(widths, 0.95);
    let p99 = percentile_sorted(widths, 0.99);

    let sum: f64 = widths.iter().sum();
    let mean = sum / n as f64;

    let variance: f64 =
        widths.iter().map(|w| (w - mean).powi(2)).sum::<f64>() / n as f64;
    let std_dev = variance.sqrt();
    let cv = if mean > 0.0 { std_dev / mean } else { 0.0 };

    let p99_median_ratio = if median > 0.0 { p99 / median } else { 0.0 };

    WidthStats {
        median,
        mean,
        p95,
        p99,
        cv,
        p99_median_ratio,
    }
}

/// Compute a percentile from a sorted slice using linear interpolation.
fn percentile_sorted(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    if sorted.len() == 1 {
        return sorted[0];
    }

    let rank = p * (sorted.len() - 1) as f64;
    let lower = rank.floor() as usize;
    let upper = rank.ceil() as usize;
    let frac = rank - lower as f64;

    if lower == upper {
        sorted[lower]
    } else {
        sorted[lower] * (1.0 - frac) + sorted[upper] * frac
    }
}

/// Extract an i64 column from a record batch, handling both Int32 and
/// Int64 physical types.
fn extract_i64_column(
    batch: &arrow::record_batch::RecordBatch,
    col_idx: usize,
    col_name: &str,
) -> Result<Int64Array> {
    let col = batch.column(col_idx);

    if let Some(arr) = col.as_any().downcast_ref::<Int64Array>() {
        return Ok(arr.clone());
    }

    if let Some(arr) =
        col.as_any().downcast_ref::<arrow::array::Int32Array>()
    {
        let converted: Int64Array = arr
            .iter()
            .map(|v| v.map(|x| x as i64))
            .collect();
        return Ok(converted);
    }

    Err(datafusion::error::DataFusionError::Plan(format!(
        "Column '{col_name}' is not Int32 or Int64"
    )))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_width_stats_uniform() {
        let mut widths = vec![100.0; 1000];
        let stats = compute_width_stats(&mut widths);

        assert!((stats.median - 100.0).abs() < 1e-6);
        assert!((stats.mean - 100.0).abs() < 1e-6);
        assert!((stats.cv).abs() < 1e-6);
        assert!((stats.p99_median_ratio - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_compute_width_stats_heavy_tailed() {
        // 950 intervals of width 100, 50 intervals of width 10000.
        // With 5% outliers, p99 lands squarely in the 10000 range.
        let mut widths: Vec<f64> = vec![100.0; 950];
        widths.extend(vec![10000.0; 50]);
        let stats = compute_width_stats(&mut widths);

        assert!(stats.p99_median_ratio > 10.0);
        assert!(stats.cv > 1.0);
    }

    #[test]
    fn test_percentile_sorted() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!((percentile_sorted(&data, 0.0) - 1.0).abs() < 1e-6);
        assert!((percentile_sorted(&data, 0.5) - 3.0).abs() < 1e-6);
        assert!((percentile_sorted(&data, 1.0) - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_percentile_sorted_single() {
        assert!((percentile_sorted(&[42.0], 0.5) - 42.0).abs() < 1e-6);
    }

    #[test]
    fn test_sample_widths_uniform_parquet() {
        use arrow::array::{Int64Array, StringArray};
        use arrow::datatypes::{DataType, Field, Schema};
        use arrow::record_batch::RecordBatch;
        use parquet::arrow::ArrowWriter;
        use std::sync::Arc;
        use tempfile::NamedTempFile;

        // Write a Parquet file with uniform 100bp intervals
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let starts: Vec<i64> = (0..100).map(|i| i * 200).collect();
        let ends: Vec<i64> = starts.iter().map(|s| s + 100).collect();
        let chroms: Vec<&str> = vec!["chr1"; 100];
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(chroms)),
                Arc::new(Int64Array::from(starts)),
                Arc::new(Int64Array::from(ends)),
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
            sample_widths(file.path(), "start", "end", &[0]).unwrap();
        assert!((stats.median - 100.0).abs() < 1e-6);
        assert!(stats.cv < 0.01);
    }

    #[test]
    fn test_sample_widths_missing_column() {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["chr1"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![200])),
            ],
        )
        .unwrap();
        let file = NamedTempFile::new().unwrap();
        let mut writer =
            ArrowWriter::try_new(file.reopen().unwrap(), schema, None)
                .unwrap();
        writer.write(&batch).unwrap();
        writer.close().unwrap();

        let result =
            sample_widths(file.path(), "nonexistent", "end", &[0]);
        assert!(result.is_err());
    }

    #[test]
    fn test_sample_widths_empty_row_groups() {
        let result = sample_widths(
            std::path::Path::new("/tmp/fake.parquet"),
            "start",
            "end",
            &[],
        );
        assert!(result.is_err());
    }
}
