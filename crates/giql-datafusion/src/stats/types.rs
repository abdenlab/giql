/// Per-row-group statistics derived from Parquet column chunk metadata.
#[derive(Debug, Clone)]
pub struct RowGroupBounds {
    /// Minimum value of the start column in this row group.
    pub min_start: i64,
    /// Maximum value of the start column in this row group.
    pub max_start: i64,
    /// Minimum value of the end column in this row group.
    pub min_end: i64,
    /// Maximum value of the end column in this row group.
    pub max_end: i64,
    /// Number of rows in this row group.
    pub row_count: usize,
}

impl RowGroupBounds {
    /// Loose upper bound on interval width within this row group.
    ///
    /// No interval can be wider than `max(end) - min(start)`.
    pub fn width_upper_bound(&self) -> i64 {
        self.max_end - self.min_start
    }

    /// Width signal from the high end of the sort order.
    ///
    /// If `max(end) - max(start)` is small, the widest intervals at
    /// the end of the row group are narrow.
    pub fn width_at_max(&self) -> i64 {
        self.max_end - self.max_start
    }

    /// Width signal from the low end of the sort order.
    ///
    /// If `min(end) - min(start)` is small, the earliest intervals
    /// are narrow.
    pub fn width_at_min(&self) -> i64 {
        self.min_end - self.min_start
    }
}

/// Width distribution statistics computed from sampled intervals.
#[derive(Debug, Clone)]
pub struct WidthStats {
    /// Median interval width.
    pub median: f64,
    /// Mean interval width.
    pub mean: f64,
    /// 95th percentile width.
    pub p95: f64,
    /// 99th percentile width.
    pub p99: f64,
    /// Coefficient of variation (std_dev / mean).
    pub cv: f64,
    /// Ratio of p99 to median width.
    pub p99_median_ratio: f64,
}

/// Aggregate interval statistics for a Parquet file, combining
/// metadata-derived bounds with sampled width distribution.
#[derive(Debug, Clone)]
pub struct IntervalStats {
    /// Total row count across all row groups.
    pub row_count: usize,
    /// Global minimum start position.
    pub domain_min: i64,
    /// Global maximum end position.
    pub domain_max: i64,
    /// Whether the file is sorted by the start column
    /// (from `sorting_columns` metadata).
    pub is_sorted_by_start: bool,
    /// Per-row-group bounds from column chunk statistics.
    pub row_group_bounds: Vec<RowGroupBounds>,
    /// Width distribution from sampling.
    pub width: WidthStats,
}

impl IntervalStats {
    /// Domain span: total coordinate range covered by the file.
    pub fn domain_span(&self) -> i64 {
        self.domain_max - self.domain_min
    }
}
