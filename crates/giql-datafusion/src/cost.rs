use crate::stats::IntervalStats;
use crate::IntersectsOptimizerConfig;

/// Join algorithm selected by the cost model.
#[derive(Debug, Clone, PartialEq)]
pub enum JoinStrategy {
    /// Default nested-loop join (no plan rewrite).
    NestedLoop,
    /// Sweep-line join: sort both sides by start, sweep with an active
    /// set. O((n+m) log(n+m) + k).
    SweepLine {
        /// True if the input is already sorted and the sort step can
        /// be skipped.
        skip_sort: bool,
    },
    /// Binned equi-join: expand intervals into genome bins, hash-join
    /// on bin ID. O(n+m+k) amortized for uniform widths.
    BinnedJoin {
        /// Bin width in base pairs.
        bin_size: usize,
    },
}

/// Cost model for choosing the optimal INTERSECTS join algorithm.
///
/// Uses two fast short-circuit signals (p99/median ratio and CV) to
/// detect distributions where binning degrades, and falls back to a
/// cost comparison for ambiguous cases.
pub struct CostModel {
    p99_median_threshold: f64,
    cv_threshold: f64,
}

/// Relative cost constants for the cost comparison branch.
/// These are unitless scaling factors, not absolute times.
const HASH_COST: f64 = 1.0;
const COMPARE_COST: f64 = 2.0;

/// Minimum and maximum bin sizes to clamp the adaptive bin size.
const MIN_BIN_SIZE: usize = 1_000;
const MAX_BIN_SIZE: usize = 1_000_000;

impl CostModel {
    pub fn new(config: &IntersectsOptimizerConfig) -> Self {
        Self {
            p99_median_threshold: config.p99_median_threshold,
            cv_threshold: config.cv_threshold,
        }
    }

    /// Choose the optimal join strategy based on interval statistics
    /// from both sides of the join.
    pub fn decide(
        &self,
        left: &IntervalStats,
        right: &IntervalStats,
    ) -> JoinStrategy {
        // Short-circuit 1: heavy-tailed distribution.
        // If p99/median > threshold on either side, wide outliers will
        // replicate across many bins, destroying binning performance.
        if left.width.p99_median_ratio > self.p99_median_threshold
            || right.width.p99_median_ratio > self.p99_median_threshold
        {
            let skip_sort = left.is_sorted_by_start
                || right.is_sorted_by_start;
            return JoinStrategy::SweepLine { skip_sort };
        }

        // Short-circuit 2: high width variance.
        // No single bin size works well when CV is high.
        if left.width.cv > self.cv_threshold
            || right.width.cv > self.cv_threshold
        {
            let skip_sort = left.is_sorted_by_start
                || right.is_sorted_by_start;
            return JoinStrategy::SweepLine { skip_sort };
        }

        // Cost comparison: estimate binned vs sweep costs.
        let bin_size = self.estimate_optimal_bin_size(left, right);
        let binned_cost = self.estimate_binned_cost(left, right, bin_size);
        let sweep_cost = self.estimate_sweep_cost(left, right);

        if binned_cost < sweep_cost {
            JoinStrategy::BinnedJoin { bin_size }
        } else {
            let skip_sort = left.is_sorted_by_start
                || right.is_sorted_by_start;
            JoinStrategy::SweepLine { skip_sort }
        }
    }

    /// Estimate optimal bin size from the width distributions.
    ///
    /// Uses p95 as the bin width so that ~95% of intervals fit in a
    /// single bin (replication factor ≈ 1).
    fn estimate_optimal_bin_size(
        &self,
        left: &IntervalStats,
        right: &IntervalStats,
    ) -> usize {
        // Use the larger p95 so both sides have low replication.
        let p95 = left.width.p95.max(right.width.p95);
        let raw = p95.ceil() as usize;
        raw.clamp(MIN_BIN_SIZE, MAX_BIN_SIZE)
    }

    /// Estimate cost of binned equi-join.
    ///
    /// Each interval is replicated into `mean_width / bin_size + 1`
    /// bin entries, then hash-joined. Cost scales with total expanded
    /// row count.
    fn estimate_binned_cost(
        &self,
        left: &IntervalStats,
        right: &IntervalStats,
        bin_size: usize,
    ) -> f64 {
        let left_replication =
            left.width.mean / bin_size as f64 + 1.0;
        let right_replication =
            right.width.mean / bin_size as f64 + 1.0;

        let expanded_left =
            left.row_count as f64 * left_replication;
        let expanded_right =
            right.row_count as f64 * right_replication;

        (expanded_left + expanded_right) * HASH_COST
    }

    /// Estimate cost of sweep-line join.
    ///
    /// Dominated by sorting both sides: O((n+m) log(n+m)). If either
    /// side is already sorted, the cost drops by ~50%.
    fn estimate_sweep_cost(
        &self,
        left: &IntervalStats,
        right: &IntervalStats,
    ) -> f64 {
        let n = left.row_count as f64;
        let m = right.row_count as f64;
        let total = n + m;

        let mut cost = total * total.log2() * COMPARE_COST;

        // If either side is sorted, we skip one of the two sorts.
        if left.is_sorted_by_start || right.is_sorted_by_start {
            cost *= 0.5;
        }

        cost
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::{RowGroupBounds, WidthStats};

    fn make_stats(
        row_count: usize,
        median: f64,
        mean: f64,
        p95: f64,
        p99: f64,
        cv: f64,
        sorted: bool,
    ) -> IntervalStats {
        IntervalStats {
            row_count,
            domain_min: 0,
            domain_max: 1_000_000,
            is_sorted_by_start: sorted,
            row_group_bounds: vec![RowGroupBounds {
                min_start: 0,
                max_start: 900_000,
                min_end: 100,
                max_end: 1_000_000,
                row_count,
            }],
            width: WidthStats {
                median,
                mean,
                p95,
                p99,
                cv,
                p99_median_ratio: if median > 0.0 {
                    p99 / median
                } else {
                    0.0
                },
            },
        }
    }

    fn default_config() -> IntersectsOptimizerConfig {
        IntersectsOptimizerConfig::default()
    }

    #[test]
    fn test_heavy_tailed_selects_sweep_line() {
        let model = CostModel::new(&default_config());
        // p99/median = 100/5 = 20 > 10
        let left = make_stats(100_000, 5.0, 10.0, 50.0, 100.0, 0.8, false);
        let right = make_stats(100_000, 100.0, 100.0, 100.0, 100.0, 0.0, false);

        match model.decide(&left, &right) {
            JoinStrategy::SweepLine { .. } => {}
            other => panic!("Expected SweepLine, got {:?}", other),
        }
    }

    #[test]
    fn test_high_cv_selects_sweep_line() {
        let model = CostModel::new(&default_config());
        // CV = 2.0 > 1.5
        let left = make_stats(100_000, 100.0, 100.0, 100.0, 100.0, 2.0, false);
        let right = make_stats(100_000, 100.0, 100.0, 100.0, 100.0, 0.5, false);

        match model.decide(&left, &right) {
            JoinStrategy::SweepLine { .. } => {}
            other => panic!("Expected SweepLine, got {:?}", other),
        }
    }

    #[test]
    fn test_uniform_selects_binned() {
        let model = CostModel::new(&default_config());
        // Uniform 100bp intervals, low CV, p99/median ≈ 1
        let left = make_stats(1_000_000, 100.0, 100.0, 100.0, 100.0, 0.0, false);
        let right = make_stats(1_000_000, 100.0, 100.0, 100.0, 100.0, 0.0, false);

        match model.decide(&left, &right) {
            JoinStrategy::BinnedJoin { bin_size } => {
                assert!(bin_size >= MIN_BIN_SIZE);
            }
            other => panic!("Expected BinnedJoin, got {:?}", other),
        }
    }

    #[test]
    fn test_sorted_input_sets_skip_sort() {
        let model = CostModel::new(&default_config());
        // High CV triggers sweep line; sorted input should set skip_sort
        let left = make_stats(1_000_000, 100.0, 500.0, 1000.0, 5000.0, 2.0, true);
        let right = make_stats(1_000_000, 100.0, 500.0, 1000.0, 5000.0, 0.5, false);

        match model.decide(&left, &right) {
            JoinStrategy::SweepLine { skip_sort } => {
                assert!(skip_sort);
            }
            other => panic!("Expected SweepLine with skip_sort, got {:?}", other),
        }
    }

    #[test]
    fn test_bin_size_clamped() {
        let model = CostModel::new(&default_config());
        // Very small p95 — bin size should clamp to MIN_BIN_SIZE
        let left = make_stats(100_000, 10.0, 10.0, 10.0, 10.0, 0.1, false);
        let right = make_stats(100_000, 10.0, 10.0, 10.0, 10.0, 0.1, false);

        let bin_size = model.estimate_optimal_bin_size(&left, &right);
        assert_eq!(bin_size, MIN_BIN_SIZE);
    }
}
