"""Extended correctness tests for GIQL MERGE operator vs bedtools merge.

These tests cover transitive chains, topology variations, and scale scenarios
to ensure comprehensive GIQL/bedtools equivalence for merge operations.
"""

import pytest

from giql import transpile

from .utils.bedtools_wrapper import merge
from .utils.comparison import compare_results
from .utils.data_models import GenomicInterval
from .utils.duckdb_loader import load_intervals
from .utils.random_intervals import generate_random_intervals

pytestmark = pytest.mark.integration


def _run_merge_comparison(duckdb_connection, intervals, strand_mode=None):
    """Run GIQL MERGE and bedtools merge, return ComparisonResult."""
    load_intervals(
        duckdb_connection,
        "intervals",
        [i.to_tuple() for i in intervals],
    )

    bedtools_result = merge(
        [i.to_tuple() for i in intervals],
        strand_mode=strand_mode,
    )

    if strand_mode == "same":
        giql_sql = "SELECT MERGE(interval, stranded := true) FROM intervals"
    else:
        giql_sql = "SELECT MERGE(interval) FROM intervals"

    sql = transpile(giql_sql, tables=["intervals"])
    giql_result = duckdb_connection.execute(sql).fetchall()

    return compare_results(giql_result, bedtools_result)


def test_merge_should_combine_transitive_chain_into_single_interval(duckdb_connection):
    """Test MERGE collapses a transitive overlap chain.

    Given:
        A chain A overlaps B, B overlaps C (but A does not overlap C
        directly)
    When:
        GIQL MERGE is compared to bedtools merge
    Then:
        It should merge the entire chain into a single interval
    """
    # Arrange
    intervals = [
        GenomicInterval("chr1", 100, 200, "i1", 0, "+"),
        GenomicInterval("chr1", 180, 300, "i2", 0, "+"),
        GenomicInterval("chr1", 280, 400, "i3", 0, "+"),
    ]

    # Act
    comparison = _run_merge_comparison(duckdb_connection, intervals)

    # Assert
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_merge_should_return_interval_unchanged_when_input_is_single_interval(
    duckdb_connection,
):
    """Test MERGE is a no-op for a single-interval input.

    Given:
        A single interval
    When:
        GIQL MERGE is compared to bedtools merge
    Then:
        It should return the single interval unchanged
    """
    # Arrange
    intervals = [GenomicInterval("chr1", 100, 200, "i1", 0, "+")]

    # Act
    comparison = _run_merge_comparison(duckdb_connection, intervals)

    # Assert
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_merge_should_produce_one_region_when_all_intervals_overlap(
    duckdb_connection,
):
    """Test MERGE collapses fully overlapping intervals into one region.

    Given:
        All intervals on a chromosome overlap forming one big region
    When:
        GIQL MERGE is compared to bedtools merge
    Then:
        It should return a single merged interval
    """
    # Arrange
    intervals = [
        GenomicInterval("chr1", 100, 500, "i1", 0, "+"),
        GenomicInterval("chr1", 200, 400, "i2", 0, "+"),
        GenomicInterval("chr1", 300, 600, "i3", 0, "+"),
        GenomicInterval("chr1", 150, 550, "i4", 0, "+"),
    ]

    # Act
    comparison = _run_merge_comparison(duckdb_connection, intervals)

    # Assert
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_merge_should_return_correct_region_count_when_topology_is_mixed(
    duckdb_connection,
):
    """Test MERGE handles a mix of overlapping clusters and isolated intervals.

    Given:
        A mix of overlapping clusters and isolated intervals
    When:
        GIQL MERGE is compared to bedtools merge
    Then:
        It should produce the correct number of merged regions
    """
    # Arrange
    intervals = [
        # Cluster 1: overlapping
        GenomicInterval("chr1", 100, 200, "c1a", 0, "+"),
        GenomicInterval("chr1", 150, 300, "c1b", 0, "+"),
        # Isolated
        GenomicInterval("chr1", 500, 600, "iso", 0, "+"),
        # Cluster 2: overlapping
        GenomicInterval("chr1", 800, 900, "c2a", 0, "+"),
        GenomicInterval("chr1", 850, 1000, "c2b", 0, "+"),
    ]

    # Act
    comparison = _run_merge_comparison(duckdb_connection, intervals)

    # Assert
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 3


def test_merge_should_combine_intervals_when_overlap_is_one_base(duckdb_connection):
    """Test MERGE triggers on a single-base overlap.

    Given:
        Intervals with exactly 1bp overlap
    When:
        GIQL MERGE is compared to bedtools merge
    Then:
        It should treat the 1bp overlap as sufficient to merge
    """
    # Arrange
    intervals = [
        GenomicInterval("chr1", 100, 200, "i1", 0, "+"),
        GenomicInterval("chr1", 199, 300, "i2", 0, "+"),
    ]

    # Act
    comparison = _run_merge_comparison(duckdb_connection, intervals)

    # Assert
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_merge_should_match_bedtools_when_input_is_unsorted(duckdb_connection):
    """Test MERGE is insensitive to input ordering.

    Given:
        Intervals inserted in non-sorted order
    When:
        GIQL MERGE is compared to bedtools merge
    Then:
        It should produce the same results regardless of input order
    """
    # Arrange
    intervals = [
        GenomicInterval("chr1", 400, 500, "i3", 0, "+"),
        GenomicInterval("chr1", 100, 200, "i1", 0, "+"),
        GenomicInterval("chr1", 150, 250, "i2", 0, "+"),
    ]

    # Act
    comparison = _run_merge_comparison(duckdb_connection, intervals)

    # Assert
    assert comparison.match, comparison.failure_message()


def test_merge_should_operate_per_chromosome_when_input_spans_multiple_chromosomes(
    duckdb_connection,
):
    """Test MERGE groups merges per chromosome.

    Given:
        Overlapping intervals on separate chromosomes
    When:
        GIQL MERGE is compared to bedtools merge
    Then:
        It should merge per-chromosome independently
    """
    # Arrange
    intervals = [
        GenomicInterval("chr1", 100, 200, "c1a", 0, "+"),
        GenomicInterval("chr1", 150, 300, "c1b", 0, "+"),
        GenomicInterval("chr2", 100, 200, "c2a", 0, "+"),
        GenomicInterval("chr2", 150, 300, "c2b", 0, "+"),
        GenomicInterval("chr3", 100, 200, "c3", 0, "+"),  # no overlap
    ]

    # Act
    comparison = _run_merge_comparison(duckdb_connection, intervals)

    # Assert
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 3  # 1 per chrom


def test_merge_should_preserve_strand_when_stranded_true(duckdb_connection):
    """Test MERGE with stranded=true matches bedtools merge -s.

    Given:
        Overlapping intervals on different strands
    When:
        GIQL MERGE(stranded=true) is compared to bedtools merge -s
    Then:
        It should produce the same per-strand merge count as bedtools
    """
    # Arrange
    intervals = [
        GenomicInterval("chr1", 100, 200, "i1", 0, "+"),
        GenomicInterval("chr1", 150, 250, "i2", 0, "+"),
        GenomicInterval("chr1", 120, 220, "i3", 0, "-"),
        GenomicInterval("chr1", 180, 280, "i4", 0, "-"),
    ]
    load_intervals(
        duckdb_connection,
        "intervals",
        [i.to_tuple() for i in intervals],
    )

    bedtools_result = merge(
        [i.to_tuple() for i in intervals],
        strand_mode="same",
    )

    sql = transpile(
        "SELECT MERGE(interval, stranded := true) FROM intervals",
        tables=["intervals"],
    )

    # Act
    giql_result = duckdb_connection.execute(sql).fetchall()

    # Assert
    # Both should have 2 merged intervals (one per strand)
    assert len(giql_result) == len(bedtools_result)


def test_merge_should_match_bedtools_when_dataset_is_large(duckdb_connection):
    """Test MERGE agrees with bedtools on a large synthetic dataset.

    Given:
        100+ intervals across 3 chromosomes
    When:
        GIQL MERGE is compared to bedtools merge
    Then:
        It should produce results matching bedtools on the full dataset
    """
    # Arrange
    intervals = generate_random_intervals(
        seed=42,
        prefix="chr",
        count_per_chrom=100,
        n_chroms=3,
        start_max=500_000,
        max_size=2000,
    )

    # Act
    comparison = _run_merge_comparison(duckdb_connection, intervals)

    # Assert
    assert comparison.match, comparison.failure_message()
