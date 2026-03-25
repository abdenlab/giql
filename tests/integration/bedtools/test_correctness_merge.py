"""Extended correctness tests for GIQL MERGE operator vs bedtools merge.

These tests cover transitive chains, topology variations, and scale scenarios
to ensure comprehensive GIQL/bedtools equivalence for merge operations.
"""

from giql import transpile

from .utils.bedtools_wrapper import merge
from .utils.comparison import compare_results
from .utils.data_models import GenomicInterval
from .utils.duckdb_loader import load_intervals


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


def test_merge_transitive_chain(duckdb_connection):
    """
    GIVEN a chain A overlaps B, B overlaps C (but A doesn't overlap C directly)
    WHEN GIQL MERGE is compared to bedtools merge
    THEN entire chain merged into single interval
    """
    intervals = [
        GenomicInterval("chr1", 100, 200, "i1", 0, "+"),
        GenomicInterval("chr1", 180, 300, "i2", 0, "+"),
        GenomicInterval("chr1", 280, 400, "i3", 0, "+"),
    ]
    comparison = _run_merge_comparison(duckdb_connection, intervals)
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_merge_single_interval(duckdb_connection):
    """
    GIVEN a single interval
    WHEN GIQL MERGE is compared to bedtools merge
    THEN single interval returned unchanged
    """
    intervals = [GenomicInterval("chr1", 100, 200, "i1", 0, "+")]
    comparison = _run_merge_comparison(duckdb_connection, intervals)
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_merge_complete_overlap(duckdb_connection):
    """
    GIVEN all intervals on chromosome overlap (one big region)
    WHEN GIQL MERGE is compared to bedtools merge
    THEN single merged interval
    """
    intervals = [
        GenomicInterval("chr1", 100, 500, "i1", 0, "+"),
        GenomicInterval("chr1", 200, 400, "i2", 0, "+"),
        GenomicInterval("chr1", 300, 600, "i3", 0, "+"),
        GenomicInterval("chr1", 150, 550, "i4", 0, "+"),
    ]
    comparison = _run_merge_comparison(duckdb_connection, intervals)
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_merge_mixed_topology(duckdb_connection):
    """
    GIVEN a mix of overlapping clusters and isolated intervals
    WHEN GIQL MERGE is compared to bedtools merge
    THEN correct number of merged regions
    """
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
    comparison = _run_merge_comparison(duckdb_connection, intervals)
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 3


def test_merge_minimal_overlap(duckdb_connection):
    """
    GIVEN intervals with exactly 1bp overlap
    WHEN GIQL MERGE is compared to bedtools merge
    THEN 1bp overlap triggers merge
    """
    intervals = [
        GenomicInterval("chr1", 100, 200, "i1", 0, "+"),
        GenomicInterval("chr1", 199, 300, "i2", 0, "+"),
    ]
    comparison = _run_merge_comparison(duckdb_connection, intervals)
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 1


def test_merge_unsorted_input(duckdb_connection):
    """
    GIVEN intervals inserted in non-sorted order
    WHEN GIQL MERGE is compared to bedtools merge
    THEN results match regardless of input order
    """
    intervals = [
        GenomicInterval("chr1", 400, 500, "i3", 0, "+"),
        GenomicInterval("chr1", 100, 200, "i1", 0, "+"),
        GenomicInterval("chr1", 150, 250, "i2", 0, "+"),
    ]
    comparison = _run_merge_comparison(duckdb_connection, intervals)
    assert comparison.match, comparison.failure_message()


def test_merge_per_chromosome(duckdb_connection):
    """
    GIVEN overlapping intervals on separate chromosomes
    WHEN GIQL MERGE is compared to bedtools merge
    THEN merging occurs per-chromosome independently
    """
    intervals = [
        GenomicInterval("chr1", 100, 200, "c1a", 0, "+"),
        GenomicInterval("chr1", 150, 300, "c1b", 0, "+"),
        GenomicInterval("chr2", 100, 200, "c2a", 0, "+"),
        GenomicInterval("chr2", 150, 300, "c2b", 0, "+"),
        GenomicInterval("chr3", 100, 200, "c3", 0, "+"),  # no overlap
    ]
    comparison = _run_merge_comparison(duckdb_connection, intervals)
    assert comparison.match, comparison.failure_message()
    assert comparison.giql_row_count == 3  # 1 per chrom


def test_merge_strand_specific_correctness(duckdb_connection):
    """
    GIVEN overlapping intervals on different strands
    WHEN GIQL MERGE(stranded=true) is compared to bedtools merge -s
    THEN per-strand merge count matches
    """
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
    giql_result = duckdb_connection.execute(sql).fetchall()

    # Both should have 2 merged intervals (one per strand)
    assert len(giql_result) == len(bedtools_result)


def test_merge_large_scale(duckdb_connection):
    """
    GIVEN 100+ intervals across 3 chromosomes
    WHEN GIQL MERGE is compared to bedtools merge
    THEN results match on the full dataset
    """
    import random

    rng = random.Random(42)
    intervals = []

    for chrom_num in range(1, 4):
        chrom = f"chr{chrom_num}"
        for i in range(100):
            start = rng.randint(0, 500_000)
            size = rng.randint(100, 2000)
            intervals.append(
                GenomicInterval(chrom, start, start + size, f"{chrom}_{i}", 0, "+")
            )

    comparison = _run_merge_comparison(duckdb_connection, intervals)
    assert comparison.match, comparison.failure_message()
