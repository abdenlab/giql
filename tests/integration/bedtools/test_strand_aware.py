"""Integration tests for GIQL strand-aware operations.

These tests validate that GIQL correctly handles strand-specific interval
operations, matching bedtools behavior with -s and -S flags.
"""

from giql import transpile

from .utils.bed_export import load_intervals
from .utils.bedtools_wrapper import closest
from .utils.bedtools_wrapper import intersect
from .utils.bedtools_wrapper import merge
from .utils.comparison import compare_results
from .utils.data_models import GenomicInterval


def test_intersect_same_strand(duckdb_connection):
    """
    Given:
        Intervals on both same and opposite strands
    When:
        INTERSECTS with same-strand filter is applied
    Then:
        Only same-strand overlaps are reported
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
        GenomicInterval("chr1", 300, 400, "a2", 150, "-"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 150, 250, "b1", 100, "+"),
        GenomicInterval("chr1", 350, 450, "b2", 150, "-"),
        GenomicInterval("chr1", 150, 250, "b3", 200, "-"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals_a",
        [i.to_tuple() for i in intervals_a],
    )
    load_intervals(
        duckdb_connection,
        "intervals_b",
        [i.to_tuple() for i in intervals_b],
    )

    bedtools_result = intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
        strand_mode="same",
    )

    sql = transpile(
        """
        SELECT DISTINCT a.*
        FROM intervals_a a, intervals_b b
        WHERE a.interval INTERSECTS b.interval
          AND a.strand = b.strand
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, comparison.failure_message()


def test_intersect_opposite_strand(duckdb_connection):
    """
    Given:
        Intervals on both same and opposite strands
    When:
        INTERSECTS with opposite-strand filter is applied
    Then:
        Only opposite-strand overlaps are reported
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
        GenomicInterval("chr1", 300, 400, "a2", 150, "-"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 150, 250, "b1", 100, "-"),
        GenomicInterval("chr1", 350, 450, "b2", 150, "+"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals_a",
        [i.to_tuple() for i in intervals_a],
    )
    load_intervals(
        duckdb_connection,
        "intervals_b",
        [i.to_tuple() for i in intervals_b],
    )

    bedtools_result = intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
        strand_mode="opposite",
    )

    sql = transpile(
        """
        SELECT DISTINCT a.*
        FROM intervals_a a, intervals_b b
        WHERE a.interval INTERSECTS b.interval
          AND a.strand != b.strand
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, comparison.failure_message()


def test_intersect_ignore_strand(duckdb_connection):
    """
    Given:
        Intervals with various strand combinations
    When:
        INTERSECTS without strand requirements is applied
    Then:
        All overlaps are reported regardless of strand
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 150, 250, "b1", 100, "+"),
        GenomicInterval("chr1", 150, 250, "b2", 150, "-"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals_a",
        [i.to_tuple() for i in intervals_a],
    )
    load_intervals(
        duckdb_connection,
        "intervals_b",
        [i.to_tuple() for i in intervals_b],
    )

    bedtools_result = intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
    )

    sql = transpile(
        """
        SELECT DISTINCT a.*
        FROM intervals_a a, intervals_b b
        WHERE a.interval INTERSECTS b.interval
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, comparison.failure_message()


def test_intersect_mixed_strands(duckdb_connection):
    """
    Given:
        Complex scenario with +, -, and unstranded intervals
    When:
        INTERSECTS with same-strand requirement is applied
    Then:
        Results correctly handle strand matching logic
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
        GenomicInterval("chr1", 300, 400, "a2", 150, "-"),
        GenomicInterval("chr1", 500, 600, "a3", 200, "."),
    ]
    intervals_b = [
        GenomicInterval("chr1", 150, 250, "b1", 100, "+"),
        GenomicInterval("chr1", 350, 450, "b2", 150, "-"),
        GenomicInterval("chr1", 550, 650, "b3", 200, "."),
    ]

    load_intervals(
        duckdb_connection,
        "intervals_a",
        [i.to_tuple() for i in intervals_a],
    )
    load_intervals(
        duckdb_connection,
        "intervals_b",
        [i.to_tuple() for i in intervals_b],
    )

    bedtools_result = intersect(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
        strand_mode="same",
    )

    sql = transpile(
        """
        SELECT DISTINCT a.*
        FROM intervals_a a, intervals_b b
        WHERE a.interval INTERSECTS b.interval
            AND a.strand = b.strand
            AND a.strand != '.'
            AND b.strand != '.'
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, comparison.failure_message()


def test_nearest_same_strand(duckdb_connection):
    """
    Given:
        Intervals with candidates on same and opposite strands
    When:
        NEAREST with stranded := true is applied
    Then:
        Only same-strand nearest intervals are reported
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 250, 300, "b1", 100, "+"),
        GenomicInterval("chr1", 220, 240, "b2", 150, "-"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals_a",
        [i.to_tuple() for i in intervals_a],
    )
    load_intervals(
        duckdb_connection,
        "intervals_b",
        [i.to_tuple() for i in intervals_b],
    )

    bedtools_result = closest(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
        strand_mode="same",
    )

    sql = transpile(
        """
        SELECT a.*, b.*
        FROM intervals_a a
        CROSS JOIN LATERAL NEAREST(
            intervals_b,
            reference := a.interval,
            k := 1,
            stranded := true
        ) b
        ORDER BY a.chrom, a.start
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    assert len(giql_result) == len(bedtools_result)
    # Should find b1 (same strand +), not b2 (opposite strand -)
    assert giql_result[0][9] == "b1"


def test_nearest_opposite_strand(duckdb_connection):
    """
    Given:
        Intervals with candidates on same and opposite strands
    When:
        bedtools closest with opposite-strand requirement is applied
    Then:
        Only opposite-strand nearest intervals are reported
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 250, 300, "b1", 100, "-"),
        GenomicInterval("chr1", 220, 240, "b2", 150, "+"),
    ]

    bedtools_result = closest(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
        strand_mode="opposite",
    )

    # Verify bedtools returns the opposite-strand interval
    assert len(bedtools_result) == 1
    assert bedtools_result[0][3] == "a1"
    assert bedtools_result[0][9] == "b1"


def test_nearest_ignore_strand(duckdb_connection):
    """
    Given:
        Intervals on different strands
    When:
        NEAREST without strand requirements is applied
    Then:
        Closest interval is found regardless of strand
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 250, 300, "b1", 100, "+"),
        GenomicInterval("chr1", 220, 240, "b2", 150, "-"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals_a",
        [i.to_tuple() for i in intervals_a],
    )
    load_intervals(
        duckdb_connection,
        "intervals_b",
        [i.to_tuple() for i in intervals_b],
    )

    bedtools_result = closest(
        [i.to_tuple() for i in intervals_a],
        [i.to_tuple() for i in intervals_b],
    )

    sql = transpile(
        """
        SELECT a.*, b.*
        FROM intervals_a a
        CROSS JOIN LATERAL NEAREST(
            intervals_b,
            reference := a.interval,
            k := 1
        ) b
        ORDER BY a.chrom, a.start
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    assert len(giql_result) == len(bedtools_result)
    # b2 is closer (gap of 20bp) regardless of strand
    assert giql_result[0][9] == "b2"


def test_merge_strand_specific(duckdb_connection):
    """
    Given:
        Overlapping intervals on different strands
    When:
        bedtools merge with strand-specific flag is applied
    Then:
        Intervals are merged per-strand
    """
    intervals = [
        GenomicInterval("chr1", 100, 200, "i1", 100, "+"),
        GenomicInterval("chr1", 150, 250, "i2", 150, "+"),
        GenomicInterval("chr1", 120, 180, "i3", 200, "-"),
        GenomicInterval("chr1", 160, 240, "i4", 100, "-"),
    ]

    bedtools_result = merge(
        [i.to_tuple() for i in intervals],
        strand_mode="same",
    )

    # Should produce at least 2 merged intervals (one per strand)
    assert len(bedtools_result) >= 2
