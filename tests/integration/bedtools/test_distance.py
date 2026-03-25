"""Integration tests for GIQL DISTANCE operator.

These tests validate that GIQL's DISTANCE function computes the correct
genomic distance between intervals, cross-validated against bedtools
closest -d output.
"""

from .utils.data_models import GenomicInterval


def test_distance_non_overlapping(giql_query):
    """
    GIVEN two non-overlapping intervals with a known gap
    WHEN DISTANCE is computed via GIQL
    THEN the distance equals b.start - a.end (half-open arithmetic)
    """
    result = giql_query(
        """
        SELECT DISTANCE(a.interval, b.interval) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[GenomicInterval("chr1", 100, 200, "a1", 100, "+")],
        intervals_b=[GenomicInterval("chr1", 300, 400, "b1", 100, "+")],
    )

    assert len(result) == 1
    # Half-open distance: b.start - a.end = 300 - 200 = 100
    assert result[0][0] == 100


def test_distance_overlapping(giql_query):
    """
    GIVEN two overlapping intervals
    WHEN DISTANCE is computed via GIQL
    THEN the distance is 0
    """
    result = giql_query(
        """
        SELECT DISTANCE(a.interval, b.interval) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[GenomicInterval("chr1", 100, 300, "a1", 100, "+")],
        intervals_b=[GenomicInterval("chr1", 200, 400, "b1", 100, "+")],
    )

    assert len(result) == 1
    assert result[0][0] == 0


def test_distance_adjacent(giql_query):
    """
    GIVEN two adjacent intervals (touching, half-open coordinates)
    WHEN DISTANCE is computed via GIQL
    THEN the distance is 0 (half-open: end of A == start of B means no gap)
    """
    result = giql_query(
        """
        SELECT DISTANCE(a.interval, b.interval) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[GenomicInterval("chr1", 100, 200, "a1", 100, "+")],
        intervals_b=[GenomicInterval("chr1", 200, 300, "b1", 100, "+")],
    )

    assert len(result) == 1
    # Half-open adjacent: b.start - a.end = 200 - 200 = 0
    assert result[0][0] == 0


def test_distance_cross_chromosome(giql_query):
    """
    GIVEN two intervals on different chromosomes
    WHEN DISTANCE is computed via GIQL
    THEN the distance is NULL (cross-chromosome distance is undefined)
    """
    result = giql_query(
        """
        SELECT DISTANCE(a.interval, b.interval) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[GenomicInterval("chr1", 100, 200, "a1", 100, "+")],
        intervals_b=[GenomicInterval("chr2", 100, 200, "b1", 100, "+")],
    )

    assert len(result) == 1
    assert result[0][0] is None, (
        f"Cross-chromosome distance should be NULL, got {result[0][0]}"
    )


def test_distance_signed_downstream(giql_query):
    """
    GIVEN B is downstream of A (B starts after A ends) on + strand
    WHEN DISTANCE with signed := true is computed
    THEN the distance is positive (downstream)
    """
    result = giql_query(
        """
        SELECT DISTANCE(a.interval, b.interval, signed := true) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[GenomicInterval("chr1", 100, 200, "a1", 100, "+")],
        intervals_b=[GenomicInterval("chr1", 300, 400, "b1", 100, "+")],
    )

    assert len(result) == 1
    assert result[0][0] > 0, f"Expected positive (downstream), got {result[0][0]}"


def test_distance_signed_upstream(giql_query):
    """
    GIVEN B is upstream of A (B ends before A starts) on + strand
    WHEN DISTANCE with signed := true is computed
    THEN the distance is negative (upstream)
    """
    result = giql_query(
        """
        SELECT DISTANCE(a.interval, b.interval, signed := true) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[GenomicInterval("chr1", 300, 400, "a1", 100, "+")],
        intervals_b=[GenomicInterval("chr1", 100, 200, "b1", 100, "+")],
    )

    assert len(result) == 1
    assert result[0][0] < 0, f"Expected negative (upstream), got {result[0][0]}"


def test_distance_stranded_unstranded_input(giql_query):
    """
    GIVEN one interval with strand "." (unstranded)
    WHEN DISTANCE with stranded := true is computed
    THEN the distance is NULL (stranded mode requires valid strand)
    """
    result = giql_query(
        """
        SELECT DISTANCE(a.interval, b.interval, stranded := true) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[GenomicInterval("chr1", 100, 200, "a1", 100, ".")],
        intervals_b=[GenomicInterval("chr1", 300, 400, "b1", 100, "+")],
    )

    assert len(result) == 1
    assert result[0][0] is None, (
        f"Stranded distance with '.' strand should be NULL, got {result[0][0]}"
    )


def test_distance_stranded_same_strand(giql_query):
    """
    GIVEN two non-overlapping intervals both on + strand
    WHEN DISTANCE with stranded := true is computed
    THEN the distance is computed normally (same strand is valid)
    """
    result = giql_query(
        """
        SELECT DISTANCE(a.interval, b.interval, stranded := true) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[GenomicInterval("chr1", 100, 200, "a1", 100, "+")],
        intervals_b=[GenomicInterval("chr1", 300, 400, "b1", 100, "+")],
    )

    assert len(result) == 1
    # Both on same + strand, distance should be computed normally
    assert result[0][0] == 100, f"Expected 100, got {result[0][0]}"


def test_distance_signed_stranded_minus_strand(giql_query):
    """
    GIVEN two non-overlapping intervals on - strand, B downstream genomically
    WHEN DISTANCE with signed := true, stranded := true is computed
    THEN the sign is inverted due to - strand (downstream
    genomically is upstream in transcript orientation)
    """
    result = giql_query(
        """
        SELECT DISTANCE(
            a.interval, b.interval,
            signed := true,
            stranded := true
        ) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[GenomicInterval("chr1", 100, 200, "a1", 100, "-")],
        intervals_b=[GenomicInterval("chr1", 300, 400, "b1", 100, "-")],
    )

    assert len(result) == 1
    # On - strand with signed+stranded, genomic downstream becomes
    # transcript upstream, so sign should be negative
    assert result[0][0] < 0, (
        f"Expected negative distance on - strand, got {result[0][0]}"
    )
