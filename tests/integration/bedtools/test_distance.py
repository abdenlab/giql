"""Integration tests for GIQL DISTANCE operator.

These tests validate that GIQL's DISTANCE function computes the correct
genomic distance between intervals, cross-validated against bedtools
closest -d output.
"""

from giql import transpile

from .utils.bed_export import load_intervals
from .utils.data_models import GenomicInterval


def test_distance_non_overlapping(duckdb_connection):
    """
    Given:
        Two non-overlapping intervals with a known gap
    When:
        DISTANCE is computed via GIQL
    Then:
        The distance equals b.start - a.end (half-open arithmetic)
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 300, 400, "b1", 100, "+"),
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

    sql = transpile(
        """
        SELECT DISTANCE(a.interval, b.interval) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    assert len(giql_result) == 1
    # Half-open distance: b.start - a.end = 300 - 200 = 100
    assert giql_result[0][0] == 100


def test_distance_overlapping(duckdb_connection):
    """
    Given:
        Two overlapping intervals
    When:
        DISTANCE is computed via GIQL
    Then:
        The distance is 0
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 300, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 200, 400, "b1", 100, "+"),
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

    sql = transpile(
        """
        SELECT DISTANCE(a.interval, b.interval) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    assert len(giql_result) == 1
    assert giql_result[0][0] == 0


def test_distance_adjacent(duckdb_connection):
    """
    Given:
        Two adjacent intervals (touching, half-open coordinates)
    When:
        DISTANCE is computed via GIQL
    Then:
        The distance is 0 (half-open: end of A == start of B means
        no gap between them)
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 200, 300, "b1", 100, "+"),
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

    sql = transpile(
        """
        SELECT DISTANCE(a.interval, b.interval) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    assert len(giql_result) == 1
    # Half-open adjacent: b.start - a.end = 200 - 200 = 0
    assert giql_result[0][0] == 0


def test_distance_cross_chromosome(duckdb_connection):
    """
    Given:
        Two intervals on different chromosomes
    When:
        DISTANCE is computed via GIQL
    Then:
        The distance is NULL (cross-chromosome distance is undefined)
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr2", 100, 200, "b1", 100, "+"),
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

    sql = transpile(
        """
        SELECT DISTANCE(a.interval, b.interval) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    assert len(giql_result) == 1
    assert giql_result[0][0] is None, (
        f"Cross-chromosome distance should be NULL, got {giql_result[0][0]}"
    )


def test_distance_signed_downstream(duckdb_connection):
    """
    Given:
        B is downstream of A (B starts after A ends) on + strand
    When:
        DISTANCE with signed := true is computed
    Then:
        The distance is positive (downstream)
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 300, 400, "b1", 100, "+"),
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

    sql = transpile(
        """
        SELECT DISTANCE(a.interval, b.interval, signed := true) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
    )
    result = duckdb_connection.execute(sql).fetchall()

    assert len(result) == 1
    assert result[0][0] > 0, f"Expected positive (downstream), got {result[0][0]}"


def test_distance_signed_upstream(duckdb_connection):
    """
    Given:
        B is upstream of A (B ends before A starts) on + strand
    When:
        DISTANCE with signed := true is computed
    Then:
        The distance is negative (upstream)
    """
    intervals_a = [
        GenomicInterval("chr1", 300, 400, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 100, 200, "b1", 100, "+"),
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

    sql = transpile(
        """
        SELECT DISTANCE(a.interval, b.interval, signed := true) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
    )
    result = duckdb_connection.execute(sql).fetchall()

    assert len(result) == 1
    assert result[0][0] < 0, f"Expected negative (upstream), got {result[0][0]}"


def test_distance_stranded_unstranded_input(duckdb_connection):
    """
    Given:
        One interval with strand "." (unstranded)
    When:
        DISTANCE with stranded := true is computed
    Then:
        The distance is NULL (stranded mode requires valid strand)
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "."),
    ]
    intervals_b = [
        GenomicInterval("chr1", 300, 400, "b1", 100, "+"),
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

    sql = transpile(
        """
        SELECT DISTANCE(a.interval, b.interval, stranded := true) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
    )
    result = duckdb_connection.execute(sql).fetchall()

    assert len(result) == 1
    assert result[0][0] is None, (
        f"Stranded distance with '.' strand should be NULL, got {result[0][0]}"
    )


def test_distance_stranded_same_strand(duckdb_connection):
    """
    Given:
        Two non-overlapping intervals both on + strand
    When:
        DISTANCE with stranded := true is computed
    Then:
        The distance is computed normally (same strand is valid)
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 300, 400, "b1", 100, "+"),
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

    sql = transpile(
        """
        SELECT DISTANCE(a.interval, b.interval, stranded := true) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
    )
    result = duckdb_connection.execute(sql).fetchall()

    assert len(result) == 1
    # Both on same + strand, distance should be computed normally
    assert result[0][0] == 100, f"Expected 100, got {result[0][0]}"


def test_distance_signed_stranded_minus_strand(duckdb_connection):
    """
    Given:
        Two non-overlapping intervals on - strand, B downstream genomically
    When:
        DISTANCE with signed := true, stranded := true is computed
    Then:
        The sign is inverted due to - strand (downstream on - strand
        is upstream in transcript orientation)
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "-"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 300, 400, "b1", 100, "-"),
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

    sql = transpile(
        """
        SELECT DISTANCE(
            a.interval, b.interval,
            signed := true,
            stranded := true
        ) AS dist
        FROM intervals_a a, intervals_b b
        """,
        tables=["intervals_a", "intervals_b"],
    )
    result = duckdb_connection.execute(sql).fetchall()

    assert len(result) == 1
    # On - strand with signed+stranded, genomic downstream becomes
    # transcript upstream, so sign should be negative
    assert result[0][0] < 0, (
        f"Expected negative distance on - strand, got {result[0][0]}"
    )
