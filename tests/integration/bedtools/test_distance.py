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
