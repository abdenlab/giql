"""Integration tests for GIQL INTERSECTS operator.

These tests validate that GIQL's INTERSECTS operator produces identical
results to bedtools intersect command.
"""

from giql import transpile

from .utils.bed_export import load_intervals
from .utils.bedtools_wrapper import intersect
from .utils.comparison import compare_results
from .utils.data_models import GenomicInterval


def test_intersect_basic_overlap(duckdb_connection):
    """
    Given:
        Two tables with genomic intervals where some intervals overlap
    When:
        A GIQL query uses INTERSECTS predicate in WHERE clause
    Then:
        Results match bedtools intersect output exactly
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
        GenomicInterval("chr1", 150, 250, "a2", 200, "+"),
        GenomicInterval("chr1", 300, 400, "a3", 150, "-"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 180, 220, "b1", 100, "+"),
        GenomicInterval("chr1", 350, 450, "b2", 200, "-"),
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


def test_intersect_partial_overlap(duckdb_connection):
    """
    Given:
        Intervals with partial overlaps
    When:
        INTERSECTS query is executed
    Then:
        Results match bedtools partial overlap behavior
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 250, "a1", 100, "+"),
        GenomicInterval("chr1", 300, 400, "a2", 200, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 200, 350, "b1", 150, "+"),
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


def test_intersect_no_overlap(duckdb_connection):
    """
    Given:
        Two sets of intervals with no overlaps
    When:
        INTERSECTS query is executed
    Then:
        No results returned (matches bedtools empty output)
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 300, 400, "b1", 150, "+"),
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


def test_intersect_adjacent_intervals(duckdb_connection):
    """
    Given:
        Intervals that touch but don't overlap (half-open coordinates)
    When:
        INTERSECTS query is executed
    Then:
        No results returned (adjacent != overlapping)
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 200, 300, "b1", 150, "+"),
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


def test_intersect_multiple_chromosomes(duckdb_connection):
    """
    Given:
        Intervals on different chromosomes
    When:
        INTERSECTS query is executed
    Then:
        Only same-chromosome overlaps are returned
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
        GenomicInterval("chr2", 150, 250, "a2", 200, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 150, 250, "b1", 150, "+"),
        GenomicInterval("chr2", 200, 300, "b2", 100, "+"),
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
