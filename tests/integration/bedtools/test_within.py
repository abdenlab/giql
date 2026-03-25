"""Integration tests for GIQL WITHIN operator.

These tests validate that GIQL's WITHIN predicate correctly identifies
intervals that fall entirely within a given range. No direct bedtools
equivalent exists, so tests validate against known expected results.
"""

from giql import transpile

from .utils.bed_export import load_intervals
from .utils.data_models import GenomicInterval


def test_within_basic(duckdb_connection):
    """
    Given:
        A table with intervals of varying sizes
    When:
        WITHIN is used with a range literal
    Then:
        Only intervals fully within the range are returned
    """
    intervals = [
        GenomicInterval("chr1", 150, 250, "inside", 100, "+"),
        GenomicInterval("chr1", 50, 150, "partial_left", 100, "+"),
        GenomicInterval("chr1", 250, 350, "partial_right", 100, "+"),
        GenomicInterval("chr1", 500, 600, "outside", 100, "+"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals",
        [i.to_tuple() for i in intervals],
    )

    sql = transpile(
        "SELECT * FROM intervals WHERE interval WITHIN 'chr1:100-300'",
        tables=["intervals"],
    )
    result = duckdb_connection.execute(sql).fetchall()

    names = {row[3] for row in result}
    # inside [150,250) is within [100,300)
    # partial_left [50,150) starts before 100
    # partial_right [250,350) ends after 300
    assert names == {"inside"}, f"Expected only inside, got {names}"


def test_within_narrow_range(duckdb_connection):
    """
    Given:
        A table with intervals of varying sizes
    When:
        WITHIN is used with a narrow range
    Then:
        Only intervals small enough to fit are returned
    """
    intervals = [
        GenomicInterval("chr1", 152, 158, "tiny", 100, "+"),
        GenomicInterval("chr1", 140, 160, "medium", 100, "+"),
        GenomicInterval("chr1", 100, 200, "large", 100, "+"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals",
        [i.to_tuple() for i in intervals],
    )

    sql = transpile(
        "SELECT * FROM intervals WHERE interval WITHIN 'chr1:150-160'",
        tables=["intervals"],
    )
    result = duckdb_connection.execute(sql).fetchall()

    names = {row[3] for row in result}
    # tiny [152,158) is within [150,160)
    # medium [140,160) starts before 150
    # large [100,200) starts before 150
    assert names == {"tiny"}, f"Expected only tiny, got {names}"


def test_within_column_to_column(duckdb_connection):
    """
    Given:
        Two tables where some intervals in A are within intervals in B
    When:
        a.interval WITHIN b.interval is used in WHERE clause
    Then:
        Only pairs where A is fully within B are returned
    """
    intervals_a = [
        GenomicInterval("chr1", 150, 250, "a_inner", 100, "+"),
        GenomicInterval("chr1", 50, 400, "a_outer", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 100, 300, "b1", 100, "+"),
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
        SELECT a.name, b.name
        FROM intervals_a a, intervals_b b
        WHERE a.interval WITHIN b.interval
        """,
        tables=["intervals_a", "intervals_b"],
    )
    result = duckdb_connection.execute(sql).fetchall()

    pairs = {(row[0], row[1]) for row in result}
    # a_inner [150,250) is within b1 [100,300)
    # a_outer [50,400) is NOT within b1 [100,300)
    assert pairs == {("a_inner", "b1")}, f"Expected one pair, got {pairs}"


def test_within_exact_boundary(duckdb_connection):
    """
    Given:
        An interval whose boundaries exactly match the query range
    When:
        WITHIN is used with that exact range
    Then:
        The interval is returned (exact match counts as within)
    """
    intervals = [
        GenomicInterval("chr1", 100, 200, "exact", 100, "+"),
        GenomicInterval("chr1", 99, 200, "start_outside", 100, "+"),
        GenomicInterval("chr1", 100, 201, "end_outside", 100, "+"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals",
        [i.to_tuple() for i in intervals],
    )

    sql = transpile(
        "SELECT * FROM intervals WHERE interval WITHIN 'chr1:100-200'",
        tables=["intervals"],
    )
    result = duckdb_connection.execute(sql).fetchall()

    names = {row[3] for row in result}
    assert names == {"exact"}, f"Expected only exact, got {names}"
