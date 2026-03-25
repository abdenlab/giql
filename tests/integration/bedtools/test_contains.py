"""Integration tests for GIQL CONTAINS operator.

These tests validate that GIQL's CONTAINS predicate correctly identifies
intervals that fully contain a point or range. No direct bedtools
equivalent exists, so tests validate against known expected results.
"""

from giql import transpile

from .utils.bed_export import load_intervals
from .utils.data_models import GenomicInterval


def test_contains_point(duckdb_connection):
    """
    Given:
        A table with intervals of varying sizes on chr1
    When:
        CONTAINS is used with a point literal
    Then:
        Only intervals that contain the point are returned
    """
    intervals = [
        GenomicInterval("chr1", 100, 300, "wide", 100, "+"),
        GenomicInterval("chr1", 140, 160, "narrow", 100, "+"),
        GenomicInterval("chr1", 200, 400, "partial", 100, "+"),
        GenomicInterval("chr1", 500, 600, "far", 100, "+"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals",
        [i.to_tuple() for i in intervals],
    )

    sql = transpile(
        "SELECT * FROM intervals WHERE interval CONTAINS 'chr1:150'",
        tables=["intervals"],
    )
    result = duckdb_connection.execute(sql).fetchall()

    names = {row[3] for row in result}
    # wide [100,300) contains 150, narrow [140,160) contains 150
    # partial [200,400) does not contain 150, far [500,600) does not
    assert names == {"wide", "narrow"}, f"Expected wide+narrow, got {names}"


def test_contains_range(duckdb_connection):
    """
    Given:
        A table with intervals of varying sizes
    When:
        CONTAINS is used with a range literal
    Then:
        Only intervals that fully contain the range are returned
    """
    intervals = [
        GenomicInterval("chr1", 100, 400, "large", 100, "+"),
        GenomicInterval("chr1", 150, 250, "medium", 100, "+"),
        GenomicInterval("chr1", 180, 220, "small", 100, "+"),
        GenomicInterval("chr1", 500, 600, "far", 100, "+"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals",
        [i.to_tuple() for i in intervals],
    )

    sql = transpile(
        "SELECT * FROM intervals WHERE interval CONTAINS 'chr1:150-250'",
        tables=["intervals"],
    )
    result = duckdb_connection.execute(sql).fetchall()

    names = {row[3] for row in result}
    # large [100,400) contains [150,250), medium [150,250) contains [150,250)
    # small [180,220) does not fully contain [150,250)
    assert names == {"large", "medium"}, f"Expected large+medium, got {names}"


def test_contains_column_to_column(duckdb_connection):
    """
    Given:
        Two tables where some intervals in A fully contain intervals in B
    When:
        a.interval CONTAINS b.interval is used in WHERE clause
    Then:
        Only pairs where A fully contains B are returned
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 400, "a_large", 100, "+"),
        GenomicInterval("chr1", 200, 250, "a_small", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 150, 300, "b1", 100, "+"),
        GenomicInterval("chr1", 210, 240, "b2", 100, "+"),
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
        WHERE a.interval CONTAINS b.interval
        """,
        tables=["intervals_a", "intervals_b"],
    )
    result = duckdb_connection.execute(sql).fetchall()

    pairs = {(row[0], row[1]) for row in result}
    # a_large [100,400) contains b1 [150,300) and b2 [210,240)
    # a_small [200,250) contains b2 [210,240) but not b1 [150,300)
    assert pairs == {("a_large", "b1"), ("a_large", "b2"), ("a_small", "b2")}, (
        f"Expected 3 containment pairs, got {pairs}"
    )


def test_contains_cross_chromosome(duckdb_connection):
    """
    Given:
        A table with intervals on multiple chromosomes
    When:
        CONTAINS is used with a chr1 point literal
    Then:
        Only chr1 intervals are considered
    """
    intervals = [
        GenomicInterval("chr1", 100, 300, "chr1_hit", 100, "+"),
        GenomicInterval("chr2", 100, 300, "chr2_miss", 100, "+"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals",
        [i.to_tuple() for i in intervals],
    )

    sql = transpile(
        "SELECT * FROM intervals WHERE interval CONTAINS 'chr1:150'",
        tables=["intervals"],
    )
    result = duckdb_connection.execute(sql).fetchall()

    names = {row[3] for row in result}
    assert names == {"chr1_hit"}


def test_contains_all_set_predicate(duckdb_connection):
    """
    Given:
        A table with intervals of varying sizes
    When:
        CONTAINS ALL is used with multiple points
    Then:
        Only intervals containing all points are returned
    """
    intervals = [
        GenomicInterval("chr1", 100, 400, "large", 100, "+"),
        GenomicInterval("chr1", 100, 200, "left", 100, "+"),
        GenomicInterval("chr1", 250, 400, "right", 100, "+"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals",
        [i.to_tuple() for i in intervals],
    )

    sql = transpile(
        """
        SELECT * FROM intervals
        WHERE interval CONTAINS ALL('chr1:150', 'chr1:300')
        """,
        tables=["intervals"],
    )
    result = duckdb_connection.execute(sql).fetchall()

    names = {row[3] for row in result}
    # Only large [100,400) contains both 150 and 300
    assert names == {"large"}, f"Expected only large, got {names}"
