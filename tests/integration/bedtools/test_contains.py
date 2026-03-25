"""Integration tests for GIQL CONTAINS operator.

These tests validate that GIQL's CONTAINS predicate correctly identifies
intervals that fully contain a point or range. No direct bedtools
equivalent exists, so tests validate against known expected results.
"""

from .utils.data_models import GenomicInterval


def test_contains_point(giql_query):
    """
    GIVEN a table with intervals of varying sizes on chr1
    WHEN CONTAINS is used with a point literal
    THEN only intervals that contain the point are returned
    """
    result = giql_query(
        "SELECT * FROM intervals WHERE interval CONTAINS 'chr1:150'",
        tables=["intervals"],
        intervals=[
            GenomicInterval("chr1", 100, 300, "wide", 100, "+"),
            GenomicInterval("chr1", 140, 160, "narrow", 100, "+"),
            GenomicInterval("chr1", 200, 400, "partial", 100, "+"),
            GenomicInterval("chr1", 500, 600, "far", 100, "+"),
        ],
    )

    names = {row[3] for row in result}
    # wide [100,300) contains 150, narrow [140,160) contains 150
    # partial [200,400) does not contain 150, far [500,600) does not
    assert names == {"wide", "narrow"}, f"Expected wide+narrow, got {names}"


def test_contains_range(giql_query):
    """
    GIVEN a table with intervals of varying sizes
    WHEN CONTAINS is used with a range literal
    THEN only intervals that fully contain the range are returned
    """
    result = giql_query(
        "SELECT * FROM intervals WHERE interval CONTAINS 'chr1:150-250'",
        tables=["intervals"],
        intervals=[
            GenomicInterval("chr1", 100, 400, "large", 100, "+"),
            GenomicInterval("chr1", 150, 250, "medium", 100, "+"),
            GenomicInterval("chr1", 180, 220, "small", 100, "+"),
            GenomicInterval("chr1", 500, 600, "far", 100, "+"),
        ],
    )

    names = {row[3] for row in result}
    # large [100,400) contains [150,250), medium [150,250) contains [150,250)
    # small [180,220) does not fully contain [150,250)
    assert names == {"large", "medium"}, f"Expected large+medium, got {names}"


def test_contains_column_to_column(giql_query):
    """
    GIVEN two tables where some intervals in A fully contain intervals in B
    WHEN a.interval CONTAINS b.interval is used in WHERE clause
    THEN only pairs where A fully contains B are returned
    """
    result = giql_query(
        """
        SELECT a.name, b.name
        FROM intervals_a a, intervals_b b
        WHERE a.interval CONTAINS b.interval
        """,
        tables=["intervals_a", "intervals_b"],
        intervals_a=[
            GenomicInterval("chr1", 100, 400, "a_large", 100, "+"),
            GenomicInterval("chr1", 200, 250, "a_small", 100, "+"),
        ],
        intervals_b=[
            GenomicInterval("chr1", 150, 300, "b1", 100, "+"),
            GenomicInterval("chr1", 210, 240, "b2", 100, "+"),
        ],
    )

    pairs = {(row[0], row[1]) for row in result}
    # a_large [100,400) contains b1 [150,300) and b2 [210,240)
    # a_small [200,250) contains b2 [210,240) but not b1 [150,300)
    assert pairs == {("a_large", "b1"), ("a_large", "b2"), ("a_small", "b2")}, (
        f"Expected 3 containment pairs, got {pairs}"
    )


def test_contains_cross_chromosome(giql_query):
    """
    GIVEN a table with intervals on multiple chromosomes
    WHEN CONTAINS is used with a chr1 point literal
    THEN only chr1 intervals are considered
    """
    result = giql_query(
        "SELECT * FROM intervals WHERE interval CONTAINS 'chr1:150'",
        tables=["intervals"],
        intervals=[
            GenomicInterval("chr1", 100, 300, "chr1_hit", 100, "+"),
            GenomicInterval("chr2", 100, 300, "chr2_miss", 100, "+"),
        ],
    )

    names = {row[3] for row in result}
    assert names == {"chr1_hit"}


def test_contains_all_set_predicate(giql_query):
    """
    GIVEN a table with intervals of varying sizes
    WHEN CONTAINS ALL is used with multiple points
    THEN only intervals containing all points are returned
    """
    result = giql_query(
        """
        SELECT * FROM intervals
        WHERE interval CONTAINS ALL('chr1:150', 'chr1:300')
        """,
        tables=["intervals"],
        intervals=[
            GenomicInterval("chr1", 100, 400, "large", 100, "+"),
            GenomicInterval("chr1", 100, 200, "left", 100, "+"),
            GenomicInterval("chr1", 250, 400, "right", 100, "+"),
        ],
    )

    names = {row[3] for row in result}
    # Only large [100,400) contains both 150 and 300
    assert names == {"large"}, f"Expected only large, got {names}"
