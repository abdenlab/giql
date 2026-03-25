"""Integration tests for GIQL INTERSECTS operator.

These tests validate that GIQL's INTERSECTS operator produces identical
results to bedtools intersect command.
"""

from giql import transpile

from .utils.bedtools_wrapper import intersect
from .utils.comparison import compare_results
from .utils.data_models import GenomicInterval
from .utils.duckdb_loader import load_intervals


def test_intersect_basic_overlap(duckdb_connection):
    """
    GIVEN two tables with genomic intervals where some intervals overlap
    WHEN a GIQL query uses INTERSECTS predicate in WHERE clause
    THEN results match bedtools intersect output exactly
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
    GIVEN intervals with partial overlaps
    WHEN INTERSECTS query is executed
    THEN results match bedtools partial overlap behavior
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
    GIVEN two sets of intervals with no overlaps
    WHEN INTERSECTS query is executed
    THEN no results returned matching bedtools empty output
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
    GIVEN intervals that touch but don't overlap using half-open coordinates
    WHEN INTERSECTS query is executed
    THEN no results returned because adjacent does not mean overlapping
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
    GIVEN intervals on different chromosomes
    WHEN INTERSECTS query is executed
    THEN only same-chromosome overlaps are returned
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


def test_intersect_literal_range(giql_query):
    """
    GIVEN a table with intervals on chr1
    WHEN INTERSECTS is used with a literal range string
    THEN only intervals overlapping the literal range are returned
    """
    intervals = [
        GenomicInterval("chr1", 100, 200, "hit1", 100, "+"),
        GenomicInterval("chr1", 180, 250, "hit2", 100, "+"),
        GenomicInterval("chr1", 300, 400, "miss", 100, "+"),
    ]

    result = giql_query(
        "SELECT * FROM intervals WHERE interval INTERSECTS 'chr1:150-220'",
        tables=["intervals"],
        intervals=intervals,
    )

    names = {row[3] for row in result}
    assert names == {"hit1", "hit2"}, f"Expected hit1+hit2, got {names}"


def test_intersect_literal_cross_chromosome(giql_query):
    """
    GIVEN a table with intervals on chr1 and chr2
    WHEN INTERSECTS is used with a chr2 literal range
    THEN only chr2 intervals overlapping the range are returned
    """
    intervals = [
        GenomicInterval("chr1", 100, 200, "chr1_int", 100, "+"),
        GenomicInterval("chr2", 100, 200, "chr2_hit", 100, "+"),
        GenomicInterval("chr2", 300, 400, "chr2_miss", 100, "+"),
    ]

    result = giql_query(
        "SELECT * FROM intervals WHERE interval INTERSECTS 'chr2:150-250'",
        tables=["intervals"],
        intervals=intervals,
    )

    names = {row[3] for row in result}
    assert names == {"chr2_hit"}, f"Expected only chr2_hit, got {names}"


def test_intersect_any_set_predicate(giql_query):
    """
    GIVEN a table with intervals across chromosomes
    WHEN INTERSECTS ANY is used with multiple ranges
    THEN intervals overlapping any of the ranges are returned
    """
    intervals = [
        GenomicInterval("chr1", 100, 200, "chr1_hit", 100, "+"),
        GenomicInterval("chr2", 300, 400, "chr2_hit", 100, "+"),
        GenomicInterval("chr3", 100, 200, "chr3_miss", 100, "+"),
    ]

    result = giql_query(
        """
        SELECT * FROM intervals
        WHERE interval INTERSECTS ANY('chr1:150-250', 'chr2:350-450')
        """,
        tables=["intervals"],
        intervals=intervals,
    )

    names = {row[3] for row in result}
    assert names == {"chr1_hit", "chr2_hit"}, f"Expected chr1_hit+chr2_hit, got {names}"


def test_intersect_all_set_predicate(giql_query):
    """
    GIVEN a table with intervals of varying sizes
    WHEN INTERSECTS ALL is used with two ranges
    THEN only intervals overlapping both ranges are returned
    """
    intervals = [
        GenomicInterval("chr1", 100, 400, "large", 100, "+"),
        GenomicInterval("chr1", 100, 200, "left_only", 100, "+"),
        GenomicInterval("chr1", 250, 400, "right_only", 100, "+"),
    ]

    result = giql_query(
        """
        SELECT * FROM intervals
        WHERE interval INTERSECTS ALL('chr1:120-180', 'chr1:280-350')
        """,
        tables=["intervals"],
        intervals=intervals,
    )

    names = {row[3] for row in result}
    assert names == {"large"}, f"Expected only large, got {names}"
