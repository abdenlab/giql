"""Integration tests for GIQL MERGE operator.

These tests validate that GIQL's MERGE operator produces identical
results to bedtools merge command.
"""

from giql import transpile

from .utils.bed_export import load_intervals
from .utils.bedtools_wrapper import merge
from .utils.comparison import compare_results
from .utils.data_models import GenomicInterval


def test_merge_adjacent_intervals(duckdb_connection):
    """
    Given:
        A set of adjacent intervals (bookended, half-open)
    When:
        MERGE operator is applied
    Then:
        Adjacent intervals are merged into single intervals
    """
    intervals = [
        GenomicInterval("chr1", 100, 200, "i1", 100, "+"),
        GenomicInterval("chr1", 200, 300, "i2", 150, "+"),
        GenomicInterval("chr1", 300, 400, "i3", 200, "+"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals",
        [i.to_tuple() for i in intervals],
    )

    bedtools_result = merge([i.to_tuple() for i in intervals])

    sql = transpile(
        "SELECT MERGE(interval) FROM intervals",
        tables=["intervals"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, comparison.failure_message()


def test_merge_overlapping_intervals(duckdb_connection):
    """
    Given:
        A set of overlapping intervals
    When:
        MERGE operator is applied
    Then:
        Overlapping intervals are merged
    """
    intervals = [
        GenomicInterval("chr1", 100, 250, "i1", 100, "+"),
        GenomicInterval("chr1", 200, 350, "i2", 150, "+"),
        GenomicInterval("chr1", 300, 400, "i3", 200, "+"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals",
        [i.to_tuple() for i in intervals],
    )

    bedtools_result = merge([i.to_tuple() for i in intervals])

    sql = transpile(
        "SELECT MERGE(interval) FROM intervals",
        tables=["intervals"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, comparison.failure_message()


def test_merge_separated_intervals(duckdb_connection):
    """
    Given:
        Intervals with gaps between them
    When:
        MERGE operator is applied
    Then:
        Separated intervals remain separate
    """
    intervals = [
        GenomicInterval("chr1", 100, 200, "i1", 100, "+"),
        GenomicInterval("chr1", 300, 400, "i2", 150, "+"),
        GenomicInterval("chr1", 500, 600, "i3", 200, "+"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals",
        [i.to_tuple() for i in intervals],
    )

    bedtools_result = merge([i.to_tuple() for i in intervals])

    sql = transpile(
        "SELECT MERGE(interval) FROM intervals",
        tables=["intervals"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, comparison.failure_message()


def test_merge_multiple_chromosomes(duckdb_connection):
    """
    Given:
        Intervals on different chromosomes with overlaps within each
    When:
        MERGE operator is applied
    Then:
        Merging occurs per chromosome independently
    """
    intervals = [
        GenomicInterval("chr1", 100, 200, "i1", 100, "+"),
        GenomicInterval("chr1", 180, 300, "i2", 150, "+"),
        GenomicInterval("chr2", 100, 200, "i3", 100, "+"),
        GenomicInterval("chr2", 180, 300, "i4", 150, "+"),
    ]

    load_intervals(
        duckdb_connection,
        "intervals",
        [i.to_tuple() for i in intervals],
    )

    bedtools_result = merge([i.to_tuple() for i in intervals])

    sql = transpile(
        "SELECT MERGE(interval) FROM intervals",
        tables=["intervals"],
    )
    giql_result = duckdb_connection.execute(sql).fetchall()

    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, comparison.failure_message()
