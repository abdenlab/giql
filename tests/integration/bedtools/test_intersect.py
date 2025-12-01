"""Integration tests for GIQL INTERSECTS operator.

These tests validate that GIQL's INTERSECTS operator produces identical
results to bedtools intersect command.
"""

from giql import GIQLEngine

from .utils.bed_export import load_intervals
from .utils.bedtools_wrapper import intersect
from .utils.comparison import compare_results
from .utils.data_models import GenomicInterval


def _setup_giql_engine(duckdb_connection):
    """Helper to set up GIQL engine with table schemas."""
    engine = GIQLEngine(target_dialect="duckdb", verbose=False)
    engine.conn = duckdb_connection

    schema = {
        "chromosome": "VARCHAR",
        "start_pos": "BIGINT",
        "end_pos": "BIGINT",
        "name": "VARCHAR",
        "score": "BIGINT",
        "strand": "VARCHAR",
    }

    engine.register_table_schema("intervals_a", schema, genomic_column="interval")
    engine.register_table_schema("intervals_b", schema, genomic_column="interval")

    return engine


def test_intersect_basic_overlap(duckdb_connection, interval_generator):
    """Test INTERSECTS predicate finds overlapping intervals.

    Given:
        Two tables with genomic intervals where some intervals overlap
    When:
        A GIQL query uses INTERSECTS predicate in WHERE clause
    Then:
        Results match bedtools intersect output exactly
    """
    # Arrange: Create overlapping intervals
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
        GenomicInterval("chr1", 150, 250, "a2", 200, "+"),
        GenomicInterval("chr1", 300, 400, "a3", 150, "-"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 180, 220, "b1", 100, "+"),
        GenomicInterval("chr1", 350, 450, "b2", 200, "-"),
    ]

    # Load into DuckDB
    load_intervals(
        duckdb_connection,
        "intervals_a",
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
    )
    load_intervals(
        duckdb_connection,
        "intervals_b",
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
    )

    # Act: Execute bedtools operation using pybedtools
    bedtools_result = intersect(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
    )

    # Act: Execute GIQL query
    engine = _setup_giql_engine(duckdb_connection)
    giql_query = """
        SELECT DISTINCT a.*
        FROM intervals_a a, intervals_b b
        WHERE a.interval INTERSECTS b.interval
    """
    sql = engine.transpile(giql_query)
    giql_result = duckdb_connection.execute(sql).fetchall()

    # Assert: Compare GIQL and bedtools results
    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, (
        f"GIQL results don't match bedtools:\n"
        f"Differences: {comparison.differences}\n"
        f"GIQL rows: {len(giql_result)}, bedtools rows: {len(bedtools_result)}"
    )


def test_intersect_partial_overlap(duckdb_connection, interval_generator):
    """Test INTERSECTS with partially overlapping intervals.

    Given:
        Intervals with partial overlaps
    When:
        INTERSECTS query is executed
    Then:
        Results match bedtools partial overlap behavior
    """
    # Arrange
    intervals_a = [
        GenomicInterval("chr1", 100, 250, "a1", 100, "+"),
        GenomicInterval("chr1", 300, 400, "a2", 200, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 200, 350, "b1", 150, "+"),
    ]

    # Load into DuckDB
    load_intervals(
        duckdb_connection,
        "intervals_a",
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
    )
    load_intervals(
        duckdb_connection,
        "intervals_b",
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
    )

    # Act: Execute bedtools operation using pybedtools
    bedtools_result = intersect(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
    )

    # Act: Execute GIQL query
    engine = _setup_giql_engine(duckdb_connection)
    giql_query = """
        SELECT DISTINCT a.*
        FROM intervals_a a, intervals_b b
        WHERE a.interval INTERSECTS b.interval
    """
    sql = engine.transpile(giql_query)
    giql_result = duckdb_connection.execute(sql).fetchall()

    # Assert: Compare GIQL and bedtools results
    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, (
        f"GIQL results don't match bedtools:\n"
        f"Differences: {comparison.differences}\n"
        f"GIQL rows: {len(giql_result)}, bedtools rows: {len(bedtools_result)}"
    )


def test_intersect_no_overlap(duckdb_connection, interval_generator):
    """Test INTERSECTS with non-overlapping intervals.

    Given:
        Two sets of intervals with no overlaps
    When:
        INTERSECTS query is executed
    Then:
        No results returned (matches bedtools empty output)
    """
    # Arrange
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 300, 400, "b1", 150, "+"),
    ]

    # Load into DuckDB
    load_intervals(
        duckdb_connection,
        "intervals_a",
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
    )
    load_intervals(
        duckdb_connection,
        "intervals_b",
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
    )

    # Act: Execute bedtools operation using pybedtools
    bedtools_result = intersect(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
    )

    # Act: Execute GIQL query
    engine = _setup_giql_engine(duckdb_connection)
    giql_query = """
        SELECT DISTINCT a.*
        FROM intervals_a a, intervals_b b
        WHERE a.interval INTERSECTS b.interval
    """
    sql = engine.transpile(giql_query)
    giql_result = duckdb_connection.execute(sql).fetchall()

    # Assert: Compare GIQL and bedtools results
    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, (
        f"GIQL results don't match bedtools:\n"
        f"Differences: {comparison.differences}\n"
        f"GIQL rows: {len(giql_result)}, bedtools rows: {len(bedtools_result)}"
    )


def test_intersect_adjacent_intervals(duckdb_connection, interval_generator):
    """Test INTERSECTS with adjacent (touching) intervals.

    Given:
        Intervals that touch but don't overlap
    When:
        INTERSECTS query is executed
    Then:
        No results returned (adjacent != overlapping)
    """
    # Arrange: Adjacent intervals (end of a1 == start of b1)
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 200, 300, "b1", 150, "+"),
    ]

    # Load into DuckDB
    load_intervals(
        duckdb_connection,
        "intervals_a",
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
    )
    load_intervals(
        duckdb_connection,
        "intervals_b",
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
    )

    # Act: Execute bedtools operation using pybedtools
    bedtools_result = intersect(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
    )

    # Act: Execute GIQL query
    engine = _setup_giql_engine(duckdb_connection)
    giql_query = """
        SELECT DISTINCT a.*
        FROM intervals_a a, intervals_b b
        WHERE a.interval INTERSECTS b.interval
    """
    sql = engine.transpile(giql_query)
    giql_result = duckdb_connection.execute(sql).fetchall()

    # Assert: Compare GIQL and bedtools results
    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, (
        f"GIQL results don't match bedtools:\n"
        f"Differences: {comparison.differences}\n"
        f"GIQL rows: {len(giql_result)}, bedtools rows: {len(bedtools_result)}"
    )


def test_intersect_multiple_chromosomes(duckdb_connection, interval_generator):
    """Test INTERSECTS across multiple chromosomes.

    Given:
        Intervals on different chromosomes
    When:
        INTERSECTS query is executed
    Then:
        Only same-chromosome overlaps are returned
    """
    # Arrange
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
        GenomicInterval("chr2", 150, 250, "a2", 200, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 150, 250, "b1", 150, "+"),
        GenomicInterval("chr2", 200, 300, "b2", 100, "+"),
    ]

    # Load into DuckDB
    load_intervals(
        duckdb_connection,
        "intervals_a",
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
    )
    load_intervals(
        duckdb_connection,
        "intervals_b",
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
    )

    # Act: Execute bedtools operation using pybedtools
    bedtools_result = intersect(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
    )

    # Act: Execute GIQL query
    engine = _setup_giql_engine(duckdb_connection)
    giql_query = """
        SELECT DISTINCT a.*
        FROM intervals_a a, intervals_b b
        WHERE a.interval INTERSECTS b.interval
    """
    sql = engine.transpile(giql_query)
    giql_result = duckdb_connection.execute(sql).fetchall()

    # Assert: Compare GIQL and bedtools results
    comparison = compare_results(giql_result, bedtools_result)
    assert comparison.match, (
        f"GIQL results don't match bedtools:\n"
        f"Differences: {comparison.differences}\n"
        f"GIQL rows: {len(giql_result)}, bedtools rows: {len(bedtools_result)}"
    )
