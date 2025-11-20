"""Integration tests for GIQL MERGE operator.

These tests validate that GIQL's MERGE operator produces identical
results to bedtools merge command.
"""

from giql import GIQLEngine

from .utils.bed_export import load_intervals
from .utils.bedtools_wrapper import merge
from .utils.comparison import compare_results
from .utils.data_models import GenomicInterval


def _setup_giql_engine(duckdb_connection):
    """Helper to set up GIQL engine with table schema."""
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

    engine.register_table_schema(
        "intervals",
        schema,
        genomic_column="position",
    )

    return engine


def test_merge_adjacent_intervals(duckdb_connection):
    """Test MERGE with adjacent intervals.

    Given:
        A set of adjacent intervals
    When:
        MERGE operator is applied
    Then:
        Adjacent intervals are merged into single intervals
    """
    # Arrange
    intervals = [
        GenomicInterval("chr1", 100, 200, "i1", 100, "+"),
        GenomicInterval("chr1", 200, 300, "i2", 150, "+"),
        GenomicInterval("chr1", 300, 400, "i3", 200, "+"),
    ]

    # Load into DuckDB
    load_intervals(
        duckdb_connection,
        "intervals",
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals],
    )

    # Act: Execute bedtools operation using pybedtools
    bedtools_result = merge(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals]
    )

    # Act: Execute GIQL query
    engine = _setup_giql_engine(duckdb_connection)
    giql_query = """
        SELECT MERGE(position)
        FROM intervals
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


def test_merge_overlapping_intervals(duckdb_connection):
    """Test MERGE with overlapping intervals.

    Given:
        A set of overlapping intervals
    When:
        MERGE operator is applied
    Then:
        Overlapping intervals are merged
    """
    # Arrange
    intervals = [
        GenomicInterval("chr1", 100, 250, "i1", 100, "+"),
        GenomicInterval("chr1", 200, 350, "i2", 150, "+"),
        GenomicInterval("chr1", 300, 400, "i3", 200, "+"),
    ]

    # Load into DuckDB
    load_intervals(
        duckdb_connection,
        "intervals",
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals],
    )

    # Act: Execute bedtools operation using pybedtools
    bedtools_result = merge(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals]
    )

    # Act: Execute GIQL query
    engine = _setup_giql_engine(duckdb_connection)
    giql_query = """
        SELECT MERGE(position)
        FROM intervals
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


def test_merge_separated_intervals(duckdb_connection):
    """Test MERGE with separated intervals.

    Given:
        Intervals with gaps between them
    When:
        MERGE operator is applied
    Then:
        Separated intervals remain separate
    """
    # Arrange
    intervals = [
        GenomicInterval("chr1", 100, 200, "i1", 100, "+"),
        GenomicInterval("chr1", 300, 400, "i2", 150, "+"),
        GenomicInterval("chr1", 500, 600, "i3", 200, "+"),
    ]

    # Load into DuckDB
    load_intervals(
        duckdb_connection,
        "intervals",
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals],
    )

    # Act: Execute bedtools operation using pybedtools
    bedtools_result = merge(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals]
    )

    # Act: Execute GIQL query
    engine = _setup_giql_engine(duckdb_connection)
    giql_query = """
        SELECT MERGE(position)
        FROM intervals
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


def test_merge_multiple_chromosomes(duckdb_connection):
    """Test MERGE across multiple chromosomes.

    Given:
        Intervals on different chromosomes
    When:
        MERGE operator is applied
    Then:
        Merging occurs per chromosome
    """
    # Arrange
    intervals = [
        GenomicInterval("chr1", 100, 200, "i1", 100, "+"),
        GenomicInterval("chr1", 180, 300, "i2", 150, "+"),
        GenomicInterval("chr2", 100, 200, "i3", 100, "+"),
        GenomicInterval("chr2", 180, 300, "i4", 150, "+"),
    ]

    # Load into DuckDB
    load_intervals(
        duckdb_connection,
        "intervals",
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals],
    )

    # Act: Execute bedtools operation using pybedtools
    bedtools_result = merge(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals]
    )

    # Act: Execute GIQL query
    engine = _setup_giql_engine(duckdb_connection)
    giql_query = """
        SELECT MERGE(position)
        FROM intervals
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
