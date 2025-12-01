"""Integration tests for GIQL NEAREST operator.

These tests validate that GIQL's NEAREST operator produces identical
results to bedtools closest command.
"""

from giql import GIQLEngine

from .utils.bed_export import load_intervals
from .utils.bedtools_wrapper import closest
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

    engine.register_table_schema(
        "intervals_a",
        schema,
        genomic_column="interval",
        interval_type="closed",  # Match bedtools distance calculation
    )
    engine.register_table_schema(
        "intervals_b",
        schema,
        genomic_column="interval",
        interval_type="closed",  # Match bedtools distance calculation
    )

    return engine


def test_nearest_non_overlapping(duckdb_connection):
    """Test NEAREST with non-overlapping intervals.

    Given:
        Two sets of non-overlapping intervals
    When:
        NEAREST operator is applied
    Then:
        Each interval in A finds its closest neighbor in B
    """
    # Arrange
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
        GenomicInterval("chr1", 500, 600, "a2", 150, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 250, 300, "b1", 100, "+"),
        GenomicInterval("chr1", 350, 400, "b2", 150, "+"),
        GenomicInterval("chr1", 700, 800, "b3", 200, "+"),
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
    bedtools_result = closest(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
    )

    # Act: Execute GIQL query
    engine = _setup_giql_engine(duckdb_connection)
    giql_query = """
        SELECT a.*, b.*
        FROM intervals_a a, NEAREST(intervals_b, k=1) b
        ORDER BY a.chromosome, a.start_pos
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


def test_nearest_multiple_candidates(duckdb_connection):
    """Test NEAREST with equidistant intervals.

    Given:
        Interval in A with multiple equidistant intervals in B
    When:
        NEAREST operator is applied
    Then:
        Bedtools reports one of the equidistant intervals (tie-breaking behavior)
    """
    # Arrange: a1 is equidistant from b1 and b2
    intervals_a = [
        GenomicInterval("chr1", 300, 400, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 100, 200, "b1", 100, "+"),  # Distance: 100 bp
        GenomicInterval("chr1", 500, 600, "b2", 150, "+"),  # Distance: 100 bp
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
    bedtools_result = closest(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
    )

    # Act: Execute GIQL query
    engine = _setup_giql_engine(duckdb_connection)
    giql_query = """
        SELECT a.*, b.*
        FROM intervals_a a, NEAREST(intervals_b, k=1) b
        ORDER BY a.chromosome, a.start_pos
    """
    sql = engine.transpile(giql_query)
    giql_result = duckdb_connection.execute(sql).fetchall()

    # Assert: Compare GIQL and bedtools results (allowing tie-breaking differences)
    assert len(giql_result) == len(bedtools_result)
    # The nearest interval is either b1 or b2 (both equidistant)
    assert giql_result[0][3] == "a1"  # Interval A name
    assert giql_result[0][9] in ("b1", "b2")  # Nearest could be either


def test_nearest_cross_chromosome(duckdb_connection):
    """Test NEAREST across multiple chromosomes.

    Given:
        Intervals on different chromosomes
    When:
        NEAREST operator is applied
    Then:
        Each interval finds nearest only on same chromosome
    """
    # Arrange
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
        GenomicInterval("chr2", 100, 200, "a2", 150, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 300, 400, "b1", 100, "+"),
        GenomicInterval("chr2", 300, 400, "b2", 150, "+"),
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
    bedtools_result = closest(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
    )

    # Act: Execute GIQL query
    engine = _setup_giql_engine(duckdb_connection)
    giql_query = """
        SELECT a.*, b.*
        FROM intervals_a a, NEAREST(intervals_b, k=1) b
        ORDER BY a.chromosome, a.start_pos
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


def test_nearest_boundary_cases(duckdb_connection):
    """Test NEAREST with boundary cases.

    Given:
        Adjacent intervals (touching but not overlapping)
    When:
        NEAREST operator is applied
    Then:
        Adjacent intervals are reported as nearest (distance = 0)
    """
    # Arrange: a1 ends where b1 starts (adjacent, distance = 0)
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 200, 300, "b1", 150, "+"),  # Adjacent to a1
        GenomicInterval("chr1", 500, 600, "b2", 200, "+"),  # Far away
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
    bedtools_result = closest(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
    )

    # Act: Execute GIQL query
    engine = _setup_giql_engine(duckdb_connection)
    giql_query = """
        SELECT a.*, b.*
        FROM intervals_a a, NEAREST(intervals_b, k=1) b
        ORDER BY a.chromosome, a.start_pos
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
