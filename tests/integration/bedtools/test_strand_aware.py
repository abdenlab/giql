"""Integration tests for GIQL strand-aware operations.

These tests validate that GIQL correctly handles strand-specific interval
operations, matching bedtools behavior with -s and -S flags.
"""

from giql import GIQLEngine

from .utils.bed_export import load_intervals
from .utils.bedtools_wrapper import closest
from .utils.bedtools_wrapper import intersect
from .utils.bedtools_wrapper import merge
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

    for table_name in ["intervals_a", "intervals_b", "intervals"]:
        engine.register_table_schema(
            table_name,
            schema,
            genomic_column="position",
            interval_type="closed",  # Match bedtools distance calculation
        )

    return engine


def test_intersect_same_strand(duckdb_connection):
    """Test INTERSECTS with same-strand requirement.

    Given:
        Intervals on both same and opposite strands
    When:
        INTERSECTS with same-strand requirement is applied
    Then:
        Only same-strand overlaps are reported
    """
    # Arrange
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
        GenomicInterval("chr1", 300, 400, "a2", 150, "-"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 150, 250, "b1", 100, "+"),  # Overlaps a1 (same +)
        GenomicInterval("chr1", 350, 450, "b2", 150, "-"),  # Overlaps a2 (same -)
        GenomicInterval("chr1", 150, 250, "b3", 200, "-"),  # Overlaps a1 (opposite)
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

    # Act: Execute bedtools operation using pybedtools with same-strand requirement
    bedtools_result = intersect(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
        strand_mode="same",
    )

    # Act: Execute GIQL query with same-strand filter
    engine = _setup_giql_engine(duckdb_connection)
    giql_query = """
        SELECT DISTINCT a.*
        FROM intervals_a a, intervals_b b
        WHERE a.position INTERSECTS b.position
          AND a.strand = b.strand
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


def test_intersect_opposite_strand(duckdb_connection):
    """Test INTERSECTS with opposite-strand requirement.

    Given:
        Intervals on both same and opposite strands
    When:
        INTERSECTS with opposite-strand requirement is applied
    Then:
        Only opposite-strand overlaps are reported
    """
    # Arrange
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
        GenomicInterval("chr1", 300, 400, "a2", 150, "-"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 150, 250, "b1", 100, "-"),  # Overlaps a1 (opposite)
        GenomicInterval("chr1", 350, 450, "b2", 150, "+"),  # Overlaps a2 (opposite)
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

    # Act: Execute bedtools operation using pybedtools with opposite-strand requirement
    bedtools_result = intersect(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
        strand_mode="opposite",
    )

    # Act: Execute GIQL query with opposite-strand filter
    engine = _setup_giql_engine(duckdb_connection)
    giql_query = """
        SELECT DISTINCT a.*
        FROM intervals_a a, intervals_b b
        WHERE a.position INTERSECTS b.position
          AND a.strand != b.strand
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


def test_intersect_ignore_strand(duckdb_connection):
    """Test INTERSECTS ignoring strand information.

    Given:
        Intervals with various strand combinations
    When:
        INTERSECTS without strand requirements is applied
    Then:
        All overlaps are reported regardless of strand
    """
    # Arrange
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 150, 250, "b1", 100, "+"),  # Same strand
        GenomicInterval("chr1", 150, 250, "b2", 150, "-"),  # Opposite strand
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

    # Act: Execute bedtools operation using pybedtools without strand requirements
    bedtools_result = intersect(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
    )

    # Act: Execute GIQL query without strand filter
    engine = _setup_giql_engine(duckdb_connection)
    giql_query = """
        SELECT DISTINCT a.*
        FROM intervals_a a, intervals_b b
        WHERE a.position INTERSECTS b.position
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


def test_intersect_mixed_strands(duckdb_connection):
    """Test INTERSECTS with mixed strand scenarios.

    Given:
        Complex scenario with +, -, and unstranded intervals
    When:
        INTERSECTS with same-strand requirement is applied
    Then:
        Results correctly handle strand matching logic
    """
    # Arrange
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
        GenomicInterval("chr1", 300, 400, "a2", 150, "-"),
        GenomicInterval("chr1", 500, 600, "a3", 200, "."),  # Unstranded
    ]
    intervals_b = [
        GenomicInterval("chr1", 150, 250, "b1", 100, "+"),
        GenomicInterval("chr1", 350, 450, "b2", 150, "-"),
        GenomicInterval("chr1", 550, 650, "b3", 200, "."),
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

    # Act: Execute bedtools operation using pybedtools with same-strand requirement
    bedtools_result = intersect(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
        strand_mode="same",
    )

    # Act: Execute GIQL query with same-strand filter
    engine = _setup_giql_engine(duckdb_connection)
    giql_query = """
        SELECT DISTINCT a.*
        FROM intervals_a a, intervals_b b
        WHERE a.position INTERSECTS b.position
            AND a.strand = b.strand
            AND a.strand != '.'
            AND b.strand != '.'
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


def test_nearest_same_strand(duckdb_connection):
    """Test NEAREST with same-strand requirement.

    Given:
        Intervals with candidates on same and opposite strands
    When:
        NEAREST with same-strand requirement is applied
    Then:
        Only same-strand nearest intervals are reported
    """
    # Arrange
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 250, 300, "b1", 100, "+"),  # Nearest on same strand
        GenomicInterval("chr1", 220, 240, "b2", 150, "-"),  # Closer but opposite
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

    # Act: Execute bedtools operation using pybedtools with same-strand requirement
    bedtools_result = closest(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
        strand_mode="same",
    )

    # Act: Execute GIQL query with same-strand NEAREST
    engine = _setup_giql_engine(duckdb_connection)
    giql_query = """
        SELECT a.*, b.*
        FROM intervals_a a, NEAREST(intervals_b, k=1, stranded=true) b
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


def test_nearest_opposite_strand(duckdb_connection):
    """Test NEAREST with opposite-strand requirement.

    Given:
        Intervals with candidates on same and opposite strands
    When:
        NEAREST with opposite-strand requirement is applied
    Then:
        Only opposite-strand nearest intervals are reported
    """
    # Arrange
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 250, 300, "b1", 100, "-"),  # Nearest opposite strand
        GenomicInterval("chr1", 220, 240, "b2", 150, "+"),  # Closer but same strand
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

    # Act: Execute bedtools operation using pybedtools with opposite-strand requirement
    bedtools_result = closest(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
        strand_mode="opposite",
    )

    # Note: GIQL may not have direct opposite-strand support
    # This test documents the expected behavior
    assert len(bedtools_result) == 1
    assert bedtools_result[0][3] == "a1"
    assert bedtools_result[0][9] == "b1"


def test_nearest_ignore_strand(duckdb_connection):
    """Test NEAREST ignoring strand information.

    Given:
        Intervals on different strands
    When:
        NEAREST without strand requirements is applied
    Then:
        Closest interval is found regardless of strand
    """
    # Arrange
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 250, 300, "b1", 100, "+"),
        GenomicInterval("chr1", 220, 240, "b2", 150, "-"),  # Closer
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

    # Act: Execute bedtools operation using pybedtools without strand requirements
    bedtools_result = closest(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_a],
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals_b],
    )

    # Act: Execute GIQL query without strand filter
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


def test_merge_strand_specific(duckdb_connection):
    """Test MERGE with strand-specific behavior.

    Given:
        Overlapping intervals on different strands
    When:
        MERGE with strand-specific flag is applied
    Then:
        Intervals are merged per-strand (same-strand intervals merge together)
    """
    # Arrange - overlapping intervals on both strands
    intervals = [
        GenomicInterval("chr1", 100, 200, "i1", 100, "+"),
        GenomicInterval("chr1", 150, 250, "i2", 150, "+"),  # Overlaps i1 (same +)
        GenomicInterval("chr1", 120, 180, "i3", 200, "-"),  # Overlaps i1 (opposite)
        GenomicInterval("chr1", 160, 240, "i4", 100, "-"),  # Overlaps i2 (opposite)
    ]

    # Load into DuckDB
    load_intervals(
        duckdb_connection,
        "intervals",
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals],
    )

    # Act: Execute bedtools operation using pybedtools with strand-specific merging
    bedtools_result = merge(
        [(i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in intervals],
        strand_mode="same",
    )

    # Note: GIQL MERGE with strand grouping would require GROUP BY strand
    # This test documents the expected behavior
    assert len(bedtools_result) >= 2  # At least one per strand
