"""Unit tests for DISTANCE operator SQL generation and behavior.

Tests verify the distance calculation logic by checking the transpiled
SQL output and executed results.
"""

import duckdb
import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st
from sqlglot import exp
from sqlglot import parse_one

import giql  # noqa: F401  (ensures the built-in expanders are registered)
from giql.canonicalizer import canonicalize_coordinates
from giql.dialect import GIQLDialect
from giql.expander import ExpandOperators
from giql.generators import BaseGIQLGenerator
from giql.resolver import resolve_operator_refs
from giql.table import Tables
from giql.targets import GenericTarget

#: This module executes generated SQL against a real in-memory DuckDB, so every
#: test here is an integration test (the marker is registered in pyproject).
pytestmark = pytest.mark.integration


def _generate(sql: str) -> str:
    """Parse, run normalization passes 1-3, then generate SQL.

    DISTANCE operand resolution and coordinate canonicalization moved out of the
    emitter and into the ResolveOperatorRefs / CanonicalizeCoordinates passes
    (epic #114, issues #119 / #123), and DISTANCE generation itself moved onto
    the registry's AST-expansion pass (epic #137, issue #140). These behavioral
    tests must run all three passes before generating, exactly as
    :func:`giql.transpile.transpile` does, rather than calling ``generate`` on a
    bare parsed AST.
    """
    tables = Tables()
    ast = parse_one(sql, dialect=GIQLDialect)
    ast = resolve_operator_refs(ast, tables)
    ast = canonicalize_coordinates(ast)
    ast = ExpandOperators(GenericTarget(), tables).transform(ast)
    return BaseGIQLGenerator().generate(ast)


def _run(sql: str):
    """Generate, execute against in-memory DuckDB, and return the first column."""
    conn = duckdb.connect(":memory:")
    try:
        return conn.execute(_generate(sql)).fetchone()[0]
    finally:
        conn.close()


class TestDistanceCalculation:
    """Unit tests for basic distance calculation logic."""

    def test_overlapping_intervals_return_zero(self):
        """
        GIVEN two overlapping genomic intervals
        WHEN DISTANCE() is calculated between them
        THEN the distance should be 0
        """
        # Create a test query with DISTANCE()
        sql = """
        SELECT
            DISTANCE(a.interval, b.interval) as distance
        FROM
            (SELECT 'chr1' as chrom, 100 as start, 200 as end) a
        CROSS JOIN
            (SELECT 'chr1' as chrom, 150 as start, 250 as end) b
        """

        # Parse and generate SQL
        output_sql = _generate(sql)

        # Verify SQL contains CASE expression logic for overlaps
        assert "CASE" in output_sql
        assert "WHEN" in output_sql

        # Execute with DuckDB to verify behavior
        conn = duckdb.connect(":memory:")
        result = conn.execute(output_sql).fetchone()

        # Overlapping intervals should return distance = 0
        assert result[0] == 0, (
            f"Expected distance 0 for overlapping intervals, got {result[0]}"
        )

        conn.close()

    def test_non_overlapping_intervals_return_positive_distance(self):
        """
        GIVEN two non-overlapping genomic intervals with a gap
        WHEN DISTANCE() is calculated between them
        THEN the distance should be a positive integer (gap size)
        """
        # Arrange
        # Interval A: chr1:100-200, Interval B: chr1:300-400
        # Half-open gap: 300 - 200 = 100; bedtools closest -d parity adds 1 -> 101
        sql = """
        SELECT
            DISTANCE(a.interval, b.interval) as distance
        FROM
            (SELECT 'chr1' as chrom, 100 as start, 200 as end) a
        CROSS JOIN
            (SELECT 'chr1' as chrom, 300 as start, 400 as end) b
        """

        # Act
        result = _run(sql)

        # Assert
        assert result == 101, f"Expected distance 101, got {result}"

    def test_different_chromosomes_return_null(self):
        """
        GIVEN two intervals on different chromosomes
        WHEN DISTANCE() is calculated between them
        THEN the distance should be NULL
        """
        sql = """
        SELECT
            DISTANCE(a.interval, b.interval) as distance
        FROM
            (SELECT 'chr1' as chrom, 100 as start, 200 as end) a
        CROSS JOIN
            (SELECT 'chr2' as chrom, 150 as start, 250 as end) b
        """

        output_sql = _generate(sql)

        conn = duckdb.connect(":memory:")
        result = conn.execute(output_sql).fetchone()

        # Different chromosomes should return NULL
        assert result[0] is None, (
            f"Expected NULL for different chromosomes, got {result[0]}"
        )

        conn.close()

    def test_adjacent_bookended_intervals_return_one(self):
        """
        GIVEN two adjacent intervals where end_a == start_b (bookended)
        WHEN DISTANCE() is calculated between them
        THEN the distance should be 1 (bedtools closest -d reports 1 for
            book-ended features)
        """
        # Arrange
        # Interval A: chr1:100-200, Interval B: chr1:200-300 (starts where A ends)
        sql = """
        SELECT
            DISTANCE(a.interval, b.interval) as distance
        FROM
            (SELECT 'chr1' as chrom, 100 as start, 200 as end) a
        CROSS JOIN
            (SELECT 'chr1' as chrom, 200 as start, 300 as end) b
        """

        # Act
        result = _run(sql)

        # Assert
        assert result == 1, f"Expected distance 1 for bookended intervals, got {result}"

    def test_zero_width_intervals_point_features(self):
        """
        GIVEN a zero-width interval (point feature) and a regular interval
        WHEN DISTANCE() is calculated
        THEN the distance should be calculated correctly
        """
        # Arrange
        # Point feature at chr1:150 (start=150, end=150), interval at chr1:300-400
        # Half-open gap: 300 - 150 = 150; bedtools parity adds 1 -> 151
        sql = """
        SELECT
            DISTANCE(a.interval, b.interval) as distance
        FROM
            (SELECT 'chr1' as chrom, 150 as start, 150 as end) a
        CROSS JOIN
            (SELECT 'chr1' as chrom, 300 as start, 400 as end) b
        """

        # Act
        result = _run(sql)

        # Assert
        assert result == 151, f"Expected distance 151, got {result}"

    def test_one_base_gap_returns_two(self):
        """
        GIVEN two intervals separated by a 1 bp half-open gap (A chr1:100-200,
            B chr1:201-300)
        WHEN DISTANCE() is calculated between them
        THEN the distance should be 2 (half-open gap 1 + bedtools parity 1)
        """
        # Arrange
        sql = """
        SELECT DISTANCE(a.interval, b.interval) as distance
        FROM (SELECT 'chr1' as chrom, 100 as start, 200 as end) a
        CROSS JOIN (SELECT 'chr1' as chrom, 201 as start, 300 as end) b
        """

        # Act
        result = _run(sql)

        # Assert
        assert result == 2, f"Expected distance 2, got {result}"

    @pytest.mark.parametrize(
        ("b_start", "expected"),
        [(199, 0), (200, 1), (201, 2), (202, 3)],
    )
    def test_distance_increments_by_one_across_the_overlap_boundary(
        self, b_start, expected
    ):
        """
        GIVEN A chr1:100-200 and B starting at 199, 200, 201, or 202 (1 bp
            overlap, book-ended, 1 bp gap, 2 bp gap)
        WHEN DISTANCE() is calculated between them
        THEN the distance steps 0, 1, 2, 3 as the gap opens, pinning the
            overlap -> book-ended -> gap transition
        """
        # Arrange
        sql = f"""
        SELECT DISTANCE(a.interval, b.interval) as distance
        FROM (SELECT 'chr1' as chrom, 100 as start, 200 as end) a
        CROSS JOIN (SELECT 'chr1' as chrom, {b_start} as start, 300 as end) b
        """

        # Act
        result = _run(sql)

        # Assert
        assert result == expected, f"Expected distance {expected}, got {result}"

    def test_single_base_overlap_returns_zero(self):
        """
        GIVEN two intervals sharing exactly 1 bp (A chr1:100-200, B chr1:199-300)
        WHEN DISTANCE() is calculated between them
        THEN the distance should be 0, distinguishing a minimal overlap from a
            book-ended pair
        """
        # Arrange
        sql = """
        SELECT DISTANCE(a.interval, b.interval) as distance
        FROM (SELECT 'chr1' as chrom, 100 as start, 200 as end) a
        CROSS JOIN (SELECT 'chr1' as chrom, 199 as start, 300 as end) b
        """

        # Act
        result = _run(sql)

        # Assert
        assert result == 0, f"Expected distance 0, got {result}"

    def test_upstream_operand_order_returns_positive_distance(self):
        """
        GIVEN A downstream of B (A chr1:300-400, B chr1:100-200), unsigned
        WHEN DISTANCE() is calculated between them
        THEN the distance should be 101, proving the parity +1 is applied
            symmetrically on the upstream (ELSE) branch as well
        """
        # Arrange
        sql = """
        SELECT DISTANCE(a.interval, b.interval) as distance
        FROM (SELECT 'chr1' as chrom, 300 as start, 400 as end) a
        CROSS JOIN (SELECT 'chr1' as chrom, 100 as start, 200 as end) b
        """

        # Act
        result = _run(sql)

        # Assert
        assert result == 101, f"Expected distance 101, got {result}"


class TestSignedDistance:
    """Unit tests for signed (non-stranded) distance calculation."""

    def test_signed_book_ended_downstream_returns_positive_one(self):
        """
        GIVEN a book-ended pair with B downstream of A (A chr1:100-200,
            B chr1:200-300)
        WHEN DISTANCE() is calculated with signed := true
        THEN the distance should be +1, the issue's named book-ended sign case
            (sign applied after the parity +1)
        """
        # Arrange
        sql = """
        SELECT DISTANCE(a.interval, b.interval, signed := true) as distance
        FROM (SELECT 'chr1' as chrom, 100 as start, 200 as end) a
        CROSS JOIN (SELECT 'chr1' as chrom, 200 as start, 300 as end) b
        """

        # Act
        result = _run(sql)

        # Assert
        assert result == 1, f"Expected distance 1, got {result}"

    def test_signed_book_ended_upstream_returns_negative_one(self):
        """
        GIVEN a book-ended pair with B upstream of A (A chr1:200-300,
            B chr1:100-200)
        WHEN DISTANCE() is calculated with signed := true
        THEN the distance should be -1, proving the sign wraps the post-parity
            magnitude on the upstream branch
        """
        # Arrange
        sql = """
        SELECT DISTANCE(a.interval, b.interval, signed := true) as distance
        FROM (SELECT 'chr1' as chrom, 200 as start, 300 as end) a
        CROSS JOIN (SELECT 'chr1' as chrom, 100 as start, 200 as end) b
        """

        # Act
        result = _run(sql)

        # Assert
        assert result == -1, f"Expected distance -1, got {result}"


class TestStrandedDistance:
    """Tests for stranded distance calculation."""

    def test_stranded_same_strand_plus(self):
        """
        GIVEN two intervals on the same chromosome and same '+' strand
        WHEN DISTANCE() is calculated with stranded := true
        THEN the distance should be calculated normally (positive value)
        """
        # Arrange
        sql = """
        SELECT
            DISTANCE(a.interval, b.interval, stranded := true) as distance
        FROM
            (SELECT 'chr1' as chrom, 100 as start, 200 as end, '+' as strand) a
        CROSS JOIN
            (SELECT 'chr1' as chrom, 300 as start, 400 as end, '+' as strand) b
        """

        # Act
        result = _run(sql)

        # Assert
        # Gap distance should be 101 (positive, since strand is '+')
        assert result == 101, f"Expected distance 101, got {result}"

    def test_stranded_same_strand_minus(self):
        """
        GIVEN two intervals on the same chromosome and same '-' strand
        WHEN DISTANCE() is calculated with stranded := true
        THEN the distance should be negative (multiplied by -1)
        """
        # Arrange
        sql = """
        SELECT
            DISTANCE(a.interval, b.interval, stranded := true) as distance
        FROM
            (SELECT 'chr1' as chrom, 100 as start, 200 as end, '-' as strand) a
        CROSS JOIN
            (SELECT 'chr1' as chrom, 300 as start, 400 as end, '-' as strand) b
        """

        # Act
        result = _run(sql)

        # Assert
        # Gap distance should be -101 (negative, since first interval strand is '-')
        assert result == -101, f"Expected distance -101, got {result}"

    def test_stranded_different_strands_calculates_distance(self):
        """
        GIVEN two intervals on different strands ('+' and '-')
        WHEN DISTANCE() is calculated with stranded := true
        THEN the distance should be calculated normally (positive, since first interval is '+')
        """
        # Arrange
        sql = """
        SELECT
            DISTANCE(a.interval, b.interval, stranded := true) as distance
        FROM
            (SELECT 'chr1' as chrom, 100 as start, 200 as end, '+' as strand) a
        CROSS JOIN
            (SELECT 'chr1' as chrom, 300 as start, 400 as end, '-' as strand) b
        """

        # Act
        result = _run(sql)

        # Assert
        # Different strands still calculate distance, sign based on first interval
        assert result == 101, f"Expected distance 101, got {result}"

    def test_stranded_different_strands_minus_first(self):
        """
        GIVEN two intervals on different strands ('-' first, then '+')
        WHEN DISTANCE() is calculated with stranded := true
        THEN the distance should be negative (based on first interval's strand)
        """
        # Arrange
        sql = """
        SELECT
            DISTANCE(a.interval, b.interval, stranded := true) as distance
        FROM
            (SELECT 'chr1' as chrom, 100 as start, 200 as end, '-' as strand) a
        CROSS JOIN
            (SELECT 'chr1' as chrom, 300 as start, 400 as end, '+' as strand) b
        """

        # Act
        result = _run(sql)

        # Assert
        # Distance should be negative since first interval is '-'
        assert result == -101, f"Expected distance -101, got {result}"

    def test_stranded_dot_strand_returns_null(self):
        """
        GIVEN intervals with '.' strand (unspecified)
        WHEN DISTANCE() is calculated with stranded := true
        THEN the distance should be NULL
        """
        sql = """
        SELECT
            DISTANCE(a.interval, b.interval, stranded := true) as distance
        FROM
            (SELECT 'chr1' as chrom, 100 as start, 200 as end, '.' as strand) a
        CROSS JOIN
            (SELECT 'chr1' as chrom, 300 as start, 400 as end, '.' as strand) b
        """

        output_sql = _generate(sql)

        conn = duckdb.connect(":memory:")
        result = conn.execute(output_sql).fetchone()

        # '.' strand should return NULL
        assert result[0] is None, f"Expected NULL for '.' strand, got {result[0]}"

        conn.close()

    def test_stranded_question_mark_strand_returns_null(self):
        """
        GIVEN intervals with '?' strand (unknown)
        WHEN DISTANCE() is calculated with stranded := true
        THEN the distance should be NULL
        """
        sql = """
        SELECT
            DISTANCE(a.interval, b.interval, stranded := true) as distance
        FROM
            (SELECT 'chr1' as chrom, 100 as start, 200 as end, '?' as strand) a
        CROSS JOIN
            (SELECT 'chr1' as chrom, 300 as start, 400 as end, '+' as strand) b
        """

        output_sql = _generate(sql)

        conn = duckdb.connect(":memory:")
        result = conn.execute(output_sql).fetchone()

        # '?' strand should return NULL
        assert result[0] is None, f"Expected NULL for '?' strand, got {result[0]}"

        conn.close()

    def test_stranded_null_strand_returns_null(self):
        """
        GIVEN intervals with NULL strand
        WHEN DISTANCE() is calculated with stranded := true
        THEN the distance should be NULL
        """
        sql = """
        SELECT
            DISTANCE(a.interval, b.interval, stranded := true) as distance
        FROM
            (SELECT 'chr1' as chrom, 100 as start, 200 as end, NULL as strand) a
        CROSS JOIN
            (SELECT 'chr1' as chrom, 300 as start, 400 as end, '+' as strand) b
        """

        output_sql = _generate(sql)

        conn = duckdb.connect(":memory:")
        result = conn.execute(output_sql).fetchone()

        # NULL strand should return NULL
        assert result[0] is None, f"Expected NULL for NULL strand, got {result[0]}"

        conn.close()

    def test_stranded_overlapping_intervals_minus_strand(self):
        """
        GIVEN two overlapping intervals on '-' strand
        WHEN DISTANCE() is calculated with stranded := true
        THEN the distance should be 0 (overlaps have distance 0 regardless of strand)
        """
        sql = """
        SELECT
            DISTANCE(a.interval, b.interval, stranded := true) as distance
        FROM
            (SELECT 'chr1' as chrom, 100 as start, 200 as end, '-' as strand) a
        CROSS JOIN
            (SELECT 'chr1' as chrom, 150 as start, 250 as end, '-' as strand) b
        """

        output_sql = _generate(sql)

        conn = duckdb.connect(":memory:")
        result = conn.execute(output_sql).fetchone()

        # Overlapping intervals should return 0
        assert result[0] == 0, (
            f"Expected distance 0 for overlapping intervals, got {result[0]}"
        )

        conn.close()

    def test_stranded_book_ended_plus_strand_returns_one(self):
        """
        GIVEN a book-ended pair both on '+' strand (A chr1:100-200,
            B chr1:200-300)
        WHEN DISTANCE() is calculated with stranded := true
        THEN the distance should be 1 (parity +1 reaches the stranded branch;
            '+' applies no flip)
        """
        # Arrange
        sql = """
        SELECT DISTANCE(a.interval, b.interval, stranded := true) as distance
        FROM (SELECT 'chr1' as chrom, 100 as start, 200 as end, '+' as strand) a
        CROSS JOIN (SELECT 'chr1' as chrom, 200 as start, 300 as end, '+' as strand) b
        """

        # Act
        result = _run(sql)

        # Assert
        assert result == 1, f"Expected distance 1, got {result}"

    def test_stranded_book_ended_minus_strand_returns_negative_one(self):
        """
        GIVEN a book-ended pair both on '-' strand (A chr1:100-200,
            B chr1:200-300)
        WHEN DISTANCE() is calculated with stranded := true
        THEN the distance should be -1, the '-' strand flipping the sign of the
            post-parity magnitude
        """
        # Arrange
        sql = """
        SELECT DISTANCE(a.interval, b.interval, stranded := true) as distance
        FROM (SELECT 'chr1' as chrom, 100 as start, 200 as end, '-' as strand) a
        CROSS JOIN (SELECT 'chr1' as chrom, 200 as start, 300 as end, '-' as strand) b
        """

        # Act
        result = _run(sql)

        # Assert
        assert result == -1, f"Expected distance -1, got {result}"

    def test_stranded_upstream_plus_strand_returns_positive(self):
        """
        GIVEN B upstream of A both on '+' strand (A chr1:300-400, B chr1:100-200)
        WHEN DISTANCE() is calculated with stranded := true
        THEN the distance should be 101, exercising the stranded upstream (ELSE)
            branch that the existing stranded tests never reach
        """
        # Arrange
        sql = """
        SELECT DISTANCE(a.interval, b.interval, stranded := true) as distance
        FROM (SELECT 'chr1' as chrom, 300 as start, 400 as end, '+' as strand) a
        CROSS JOIN (SELECT 'chr1' as chrom, 100 as start, 200 as end, '+' as strand) b
        """

        # Act
        result = _run(sql)

        # Assert
        assert result == 101, f"Expected distance 101, got {result}"

    def test_stranded_upstream_minus_strand_returns_negative(self):
        """
        GIVEN B upstream of A both on '-' strand (A chr1:300-400, B chr1:100-200)
        WHEN DISTANCE() is calculated with stranded := true
        THEN the distance should be -101, the '-' flip applied to the upstream
            branch's post-parity magnitude
        """
        # Arrange
        sql = """
        SELECT DISTANCE(a.interval, b.interval, stranded := true) as distance
        FROM (SELECT 'chr1' as chrom, 300 as start, 400 as end, '-' as strand) a
        CROSS JOIN (SELECT 'chr1' as chrom, 100 as start, 200 as end, '-' as strand) b
        """

        # Act
        result = _run(sql)

        # Assert
        assert result == -101, f"Expected distance -101, got {result}"

    def test_stranded_question_strand_returns_null(self):
        """
        GIVEN one interval with '?' strand (A chr1:100-200 '?', B chr1:300-400 '+')
        WHEN DISTANCE() is calculated with stranded := true
        THEN the distance should be NULL, confirming the parity +1 did not
            perturb the '?' NULL short-circuit
        """
        # Arrange
        sql = """
        SELECT DISTANCE(a.interval, b.interval, stranded := true) as distance
        FROM (SELECT 'chr1' as chrom, 100 as start, 200 as end, '?' as strand) a
        CROSS JOIN (SELECT 'chr1' as chrom, 300 as start, 400 as end, '+' as strand) b
        """

        # Act
        result = _run(sql)

        # Assert
        assert result is None, f"Expected NULL for '?' strand, got {result}"


class TestStrandedSignedDistance:
    """Unit tests for the combined stranded + signed distance calculation."""

    def test_stranded_signed_upstream_plus_strand_returns_negative(self):
        """
        GIVEN B upstream of A both on '+' strand (A chr1:300-400, B chr1:100-200)
        WHEN DISTANCE() is calculated with stranded := true, signed := true
        THEN the distance should be -101, the upstream-branch sign that DIVERGES
            from the stranded-only variant ('+' -> -(magnitude+1))
        """
        # Arrange
        sql = """
        SELECT DISTANCE(a.interval, b.interval, stranded := true, signed := true)
            as distance
        FROM (SELECT 'chr1' as chrom, 300 as start, 400 as end, '+' as strand) a
        CROSS JOIN (SELECT 'chr1' as chrom, 100 as start, 200 as end, '+' as strand) b
        """

        # Act
        result = _run(sql)

        # Assert
        assert result == -101, f"Expected distance -101, got {result}"

    def test_stranded_signed_upstream_minus_strand_returns_positive(self):
        """
        GIVEN B upstream of A both on '-' strand (A chr1:300-400, B chr1:100-200)
        WHEN DISTANCE() is calculated with stranded := true, signed := true
        THEN the distance should be +101, the only ELSE-branch case where the
            sign is +(magnitude+1) ('-' strand double-flip)
        """
        # Arrange
        sql = """
        SELECT DISTANCE(a.interval, b.interval, stranded := true, signed := true)
            as distance
        FROM (SELECT 'chr1' as chrom, 300 as start, 400 as end, '-' as strand) a
        CROSS JOIN (SELECT 'chr1' as chrom, 100 as start, 200 as end, '-' as strand) b
        """

        # Act
        result = _run(sql)

        # Assert
        assert result == 101, f"Expected distance 101, got {result}"

    def test_stranded_signed_book_ended_minus_strand_returns_negative_one(self):
        """
        GIVEN a book-ended pair both on '-' strand, B downstream (A chr1:100-200,
            B chr1:200-300)
        WHEN DISTANCE() is calculated with stranded := true, signed := true
        THEN the distance should be -1, pinning the sign at the smallest
            (book-ended) post-parity magnitude
        """
        # Arrange
        sql = """
        SELECT DISTANCE(a.interval, b.interval, stranded := true, signed := true)
            as distance
        FROM (SELECT 'chr1' as chrom, 100 as start, 200 as end, '-' as strand) a
        CROSS JOIN (SELECT 'chr1' as chrom, 200 as start, 300 as end, '-' as strand) b
        """

        # Act
        result = _run(sql)

        # Assert
        assert result == -1, f"Expected distance -1, got {result}"

    def test_stranded_signed_null_strand_returns_null(self):
        """
        GIVEN one interval with '.' strand (A chr1:100-200 '.', B chr1:200-300 '+')
        WHEN DISTANCE() is calculated with stranded := true, signed := true
        THEN the distance should be NULL, confirming the parity +1 left the NULL
            short-circuit intact in the combined variant
        """
        # Arrange
        sql = """
        SELECT DISTANCE(a.interval, b.interval, stranded := true, signed := true)
            as distance
        FROM (SELECT 'chr1' as chrom, 100 as start, 200 as end, '.' as strand) a
        CROSS JOIN (SELECT 'chr1' as chrom, 200 as start, 300 as end, '+' as strand) b
        """

        # Act
        result = _run(sql)

        # Assert
        assert result is None, f"Expected NULL for '.' strand, got {result}"


# --- Drift guard: expand_distance vs the legacy _generate_distance_case -------
#
# DISTANCE moved onto the AST-expansion pass (expand_distance), but
# BaseGIQLGenerator._generate_distance_case is retained because NEAREST still
# calls it. The two compute the same distance by different routes; these tests
# pin that they stay semantically equivalent until NEAREST migrates and the
# legacy method can be deleted.

#: Column expressions both routes are evaluated over. ``a``/``b`` are the two
#: operand relations supplied by the parity harness's VALUES row.
_CHROM_A, _START_A, _END_A, _STRAND_A = 'a."chrom"', 'a."start"', 'a."end"', 'a."strand"'
_CHROM_B, _START_B, _END_B, _STRAND_B = 'b."chrom"', 'b."start"', 'b."end"', 'b."strand"'

#: The four DISTANCE shapes: (id, stranded, signed).
_SHAPES = [
    ("unsigned_nonstranded", False, False),
    ("signed_nonstranded", False, True),
    ("unsigned_stranded", True, False),
    ("signed_stranded", True, True),
]

#: Rows exercised by the parity test: ordinary downstream/upstream gaps,
#: book-ended pairs, overlaps, a chrom mismatch, and strand-invalid ('.'/'?'/
#: NULL) rows so every WHEN branch of both routes is covered.
_PARITY_ROWS = [
    ("chr1", 100, 200, "+", "chr1", 300, 400, "+"),  # downstream gap
    ("chr1", 300, 400, "+", "chr1", 100, 200, "+"),  # upstream gap
    ("chr1", 100, 200, "+", "chr1", 200, 300, "+"),  # book-ended downstream
    ("chr1", 200, 300, "+", "chr1", 100, 200, "+"),  # book-ended upstream
    ("chr1", 100, 200, "+", "chr1", 150, 250, "+"),  # overlap
    ("chr1", 100, 200, "-", "chr1", 300, 400, "+"),  # '-' strand A, downstream
    ("chr1", 300, 400, "-", "chr1", 100, 200, "+"),  # '-' strand A, upstream
    ("chr1", 100, 200, "+", "chr2", 300, 400, "+"),  # chrom mismatch -> NULL
    ("chr1", 100, 200, ".", "chr1", 300, 400, "+"),  # strand A '.' -> NULL
    ("chr1", 100, 200, "?", "chr1", 300, 400, "+"),  # strand A '?' -> NULL
    ("chr1", 100, 200, "+", "chr1", 300, 400, "."),  # strand B '.' -> NULL
    ("chr1", 100, 200, None, "chr1", 300, 400, "+"),  # strand A NULL -> NULL
]


def _row_values_cte(row) -> str:
    """Render one parity row as ``a``/``b`` relations via SELECT subqueries.

    *row* is ``(chrom_a, start_a, end_a, strand_a, chrom_b, ...)``; strands may
    be ``None`` (rendered as SQL ``NULL``).
    """
    ca, sa, ea, ta, cb, sb, eb, tb = row

    def _strand(value):
        return "NULL" if value is None else f"'{value}'"

    a = (
        f"(SELECT '{ca}' AS chrom, {sa} AS start, {ea} AS \"end\", "
        f"{_strand(ta)} AS strand) a"
    )
    b = (
        f"(SELECT '{cb}' AS chrom, {sb} AS start, {eb} AS \"end\", "
        f"{_strand(tb)} AS strand) b"
    )
    return f"{a} CROSS JOIN {b}"


def _expander_distance_case(stranded: bool, signed: bool) -> str:
    """Return the DISTANCE CASE the expander builds, isolated from its SELECT.

    Runs the real transpile passes over a DISTANCE query whose operands resolve
    to ``a``/``b`` default columns, then lifts the generated CASE expression so
    it can be re-embedded over arbitrary VALUES rows.
    """
    args = ["a.interval", "b.interval"]
    if stranded:
        args.append("stranded := true")
    if signed:
        args.append("signed := true")
    sql = (
        f"SELECT DISTANCE({', '.join(args)}) AS d "
        "FROM (SELECT 'x' AS chrom, 0 AS start, 0 AS \"end\", '+' AS strand) a "
        "CROSS JOIN (SELECT 'x' AS chrom, 0 AS start, 0 AS \"end\", '+' AS strand) b"
    )
    tables = Tables()
    ast = parse_one(sql, dialect=GIQLDialect)
    ast = resolve_operator_refs(ast, tables)
    ast = canonicalize_coordinates(ast)
    ast = ExpandOperators(GenericTarget(), tables).transform(ast)
    # The single projected expression is the expander's CASE.
    return ast.find(exp.Select).expressions[0].this.sql()


def _legacy_distance_case(stranded: bool, signed: bool) -> str:
    """Return the CASE the legacy _generate_distance_case builds for ``a``/``b``."""
    return BaseGIQLGenerator()._generate_distance_case(
        _CHROM_A,
        _START_A,
        _END_A,
        _STRAND_A if stranded else None,
        _CHROM_B,
        _START_B,
        _END_B,
        _STRAND_B if stranded else None,
        stranded=stranded,
        signed=signed,
    )


def _eval_case(conn, case_sql: str, row) -> object:
    """Execute one distance CASE over one parity row, returning the scalar.

    Reuses the caller-supplied *conn* (a module-scoped in-memory DuckDB) rather
    than opening a fresh connection per row — every row is a standalone
    ``SELECT ... FROM (VALUES)`` against no persistent state, so one connection
    serves the whole parity sweep.
    """
    query = f"SELECT {case_sql} AS d FROM {_row_values_cte(row)}"
    return conn.execute(query).fetchone()[0]


@pytest.fixture(scope="module")
def parity_conn():
    """A module-scoped in-memory DuckDB connection for the parity/property sweep.

    Opened once and shared across every parity row and Hypothesis example
    instead of reconnecting per row; each evaluation is a self-contained
    ``SELECT`` over inline ``VALUES``, so no per-row isolation is needed.
    """
    conn = duckdb.connect(":memory:")
    try:
        yield conn
    finally:
        conn.close()


class TestDistanceExpanderLegacyParity:
    """expand_distance and the retained _generate_distance_case agree row-for-row."""

    @pytest.mark.parametrize(
        "shape_id, stranded, signed", _SHAPES, ids=[s[0] for s in _SHAPES]
    )
    def test_expander_matches_legacy_distance_case(
        self, parity_conn, shape_id, stranded, signed
    ):
        """
        GIVEN the four DISTANCE shapes (unsigned/signed x non-stranded/stranded)
            plus overlap, chrom-mismatch, and strand-invalid input rows
        WHEN the same inputs run through expand_distance and the retained
            _generate_distance_case
        THEN both routes return the identical scalar for every row, pinning the
            two distance implementations against drift until NEAREST migrates.
        """
        # Arrange
        expander_case = _expander_distance_case(stranded, signed)
        legacy_case = _legacy_distance_case(stranded, signed)

        # Act & assert
        for row in _PARITY_ROWS:
            expander_result = _eval_case(parity_conn, expander_case, row)
            legacy_result = _eval_case(parity_conn, legacy_case, row)
            assert expander_result == legacy_result, (
                f"{shape_id}: expander {expander_result!r} != "
                f"legacy {legacy_result!r} for row {row}"
            )


class TestDistanceExpanderProperties:
    """Property-based invariants of the expander's distance CASE."""

    @settings(max_examples=200, deadline=None)
    @given(
        start_a=st.integers(min_value=0, max_value=10_000),
        len_a=st.integers(min_value=1, max_value=5_000),
        start_b=st.integers(min_value=0, max_value=10_000),
        len_b=st.integers(min_value=1, max_value=5_000),
    )
    def test_distance_invariants_hold(
        self, parity_conn, start_a, len_a, start_b, len_b
    ):
        """
        GIVEN random A and B intervals (start + positive length)
        WHEN DISTANCE is evaluated unsigned, signed, and cross-chromosome
        THEN unsigned == abs(signed), overlapping intervals report 0, a
            non-overlapping same-chrom pair reports the half-open gap + 1
            (bedtools parity ground truth), and a cross-chromosome pair reports
            NULL.
        """
        # Arrange
        end_a = start_a + len_a
        end_b = start_b + len_b
        same_chrom = ("chr1", start_a, end_a, "+", "chr1", start_b, end_b, "+")
        cross_chrom = ("chr1", start_a, end_a, "+", "chr2", start_b, end_b, "+")
        unsigned_case = _expander_distance_case(stranded=False, signed=False)
        signed_case = _expander_distance_case(stranded=False, signed=True)

        # Act
        unsigned = _eval_case(parity_conn, unsigned_case, same_chrom)
        signed = _eval_case(parity_conn, signed_case, same_chrom)
        cross = _eval_case(parity_conn, unsigned_case, cross_chrom)

        # Assert
        assert unsigned == abs(signed)
        overlaps = start_a < end_b and end_a > start_b
        if overlaps:
            assert unsigned == 0
        else:
            # Ground truth (bedtools closest -d): the unsigned distance of two
            # non-overlapping same-chrom intervals is the raw half-open gap plus
            # one. This ties the property test to the +1 offset directly, so
            # dropping the +1 from the expander would fail here (not only in the
            # legacy-parity test). The gap is end_a..start_b downstream or
            # end_b..start_a upstream, whichever is non-negative.
            gap_in_bases = max(start_b - end_a, start_a - end_b)
            assert unsigned == gap_in_bases + 1
        assert cross is None
