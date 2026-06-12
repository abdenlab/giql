"""Unit tests for DISTANCE operator SQL generation and behavior.

Tests verify the distance calculation logic by checking the transpiled
SQL output and executed results.
"""

import duckdb
import pytest
from sqlglot import parse_one

from giql.canonicalizer import canonicalize_coordinates
from giql.dialect import GIQLDialect
from giql.generators import BaseGIQLGenerator
from giql.resolver import resolve_operator_refs
from giql.table import Tables


def _generate(sql: str) -> str:
    """Parse, run normalization passes 1 and 2, then generate SQL.

    DISTANCE operand resolution and coordinate canonicalization moved out of the
    emitter and into the ResolveOperatorRefs / CanonicalizeCoordinates passes
    (epic #114, issues #119 / #123). These behavioral tests must run both passes
    before generating, exactly as :func:`giql.transpile.transpile` does, rather
    than calling ``generate`` on a bare parsed AST.
    """
    ast = parse_one(sql, dialect=GIQLDialect)
    ast = resolve_operator_refs(ast, Tables())
    ast = canonicalize_coordinates(ast)
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
