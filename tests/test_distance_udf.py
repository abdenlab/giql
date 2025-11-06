"""Unit tests for DISTANCE operator SQL generation and behavior.

Tests verify the distance calculation logic by checking the transpiled
SQL output and executed results.
"""

import duckdb
import pytest
from sqlglot import parse_one

from giql.dialect import GIQLDialect
from giql.generators import BaseGIQLGenerator


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
            DISTANCE(a.position, b.position) as distance
        FROM
            (SELECT 'chr1' as chromosome, 100 as start_pos, 200 as end_pos) a
        CROSS JOIN
            (SELECT 'chr1' as chromosome, 150 as start_pos, 250 as end_pos) b
        """

        # Parse and generate SQL
        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator()
        output_sql = generator.generate(ast)

        # Verify SQL contains CASE expression logic for overlaps
        assert "CASE" in output_sql
        assert "WHEN" in output_sql

        # Execute with DuckDB to verify behavior
        conn = duckdb.connect(":memory:")
        result = conn.execute(output_sql).fetchone()

        # Overlapping intervals should return distance = 0
        assert result[0] == 0, f"Expected distance 0 for overlapping intervals, got {result[0]}"

        conn.close()

    def test_non_overlapping_intervals_return_positive_distance(self):
        """
        GIVEN two non-overlapping genomic intervals with a gap
        WHEN DISTANCE() is calculated between them
        THEN the distance should be a positive integer (gap size)
        """
        # Interval A: chr1:100-200
        # Interval B: chr1:300-400
        # Gap: 300 - 200 = 100 base pairs
        sql = """
        SELECT
            DISTANCE(a.position, b.position) as distance
        FROM
            (SELECT 'chr1' as chromosome, 100 as start_pos, 200 as end_pos) a
        CROSS JOIN
            (SELECT 'chr1' as chromosome, 300 as start_pos, 400 as end_pos) b
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator()
        output_sql = generator.generate(ast)

        conn = duckdb.connect(":memory:")
        result = conn.execute(output_sql).fetchone()

        # Gap distance should be 100
        assert result[0] == 100, f"Expected distance 100, got {result[0]}"

        conn.close()

    def test_different_chromosomes_return_null(self):
        """
        GIVEN two intervals on different chromosomes
        WHEN DISTANCE() is calculated between them
        THEN the distance should be NULL
        """
        sql = """
        SELECT
            DISTANCE(a.position, b.position) as distance
        FROM
            (SELECT 'chr1' as chromosome, 100 as start_pos, 200 as end_pos) a
        CROSS JOIN
            (SELECT 'chr2' as chromosome, 150 as start_pos, 250 as end_pos) b
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator()
        output_sql = generator.generate(ast)

        conn = duckdb.connect(":memory:")
        result = conn.execute(output_sql).fetchone()

        # Different chromosomes should return NULL
        assert result[0] is None, f"Expected NULL for different chromosomes, got {result[0]}"

        conn.close()

    def test_adjacent_bookended_intervals_return_zero(self):
        """
        GIVEN two adjacent intervals where end_a == start_b (bookended)
        WHEN DISTANCE() is calculated between them
        THEN the distance should be 0 (following bedtools convention)
        """
        # Interval A: chr1:100-200
        # Interval B: chr1:200-300 (starts exactly where A ends)
        sql = """
        SELECT
            DISTANCE(a.position, b.position) as distance
        FROM
            (SELECT 'chr1' as chromosome, 100 as start_pos, 200 as end_pos) a
        CROSS JOIN
            (SELECT 'chr1' as chromosome, 200 as start_pos, 300 as end_pos) b
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator()
        output_sql = generator.generate(ast)

        conn = duckdb.connect(":memory:")
        result = conn.execute(output_sql).fetchone()

        # Bookended intervals should return distance = 0
        assert result[0] == 0, f"Expected distance 0 for bookended intervals, got {result[0]}"

        conn.close()

    def test_zero_width_intervals_point_features(self):
        """
        GIVEN a zero-width interval (point feature) and a regular interval
        WHEN DISTANCE() is calculated
        THEN the distance should be calculated correctly
        """
        # Point feature at chr1:150 (start=150, end=150)
        # Interval at chr1:300-400
        # Distance: 300 - 150 = 150
        sql = """
        SELECT
            DISTANCE(a.position, b.position) as distance
        FROM
            (SELECT 'chr1' as chromosome, 150 as start_pos, 150 as end_pos) a
        CROSS JOIN
            (SELECT 'chr1' as chromosome, 300 as start_pos, 400 as end_pos) b
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator()
        output_sql = generator.generate(ast)

        conn = duckdb.connect(":memory:")
        result = conn.execute(output_sql).fetchone()

        # Distance should be 150
        assert result[0] == 150, f"Expected distance 150, got {result[0]}"

        conn.close()
