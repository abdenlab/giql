"""Transpilation tests for NEAREST operator SQL generation.

Tests verify that NEAREST() is correctly transpiled to SQL
(LATERAL joins for correlated queries, ORDER BY + LIMIT for standalone).
"""

import pytest
from sqlglot import parse_one

from giql import Table
from giql.dialect import GIQLDialect
from giql.generators import BaseGIQLGenerator
from giql.table import Tables


@pytest.fixture
def tables_with_peaks_and_genes():
    """Tables container with peaks and genes tables."""
    tables = Tables()
    tables.register("peaks", Table())
    tables.register("genes", Table())
    return tables


class TestNearestTranspilation:
    """Tests for NEAREST transpilation to SQL."""

    def test_nearest_basic_k3(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k=3)
        WHEN transpiling to SQL
        THEN should generate LATERAL join with DISTANCE and LIMIT 3
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=3)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.generate(ast)

        # Expectations:
        # - LATERAL subquery
        # - DISTANCE(...) AS distance in SELECT
        # - WHERE peaks.chromosome = genes.chromosome (pre-filter)
        # - ORDER BY distance
        # - LIMIT 3
        assert "LATERAL" in output.upper()
        assert "CASE" in output or "DISTANCE" in output  # Distance calculation
        assert " AS distance" in output or " as distance" in output.lower()
        assert "LIMIT 3" in output
        assert "ORDER BY" in output

    def test_nearest_with_max_distance(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k=5, max_distance=100000)
        WHEN transpiling to SQL
        THEN should generate LATERAL join with distance filter
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=5, max_distance=100000)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.generate(ast)

        # Expectations:
        # - LATERAL subquery
        # - Distance filter: <= 100000
        # - LIMIT 5
        assert "LATERAL" in output.upper()
        assert "100000" in output
        assert "LIMIT 5" in output

    def test_nearest_standalone_literal(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQL query with literal reference NEAREST(genes, reference='chr1:1000-2000', k=3)
        WHEN transpiling to SQL
        THEN should generate standalone query without LATERAL
        """
        sql = """
        SELECT *
        FROM NEAREST(genes, reference='chr1:1000-2000', k=3)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.generate(ast)

        # Expectations:
        # - No LATERAL (standalone mode)
        # - Distance calculation with literal 'chr1', 1000, 2000
        # - ORDER BY distance
        # - LIMIT 3
        assert "LATERAL" not in output.upper()
        assert "chr1" in output.lower()
        assert "LIMIT 3" in output

    def test_nearest_with_stranded(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k=3, stranded=true)
        WHEN transpiling to SQL
        THEN should generate SQL with strand filtering
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=3, stranded=true)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.generate(ast)

        # Expectations:
        # - LATERAL subquery
        # - Strand filtering in WHERE clause
        # - LIMIT 3
        assert "LATERAL" in output.upper()
        assert "strand" in output.lower()
        assert "LIMIT 3" in output

    def test_nearest_with_signed(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k=3, signed=true)
        WHEN transpiling to SQL
        THEN should generate SQL with signed distance column
            (negative for upstream, positive for downstream)
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=3, signed=true)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.generate(ast)

        # Expectations:
        # - LATERAL subquery
        # - Signed distance calculation (includes negation for upstream)
        # - LIMIT 3
        assert "LATERAL" in output.upper()
        assert "LIMIT 3" in output
        # Check for signed distance: the ELSE branch should have a negation
        # for upstream features (B before A)
        assert "ELSE -(" in output, (
            f"Expected signed distance with negation for upstream, got:\n{output}"
        )
