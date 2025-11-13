"""Transpilation tests for NEAREST operator SQL generation.

Tests verify that NEAREST() is correctly transpiled to dialect-specific SQL
(LATERAL joins for PostgreSQL/DuckDB, window functions for SQLite).
"""

import pytest
from sqlglot import parse_one

from giql.dialect import GIQLDialect
from giql.generators import BaseGIQLGenerator
from giql.generators import GIQLDuckDBGenerator
from giql.schema import ColumnInfo
from giql.schema import SchemaInfo
from giql.schema import TableSchema


@pytest.fixture
def schema_with_peaks_and_genes():
    """Schema info with peaks and genes tables."""
    schema = SchemaInfo()

    # Register peaks table
    peaks_table = TableSchema(name="peaks", columns={})
    peaks_table.columns["peak_id"] = ColumnInfo(name="peak_id", type="INTEGER")
    peaks_table.columns["position"] = ColumnInfo(
        name="position",
        type="VARCHAR",
        is_genomic=True,
        chrom_col="chromosome",
        start_col="start_pos",
        end_col="end_pos",
        strand_col="strand",
    )
    schema.tables["peaks"] = peaks_table

    # Register genes table
    genes_table = TableSchema(name="genes", columns={})
    genes_table.columns["gene_id"] = ColumnInfo(name="gene_id", type="INTEGER")
    genes_table.columns["name"] = ColumnInfo(name="name", type="VARCHAR")
    genes_table.columns["position"] = ColumnInfo(
        name="position",
        type="VARCHAR",
        is_genomic=True,
        chrom_col="chromosome",
        start_col="start_pos",
        end_col="end_pos",
        strand_col="strand",
    )
    schema.tables["genes"] = genes_table

    return schema


class TestNearestTranspilationDuckDB:
    """Tests for NEAREST transpilation to DuckDB SQL (LATERAL joins)."""

    def test_nearest_basic_k3_duckdb(self, schema_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k=3)
        WHEN transpiling to DuckDB SQL
        THEN should generate LATERAL join with DISTANCE and LIMIT 3
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=3)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = GIQLDuckDBGenerator(schema_info=schema_with_peaks_and_genes)
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

    def test_nearest_with_max_distance_duckdb(self, schema_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k=5, max_distance=100000)
        WHEN transpiling to DuckDB SQL
        THEN should generate LATERAL join with distance filter
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=5, max_distance=100000)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = GIQLDuckDBGenerator(schema_info=schema_with_peaks_and_genes)
        output = generator.generate(ast)

        # Expectations:
        # - LATERAL subquery
        # - Distance filter: <= 100000
        # - LIMIT 5
        assert "LATERAL" in output.upper()
        assert "100000" in output
        assert "LIMIT 5" in output

    def test_nearest_standalone_literal_duckdb(self, schema_with_peaks_and_genes):
        """
        GIVEN a GIQL query with literal reference NEAREST(genes, reference='chr1:1000-2000', k=3)
        WHEN transpiling to DuckDB SQL
        THEN should generate standalone query without LATERAL
        """
        sql = """
        SELECT *
        FROM NEAREST(genes, reference='chr1:1000-2000', k=3)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = GIQLDuckDBGenerator(schema_info=schema_with_peaks_and_genes)
        output = generator.generate(ast)

        # Expectations:
        # - No LATERAL (standalone mode)
        # - Distance calculation with literal 'chr1', 1000, 2000
        # - ORDER BY distance
        # - LIMIT 3
        assert "LATERAL" not in output.upper()
        assert "chr1" in output.lower()
        assert "LIMIT 3" in output

    def test_nearest_with_stranded_duckdb(self, schema_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k=3, stranded=true)
        WHEN transpiling to DuckDB SQL
        THEN should generate SQL with strand filtering
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=3, stranded=true)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = GIQLDuckDBGenerator(schema_info=schema_with_peaks_and_genes)
        output = generator.generate(ast)

        # Expectations:
        # - LATERAL subquery
        # - Strand filtering in WHERE clause
        # - LIMIT 3
        assert "LATERAL" in output.upper()
        assert "strand" in output.lower()
        assert "LIMIT 3" in output

    def test_nearest_with_signed_duckdb(self, schema_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k=3, signed=true)
        WHEN transpiling to DuckDB SQL
        THEN should generate SQL with signed distance calculation
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=3, signed=true)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = GIQLDuckDBGenerator(schema_info=schema_with_peaks_and_genes)
        output = generator.generate(ast)

        # Expectations:
        # - LATERAL subquery
        # - Signed distance calculation (negative for upstream)
        # - LIMIT 3
        assert "LATERAL" in output.upper()
        # Check for negative distance calculation in ELSE clause
        # The signed version should have -( in the ELSE branch
        assert "-(" in output or "- (" in output
        assert "LIMIT 3" in output


# PostgreSQL uses same generator as base for now
# class TestNearestTranspilationPostgreSQL:
#     """Tests for NEAREST transpilation to PostgreSQL SQL (LATERAL joins)."""
#     (Skipped - uses BaseGIQLGenerator for now)


class TestNearestTranspilationSQLite:
    """Tests for NEAREST transpilation to SQLite SQL (using LATERAL for MVP)."""

    def test_nearest_basic_k3_sqlite(self, schema_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k=3)
        WHEN transpiling to SQLite SQL
        THEN should generate LATERAL subquery with ORDER BY and LIMIT
        (Note: Using LATERAL for MVP - window function optimization to be added later)
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=3)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(schema_info=schema_with_peaks_and_genes)
        output = generator.generate(ast)

        # MVP expectations (LATERAL syntax):
        # - LATERAL subquery
        # - Distance calculation (CASE WHEN)
        # - ORDER BY distance
        # - LIMIT 3
        assert "LATERAL" in output.upper()
        assert "CASE" in output.upper()
        assert " AS distance" in output or " AS DISTANCE" in output
        assert "ORDER BY" in output.upper()
        assert "LIMIT 3" in output

    def test_nearest_with_max_distance_sqlite(self, schema_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k=5, max_distance=100000)
        WHEN transpiling to SQLite SQL
        THEN should generate LATERAL with distance filter
        (Note: Using LATERAL for MVP - window function optimization to be added later)
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=5, max_distance=100000)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(schema_info=schema_with_peaks_and_genes)
        output = generator.generate(ast)

        # MVP expectations (LATERAL syntax):
        # - LATERAL subquery
        # - Distance filter: <= 100000
        # - LIMIT 5
        assert "LATERAL" in output.upper()
        assert "100000" in output
        assert "LIMIT 5" in output

    def test_nearest_standalone_literal_sqlite(self, schema_with_peaks_and_genes):
        """
        GIVEN a GIQL query with literal reference NEAREST(genes, reference='chr1:1000-2000', k=3)
        WHEN transpiling to SQLite SQL
        THEN should generate standalone query without window functions
        """
        sql = """
        SELECT *
        FROM NEAREST(genes, reference='chr1:1000-2000', k=3)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(schema_info=schema_with_peaks_and_genes)
        output = generator.generate(ast)

        # Expectations:
        # - No CTE needed (standalone mode)
        # - Distance calculation with literal 'chr1', 1000, 2000
        # - ORDER BY distance
        # - LIMIT 3
        assert "chr1" in output.lower()
        assert "ORDER BY" in output.upper()
        assert "LIMIT 3" in output

    def test_nearest_with_stranded_sqlite(self, schema_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k=3, stranded=true)
        WHEN transpiling to SQLite SQL
        THEN should generate SQL with strand filtering
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=3, stranded=true)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(schema_info=schema_with_peaks_and_genes)
        output = generator.generate(ast)

        # Expectations:
        # - LATERAL subquery
        # - Strand filtering in WHERE clause
        # - LIMIT 3
        assert "LATERAL" in output.upper()
        assert "strand" in output.lower()
        assert "LIMIT 3" in output

    def test_nearest_with_signed_sqlite(self, schema_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k=3, signed=true)
        WHEN transpiling to SQLite SQL
        THEN should generate SQL with signed distance calculation
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=3, signed=true)
        """

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator(schema_info=schema_with_peaks_and_genes)
        output = generator.generate(ast)

        # Expectations:
        # - LATERAL subquery
        # - Signed distance calculation (negative for upstream)
        # - LIMIT 3
        assert "LATERAL" in output.upper()
        # Check for negative distance calculation in ELSE clause
        assert "-(" in output or "- (" in output
        assert "LIMIT 3" in output
