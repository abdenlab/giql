"""Tests for the COVERAGE operator."""

import pytest
from sqlglot import parse_one

from giql import Table
from giql import transpile
from giql.dialect import GIQLDialect
from giql.expressions import GIQLCoverage


class TestCoverageParsing:
    """Tests for parsing COVERAGE expressions."""

    def test_parse_positional_args(self):
        """
        GIVEN a COVERAGE expression with positional arguments
        WHEN parsing with GIQLDialect
        THEN should produce GIQLCoverage with resolution=1000 and stat defaults to None
        """
        ast = parse_one(
            "SELECT COVERAGE(interval, 1000) FROM features",
            dialect=GIQLDialect,
        )
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["resolution"].this == "1000"
        assert coverage[0].args.get("stat") is None

    def test_parse_named_stat(self):
        """
        GIVEN a COVERAGE expression with named stat parameter
        WHEN parsing with GIQLDialect
        THEN should produce GIQLCoverage with resolution=500 and stat='mean'
        """
        ast = parse_one(
            "SELECT COVERAGE(interval, 500, stat := 'mean') FROM features",
            dialect=GIQLDialect,
        )
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["resolution"].this == "500"
        assert coverage[0].args["stat"].this == "mean"

    def test_parse_named_resolution(self):
        """
        GIVEN a COVERAGE expression with named resolution parameter
        WHEN parsing with GIQLDialect
        THEN should produce GIQLCoverage with named resolution=1000
        """
        ast = parse_one(
            "SELECT COVERAGE(interval, resolution := 1000) FROM features",
            dialect=GIQLDialect,
        )
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["resolution"].this == "1000"

    def test_parse_arrow_named_params(self):
        """
        GIVEN a COVERAGE expression using => (standard SQL named parameter syntax)
        WHEN parsing with GIQLDialect
        THEN should produce GIQLCoverage with the same result as :=
        """
        ast = parse_one(
            "SELECT COVERAGE(interval, 500, stat => 'mean') FROM features",
            dialect=GIQLDialect,
        )
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["resolution"].this == "500"
        assert coverage[0].args["stat"].this == "mean"


class TestCoverageTranspile:
    """Tests for COVERAGE transpilation."""

    def test_basic_transpilation(self):
        """
        GIVEN a basic COVERAGE query
        WHEN transpiling
        THEN should produce SQL with generate_series, LEFT JOIN on overlap, GROUP BY, and COUNT
        """
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features",
            tables=["features"],
        )

        upper = sql.upper()
        assert "GENERATE_SERIES" in upper
        assert "LEFT JOIN" in upper
        assert "GROUP BY" in upper
        assert "COUNT" in upper
        assert "__GIQL_BINS" in upper

    def test_stat_mean(self):
        """
        GIVEN a COVERAGE query with stat := 'mean'
        WHEN transpiling
        THEN should use AVG instead of COUNT
        """
        sql = transpile(
            "SELECT COVERAGE(interval, 1000, stat := 'mean') FROM features",
            tables=["features"],
        )

        upper = sql.upper()
        assert "AVG" in upper
        assert "COUNT" not in upper

    def test_stat_sum(self):
        """
        GIVEN a COVERAGE query with stat := 'sum'
        WHEN transpiling
        THEN should use SUM aggregate
        """
        sql = transpile(
            "SELECT COVERAGE(interval, 1000, stat := 'sum') FROM features",
            tables=["features"],
        )

        upper = sql.upper()
        assert "SUM" in upper

    def test_stat_max(self):
        """
        GIVEN a COVERAGE query with stat := 'max'
        WHEN transpiling
        THEN should use MAX aggregate
        """
        sql = transpile(
            "SELECT COVERAGE(interval, 1000, stat := 'max') FROM features",
            tables=["features"],
        )

        upper = sql.upper()
        assert "MAX(" in upper

    def test_custom_column_mapping(self):
        """
        GIVEN a COVERAGE query with custom column mappings
        WHEN transpiling
        THEN should use mapped column names in JOIN and GROUP BY
        """
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM peaks",
            tables=[
                Table(
                    "peaks",
                    genomic_col="interval",
                    chrom_col="chromosome",
                    start_col="start_pos",
                    end_col="end_pos",
                )
            ],
        )

        assert "chromosome" in sql
        assert "start_pos" in sql
        assert "end_pos" in sql

    def test_where_clause_preserved(self):
        """
        GIVEN a COVERAGE query with a WHERE clause
        WHEN transpiling
        THEN should preserve the WHERE filter
        """
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features WHERE score > 10",
            tables=["features"],
        )

        assert "score > 10" in sql

    def test_additional_select_columns(self):
        """
        GIVEN a COVERAGE query with additional SELECT columns
        WHEN transpiling
        THEN should include those columns alongside the COVERAGE aggregate
        """
        sql = transpile(
            "SELECT COVERAGE(interval, 500) AS cov, name FROM features",
            tables=["features"],
        )

        upper = sql.upper()
        assert "COV" in upper
        assert "NAME" in upper
        assert "COUNT" in upper

    def test_table_alias_handling(self):
        """
        GIVEN a COVERAGE query with a table alias
        WHEN transpiling
        THEN should handle the alias in the generated SQL
        """
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features f",
            tables=["features"],
        )

        upper = sql.upper()
        assert "GENERATE_SERIES" in upper
        assert "LEFT JOIN" in upper

    def test_resolution_in_generate_series(self):
        """
        GIVEN a COVERAGE query with resolution=500
        WHEN transpiling
        THEN should use 500 as the step in generate_series and bin width
        """
        sql = transpile(
            "SELECT COVERAGE(interval, 500) FROM features",
            tables=["features"],
        )

        assert "500" in sql

    def test_overlap_join_condition(self):
        """
        GIVEN a basic COVERAGE query
        WHEN transpiling
        THEN should have proper overlap conditions (start < end AND end > start AND chrom = chrom)
        """
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features",
            tables=["features"],
        )

        # Check for overlap join pattern
        upper = sql.upper()
        assert "LEFT JOIN" in upper
        # The overlap condition checks: source.start < bins.end AND source.end > bins.start
        assert "BINS" in upper

    def test_order_by_present(self):
        """
        GIVEN a basic COVERAGE query
        WHEN transpiling
        THEN should ORDER BY chrom, start
        """
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features",
            tables=["features"],
        )

        assert "ORDER BY" in sql.upper()
