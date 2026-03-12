"""Tests for the COVERAGE operator.

Test specification: specs/test_coverage.md
"""

import duckdb
import pytest
from hypothesis import HealthCheck
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st
from sqlglot import exp
from sqlglot import parse_one

from giql import Table
from giql import transpile
from giql.dialect import GIQLDialect
from giql.expressions import GIQLCoverage
from giql.table import Tables
from giql.transformer import CoverageTransformer

VALID_STATS = ["count", "mean", "sum", "min", "max"]


class TestGIQLCoverage:
    """Tests for GIQLCoverage expression node parsing."""

    # ------------------------------------------------------------------
    # Example-based parsing (COV-001 to COV-007)
    # ------------------------------------------------------------------

    def test_from_arg_list_with_positional_args(self):
        """Test positional interval and resolution mapping.

        Given:
            A COVERAGE expression with positional interval and resolution
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with resolution set and
            stat/target both None
        """
        # Act
        ast = parse_one(
            "SELECT COVERAGE(interval, 1000) FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["resolution"].this == "1000"
        assert coverage[0].args.get("stat") is None
        assert coverage[0].args.get("target") is None

    def test_from_arg_list_with_walrus_named_stat(self):
        """Test named stat parameter via := syntax.

        Given:
            A COVERAGE expression with := named stat parameter
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with stat set to the given value
        """
        # Act
        ast = parse_one(
            "SELECT COVERAGE(interval, 500, stat := 'mean') FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["stat"].this == "mean"

    def test_from_arg_list_with_arrow_named_stat(self):
        """Test named stat parameter via => syntax.

        Given:
            A COVERAGE expression with => named stat parameter
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with stat set to the given value
        """
        # Act
        ast = parse_one(
            "SELECT COVERAGE(interval, 500, stat => 'mean') FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["stat"].this == "mean"

    def test_from_arg_list_with_named_resolution(self):
        """Test named resolution parameter.

        Given:
            A COVERAGE expression with named resolution parameter
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with resolution set via named param
        """
        # Act
        ast = parse_one(
            "SELECT COVERAGE(interval, resolution := 1000) FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["resolution"].this == "1000"

    def test_from_arg_list_with_walrus_named_target(self):
        """Test target parameter via := syntax.

        Given:
            A COVERAGE expression with := named target parameter
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with target set
        """
        # Act
        ast = parse_one(
            "SELECT COVERAGE(interval, 1000, target := 'score') FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["target"].this == "score"

    def test_from_arg_list_with_arrow_named_target(self):
        """Test target parameter via => syntax.

        Given:
            A COVERAGE expression with => named target parameter
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with target set
        """
        # Act
        ast = parse_one(
            "SELECT COVERAGE(interval, 1000, target => 'score') FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["target"].this == "score"

    def test_from_arg_list_with_all_named_params(self):
        """Test all parameters provided as named arguments.

        Given:
            A COVERAGE expression with stat, target, and resolution all named
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with all three params set
        """
        # Act
        ast = parse_one(
            "SELECT COVERAGE(interval, resolution := 500, "
            "stat := 'mean', target := 'score') FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["resolution"].this == "500"
        assert coverage[0].args["stat"].this == "mean"
        assert coverage[0].args["target"].this == "score"

    # ------------------------------------------------------------------
    # Property-based parsing (PBT-001 to PBT-003)
    # ------------------------------------------------------------------

    @given(
        resolution=st.integers(min_value=1, max_value=10_000_000),
        stat=st.sampled_from(VALID_STATS),
        syntax=st.sampled_from([":=", "=>"]),
    )
    @settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_from_arg_list_with_varying_stat_and_resolution(
        self, resolution, stat, syntax
    ):
        """Test stat and resolution parse correctly across input space.

        Given:
            Any valid resolution (1-10M), stat (sampled from valid values),
            and syntax (:= or =>)
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with correct resolution and stat
        """
        # Act
        sql = (
            f"SELECT COVERAGE(interval, {resolution}, "
            f"stat {syntax} '{stat}') FROM features"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["resolution"].this == str(resolution)
        assert coverage[0].args["stat"].this == stat

    @given(resolution=st.integers(min_value=1, max_value=10_000_000))
    @settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_from_arg_list_with_varying_positional_only(self, resolution):
        """Test positional-only parsing across resolution range.

        Given:
            Any valid resolution (1-10M) with no stat or target
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with resolution set and
            stat/target None
        """
        # Act
        ast = parse_one(
            f"SELECT COVERAGE(interval, {resolution}) FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["resolution"].this == str(resolution)
        assert coverage[0].args.get("stat") is None
        assert coverage[0].args.get("target") is None

    @given(syntax=st.sampled_from([":=", "=>"]))
    @settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_from_arg_list_with_varying_target_syntax(self, syntax):
        """Test target parameter parsing across syntax variants.

        Given:
            Either := or => syntax for target parameter
        When:
            Parsed with GIQLDialect
        Then:
            It should produce a GIQLCoverage node with target set
        """
        # Act
        ast = parse_one(
            f"SELECT COVERAGE(interval, 1000, target {syntax} 'score') FROM features",
            dialect=GIQLDialect,
        )

        # Assert
        coverage = list(ast.find_all(GIQLCoverage))
        assert len(coverage) == 1
        assert coverage[0].args["target"].this == "score"


class TestCoverageTransformer:
    """Tests for CoverageTransformer.transform via transpile()."""

    # ------------------------------------------------------------------
    # Instantiation (CT-001)
    # ------------------------------------------------------------------

    def test___init___with_tables(self):
        """Test CoverageTransformer stores its tables reference.

        Given:
            A Tables container with registered tables
        When:
            CoverageTransformer is instantiated
        Then:
            It should store the tables reference
        """
        # Arrange
        tables = Tables()
        tables.register("features", Table("features"))

        # Act
        transformer = CoverageTransformer(tables)

        # Assert
        assert transformer.tables is tables

    # ------------------------------------------------------------------
    # Basic transpilation (CT-002, CT-003)
    # ------------------------------------------------------------------

    def test_transform_with_basic_count(self):
        """Test basic COVERAGE produces correct SQL structure.

        Given:
            A basic COVERAGE query with count (default stat)
        When:
            Transpiled
        Then:
            It should produce SQL with __giql_bins CTE, GENERATE_SERIES,
            LEFT JOIN, GROUP BY, COUNT, and ORDER BY
        """
        # Act
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        assert "__GIQL_BINS" in upper
        assert "GENERATE_SERIES" in upper
        assert "LEFT JOIN" in upper
        assert "GROUP BY" in upper
        assert "COUNT" in upper
        assert "ORDER BY" in upper

    def test_transform_without_coverage_expression(self):
        """Test non-COVERAGE query passes through unchanged.

        Given:
            A query with no COVERAGE expression
        When:
            Transformed by CoverageTransformer
        Then:
            It should return the query unchanged
        """
        # Arrange
        tables = Tables()
        tables.register("features", Table("features"))
        transformer = CoverageTransformer(tables)
        ast = parse_one("SELECT * FROM features", dialect=GIQLDialect)

        # Act
        result = transformer.transform(ast)

        # Assert
        assert result is ast

    # ------------------------------------------------------------------
    # Stat parameter (CT-004 to CT-007)
    # ------------------------------------------------------------------

    def test_transform_with_stat_mean(self):
        """Test stat='mean' maps to AVG aggregate.

        Given:
            A COVERAGE query with stat := 'mean'
        When:
            Transpiled
        Then:
            It should use AVG aggregate, not COUNT
        """
        # Act
        sql = transpile(
            "SELECT COVERAGE(interval, 1000, stat := 'mean') FROM features",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        assert "AVG" in upper
        assert "COUNT" not in upper

    def test_transform_with_stat_sum(self):
        """Test stat='sum' maps to SUM aggregate.

        Given:
            A COVERAGE query with stat := 'sum'
        When:
            Transpiled
        Then:
            It should use SUM aggregate
        """
        # Act
        sql = transpile(
            "SELECT COVERAGE(interval, 1000, stat := 'sum') FROM features",
            tables=["features"],
        )

        # Assert
        assert "SUM" in sql.upper()

    def test_transform_with_stat_min(self):
        """Test stat='min' maps to MIN aggregate.

        Given:
            A COVERAGE query with stat := 'min'
        When:
            Transpiled
        Then:
            It should use MIN aggregate
        """
        # Act
        sql = transpile(
            "SELECT COVERAGE(interval, 1000, stat := 'min') FROM features",
            tables=["features"],
        )

        # Assert
        assert "MIN(" in sql.upper()

    def test_transform_with_stat_max(self):
        """Test stat='max' maps to MAX aggregate.

        Given:
            A COVERAGE query with stat := 'max'
        When:
            Transpiled
        Then:
            It should use MAX aggregate
        """
        # Act
        sql = transpile(
            "SELECT COVERAGE(interval, 1000, stat := 'max') FROM features",
            tables=["features"],
        )

        # Assert
        assert "MAX(" in sql.upper()

    # ------------------------------------------------------------------
    # Target parameter (CT-008, CT-009)
    # ------------------------------------------------------------------

    def test_transform_with_target_and_mean(self):
        """Test target column used with mean stat.

        Given:
            A COVERAGE query with stat := 'mean' and target := 'score'
        When:
            Transpiled
        Then:
            It should use AVG on the score column
        """
        # Act
        sql = transpile(
            "SELECT COVERAGE(interval, 1000, stat := 'mean', "
            "target := 'score') FROM features",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        assert "AVG" in upper
        assert "SCORE" in upper

    def test_transform_with_target_and_count(self):
        """Test target column used with default count stat.

        Given:
            A COVERAGE query with target := 'score' (default count)
        When:
            Transpiled
        Then:
            It should use COUNT on the score column, not COUNT(*)
        """
        # Act
        sql = transpile(
            "SELECT COVERAGE(interval, 1000, target := 'score') FROM features",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        assert "COUNT" in upper
        assert "SCORE" in upper
        assert ".*)" not in sql

    # ------------------------------------------------------------------
    # Default alias (CT-010, CT-011)
    # ------------------------------------------------------------------

    def test_transform_with_default_alias(self):
        """Test bare COVERAGE gets default 'value' alias.

        Given:
            A COVERAGE query without an explicit AS alias
        When:
            Transpiled
        Then:
            It should alias the aggregate as "value"
        """
        # Act
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features",
            tables=["features"],
        )

        # Assert
        assert "AS value" in sql

    def test_transform_with_explicit_alias(self):
        """Test explicit AS alias overrides default.

        Given:
            A COVERAGE query with explicit AS alias
        When:
            Transpiled
        Then:
            It should use the explicit alias, not "value"
        """
        # Act
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) AS depth FROM features",
            tables=["features"],
        )

        # Assert
        assert "AS depth" in sql
        assert "AS value" not in sql

    # ------------------------------------------------------------------
    # WHERE clause semantics (CT-012, CT-013, CT-014)
    # ------------------------------------------------------------------

    def test_transform_where_moves_to_join_on(self):
        """Test WHERE migrates into LEFT JOIN ON clause.

        Given:
            A COVERAGE query with a WHERE clause
        When:
            Transpiled
        Then:
            It should move the WHERE condition into the LEFT JOIN ON clause,
            not the outer WHERE
        """
        # Act
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features WHERE score > 10",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        assert "ON" in upper
        assert "SCORE > 10" in upper
        # The condition should be in the ON clause (between LEFT JOIN and GROUP BY)
        after_join = sql.split("LEFT JOIN")[1]
        on_clause = after_join.split("GROUP BY")[0]
        assert "score > 10" in on_clause

    def test_transform_where_qualifies_columns_in_on(self):
        """Test WHERE column references are qualified with source table in ON.

        Given:
            A COVERAGE query with a WHERE clause
        When:
            Transpiled
        Then:
            It should qualify unqualified column references in the JOIN ON
            with the source table
        """
        # Act
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features WHERE score > 10",
            tables=["features"],
        )

        # Assert
        after_join = sql.split("LEFT JOIN")[1]
        on_clause = after_join.split("GROUP BY")[0]
        assert "features.score" in on_clause

    def test_transform_where_applied_to_chroms_subquery(self):
        """Test WHERE is also applied to the chroms subquery.

        Given:
            A COVERAGE query with a WHERE clause
        When:
            Transpiled
        Then:
            It should also apply the WHERE to the chroms subquery with
            table-qualified columns
        """
        # Act
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features WHERE score > 10",
            tables=["features"],
        )

        # Assert
        # The chroms subquery is inside the CTE, before the outer SELECT
        cte_part = sql.split(") SELECT")[0]
        assert "features.score > 10" in cte_part

    # ------------------------------------------------------------------
    # Column mapping (CT-015)
    # ------------------------------------------------------------------

    def test_transform_with_custom_column_mapping(self):
        """Test custom column names are used throughout.

        Given:
            A COVERAGE query with custom column mappings
            (chromosome, start_pos, end_pos)
        When:
            Transpiled
        Then:
            It should use the mapped column names throughout
        """
        # Act
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

        # Assert
        assert "chromosome" in sql
        assert "start_pos" in sql
        assert "end_pos" in sql

    # ------------------------------------------------------------------
    # Additional SELECT columns (CT-016)
    # ------------------------------------------------------------------

    def test_transform_with_additional_select_columns(self):
        """Test extra SELECT columns pass through alongside COVERAGE.

        Given:
            A COVERAGE query with additional columns alongside COVERAGE
        When:
            Transpiled
        Then:
            It should include the extra columns in the output
        """
        # Act
        sql = transpile(
            "SELECT COVERAGE(interval, 500) AS cov, name FROM features",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        assert "COV" in upper
        assert "NAME" in upper
        assert "COUNT" in upper

    # ------------------------------------------------------------------
    # Table alias (CT-017)
    # ------------------------------------------------------------------

    def test_transform_with_table_alias(self):
        """Test table alias is used as source reference in JOIN.

        Given:
            A COVERAGE query with a table alias (FROM features f)
        When:
            Transpiled
        Then:
            It should use the alias as the source reference in JOIN
        """
        # Act
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features f",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        assert "GENERATE_SERIES" in upper
        assert "LEFT JOIN" in upper

    # ------------------------------------------------------------------
    # Resolution (CT-018)
    # ------------------------------------------------------------------

    def test_transform_with_resolution_propagation(self):
        """Test resolution value propagates to generate_series and bin width.

        Given:
            A COVERAGE query with resolution=500
        When:
            Transpiled
        Then:
            It should use 500 as the step in generate_series and bin width
        """
        # Act
        sql = transpile(
            "SELECT COVERAGE(interval, 500) FROM features",
            tables=["features"],
        )

        # Assert
        assert "500" in sql

    # ------------------------------------------------------------------
    # CTE nesting (CT-019)
    # ------------------------------------------------------------------

    def test_transform_with_coverage_in_cte(self):
        """Test COVERAGE inside a WITH clause is transformed correctly.

        Given:
            A COVERAGE expression inside a WITH clause
        When:
            Transpiled
        Then:
            It should correctly transform the CTE containing COVERAGE
        """
        # Act
        sql = transpile(
            "WITH cov AS (SELECT COVERAGE(interval, 1000) FROM features) "
            "SELECT * FROM cov",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        assert "GENERATE_SERIES" in upper
        assert "LEFT JOIN" in upper
        assert "COUNT" in upper

    # ------------------------------------------------------------------
    # Error handling (CT-020, CT-021)
    # ------------------------------------------------------------------

    def test_transform_with_invalid_stat(self):
        """Test invalid stat raises descriptive error.

        Given:
            A COVERAGE query with an invalid stat value
        When:
            Transpiled
        Then:
            It should raise ValueError matching "Unknown COVERAGE stat"
        """
        # Act & Assert
        with pytest.raises(ValueError, match="Unknown COVERAGE stat"):
            transpile(
                "SELECT COVERAGE(interval, 1000, stat := 'median') FROM features",
                tables=["features"],
            )

    def test_transform_with_multiple_coverage(self):
        """Test multiple COVERAGE expressions raise error.

        Given:
            A query with two COVERAGE expressions
        When:
            Transpiled
        Then:
            It should raise ValueError matching "Multiple COVERAGE"
        """
        # Act & Assert
        with pytest.raises(ValueError, match="Multiple COVERAGE"):
            transpile(
                "SELECT COVERAGE(interval, 1000), "
                "COVERAGE(interval, 500) FROM features",
                tables=["features"],
            )

    # ------------------------------------------------------------------
    # Functional / DuckDB end-to-end (CT-022 to CT-026)
    # ------------------------------------------------------------------

    def test_transform_end_to_end_basic_count(self, to_df):
        """Test count correctness with two intervals in one bin.

        Given:
            A DuckDB table with two intervals in the same 1000bp bin
        When:
            COVERAGE count is transpiled and executed
        Then:
            It should return count=2 for that bin
        """
        # Arrange
        giql_sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features",
            tables=["features"],
        )
        conn = duckdb.connect(":memory:")
        conn.execute(
            "CREATE TABLE features AS "
            "SELECT 'chr1' AS chrom, 100 AS start, 200 AS \"end\" "
            "UNION ALL SELECT 'chr1', 300, 400"
        )

        # Act
        df = to_df(conn.execute(giql_sql))
        conn.close()

        # Assert
        row = df[df["start"] == 0].iloc[0]
        assert row["value"] == 2

    def test_transform_end_to_end_zero_coverage_bins(self, to_df):
        """Test zero-coverage bins are present via LEFT JOIN.

        Given:
            A DuckDB table with intervals covering only some bins
        When:
            COVERAGE count is transpiled and executed
        Then:
            Bins beyond intervals should appear with count=0
        """
        # Arrange
        giql_sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features",
            tables=["features"],
        )
        conn = duckdb.connect(":memory:")
        conn.execute(
            "CREATE TABLE features AS "
            "SELECT 'chr1' AS chrom, 100 AS start, 200 AS \"end\" "
            "UNION ALL SELECT 'chr1', 1500, 2500"
        )

        # Act
        df = to_df(conn.execute(giql_sql))
        conn.close()

        # Assert
        assert len(df) >= 3
        assert df[df["start"] == 0].iloc[0]["value"] == 1

    def test_transform_end_to_end_where_preserves_zero_bins(self, to_df):
        """Test WHERE in ON preserves bins without matching intervals.

        Given:
            A DuckDB table with high-scoring intervals in bin [0,1000) and
            bin [2000,3000), plus a low-scoring interval in bin [1000,2000)
        When:
            COVERAGE count with WHERE score > 50 is transpiled and executed
        Then:
            All three bins should be present (the WHERE is in the ON clause
            so bins are not dropped even when no source rows match)
        """
        # Arrange
        giql_sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features WHERE score > 50",
            tables=["features"],
        )
        conn = duckdb.connect(":memory:")
        conn.execute(
            "CREATE TABLE features AS "
            "SELECT 'chr1' AS chrom, 100 AS start, 200 AS \"end\", 100 AS score "
            "UNION ALL SELECT 'chr1', 1500, 1600, 10 "
            "UNION ALL SELECT 'chr1', 2100, 2200, 80"
        )

        # Act
        df = to_df(conn.execute(giql_sql))
        conn.close()

        # Assert — all three bins are present (not filtered by WHERE)
        assert len(df) == 3
        assert set(df["start"].tolist()) == {0, 1000, 2000}

    def test_transform_end_to_end_mean_with_target(self, to_df):
        """Test mean stat with target column produces correct average.

        Given:
            A DuckDB table with a score column and two intervals in one bin
        When:
            COVERAGE with stat='mean' and target='score' is transpiled
            and executed
        Then:
            It should return the average of the score values
        """
        # Arrange
        giql_sql = transpile(
            "SELECT COVERAGE(interval, 1000, stat := 'mean', "
            "target := 'score') FROM features",
            tables=["features"],
        )
        conn = duckdb.connect(":memory:")
        conn.execute(
            "CREATE TABLE features AS "
            "SELECT 'chr1' AS chrom, 100 AS start, 200 AS \"end\", "
            "10.0 AS score "
            "UNION ALL SELECT 'chr1', 300, 400, 20.0"
        )

        # Act
        df = to_df(conn.execute(giql_sql))
        conn.close()

        # Assert
        row = df[df["start"] == 0].iloc[0]
        assert row["value"] == pytest.approx(15.0)

    def test_transform_end_to_end_min_stat(self, to_df):
        """Test min stat returns minimum interval length.

        Given:
            A DuckDB table with intervals of different lengths in one bin
        When:
            COVERAGE with stat='min' is transpiled and executed
        Then:
            It should return the minimum interval length
        """
        # Arrange
        giql_sql = transpile(
            "SELECT COVERAGE(interval, 1000, stat := 'min') FROM features",
            tables=["features"],
        )
        conn = duckdb.connect(":memory:")
        conn.execute(
            "CREATE TABLE features AS "
            "SELECT 'chr1' AS chrom, 100 AS start, 200 AS \"end\" "
            "UNION ALL SELECT 'chr1', 300, 600"
        )

        # Act
        df = to_df(conn.execute(giql_sql))
        conn.close()

        # Assert
        row = df[df["start"] == 0].iloc[0]
        assert row["value"] == 100
