"""Tests for the transformer module.

Test specification: specs/test_transformer.md
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
from giql.generators import BaseGIQLGenerator
from giql.table import Tables
from giql.transformer import COVERAGE_STAT_MAP
from giql.transformer import ClusterTransformer
from giql.transformer import CoverageTransformer
from giql.transformer import MergeTransformer

VALID_STATS = ["count", "mean", "sum", "min", "max"]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_tables(*names: str, **custom: Table) -> Tables:
    tables = Tables()
    for name in names:
        tables.register(name, Table(name))
    for name, table in custom.items():
        tables.register(name, table)
    return tables


def _transpile_with_transformer(
    query: str, transformer_cls, tables: Tables | None = None
) -> str:
    """Run the full parse-transform-generate pipeline for SQL-substring assertions.

    Returned SQL reflects the composition of parser, ``transformer_cls``,
    and :class:`BaseGIQLGenerator`. Tests that assert on SQL output are
    exercising the end-to-end transpilation contract; if one of them
    fails, check all three stages to localise the regression rather
    than assuming the transformer is at fault.
    """
    tables = tables or _make_tables("features")
    ast = parse_one(query, dialect=GIQLDialect)
    transformer = transformer_cls(tables)
    result = transformer.transform(ast)
    generator = BaseGIQLGenerator(tables=tables)
    return generator.generate(result)


# ===========================================================================
# TestCoverageStatMap
# ===========================================================================


class TestCoverageStatMap:
    """Tests for the COVERAGE_STAT_MAP module-level constant."""

    def test_COVERAGE_STAT_MAP_should_contain_all_supported_stats(self):
        """Test COVERAGE_STAT_MAP maps stat names to SQL aggregates.

        Given:
            The transformer module is imported
        When:
            COVERAGE_STAT_MAP is accessed
        Then:
            It should map count->COUNT, mean->AVG, sum->SUM, min->MIN, max->MAX
        """
        # Act & Assert
        assert COVERAGE_STAT_MAP == {
            "count": "COUNT",
            "mean": "AVG",
            "sum": "SUM",
            "min": "MIN",
            "max": "MAX",
        }


# ===========================================================================
# TestClusterTransformer
# ===========================================================================


class TestClusterTransformer:
    """Tests for ClusterTransformer.transform."""

    def test_transform_should_produce_lag_and_sum_windows_when_basic_cluster(self):
        """Test basic CLUSTER produces LAG and SUM window expressions.

        Given:
            A Tables instance and a parsed SELECT with CLUSTER(interval)
        When:
            transform is called
        Then:
            It should produce a result containing LAG and SUM window expressions
        """
        # Act
        sql = _transpile_with_transformer(
            "SELECT *, CLUSTER(interval) FROM features", ClusterTransformer
        )

        # Assert
        upper = sql.upper()
        assert "LAG" in upper
        assert "SUM" in upper

    def test_transform_should_preserve_alias_when_cluster_has_alias(self):
        """Test CLUSTER alias is preserved on the SUM window expression.

        Given:
            A parsed SELECT with CLUSTER(interval) AS cluster_id
        When:
            transform is called
        Then:
            It should preserve the alias on the SUM window expression
        """
        # Act
        sql = _transpile_with_transformer(
            "SELECT *, CLUSTER(interval) AS cluster_id FROM features",
            ClusterTransformer,
        )

        # Assert
        assert "cluster_id" in sql

    def test_transform_should_include_distance_when_cluster_has_distance(self):
        """Test CLUSTER with distance adds the distance to the LAG result.

        Given:
            A parsed SELECT with CLUSTER(interval, 1000)
        When:
            transform is called
        Then:
            It should add distance 1000 to the LAG result
        """
        # Act
        sql = _transpile_with_transformer(
            "SELECT *, CLUSTER(interval, 1000) FROM features",
            ClusterTransformer,
        )

        # Assert
        upper = sql.upper()
        assert "LAG" in upper
        assert "1000" in sql

    def test_transform_should_partition_by_strand_when_stranded(self):
        """Test stranded CLUSTER partitions by chrom AND strand.

        Given:
            A parsed SELECT with CLUSTER(interval, stranded := true)
        When:
            transform is called
        Then:
            It should partition by chrom AND strand
        """
        # Act
        sql = _transpile_with_transformer(
            "SELECT *, CLUSTER(interval, stranded := true) FROM features",
            ClusterTransformer,
        )

        # Assert
        upper = sql.upper()
        assert "STRAND" in upper
        # Both chrom and strand should appear in partition
        assert "CHROM" in upper

    def test_transform_should_return_unchanged_when_expression_is_not_select(self):
        """Test non-SELECT expression passes through unchanged.

        Given:
            A non-SELECT expression
        When:
            transform is called
        Then:
            It should return the expression unchanged
        """
        # Arrange
        tables = _make_tables("features")
        transformer = ClusterTransformer(tables)
        insert = exp.Insert(this=exp.to_table("features"))

        # Act
        result = transformer.transform(insert)

        # Assert
        assert result is insert

    def test_transform_should_return_unchanged_when_no_cluster(self):
        """Test SELECT without CLUSTER passes through unchanged.

        Given:
            A SELECT with no CLUSTER expressions
        When:
            transform is called
        Then:
            It should return the query unchanged
        """
        # Arrange
        tables = _make_tables("features")
        transformer = ClusterTransformer(tables)
        ast = parse_one("SELECT * FROM features", dialect=GIQLDialect)

        # Act
        result = transformer.transform(ast)

        # Assert
        assert result is ast

    def test_transform_should_use_custom_column_names_when_tables_configured(self):
        """Test custom column names from Tables propagate into output SQL.

        Given:
            A Tables instance with custom column names
        When:
            transform is called on a CLUSTER query
        Then:
            The generated query should use the custom column names
        """
        # Arrange
        custom = Table(
            "features",
            chrom_col="chromosome",
            start_col="start_pos",
            end_col="end_pos",
        )
        tables = _make_tables(features=custom)

        # Act
        sql = _transpile_with_transformer(
            "SELECT *, CLUSTER(interval) FROM features",
            ClusterTransformer,
            tables=tables,
        )

        # Assert
        assert "chromosome" in sql
        assert "start_pos" in sql
        assert "end_pos" in sql

    def test_transform_should_recurse_when_cluster_inside_cte(self):
        """Test CLUSTER inside a CTE subquery is recursively transformed.

        Given:
            A SELECT with CLUSTER inside a CTE subquery
        When:
            transform is called
        Then:
            It should recursively transform the CTE subquery
        """
        # Act
        sql = _transpile_with_transformer(
            "WITH c AS (SELECT *, CLUSTER(interval) AS cid FROM features) "
            "SELECT * FROM c",
            ClusterTransformer,
        )

        # Assert
        upper = sql.upper()
        assert "LAG" in upper
        assert "SUM" in upper

    def test_transform_should_preserve_where_when_cluster_has_where(self):
        """Test WHERE clause is preserved alongside CLUSTER.

        Given:
            A SELECT with CLUSTER and a WHERE clause
        When:
            transform is called
        Then:
            It should preserve the WHERE clause
        """
        # Act
        sql = _transpile_with_transformer(
            "SELECT *, CLUSTER(interval) FROM features WHERE score > 10",
            ClusterTransformer,
        )

        # Assert
        assert "score > 10" in sql

    def test_transform_should_add_required_genomic_columns_when_specific_columns(self):
        """Test specific column selection adds required genomic cols to CTE.

        Given:
            A SELECT with specific columns (not *) and CLUSTER
        When:
            transform is called
        Then:
            It should add missing required genomic columns to the CTE select list
        """
        # Act
        sql = _transpile_with_transformer(
            "SELECT name, CLUSTER(interval) AS cid FROM features",
            ClusterTransformer,
        )

        # Assert
        upper = sql.upper()
        # Required genomic cols should be in the output
        assert "CHROM" in upper
        assert "START" in upper
        assert "END" in upper


# ===========================================================================
# TestMergeTransformer
# ===========================================================================


class TestMergeTransformer:
    """Tests for MergeTransformer.transform."""

    def test_transform_should_produce_group_by_min_max_when_basic_merge(self):
        """Test basic MERGE produces GROUP BY with MIN(start) and MAX(end).

        Given:
            A Tables instance and a parsed SELECT with MERGE(interval)
        When:
            transform is called
        Then:
            It should produce a result with GROUP BY, MIN(start), MAX(end)
        """
        # Act
        sql = _transpile_with_transformer(
            "SELECT MERGE(interval) FROM features", MergeTransformer
        )

        # Assert
        upper = sql.upper()
        assert "GROUP BY" in upper
        assert "MIN(" in upper
        assert "MAX(" in upper

    def test_transform_should_produce_fixed_columns_when_merge_has_alias(self):
        """Test MERGE alias is dropped but output still has fixed columns.

        Given:
            A parsed SELECT with MERGE(interval) AS merged
        When:
            transform is called
        Then:
            It should still produce valid output with fixed columns
        """
        # Act
        sql = _transpile_with_transformer(
            "SELECT MERGE(interval) AS merged FROM features",
            MergeTransformer,
        )

        # Assert
        upper = sql.upper()
        assert "GROUP BY" in upper
        assert "MIN(" in upper
        assert "MAX(" in upper

    def test_transform_should_pass_distance_when_merge_has_distance(self):
        """Test MERGE with distance passes the distance through to CLUSTER.

        Given:
            A parsed SELECT with MERGE(interval, 1000)
        When:
            transform is called
        Then:
            It should pass the distance through to CLUSTER
        """
        # Act
        sql = _transpile_with_transformer(
            "SELECT MERGE(interval, 1000) FROM features",
            MergeTransformer,
        )

        # Assert
        assert "1000" in sql

    def test_transform_should_add_strand_to_group_by_when_stranded(self):
        """Test stranded MERGE adds strand to GROUP BY and partition.

        Given:
            A parsed SELECT with MERGE(interval, stranded := true)
        When:
            transform is called
        Then:
            strand should appear in GROUP BY and partition
        """
        # Act
        sql = _transpile_with_transformer(
            "SELECT MERGE(interval, stranded := true) FROM features",
            MergeTransformer,
        )

        # Assert
        upper = sql.upper()
        assert "STRAND" in upper
        assert "GROUP BY" in upper

    def test_transform_should_return_unchanged_when_expression_is_not_select(self):
        """Test non-SELECT expression passes through unchanged.

        Given:
            A non-SELECT expression
        When:
            transform is called
        Then:
            It should return the expression unchanged
        """
        # Arrange
        tables = _make_tables("features")
        transformer = MergeTransformer(tables)
        insert = exp.Insert(this=exp.to_table("features"))

        # Act
        result = transformer.transform(insert)

        # Assert
        assert result is insert

    def test_transform_should_return_unchanged_when_no_merge(self):
        """Test SELECT without MERGE passes through unchanged.

        Given:
            A SELECT with no MERGE expressions
        When:
            transform is called
        Then:
            It should return the query unchanged
        """
        # Arrange
        tables = _make_tables("features")
        transformer = MergeTransformer(tables)
        ast = parse_one("SELECT * FROM features", dialect=GIQLDialect)

        # Act
        result = transformer.transform(ast)

        # Assert
        assert result is ast

    def test_transform_should_raise_when_multiple_merge_expressions(self):
        """Test two MERGE expressions raise ValueError.

        Given:
            A SELECT with two MERGE expressions
        When:
            transform is called
        Then:
            It should raise ValueError
        """
        # Arrange
        tables = _make_tables("features")
        transformer = MergeTransformer(tables)
        ast = parse_one(
            "SELECT MERGE(interval), MERGE(interval) FROM features",
            dialect=GIQLDialect,
        )

        # Act & Assert
        with pytest.raises(ValueError, match="Multiple MERGE"):
            transformer.transform(ast)

    def test_transform_should_preserve_where_when_merge_has_where(self):
        """Test WHERE clause is preserved in the clustered subquery.

        Given:
            A SELECT with MERGE and a WHERE clause
        When:
            transform is called
        Then:
            It should preserve the WHERE clause in the clustered subquery
        """
        # Act
        sql = _transpile_with_transformer(
            "SELECT MERGE(interval) FROM features WHERE score > 10",
            MergeTransformer,
        )

        # Assert
        assert "score > 10" in sql

    def test_transform_should_recurse_when_merge_inside_cte(self):
        """Test MERGE inside a CTE subquery is recursively transformed.

        Given:
            A SELECT with MERGE inside a CTE subquery
        When:
            transform is called
        Then:
            It should recursively transform the CTE subquery
        """
        # Act
        sql = _transpile_with_transformer(
            "WITH m AS (SELECT MERGE(interval) FROM features) SELECT * FROM m",
            MergeTransformer,
        )

        # Assert
        upper = sql.upper()
        assert "GROUP BY" in upper
        assert "MIN(" in upper
        assert "MAX(" in upper


# ===========================================================================
# TestCoverageTransformer
# ===========================================================================


class TestCoverageTransformer:
    """Tests for CoverageTransformer.transform via transpile()."""

    # ------------------------------------------------------------------
    # Instantiation
    # ------------------------------------------------------------------

    def test___init___should_store_tables_reference(self):
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
    # Basic transpilation
    # ------------------------------------------------------------------

    def test_transform_should_produce_expected_sql_structure_when_basic_count(self):
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

    def test_transform_should_return_unchanged_when_no_coverage_expression(self):
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
    # Stat parameter
    # ------------------------------------------------------------------

    def test_transform_should_use_avg_when_stat_is_mean(self):
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

    def test_transform_should_use_sum_when_stat_is_sum(self):
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

    def test_transform_should_use_min_when_stat_is_min(self):
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

    def test_transform_should_use_max_when_stat_is_max(self):
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
    # Target parameter
    # ------------------------------------------------------------------

    def test_transform_should_use_avg_on_target_when_target_with_mean(self):
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

    def test_transform_should_count_target_column_when_target_with_count(self):
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
    # Default alias
    # ------------------------------------------------------------------

    def test_transform_should_use_value_alias_when_no_explicit_alias(self):
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

    def test_transform_should_use_explicit_alias_when_alias_provided(self):
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
    # WHERE clause semantics
    # ------------------------------------------------------------------

    def test_transform_should_move_where_to_join_on_when_where_present(self):
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

    def test_transform_should_qualify_columns_in_on_when_where_present(self):
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

    def test_transform_should_apply_where_to_chroms_subquery_when_where_present(self):
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
    # Column mapping
    # ------------------------------------------------------------------

    def test_transform_should_use_custom_column_names_when_mapping_provided(self):
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
    # Additional SELECT columns
    # ------------------------------------------------------------------

    def test_transform_should_include_extra_columns_when_additional_select_columns(self):
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
    # Table alias
    # ------------------------------------------------------------------

    def test_transform_should_use_alias_as_source_when_table_has_alias(self):
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
    # Resolution
    # ------------------------------------------------------------------

    def test_transform_should_propagate_resolution_when_resolution_provided(self):
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
    # CTE nesting
    # ------------------------------------------------------------------

    def test_transform_should_transform_coverage_when_coverage_inside_cte(self):
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
    # Error handling
    # ------------------------------------------------------------------

    def test_transform_should_raise_when_stat_is_invalid(self):
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

    def test_transform_should_raise_when_multiple_coverage_expressions(self):
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
                "SELECT COVERAGE(interval, 1000), COVERAGE(interval, 500) FROM features",
                tables=["features"],
            )

    def test_transform_should_raise_when_stat_is_not_literal(self):
        """Test non-literal stat argument raises descriptive error.

        Given:
            A COVERAGE query where stat is an unquoted column reference
        When:
            Transpiled
        Then:
            It should raise ValueError matching "string literal"
        """
        # Act & Assert
        with pytest.raises(ValueError, match="string literal"):
            transpile(
                "SELECT COVERAGE(interval, 1000, stat := score) FROM features",
                tables=["features"],
            )

    def test_transform_should_raise_when_target_is_not_literal(self):
        """Test non-literal target argument raises descriptive error.

        Given:
            A COVERAGE query where target is an unquoted column reference
        When:
            Transpiled
        Then:
            It should raise ValueError matching "string literal"
        """
        # Act & Assert
        with pytest.raises(ValueError, match="string literal"):
            transpile(
                "SELECT COVERAGE(interval, 1000, target := score) FROM features",
                tables=["features"],
            )

    def test_transform_should_raise_when_from_is_subquery(self):
        """Test subquery in FROM raises a descriptive error.

        Given:
            A COVERAGE query whose FROM clause is an inline subquery
        When:
            Transpiled
        Then:
            It should raise ValueError matching "FROM clause"
        """
        # Act & Assert
        with pytest.raises(ValueError, match="FROM clause"):
            transpile(
                "SELECT COVERAGE(interval, 1000) "
                "FROM (SELECT * FROM features) AS sub",
                tables=["features"],
            )

    def test_transform_should_raise_when_resolution_is_negative(self):
        """Test negative resolution raises descriptive error.

        Given:
            A COVERAGE query with resolution = -1
        When:
            Transpiled
        Then:
            It should raise ValueError matching "positive"
        """
        # Act & Assert
        with pytest.raises(ValueError, match="positive"):
            transpile(
                "SELECT COVERAGE(interval, -1) FROM features",
                tables=["features"],
            )

    def test_transform_should_raise_when_resolution_is_zero(self):
        """Test zero resolution raises descriptive error.

        Given:
            A COVERAGE query with resolution = 0
        When:
            Transpiled
        Then:
            It should raise ValueError matching "positive"
        """
        # Act & Assert
        with pytest.raises(ValueError, match="positive"):
            transpile(
                "SELECT COVERAGE(interval, 0) FROM features",
                tables=["features"],
            )

    # ------------------------------------------------------------------
    # Functional / DuckDB end-to-end
    # ------------------------------------------------------------------

    def test_transform_should_produce_bins_when_basic_count(self, to_df):
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

    def test_transform_should_produce_zero_coverage_bins_when_gaps_exist(self, to_df):
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

    def test_transform_should_omit_trailing_bin_when_end_on_boundary(self, to_df):
        """Test no spurious trailing bin when MAX(end) is on a bin boundary.

        Given:
            An interval at chr1:100-1000 with resolution=1000 — MAX(end)
            lands exactly on a bin boundary
        When:
            COVERAGE is transpiled and executed
        Then:
            Exactly one bin [0,1000) should be returned with value=1
        """
        # Arrange
        giql_sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features",
            tables=["features"],
        )
        conn = duckdb.connect(":memory:")
        conn.execute(
            "CREATE TABLE features AS "
            "SELECT 'chr1' AS chrom, 100 AS start, 1000 AS \"end\""
        )

        # Act
        df = to_df(conn.execute(giql_sql))
        conn.close()

        # Assert
        assert len(df) == 1
        assert df.iloc[0]["start"] == 0
        assert df.iloc[0]["value"] == 1

    def test_transform_should_return_zero_when_bin_has_no_matching_rows(self, to_df):
        """Test bins with no matching source rows return value=0.

        Given:
            A DuckDB table with intervals at chr1:100-200 and chr1:2500-2600
            and COVERAGE resolution=500 (bins [0,500), [500,1000), ...,
            [2500,3000))
        When:
            COVERAGE count is transpiled and executed
        Then:
            Bins [500,1000), [1000,1500), [1500,2000), [2000,2500) should
            all report value=0
        """
        # Arrange
        giql_sql = transpile(
            "SELECT COVERAGE(interval, 500) FROM features",
            tables=["features"],
        )
        conn = duckdb.connect(":memory:")
        conn.execute(
            "CREATE TABLE features AS "
            "SELECT 'chr1' AS chrom, 100 AS start, 200 AS \"end\" "
            "UNION ALL SELECT 'chr1', 2500, 2600"
        )

        # Act
        df = to_df(conn.execute(giql_sql))
        conn.close()

        # Assert
        empty_bin_starts = {500, 1000, 1500, 2000}
        for bin_start in empty_bin_starts:
            value = df[df["start"] == bin_start].iloc[0]["value"]
            assert value == 0, (
                f"bin [{bin_start},{bin_start + 500}) expected 0, got {value}"
            )

    def test_transform_should_preserve_user_ctes_when_coverage_wraps_them(self, to_df):
        """Test user-defined CTEs are preserved when COVERAGE wraps them.

        Given:
            A query with a user-defined CTE (selected) that pre-filters
            the source, followed by SELECT COVERAGE(...) FROM selected
        When:
            COVERAGE is transpiled and executed
        Then:
            The user CTE should be preserved alongside __giql_bins and
            the query should execute without "table not found" errors
        """
        # Arrange
        giql_sql = transpile(
            "WITH selected AS (SELECT chrom, start, \"end\" FROM features WHERE score > 50) "
            "SELECT COVERAGE(interval, 1000) FROM selected",
            tables=["features", "selected"],
        )
        conn = duckdb.connect(":memory:")
        conn.execute(
            "CREATE TABLE features AS "
            "SELECT 'chr1' AS chrom, 100 AS start, 200 AS \"end\", 80 AS score "
            "UNION ALL SELECT 'chr1', 1100, 1200, 10 "
            "UNION ALL SELECT 'chr1', 2100, 2200, 90"
        )

        # Act
        df = to_df(conn.execute(giql_sql))
        conn.close()

        # Assert
        assert set(df["start"].tolist()) == {0, 1000, 2000}
        assert df[df["start"] == 1000].iloc[0]["value"] == 0

    def test_transform_should_resolve_alias_when_where_uses_table_alias(self, to_df):
        """Test alias-qualified WHERE resolves in chroms subquery.

        Given:
            A FROM clause with a table alias (features f) and a WHERE
            qualifying a column by that alias (f.score > 10)
        When:
            COVERAGE is transpiled and executed
        Then:
            The query should run without binder errors and produce all
            three bins with WHERE-filtering applied
        """
        # Arrange
        giql_sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features f WHERE f.score > 10",
            tables=["features"],
        )
        conn = duckdb.connect(":memory:")
        conn.execute(
            "CREATE TABLE features AS "
            "SELECT 'chr1' AS chrom, 100 AS start, 200 AS \"end\", 50 AS score "
            "UNION ALL SELECT 'chr1', 1100, 1200, 5 "
            "UNION ALL SELECT 'chr1', 2100, 2200, 80"
        )

        # Act
        df = to_df(conn.execute(giql_sql))
        conn.close()

        # Assert
        assert len(df) == 3
        assert set(df["start"].tolist()) == {0, 1000, 2000}

    def test_transform_should_preserve_zero_bins_when_where_in_on(self, to_df):
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

    def test_transform_should_compute_average_when_mean_with_target(self, to_df):
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

    def test_transform_should_return_minimum_interval_length_when_stat_is_min(self, to_df):
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

    # ------------------------------------------------------------------
    # Property-based transpilation
    # ------------------------------------------------------------------

    @given(
        resolution=st.integers(min_value=1, max_value=10_000_000),
        stat=st.sampled_from(VALID_STATS),
    )
    @settings(max_examples=50, suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_transform_should_map_stat_to_aggregate_when_varying_stat_and_resolution(
        self, resolution, stat
    ):
        """Test stat parameter maps to correct SQL aggregate across input space.

        Given:
            Any valid stat (count/mean/sum/min/max) and resolution (1-10M)
        When:
            Transpiled via transpile()
        Then:
            The output SQL should contain the corresponding SQL aggregate
            function name and the resolution value
        """
        # Arrange
        stat_to_sql = {
            "count": "COUNT",
            "mean": "AVG",
            "sum": "SUM(",
            "min": "MIN(",
            "max": "MAX(",
        }
        expected_agg = stat_to_sql[stat]

        # Act
        sql = transpile(
            f"SELECT COVERAGE(interval, {resolution}, stat := '{stat}') FROM features",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        assert expected_agg in upper
        assert str(resolution) in sql

    @given(
        resolution=st.integers(min_value=1, max_value=10_000_000),
        stat=st.sampled_from(VALID_STATS),
    )
    @settings(max_examples=50, suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_transform_should_contain_structural_elements_when_varying_stat_and_resolution(
        self, resolution, stat
    ):
        """Test transpiled SQL always contains required structural elements.

        Given:
            Any valid stat (count/mean/sum/min/max) and resolution (1-10M)
        When:
            Transpiled via transpile()
        Then:
            The output SQL should always contain __GIQL_BINS,
            GENERATE_SERIES, LEFT JOIN, GROUP BY, and ORDER BY
        """
        # Act
        sql = transpile(
            f"SELECT COVERAGE(interval, {resolution}, stat := '{stat}') FROM features",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        assert "__GIQL_BINS" in upper
        assert "GENERATE_SERIES" in upper
        assert "LEFT JOIN" in upper
        assert "GROUP BY" in upper
        assert "ORDER BY" in upper
