"""Tests for the transformer module.

Test specification: specs/test_transformer.md
"""

import pytest
from sqlglot import exp
from sqlglot import parse_one

from giql import transpile
from giql.dialect import GIQLDialect
from giql.generators import BaseGIQLGenerator
from giql.table import Table
from giql.table import Tables
from giql.transformer import COVERAGE_STAT_MAP
from giql.transformer import ClusterTransformer
from giql.transformer import CoverageTransformer
from giql.transformer import MergeTransformer

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


def _transform_and_sql(query: str, transformer_cls, tables: Tables | None = None) -> str:
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

    def test_csm_001_coverage_stat_map_has_correct_mappings(self):
        """GIVEN the module is imported WHEN COVERAGE_STAT_MAP is accessed THEN it maps count->COUNT, mean->AVG, sum->SUM, min->MIN, max->MAX."""
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

    def test_ct_001_basic_cluster_has_lag_and_sum_windows(self):
        """GIVEN a Tables instance and a parsed SELECT with CLUSTER(interval) WHEN transform is called THEN the result contains LAG and SUM window expressions."""
        sql = _transform_and_sql(
            "SELECT *, CLUSTER(interval) FROM features", ClusterTransformer
        )
        upper = sql.upper()
        assert "LAG" in upper
        assert "SUM" in upper

    def test_ct_002_cluster_alias_preserved(self):
        """GIVEN a parsed SELECT with CLUSTER(interval) AS cluster_id WHEN transform is called THEN the alias is preserved on the SUM window expression."""
        sql = _transform_and_sql(
            "SELECT *, CLUSTER(interval) AS cluster_id FROM features",
            ClusterTransformer,
        )
        assert "cluster_id" in sql

    def test_ct_003_cluster_with_distance(self):
        """GIVEN a parsed SELECT with CLUSTER(interval, 1000) WHEN transform is called THEN the LAG result has distance 1000 added."""
        sql = _transform_and_sql(
            "SELECT *, CLUSTER(interval, 1000) FROM features",
            ClusterTransformer,
        )
        upper = sql.upper()
        assert "LAG" in upper
        assert "1000" in sql

    def test_ct_004_cluster_stranded_partitions_by_strand(self):
        """GIVEN a parsed SELECT with CLUSTER(interval, stranded := true) WHEN transform is called THEN the result partitions by chrom AND strand."""
        sql = _transform_and_sql(
            "SELECT *, CLUSTER(interval, stranded := true) FROM features",
            ClusterTransformer,
        )
        upper = sql.upper()
        assert "STRAND" in upper
        # Both chrom and strand should appear in partition
        assert "CHROM" in upper

    def test_ct_005_non_select_returns_unchanged(self):
        """GIVEN a non-SELECT expression WHEN transform is called THEN the expression is returned unchanged."""
        tables = _make_tables("features")
        transformer = ClusterTransformer(tables)
        insert = exp.Insert(this=exp.to_table("features"))
        result = transformer.transform(insert)
        assert result is insert

    def test_ct_006_no_cluster_returns_unchanged(self):
        """GIVEN a SELECT with no CLUSTER expressions WHEN transform is called THEN the query is returned unchanged."""
        tables = _make_tables("features")
        transformer = ClusterTransformer(tables)
        ast = parse_one("SELECT * FROM features", dialect=GIQLDialect)
        result = transformer.transform(ast)
        assert result is ast

    def test_ct_007_custom_column_names_via_tables(self):
        """GIVEN a Tables instance with custom column names WHEN transform is called on a CLUSTER query THEN the generated query uses custom column names."""
        custom = Table(
            "features",
            chrom_col="chromosome",
            start_col="start_pos",
            end_col="end_pos",
        )
        tables = _make_tables(features=custom)
        sql = _transform_and_sql(
            "SELECT *, CLUSTER(interval) FROM features",
            ClusterTransformer,
            tables=tables,
        )
        assert "chromosome" in sql
        assert "start_pos" in sql
        assert "end_pos" in sql

    def test_ct_008_cluster_inside_cte_recursive_transformation(self):
        """GIVEN a SELECT with CLUSTER inside a CTE subquery WHEN transform is called THEN the CTE subquery is recursively transformed."""
        sql = _transform_and_sql(
            "WITH c AS (SELECT *, CLUSTER(interval) AS cid FROM features) "
            "SELECT * FROM c",
            ClusterTransformer,
        )
        upper = sql.upper()
        assert "LAG" in upper
        assert "SUM" in upper

    def test_ct_009_cluster_with_where_preserved(self):
        """GIVEN a SELECT with CLUSTER and a WHERE clause WHEN transform is called THEN the WHERE clause is preserved."""
        sql = _transform_and_sql(
            "SELECT *, CLUSTER(interval) FROM features WHERE score > 10",
            ClusterTransformer,
        )
        assert "score > 10" in sql

    def test_ct_010_specific_columns_with_cluster_adds_required_cols(self):
        """GIVEN a SELECT with specific columns (not *) and CLUSTER WHEN transform is called THEN missing required genomic columns are added to the CTE select list."""
        sql = _transform_and_sql(
            "SELECT name, CLUSTER(interval) AS cid FROM features",
            ClusterTransformer,
        )
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

    def test_mt_001_basic_merge_has_group_by_min_max(self):
        """GIVEN a Tables instance and a parsed SELECT with MERGE(interval) WHEN transform is called THEN the result has GROUP BY, MIN(start), MAX(end)."""
        sql = _transform_and_sql(
            "SELECT MERGE(interval) FROM features", MergeTransformer
        )
        upper = sql.upper()
        assert "GROUP BY" in upper
        assert "MIN(" in upper
        assert "MAX(" in upper

    def test_mt_002_merge_alias_dropped_output_fixed(self):
        """GIVEN a parsed SELECT with MERGE(interval) AS merged WHEN transform is called THEN the query still produces valid output with fixed columns."""
        sql = _transform_and_sql(
            "SELECT MERGE(interval) AS merged FROM features",
            MergeTransformer,
        )
        upper = sql.upper()
        assert "GROUP BY" in upper
        assert "MIN(" in upper
        assert "MAX(" in upper

    def test_mt_003_merge_with_distance(self):
        """GIVEN a parsed SELECT with MERGE(interval, 1000) WHEN transform is called THEN the distance is passed through to CLUSTER."""
        sql = _transform_and_sql(
            "SELECT MERGE(interval, 1000) FROM features",
            MergeTransformer,
        )
        assert "1000" in sql

    def test_mt_004_merge_stranded_adds_strand_to_group_by(self):
        """GIVEN a parsed SELECT with MERGE(interval, stranded := true) WHEN transform is called THEN strand appears in GROUP BY and partition."""
        sql = _transform_and_sql(
            "SELECT MERGE(interval, stranded := true) FROM features",
            MergeTransformer,
        )
        upper = sql.upper()
        assert "STRAND" in upper
        assert "GROUP BY" in upper

    def test_mt_005_non_select_returns_unchanged(self):
        """GIVEN a non-SELECT expression WHEN transform is called THEN the expression is returned unchanged."""
        tables = _make_tables("features")
        transformer = MergeTransformer(tables)
        insert = exp.Insert(this=exp.to_table("features"))
        result = transformer.transform(insert)
        assert result is insert

    def test_mt_006_no_merge_returns_unchanged(self):
        """GIVEN a SELECT with no MERGE expressions WHEN transform is called THEN the query is returned unchanged."""
        tables = _make_tables("features")
        transformer = MergeTransformer(tables)
        ast = parse_one("SELECT * FROM features", dialect=GIQLDialect)
        result = transformer.transform(ast)
        assert result is ast

    def test_mt_007_two_merge_expressions_raises_value_error(self):
        """GIVEN a SELECT with two MERGE expressions WHEN transform is called THEN it raises ValueError."""
        tables = _make_tables("features")
        transformer = MergeTransformer(tables)
        ast = parse_one(
            "SELECT MERGE(interval), MERGE(interval) FROM features",
            dialect=GIQLDialect,
        )
        with pytest.raises(ValueError, match="Multiple MERGE"):
            transformer.transform(ast)

    def test_mt_008_merge_with_where_preserved(self):
        """GIVEN a SELECT with MERGE and a WHERE clause WHEN transform is called THEN the WHERE clause is preserved in the clustered subquery."""
        sql = _transform_and_sql(
            "SELECT MERGE(interval) FROM features WHERE score > 10",
            MergeTransformer,
        )
        assert "score > 10" in sql

    def test_mt_009_merge_inside_cte_recursive_transformation(self):
        """GIVEN a SELECT with MERGE inside a CTE subquery WHEN transform is called THEN the CTE subquery is recursively transformed."""
        sql = _transform_and_sql(
            "WITH m AS (SELECT MERGE(interval) FROM features) SELECT * FROM m",
            MergeTransformer,
        )
        upper = sql.upper()
        assert "GROUP BY" in upper
        assert "MIN(" in upper
        assert "MAX(" in upper


# ===========================================================================
# TestCoverageTransformer
# ===========================================================================


class TestCoverageTransformer:
    """Tests for CoverageTransformer.transform."""

    def test_cvt_001_basic_coverage_structure(self):
        """GIVEN a Tables instance and a parsed SELECT with COVERAGE(interval, 1000) WHEN transform is called THEN the result has __giql_bins CTE, LEFT JOIN, COUNT, and GROUP BY."""
        sql = _transform_and_sql(
            "SELECT COVERAGE(interval, 1000) FROM features",
            CoverageTransformer,
        )
        upper = sql.upper()
        assert "__GIQL_BINS" in upper
        assert "LEFT JOIN" in upper
        assert "COUNT" in upper
        assert "GROUP BY" in upper

    def test_cvt_002_stat_mean_uses_avg(self):
        """GIVEN a parsed SELECT with COVERAGE(interval, 500, stat := 'mean') WHEN transform is called THEN the result uses AVG over (end - start)."""
        sql = _transform_and_sql(
            "SELECT COVERAGE(interval, 500, stat := 'mean') FROM features",
            CoverageTransformer,
        )
        upper = sql.upper()
        assert "AVG" in upper
        assert "COUNT" not in upper

    def test_cvt_003_stat_sum(self):
        """GIVEN a parsed SELECT with COVERAGE(interval, 500, stat := 'sum') WHEN transform is called THEN the result uses SUM."""
        sql = _transform_and_sql(
            "SELECT COVERAGE(interval, 500, stat := 'sum') FROM features",
            CoverageTransformer,
        )
        assert "SUM" in sql.upper()

    def test_cvt_004_stat_min(self):
        """GIVEN a parsed SELECT with COVERAGE(interval, 500, stat := 'min') WHEN transform is called THEN the result uses MIN."""
        sql = _transform_and_sql(
            "SELECT COVERAGE(interval, 500, stat := 'min') FROM features",
            CoverageTransformer,
        )
        assert "MIN(" in sql.upper()

    def test_cvt_005_stat_max(self):
        """GIVEN a parsed SELECT with COVERAGE(interval, 500, stat := 'max') WHEN transform is called THEN the result uses MAX."""
        sql = _transform_and_sql(
            "SELECT COVERAGE(interval, 500, stat := 'max') FROM features",
            CoverageTransformer,
        )
        assert "MAX(" in sql.upper()

    def test_cvt_006_stat_mean_with_target_score(self):
        """GIVEN a parsed SELECT with COVERAGE(interval, 1000, stat := 'mean', target := 'score') WHEN transform is called THEN the result uses AVG over the score column."""
        sql = _transform_and_sql(
            "SELECT COVERAGE(interval, 1000, stat := 'mean', target := 'score') FROM features",
            CoverageTransformer,
        )
        upper = sql.upper()
        assert "AVG" in upper
        assert "SCORE" in upper

    def test_cvt_007_target_score_with_default_count(self):
        """GIVEN a parsed SELECT with COVERAGE(interval, 1000, target := 'score') and default count stat WHEN transform is called THEN the result uses COUNT over the score column."""
        sql = _transform_and_sql(
            "SELECT COVERAGE(interval, 1000, target := 'score') FROM features",
            CoverageTransformer,
        )
        upper = sql.upper()
        assert "COUNT" in upper
        assert "SCORE" in upper
        # Should NOT have COUNT(source.*)
        assert ".*)" not in sql

    def test_cvt_008_coverage_alias_preserved(self):
        """GIVEN a parsed SELECT with COVERAGE(interval, 1000) AS cov WHEN transform is called THEN the aggregate column uses the alias 'cov'."""
        sql = _transform_and_sql(
            "SELECT COVERAGE(interval, 1000) AS cov FROM features",
            CoverageTransformer,
        )
        assert "AS cov" in sql
        assert "AS value" not in sql

    def test_cvt_009_bare_coverage_default_alias_value(self):
        """GIVEN a parsed SELECT with bare COVERAGE(interval, 1000) (no alias) WHEN transform is called THEN the aggregate column is aliased as 'value'."""
        sql = _transform_and_sql(
            "SELECT COVERAGE(interval, 1000) FROM features",
            CoverageTransformer,
        )
        assert "AS value" in sql

    def test_cvt_010_non_select_returns_unchanged(self):
        """GIVEN a non-SELECT expression WHEN transform is called THEN the expression is returned unchanged."""
        tables = _make_tables("features")
        transformer = CoverageTransformer(tables)
        insert = exp.Insert(this=exp.to_table("features"))
        result = transformer.transform(insert)
        assert result is insert

    def test_cvt_011_no_coverage_returns_unchanged(self):
        """GIVEN a SELECT with no COVERAGE expressions WHEN transform is called THEN the query is returned unchanged."""
        tables = _make_tables("features")
        transformer = CoverageTransformer(tables)
        ast = parse_one("SELECT * FROM features", dialect=GIQLDialect)
        result = transformer.transform(ast)
        assert result is ast

    def test_cvt_012_two_coverage_raises_value_error(self):
        """GIVEN a SELECT with two COVERAGE expressions WHEN transform is called THEN it raises ValueError."""
        tables = _make_tables("features")
        transformer = CoverageTransformer(tables)
        ast = parse_one(
            "SELECT COVERAGE(interval, 1000), COVERAGE(interval, 500) FROM features",
            dialect=GIQLDialect,
        )
        with pytest.raises(ValueError, match="Multiple COVERAGE"):
            transformer.transform(ast)

    def test_cvt_013_where_in_join_on_and_chroms_subquery(self):
        """GIVEN a parsed SELECT with COVERAGE and a WHERE clause WHEN transform is called THEN the WHERE is merged into the LEFT JOIN ON condition AND applied to the chroms subquery."""
        sql = _transform_and_sql(
            "SELECT COVERAGE(interval, 1000) FROM features WHERE score > 10",
            CoverageTransformer,
        )
        upper = sql.upper()
        # WHERE should be in the ON clause
        after_join = sql.split("LEFT JOIN")[1]
        on_clause = after_join.split("GROUP BY")[0]
        assert "score > 10" in on_clause
        # WHERE should also be in the chroms subquery (the CTE part)
        cte_part = sql.split(") SELECT")[0]
        assert "score > 10" in cte_part

    def test_cvt_014_custom_column_names(self):
        """GIVEN a Tables instance with custom column names WHEN transform is called on a COVERAGE query THEN the generated query uses custom column names."""
        custom = Table(
            "peaks",
            chrom_col="chromosome",
            start_col="start_pos",
            end_col="end_pos",
        )
        tables = _make_tables(peaks=custom)
        sql = _transform_and_sql(
            "SELECT COVERAGE(interval, 1000) FROM peaks",
            CoverageTransformer,
            tables=tables,
        )
        assert "chromosome" in sql
        assert "start_pos" in sql
        assert "end_pos" in sql

    def test_cvt_015_non_integer_resolution_raises_value_error(self):
        """GIVEN a parsed SELECT with COVERAGE where resolution is not an integer literal WHEN transform is called THEN it raises ValueError about resolution."""
        tables = _make_tables("features")
        transformer = CoverageTransformer(tables)
        # Construct an AST manually with a non-integer resolution
        from giql.expressions import GIQLCoverage

        coverage = GIQLCoverage(
            this=exp.column("interval"),
            resolution=exp.column("some_col"),
        )
        ast = exp.Select().select(coverage).from_("features")
        with pytest.raises(ValueError, match="resolution"):
            transformer.transform(ast)

    def test_cvt_016_invalid_stat_raises_value_error(self):
        """GIVEN a parsed SELECT with COVERAGE(interval, 1000, stat := 'invalid') WHEN transform is called THEN it raises ValueError about unknown stat."""
        tables = _make_tables("features")
        transformer = CoverageTransformer(tables)
        ast = parse_one(
            "SELECT COVERAGE(interval, 1000, stat := 'invalid') FROM features",
            dialect=GIQLDialect,
        )
        with pytest.raises(ValueError, match="Unknown COVERAGE stat"):
            transformer.transform(ast)

    def test_cvt_017_coverage_inside_cte_recursive_transformation(self):
        """GIVEN a parsed SELECT with COVERAGE inside a CTE subquery WHEN transform is called THEN the CTE subquery is recursively transformed."""
        sql = _transform_and_sql(
            "WITH cov AS (SELECT COVERAGE(interval, 1000) FROM features) "
            "SELECT * FROM cov",
            CoverageTransformer,
        )
        upper = sql.upper()
        assert "__GIQL_BINS" in upper
        assert "LEFT JOIN" in upper
        assert "COUNT" in upper

    def test_cvt_018_table_alias_used_as_source_ref(self):
        """GIVEN a query FROM a table with an alias (FROM features AS f) WHEN transform is called THEN the source_ref in the generated SQL uses the alias."""
        sql = _transform_and_sql(
            "SELECT COVERAGE(interval, 1000) FROM features AS f",
            CoverageTransformer,
        )
        upper = sql.upper()
        assert "LEFT JOIN" in upper
        # The alias 'f' should appear as the source reference in the join
        assert "f." in sql or "AS f" in sql

    def test_cvt_019_bins_cte_has_generate_series_with_cross_join_lateral(self):
        """GIVEN the bins CTE in a basic COVERAGE transformation WHEN the SQL is inspected THEN it contains generate_series with CROSS JOIN LATERAL."""
        sql = _transform_and_sql(
            "SELECT COVERAGE(interval, 1000) FROM features",
            CoverageTransformer,
        )
        upper = sql.upper()
        assert "GENERATE_SERIES" in upper
        assert "CROSS JOIN" in upper
        assert "LATERAL" in upper

    def test_cvt_020_output_ordered_by_bins_chrom_bins_start(self):
        """GIVEN a COVERAGE transformation output WHEN the ORDER BY clause is inspected THEN the output is ordered by bins.chrom, bins.start."""
        sql = _transform_and_sql(
            "SELECT COVERAGE(interval, 1000) FROM features",
            CoverageTransformer,
        )
        upper = sql.upper()
        assert "ORDER BY" in upper
        # Extract ORDER BY clause
        order_by_part = sql.split("ORDER BY")[1]
        order_upper = order_by_part.upper()
        assert "BINS" in order_upper
        assert "CHROM" in order_upper
        assert "START" in order_upper
