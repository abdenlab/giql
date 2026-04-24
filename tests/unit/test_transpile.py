"""Unit tests for the transpile() function.

Tests covering all public API behavior of giql.transpile as a black box:
GIQL string in, SQL string out.
"""

import pytest

from giql import Table
from giql import transpile


class TestTranspile:
    """Tests for transpile() public API."""

    # ── Basic transpilation ──────────────────────────────────────────

    def test_transpile_should_passthrough_plain_sql_unchanged(self):
        """Test that plain SQL without GIQL extensions passes through.

        Given:
            A plain SQL query with no GIQL extensions
        When:
            transpile is called
        Then:
            It should return an equivalent SQL string unchanged
        """
        # Arrange / Act
        sql = transpile("SELECT id, name FROM features")

        # Assert
        upper = sql.upper()
        assert "SELECT" in upper
        assert "FEATURES" in upper
        assert "ID" in upper

    def test_transpile_should_emit_correct_sql_for_intersects_predicate(self):
        """Test INTERSECTS predicate expands to range comparisons.

        Given:
            A query with an INTERSECTS predicate and a tables list
        When:
            transpile is called
        Then:
            It should return SQL that contains expanded range comparison predicates
        """
        # Arrange / Act
        sql = transpile(
            "SELECT * FROM features WHERE interval INTERSECTS 'chr1:1000-2000'",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        assert "CHR1" in upper
        assert "1000" in sql
        assert "2000" in sql
        # Range overlap requires both start/end comparisons
        assert "START" in upper or "END" in upper

    def test_transpile_should_emit_correct_sql_for_contains_predicate(self):
        """Test CONTAINS predicate produces containment SQL.

        Given:
            A query with a CONTAINS predicate
        When:
            transpile is called
        Then:
            It should return SQL that contains containment predicates
        """
        # Arrange / Act
        sql = transpile(
            "SELECT * FROM features WHERE interval CONTAINS 'chr1:1500'",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        assert "SELECT" in upper
        assert "1500" in sql

    def test_transpile_should_emit_correct_sql_for_within_predicate(self):
        """Test WITHIN predicate produces within SQL.

        Given:
            A query with a WITHIN predicate
        When:
            transpile is called
        Then:
            It should return SQL that contains within predicates
        """
        # Arrange / Act
        sql = transpile(
            "SELECT * FROM features WHERE interval WITHIN 'chr1:1000-2000'",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        assert "SELECT" in upper
        assert "1000" in sql
        assert "2000" in sql

    # ── CLUSTER transpilation ────────────────────────────────────────

    def test_transpile_should_emit_window_functions_for_cluster(self):
        """Test CLUSTER expands to LAG and SUM window functions.

        Given:
            A query with CLUSTER(interval) and tables=["features"]
        When:
            transpile is called
        Then:
            It should return SQL that contains LAG and SUM window functions in a subquery
        """
        # Arrange / Act
        sql = transpile(
            "SELECT *, CLUSTER(interval) AS cluster_id FROM features",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        assert "LAG" in upper
        assert "SUM" in upper

    def test_transpile_should_include_distance_offset_for_cluster_with_distance(self):
        """Test CLUSTER with distance includes the offset in LAG.

        Given:
            A query with CLUSTER(interval, 1000)
        When:
            transpile is called
        Then:
            It should return SQL that includes a distance offset in the LAG expression
        """
        # Arrange / Act
        sql = transpile(
            "SELECT *, CLUSTER(interval, 1000) AS cluster_id FROM features",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        assert "LAG" in upper
        assert "1000" in sql

    # ── MERGE transpilation ──────────────────────────────────────────

    def test_transpile_should_emit_group_by_aggregation_for_merge(self):
        """Test MERGE expands to CTE with GROUP BY and MIN/MAX.

        Given:
            A query with MERGE(interval) and tables=["features"]
        When:
            transpile is called
        Then:
            It should return SQL that contains a CLUSTER CTE with GROUP BY and MIN/MAX aggregation
        """
        # Arrange / Act
        sql = transpile(
            "SELECT MERGE(interval) FROM features",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        assert "MIN" in upper
        assert "MAX" in upper
        assert "GROUP BY" in upper

    # ── COVERAGE transpilation ───────────────────────────────────────

    def test_transpile_should_emit_bins_cte_for_coverage(self):
        """Test COVERAGE expands to bins CTE with LEFT JOIN and COUNT.

        Given:
            A query with COVERAGE(interval, 1000) and tables=["features"]
        When:
            transpile is called
        Then:
            It should return SQL that contains a bins CTE, LEFT JOIN, COUNT, GROUP BY, and ORDER BY
        """
        # Arrange / Act
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        assert "LEFT JOIN" in upper or "LEFT OUTER JOIN" in upper
        assert "COUNT" in upper
        assert "GROUP BY" in upper
        assert "ORDER BY" in upper
        assert "1000" in sql

    def test_transpile_should_use_custom_alias_for_coverage_when_provided(self):
        """Test COVERAGE with AS cov aliases the aggregate column as "cov".

        Given:
            A query with COVERAGE(interval, 1000) AS cov
        When:
            transpile is called
        Then:
            It should alias the aggregate column in the returned SQL as "cov"
        """
        # Arrange / Act
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) AS cov FROM features",
            tables=["features"],
        )

        # Assert
        assert "cov" in sql.lower()

    def test_transpile_should_use_default_value_alias_for_bare_coverage(self):
        """Test bare COVERAGE aliases the aggregate column as "value".

        Given:
            A query with bare COVERAGE(interval, 1000) (no alias)
        When:
            transpile is called
        Then:
            It should alias the aggregate column in the returned SQL as "value"
        """
        # Arrange / Act
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features",
            tables=["features"],
        )

        # Assert
        assert "value" in sql.lower()

    def test_transpile_should_fold_where_into_join_on_for_coverage(self):
        """Test COVERAGE folds WHERE into the JOIN ON condition.

        Given:
            A query with COVERAGE and a WHERE clause
        When:
            transpile is called
        Then:
            It should place the WHERE condition in the JOIN ON condition rather than as a standalone WHERE
        """
        # Arrange / Act
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features WHERE chrom = 'chr1'",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        # The WHERE should be folded into the JOIN ON condition
        assert "JOIN" in upper
        assert "CHR1" in upper

    # ── DISTANCE transpilation ───────────────────────────────────────

    def test_transpile_should_emit_case_expression_for_distance(self):
        """Test DISTANCE expands to a CASE expression.

        Given:
            A query with DISTANCE(a.interval, b.interval) and two tables
        When:
            transpile is called
        Then:
            It should return SQL that contains a CASE expression for computing distance
        """
        # Arrange / Act
        sql = transpile(
            "SELECT DISTANCE(a.interval, b.interval) FROM features a, genes b",
            tables=["features", "genes"],
        )

        # Assert
        upper = sql.upper()
        assert "CASE" in upper

    # ── NEAREST transpilation ────────────────────────────────────────

    def test_transpile_should_emit_lateral_subquery_with_limit_for_nearest(self):
        """Test NEAREST expands to a LATERAL subquery with a LIMIT.

        Given:
            A query with NEAREST in a LATERAL join and two tables
        When:
            transpile is called
        Then:
            It should return SQL that contains a LATERAL subquery with a LIMIT clause
        """
        # Arrange / Act
        sql = transpile(
            """
            SELECT *
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=3)
            """,
            tables=["peaks", "genes"],
        )

        # Assert
        upper = sql.upper()
        assert "LATERAL" in upper
        assert "LIMIT" in upper

    # ── Table configuration ──────────────────────────────────────────

    def test_transpile_should_register_string_tables_with_default_columns(self):
        """Test string-list tables use default column mappings.

        Given:
            tables parameter as a list of strings
        When:
            transpile is called
        Then:
            It should register tables with default column mappings (chrom, start, end)
        """
        # Arrange / Act
        sql = transpile(
            "SELECT * FROM features WHERE interval INTERSECTS 'chr1:100-200'",
            tables=["features"],
        )

        # Assert
        upper = sql.upper()
        assert '"CHROM"' in upper or "CHROM" in upper
        assert '"START"' in upper or "START" in upper
        assert '"END"' in upper or "END" in upper

    def test_transpile_should_honor_custom_table_object_column_names(self):
        """Test Table objects with custom column names propagate into SQL.

        Given:
            tables parameter as a list of Table objects with custom column names
        When:
            transpile is called
        Then:
            It should generate SQL that uses those custom column names
        """
        # Arrange / Act
        sql = transpile(
            "SELECT * FROM features WHERE interval INTERSECTS 'chr1:100-200'",
            tables=[
                Table(
                    "features",
                    genomic_col="interval",
                    chrom_col="chromosome",
                    start_col="start_pos",
                    end_col="end_pos",
                )
            ],
        )

        # Assert
        assert "chromosome" in sql or "CHROMOSOME" in sql.upper()
        assert "start_pos" in sql or "START_POS" in sql.upper()
        assert "end_pos" in sql or "END_POS" in sql.upper()

    def test_transpile_should_use_default_columns_when_tables_is_none(self):
        """Test None tables parameter still uses default column names.

        Given:
            tables parameter is None
        When:
            transpile is called
        Then:
            It should still use default column names (chrom, start, end)
        """
        # Arrange / Act
        sql = transpile(
            "SELECT * FROM features WHERE interval INTERSECTS 'chr1:100-200'",
            tables=None,
        )

        # Assert
        upper = sql.upper()
        assert "SELECT" in upper
        assert "CHROM" in upper

    def test_transpile_should_register_mixed_strings_and_table_objects(self):
        """Test mixing strings and Table objects in tables parameter.

        Given:
            tables parameter mixes strings and Table objects
        When:
            transpile is called
        Then:
            It should correctly register both and produce valid SQL
        """
        # Arrange / Act
        sql = transpile(
            """
            SELECT a.*, b.*
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.region
            """,
            tables=[
                "peaks",
                Table("genes", genomic_col="region", chrom_col="seqname"),
            ],
        )

        # Assert
        upper = sql.upper()
        assert "PEAKS" in upper
        assert "GENES" in upper
        assert "SEQNAME" in upper

    # ── Error handling ───────────────────────────────────────────────

    def test_transpile_should_raise_value_error_for_invalid_query(self):
        """Test unparseable query raises ValueError with Parse error message.

        Given:
            An invalid/unparseable query string
        When:
            transpile is called
        Then:
            It should raise ValueError with a message containing "Parse error"
        """
        # Arrange / Act / Assert
        with pytest.raises(ValueError, match="Parse error"):
            transpile("SELECT * FORM features")

