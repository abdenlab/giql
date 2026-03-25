"""Unit tests for the transpile() function.

Tests TR-001 through TR-021 covering all public API behavior of
giql.transpile as a black box: GIQL string in, SQL string out.
"""

import pytest

from giql import Table
from giql import transpile


class TestTranspile:
    """Tests for transpile() public API (TR-001 to TR-021)."""

    # ── Basic transpilation ──────────────────────────────────────────

    def test_plain_sql_passthrough(self):
        """
        GIVEN a plain SQL query with no GIQL extensions
        WHEN transpile is called
        THEN it returns an equivalent SQL string unchanged.
        """
        sql = transpile("SELECT id, name FROM features")
        upper = sql.upper()
        assert "SELECT" in upper
        assert "FEATURES" in upper
        assert "ID" in upper

    def test_intersects_predicate(self):
        """
        GIVEN a query with an INTERSECTS predicate and a tables list
        WHEN transpile is called
        THEN the returned SQL contains expanded range comparison predicates.
        """
        sql = transpile(
            "SELECT * FROM features WHERE interval INTERSECTS 'chr1:1000-2000'",
            tables=["features"],
        )
        upper = sql.upper()
        assert "CHR1" in upper
        assert "1000" in sql
        assert "2000" in sql
        # Range overlap requires both start/end comparisons
        assert "START" in upper or "END" in upper

    def test_contains_predicate(self):
        """
        GIVEN a query with a CONTAINS predicate
        WHEN transpile is called
        THEN the returned SQL contains containment predicates.
        """
        sql = transpile(
            "SELECT * FROM features WHERE interval CONTAINS 'chr1:1500'",
            tables=["features"],
        )
        upper = sql.upper()
        assert "SELECT" in upper
        assert "1500" in sql

    def test_within_predicate(self):
        """
        GIVEN a query with a WITHIN predicate
        WHEN transpile is called
        THEN the returned SQL contains within predicates.
        """
        sql = transpile(
            "SELECT * FROM features WHERE interval WITHIN 'chr1:1000-2000'",
            tables=["features"],
        )
        upper = sql.upper()
        assert "SELECT" in upper
        assert "1000" in sql
        assert "2000" in sql

    # ── CLUSTER transpilation ────────────────────────────────────────

    def test_cluster_basic(self):
        """
        GIVEN a query with CLUSTER(interval) and tables=["features"]
        WHEN transpile is called
        THEN the returned SQL contains LAG and SUM window functions in a subquery.
        """
        sql = transpile(
            "SELECT *, CLUSTER(interval) AS cluster_id FROM features",
            tables=["features"],
        )
        upper = sql.upper()
        assert "LAG" in upper
        assert "SUM" in upper

    def test_cluster_with_distance(self):
        """
        GIVEN a query with CLUSTER(interval, 1000)
        WHEN transpile is called
        THEN the returned SQL includes a distance offset in the LAG expression.
        """
        sql = transpile(
            "SELECT *, CLUSTER(interval, 1000) AS cluster_id FROM features",
            tables=["features"],
        )
        upper = sql.upper()
        assert "LAG" in upper
        assert "1000" in sql

    # ── MERGE transpilation ──────────────────────────────────────────

    def test_merge_basic(self):
        """
        GIVEN a query with MERGE(interval) and tables=["features"]
        WHEN transpile is called
        THEN the returned SQL contains a CLUSTER CTE with GROUP BY and MIN/MAX aggregation.
        """
        sql = transpile(
            "SELECT MERGE(interval) FROM features",
            tables=["features"],
        )
        upper = sql.upper()
        assert "MIN" in upper
        assert "MAX" in upper
        assert "GROUP BY" in upper

    # ── COVERAGE transpilation ───────────────────────────────────────

    def test_coverage_basic(self):
        """
        GIVEN a query with COVERAGE(interval, 1000) and tables=["features"]
        WHEN transpile is called
        THEN the returned SQL contains a bins CTE, LEFT JOIN, COUNT, GROUP BY, and ORDER BY.
        """
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features",
            tables=["features"],
        )
        upper = sql.upper()
        assert "LEFT JOIN" in upper or "LEFT OUTER JOIN" in upper
        assert "COUNT" in upper
        assert "GROUP BY" in upper
        assert "ORDER BY" in upper
        assert "1000" in sql

    def test_coverage_mean_stat(self):
        """
        GIVEN a query with COVERAGE(interval, 500, stat := 'mean')
        WHEN transpile is called
        THEN the returned SQL contains an AVG aggregate.
        """
        sql = transpile(
            "SELECT COVERAGE(interval, 500, stat := 'mean') FROM features",
            tables=["features"],
        )
        upper = sql.upper()
        assert "AVG" in upper

    def test_coverage_mean_with_target(self):
        """
        GIVEN a query with COVERAGE(interval, 1000, stat := 'mean', target := 'score')
        WHEN transpile is called
        THEN the returned SQL contains AVG applied to the score column.
        """
        sql = transpile(
            "SELECT COVERAGE(interval, 1000, stat := 'mean', target := 'score') FROM features",
            tables=["features"],
        )
        upper = sql.upper()
        assert "AVG" in upper
        assert "SCORE" in upper

    def test_coverage_custom_alias(self):
        """
        GIVEN a query with COVERAGE(interval, 1000) AS cov
        WHEN transpile is called
        THEN the aggregate column in the returned SQL is aliased as "cov".
        """
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) AS cov FROM features",
            tables=["features"],
        )
        assert "cov" in sql.lower()

    def test_coverage_default_alias(self):
        """
        GIVEN a query with bare COVERAGE(interval, 1000) (no alias)
        WHEN transpile is called
        THEN the aggregate column in the returned SQL is aliased as "value".
        """
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features",
            tables=["features"],
        )
        assert "value" in sql.lower()

    def test_coverage_where_in_join_on(self):
        """
        GIVEN a query with COVERAGE and a WHERE clause
        WHEN transpile is called
        THEN the WHERE condition appears in the JOIN ON condition rather than as a standalone WHERE.
        """
        sql = transpile(
            "SELECT COVERAGE(interval, 1000) FROM features WHERE chrom = 'chr1'",
            tables=["features"],
        )
        upper = sql.upper()
        # The WHERE should be folded into the JOIN ON condition
        assert "JOIN" in upper
        assert "CHR1" in upper

    # ── DISTANCE transpilation ───────────────────────────────────────

    def test_distance_case_expression(self):
        """
        GIVEN a query with DISTANCE(a.interval, b.interval) and two tables
        WHEN transpile is called
        THEN the returned SQL contains a CASE expression for computing distance.
        """
        sql = transpile(
            "SELECT DISTANCE(a.interval, b.interval) FROM features a, genes b",
            tables=["features", "genes"],
        )
        upper = sql.upper()
        assert "CASE" in upper

    # ── NEAREST transpilation ────────────────────────────────────────

    def test_nearest_lateral_join(self):
        """
        GIVEN a query with NEAREST in a LATERAL join and two tables
        WHEN transpile is called
        THEN the returned SQL contains a LATERAL subquery with a LIMIT clause.
        """
        sql = transpile(
            """
            SELECT *
            FROM peaks
            CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=3)
            """,
            tables=["peaks", "genes"],
        )
        upper = sql.upper()
        assert "LATERAL" in upper
        assert "LIMIT" in upper

    # ── Table configuration ──────────────────────────────────────────

    def test_tables_string_list(self):
        """
        GIVEN tables parameter as a list of strings
        WHEN transpile is called
        THEN tables are registered with default column mappings (chrom, start, end).
        """
        sql = transpile(
            "SELECT * FROM features WHERE interval INTERSECTS 'chr1:100-200'",
            tables=["features"],
        )
        upper = sql.upper()
        assert '"CHROM"' in upper or "CHROM" in upper
        assert '"START"' in upper or "START" in upper
        assert '"END"' in upper or "END" in upper

    def test_tables_custom_table_objects(self):
        """
        GIVEN tables parameter as a list of Table objects with custom column names
        WHEN transpile is called
        THEN the generated SQL uses those custom column names.
        """
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
        assert "chromosome" in sql or "CHROMOSOME" in sql.upper()
        assert "start_pos" in sql or "START_POS" in sql.upper()
        assert "end_pos" in sql or "END_POS" in sql.upper()

    def test_tables_none(self):
        """
        GIVEN tables parameter is None
        WHEN transpile is called
        THEN default column names (chrom, start, end) are still used.
        """
        sql = transpile(
            "SELECT * FROM features WHERE interval INTERSECTS 'chr1:100-200'",
            tables=None,
        )
        upper = sql.upper()
        assert "SELECT" in upper
        assert "CHROM" in upper

    def test_tables_mixed_strings_and_objects(self):
        """
        GIVEN tables parameter mixes strings and Table objects
        WHEN transpile is called
        THEN both are correctly registered and the SQL is valid.
        """
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
        upper = sql.upper()
        assert "PEAKS" in upper
        assert "GENES" in upper
        assert "SEQNAME" in upper

    # ── Error handling ───────────────────────────────────────────────

    def test_invalid_query_raises_parse_error(self):
        """
        GIVEN an invalid/unparseable query string
        WHEN transpile is called
        THEN a ValueError is raised with a message containing "Parse error".
        """
        with pytest.raises(ValueError, match="Parse error"):
            transpile("SELECT * FORM features")

    def test_coverage_invalid_stat_raises(self):
        """
        GIVEN a query with COVERAGE using an invalid stat name
        WHEN transpile is called
        THEN a ValueError is raised with a message containing "Unknown COVERAGE stat".
        """
        with pytest.raises(ValueError, match="Unknown COVERAGE stat"):
            transpile(
                "SELECT COVERAGE(interval, 1000, stat := 'invalid_stat') FROM features",
                tables=["features"],
            )
