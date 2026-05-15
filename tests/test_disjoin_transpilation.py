"""Transpilation tests for DISJOIN SQL generation.

Tests verify that DISJOIN() is transpiled to the expected WITH-CTE subquery
shape and that target/reference resolution and coordinate-system
canonicalization are reflected in the generated SQL.
"""

from giql import transpile
from giql.table import Table


class TestDisjoinTranspilation:
    """Tests for DISJOIN transpilation to SQL."""

    def test_giqldisjoin_sql_should_emit_with_cte_subquery(self):
        """
        GIVEN a GIQL query selecting from DISJOIN(features)
        WHEN transpiling to SQL
        THEN the output should be a WITH-CTE subquery using UNION, LEAD and EXISTS
        """
        sql = transpile("SELECT * FROM DISJOIN(features)", tables=["features"]).upper()

        assert "WITH __GIQL_DJ_REF" in sql
        assert "UNION" in sql
        assert "LEAD(" in sql
        assert "EXISTS (" in sql

    def test_giqldisjoin_sql_should_emit_disjoin_columns(self):
        """
        GIVEN a GIQL query selecting from DISJOIN(features)
        WHEN transpiling to SQL
        THEN the sub-interval should be appended as disjoin_chrom/start/end
        """
        sql = transpile("SELECT * FROM DISJOIN(features)", tables=["features"])

        assert "AS disjoin_chrom" in sql
        assert "AS disjoin_start" in sql
        assert "AS disjoin_end" in sql

    def test_giqldisjoin_sql_should_default_reference_to_target_when_omitted(self):
        """
        GIVEN a DISJOIN call with no reference argument
        WHEN transpiling to SQL
        THEN the reference CTE should select from the target table
        """
        sql = transpile("SELECT * FROM DISJOIN(features)", tables=["features"])

        assert "__giql_dj_ref AS (SELECT * FROM features)" in sql

    def test_giqldisjoin_sql_should_use_explicit_reference_relation(self):
        """
        GIVEN a DISJOIN call with an explicit reference table
        WHEN transpiling to SQL
        THEN the reference CTE should select from that reference table
        """
        sql = transpile(
            "SELECT * FROM DISJOIN(features, reference := refs)",
            tables=["features", "refs"],
        )

        assert "__giql_dj_ref AS (SELECT * FROM refs)" in sql

    def test_giqldisjoin_sql_should_honor_custom_column_names(self):
        """
        GIVEN a target table with custom genomic column names
        WHEN transpiling a DISJOIN query
        THEN the generated SQL should reference the custom column names
        """
        sql = transpile(
            "SELECT * FROM DISJOIN(feats)",
            tables=[Table("feats", chrom_col="seqid", start_col="lo", end_col="hi")],
        )

        assert 't."seqid"' in sql
        assert 't."lo"' in sql
        assert 't."hi"' in sql

    def test_giqldisjoin_sql_should_canonicalize_one_based_closed_target(self):
        """
        GIVEN a target table stored as 1-based closed intervals
        WHEN transpiling a DISJOIN query
        THEN the start should be shifted to canonical 0-based coordinates
        """
        sql = transpile(
            "SELECT * FROM DISJOIN(features)",
            tables=[
                Table(
                    "features",
                    coordinate_system="1based",
                    interval_type="closed",
                )
            ],
        )

        assert '(t."start" - 1)' in sql

    def test_giqldisjoin_sql_should_emit_disjoin_start_in_target_encoding(self):
        """
        GIVEN a target table stored as 1-based closed intervals
        WHEN transpiling a DISJOIN query
        THEN disjoin_start should be shifted back to the target's encoding
        """
        sql = transpile(
            "SELECT * FROM DISJOIN(features)",
            tables=[
                Table(
                    "features",
                    coordinate_system="1based",
                    interval_type="closed",
                )
            ],
        )

        assert "(s.seg_start + 1) AS disjoin_start" in sql

    def test_giqldisjoin_sql_should_inline_subquery_reference(self):
        """
        GIVEN a DISJOIN call whose reference is a subquery
        WHEN transpiling to SQL
        THEN the subquery should be inlined as an aliased derived table
        """
        sql = transpile(
            "SELECT * FROM DISJOIN(features, reference := (SELECT * FROM refs))",
            tables=["features", "refs"],
        )

        assert "AS __giql_dj_rs" in sql
