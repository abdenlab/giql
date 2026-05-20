"""Transpilation tests for DISJOIN SQL generation.

Tests verify that DISJOIN() is transpiled to the expected WITH-CTE subquery
shape and that target/reference resolution and coordinate-system
canonicalization are reflected in the generated SQL.
"""

import pytest

from giql import transpile
from giql.table import Table


class TestDisjoinTranspilation:
    """Tests for DISJOIN transpilation to SQL."""

    def test_giqldisjoin_sql_should_emit_with_cte_subquery(self):
        """Test that DISJOIN transpiles to a WITH-CTE subquery.

        Given:
            A GIQL query selecting from DISJOIN(features).
        When:
            Transpiling to SQL.
        Then:
            It should be a WITH-CTE subquery using UNION and LEAD.
        """
        # Arrange & act
        sql = transpile("SELECT * FROM DISJOIN(features)", tables=["features"]).upper()

        # Assert
        assert "WITH __GIQL_DJ_REF" in sql
        assert "UNION" in sql
        assert "LEAD(" in sql

    def test_giqldisjoin_sql_should_emit_disjoin_columns(self):
        """Test that DISJOIN appends the sub-interval columns.

        Given:
            A GIQL query selecting from DISJOIN(features).
        When:
            Transpiling to SQL.
        Then:
            It should append the sub-interval as disjoin_chrom, disjoin_start,
            and disjoin_end.
        """
        # Arrange & act
        sql = transpile("SELECT * FROM DISJOIN(features)", tables=["features"])

        # Assert
        assert "AS disjoin_chrom" in sql
        assert "AS disjoin_start" in sql
        assert "AS disjoin_end" in sql

    def test_giqldisjoin_sql_should_default_reference_to_target_when_omitted(self):
        """Test that an omitted reference defaults to the target set.

        Given:
            A DISJOIN call with no reference argument.
        When:
            Transpiling to SQL.
        Then:
            It should make the reference CTE select from the target table.
        """
        # Arrange & act
        sql = transpile("SELECT * FROM DISJOIN(features)", tables=["features"])

        # Assert
        assert "__giql_dj_ref AS (SELECT * FROM features)" in sql

    def test_giqldisjoin_sql_should_use_explicit_reference_relation(self):
        """Test that an explicit reference table is used as the reference.

        Given:
            A DISJOIN call with an explicit reference table.
        When:
            Transpiling to SQL.
        Then:
            It should make the reference CTE select from that reference table.
        """
        # Arrange & act
        sql = transpile(
            "SELECT * FROM DISJOIN(features, reference := refs)",
            tables=["features", "refs"],
        )

        # Assert
        assert "__giql_dj_ref AS (SELECT * FROM refs)" in sql

    def test_giqldisjoin_sql_should_honor_custom_column_names(self):
        """Test that custom genomic column names are honored.

        Given:
            A target table with custom genomic column names.
        When:
            Transpiling a DISJOIN query.
        Then:
            It should reference the custom column names in the generated SQL.
        """
        # Arrange & act
        sql = transpile(
            "SELECT * FROM DISJOIN(feats)",
            tables=[Table("feats", chrom_col="seqid", start_col="lo", end_col="hi")],
        )

        # Assert
        assert 't."seqid"' in sql
        assert 't."lo"' in sql
        assert 't."hi"' in sql

    def test_giqldisjoin_sql_should_canonicalize_one_based_closed_target(self):
        """Test that a 1-based closed target is canonicalized.

        Given:
            A target table stored as 1-based closed intervals.
        When:
            Transpiling a DISJOIN query.
        Then:
            It should shift the start to canonical 0-based coordinates.
        """
        # Arrange & act
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

        # Assert
        assert '(t."start" - 1)' in sql

    def test_giqldisjoin_sql_should_emit_disjoin_start_in_target_encoding(self):
        """Test that disjoin_start is emitted in the target's encoding.

        Given:
            A target table stored as 1-based closed intervals.
        When:
            Transpiling a DISJOIN query.
        Then:
            It should shift disjoin_start back to the target's encoding.
        """
        # Arrange & act
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

        # Assert
        assert "(s.seg_start + 1) AS disjoin_start" in sql

    def test_giqldisjoin_sql_should_inline_subquery_reference(self):
        """Test that a subquery reference is inlined as a derived table.

        Given:
            A DISJOIN call whose reference is a subquery.
        When:
            Transpiling to SQL.
        Then:
            It should inline the subquery as an aliased derived table.
        """
        # Arrange & act
        sql = transpile(
            "SELECT * FROM DISJOIN(features, reference := (SELECT * FROM refs))",
            tables=["features", "refs"],
        )

        # Assert
        assert "AS __giql_dj_rs" in sql

    def test_giqldisjoin_sql_should_resolve_in_query_cte_reference(self):
        """Test that a reference naming an in-query CTE resolves to it.

        Given:
            A DISJOIN reference naming a CTE defined in the outer query.
        When:
            Transpiling to SQL.
        Then:
            It should make the reference CTE select from that CTE name.
        """
        # Arrange & act
        sql = transpile(
            "WITH bins AS (SELECT 1) SELECT * FROM DISJOIN(features, reference := bins)",
            tables=["features"],
        )

        # Assert
        assert "__giql_dj_ref AS (SELECT * FROM bins)" in sql

    def test_giqldisjoin_sql_should_let_cte_shadow_registered_table(self):
        """Test that an in-query CTE shadows a registered table of the same name.

        Given:
            A reference name shared by a registered 1-based table and a CTE.
        When:
            Transpiling to SQL.
        Then:
            It should let the CTE shadow the table and leave canonical columns
            unshifted.
        """
        # Arrange & act
        sql = transpile(
            "WITH refs AS (SELECT 1) SELECT * FROM DISJOIN(features, reference := refs)",
            tables=[
                "features",
                Table("refs", coordinate_system="1based", interval_type="closed"),
            ],
        )

        # Assert
        assert "__giql_dj_ref AS (SELECT * FROM refs)" in sql
        assert '("start" - 1)' not in sql

    def test_giqldisjoin_sql_should_reject_literal_range_reference(self):
        """Test that a literal range reference is rejected.

        Given:
            A DISJOIN call whose reference is a literal genomic range.
        When:
            Transpiling to SQL.
        Then:
            It should raise a ValueError rejecting the reference form.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError, match="DISJOIN reference must be"):
            transpile(
                "SELECT * FROM DISJOIN(features, reference := 'chr1:1-9')",
                tables=["features"],
            )

    def test_giqldisjoin_sql_should_reject_reserved_prefix_reference(self):
        """Test that a reference name using the reserved prefix is rejected.

        Given:
            A DISJOIN reference name using the reserved __giql_dj_ prefix.
        When:
            Transpiling to SQL.
        Then:
            It should raise a ValueError naming the reserved prefix.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError, match="reserved"):
            transpile(
                "SELECT * FROM DISJOIN(features, reference := __giql_dj_ref)",
                tables=["features"],
            )

    def test_giqldisjoin_sql_should_reject_reserved_prefix_target(self):
        """Test that a target name using the reserved prefix is rejected.

        Given:
            A DISJOIN target table whose name uses the reserved __giql_dj_
            prefix.
        When:
            Transpiling to SQL.
        Then:
            It should raise a ValueError naming the reserved prefix.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError, match="reserved"):
            transpile(
                "SELECT * FROM DISJOIN(__giql_dj_tgt)",
                tables=["__giql_dj_tgt"],
            )

    def test_giqldisjoin_sql_should_reject_unknown_reference_name(self):
        """Test that a reference matching neither a table nor a CTE is rejected.

        Given:
            A DISJOIN reference that is neither a registered table nor a CTE.
        When:
            Transpiling to SQL.
        Then:
            It should raise a ValueError rejecting the unknown name.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError, match="neither a registered table"):
            transpile(
                "SELECT * FROM DISJOIN(features, reference := missing)",
                tables=["features"],
            )

    def test_giqldisjoin_sql_should_skip_exists_clause_when_reference_omitted(self):
        """Test that an omitted reference suppresses the coverage EXISTS clause.

        Given:
            A DISJOIN call with no reference argument (self-mode).
        When:
            Transpiling to SQL.
        Then:
            It should omit the coverage EXISTS subquery, which is provably
            always true when the reference equals the target set.
        """
        # Arrange & act
        sql = transpile("SELECT * FROM DISJOIN(features)", tables=["features"])

        # Assert
        assert "EXISTS (" not in sql

    def test_giqldisjoin_sql_should_skip_exists_clause_when_reference_resolves_to_target(
        self,
    ):
        """Test that a reference naming the target table suppresses EXISTS.

        Given:
            A DISJOIN call whose explicit reference is the target table name.
        When:
            Transpiling to SQL.
        Then:
            It should omit the coverage EXISTS subquery, since the reference
            resolves to the same registered relation as the target.
        """
        # Arrange & act
        sql = transpile(
            "SELECT * FROM DISJOIN(features, reference := features)",
            tables=["features"],
        )

        # Assert
        assert "EXISTS (" not in sql

    def test_giqldisjoin_sql_should_emit_exists_clause_when_reference_is_different_table(
        self,
    ):
        """Test that an explicit distinct reference keeps the EXISTS clause.

        Given:
            A DISJOIN call with an explicit reference table distinct from the
            target.
        When:
            Transpiling to SQL.
        Then:
            It should emit the coverage EXISTS subquery against the reference
            CTE.
        """
        # Arrange & act
        sql = transpile(
            "SELECT * FROM DISJOIN(features, reference := refs)",
            tables=["features", "refs"],
        )

        # Assert
        assert "EXISTS (" in sql

    def test_giqldisjoin_sql_should_emit_exists_clause_when_cte_shadows_target_name(
        self,
    ):
        """Test that a CTE shadowing the target name keeps the EXISTS clause.

        Given:
            A reference name shared by the target table and an in-query CTE.
        When:
            Transpiling to SQL.
        Then:
            It should emit the coverage EXISTS subquery, since the CTE may
            contain rows distinct from the registered target table.
        """
        # Arrange & act
        sql = transpile(
            "WITH features AS (SELECT 1) "
            "SELECT * FROM DISJOIN(features, reference := features)",
            tables=["features"],
        )

        # Assert
        assert "EXISTS (" in sql

    def test_giqldisjoin_sql_should_emit_exists_clause_when_reference_is_subquery(self):
        """Test that a subquery reference keeps the EXISTS clause.

        Given:
            A DISJOIN call whose reference is a subquery.
        When:
            Transpiling to SQL.
        Then:
            It should emit the coverage EXISTS subquery, since a subquery is
            conservatively treated as a relation distinct from the target.
        """
        # Arrange & act
        sql = transpile(
            "SELECT * FROM DISJOIN(features, reference := (SELECT * FROM features))",
            tables=["features"],
        )

        # Assert
        assert "EXISTS (" in sql

    def test_giqldisjoin_sql_should_skip_exists_clause_when_target_uses_one_based_closed_encoding(
        self,
    ):
        """Test that self-mode skips EXISTS even with non-default coordinates.

        Given:
            A self-mode DISJOIN against a target registered as 1-based closed.
        When:
            Transpiling to SQL.
        Then:
            It should omit the coverage EXISTS subquery and still emit the
            canonical 0-based shift on the target endpoints.
        """
        # Arrange & act
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

        # Assert
        assert "EXISTS (" not in sql
        assert '(t."start" - 1)' in sql

    def test_giqldisjoin_sql_should_skip_exists_clause_when_target_uses_custom_column_names(
        self,
    ):
        """Test that self-mode skips EXISTS cleanly under custom column names.

        Given:
            A self-mode DISJOIN against a target with custom genomic column
            names.
        When:
            Transpiling to SQL.
        Then:
            It should omit the coverage EXISTS subquery and leave no dangling
            reference-side qualifier (no ``r."seqid"`` / ``r."lo"`` /
            ``r."hi"``).
        """
        # Arrange & act
        sql = transpile(
            "SELECT * FROM DISJOIN(feats)",
            tables=[Table("feats", chrom_col="seqid", start_col="lo", end_col="hi")],
        )

        # Assert
        assert "EXISTS (" not in sql
        assert 'r."seqid"' not in sql
        assert 'r."lo"' not in sql
        assert 'r."hi"' not in sql

    def test_giqldisjoin_sql_should_skip_exists_clause_when_explicit_self_reference_uses_one_based_closed_encoding(
        self,
    ):
        """Test that explicit self-reference skips EXISTS under non-default config.

        Given:
            A DISJOIN call with an explicit ``reference := features`` where
            ``features`` is registered as 1-based closed.
        When:
            Transpiling to SQL.
        Then:
            It should omit the coverage EXISTS subquery and still emit the
            canonical 0-based shift on the breakpoint CTE endpoints.
        """
        # Arrange & act
        sql = transpile(
            "SELECT * FROM DISJOIN(features, reference := features)",
            tables=[
                Table(
                    "features",
                    coordinate_system="1based",
                    interval_type="closed",
                )
            ],
        )

        # Assert
        assert "EXISTS (" not in sql
        assert '("start" - 1)' in sql

    def test_giqldisjoin_sql_should_emit_one_exists_clause_when_query_mixes_self_and_distinct_reference_disjoins(
        self,
    ):
        """Test that the EXISTS skip is decided per DISJOIN call.

        Given:
            A query containing one self-mode DISJOIN and one distinct-reference
            DISJOIN in the same SELECT.
        When:
            Transpiling to SQL.
        Then:
            It should produce SQL containing exactly one ``EXISTS (``
            occurrence — only the distinct-reference call retains the
            coverage filter.
        """
        # Arrange & act
        sql = transpile(
            "SELECT * FROM DISJOIN(features) "
            "UNION ALL "
            "SELECT * FROM DISJOIN(features, reference := refs)",
            tables=["features", "refs"],
        )

        # Assert
        assert sql.count("EXISTS (") == 1

    def test_giqldisjoin_sql_should_emit_exists_clause_when_enclosing_cte_shadows_target_from_outer_scope(
        self,
    ):
        """Test that a CTE shadowing the target name from an outer scope keeps EXISTS.

        Given:
            A nested query where the outer WITH defines a CTE named ``features``
            and an inner subquery contains
            ``DISJOIN(features, reference := features)``, with ``features``
            also a registered table.
        When:
            Transpiling to SQL.
        Then:
            It should emit the coverage EXISTS subquery — the enclosing-CTE
            check follows the ancestor chain past the immediate WITH.
        """
        # Arrange & act
        sql = transpile(
            "WITH features AS (SELECT 1) "
            "SELECT * FROM ("
            "SELECT * FROM DISJOIN(features, reference := features)"
            ") AS sub",
            tables=["features"],
        )

        # Assert
        assert "EXISTS (" in sql
