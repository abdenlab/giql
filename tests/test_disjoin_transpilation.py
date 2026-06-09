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
        assert '__giql_dj_ref AS (SELECT * FROM "features")' in sql

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
        assert '__giql_dj_ref AS (SELECT * FROM "refs")' in sql

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
        assert '__giql_dj_ref AS (SELECT * FROM "bins")' in sql

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
        assert '__giql_dj_ref AS (SELECT * FROM "refs")' in sql
        assert '("start" - 1)' not in sql

    def test_giqldisjoin_sql_should_resolve_in_query_cte_target(self):
        """Test that a target naming an in-query CTE resolves to it.

        Given:
            A DISJOIN target naming a CTE defined in the outer query, with no
            registered table of the same name.
        When:
            Transpiling to SQL.
        Then:
            It should make the target CTE select from that CTE name and use
            the canonical chrom / start / end columns.
        """
        # Arrange & act
        sql = transpile(
            "WITH bins AS (SELECT 1) SELECT * FROM DISJOIN(bins)",
            tables=[],
        )

        # Assert
        assert '__giql_dj_tgt AS (SELECT * FROM "bins")' in sql
        assert 't."chrom"' in sql
        assert 't."start"' in sql
        assert 't."end"' in sql

    def test_giqldisjoin_sql_should_let_cte_target_shadow_registered_table(self):
        """Test that an in-query CTE target shadows a registered table.

        Given:
            A target name shared by a registered 1-based table and a CTE.
        When:
            Transpiling to SQL.
        Then:
            It should let the CTE shadow the table and leave canonical columns
            unshifted.
        """
        # Arrange & act
        sql = transpile(
            "WITH features AS (SELECT 1) SELECT * FROM DISJOIN(features)",
            tables=[
                Table("features", coordinate_system="1based", interval_type="closed"),
            ],
        )

        # Assert
        assert '__giql_dj_tgt AS (SELECT * FROM "features")' in sql
        assert '(t."start" - 1)' not in sql

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

    def test_giqldisjoin_sql_should_reject_unknown_target_name(self):
        """Test that a target matching neither a table nor a CTE is rejected.

        Given:
            A DISJOIN target that is neither a registered table nor a CTE.
        When:
            Transpiling to SQL.
        Then:
            It should raise a ValueError rejecting the unknown name.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError, match="neither a registered table"):
            transpile(
                "SELECT * FROM DISJOIN(missing)",
                tables=[],
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

    def test_giqldisjoin_sql_should_emit_exists_clause_when_target_is_cte_and_reference_is_omitted(
        self,
    ):
        """Test that EXISTS is preserved when an omitted reference defaults to a CTE target.

        Given:
            A self-mode DISJOIN whose target is an in-query CTE.
        When:
            Transpiling to SQL.
        Then:
            It should emit the coverage EXISTS subquery: the engine may
            inline the CTE body more than once and the two expansions are
            not provably identical (e.g. for volatile expressions), so the
            self-mode optimisation is unsafe for CTE targets.
        """
        # Arrange & act
        sql = transpile(
            "WITH bins AS (SELECT 1) SELECT * FROM DISJOIN(bins)",
            tables=[],
        )

        # Assert
        assert "EXISTS (" in sql

    def test_giqldisjoin_sql_should_emit_exists_clause_when_cte_target_has_distinct_reference(
        self,
    ):
        """Test that a distinct reference table keeps EXISTS for a CTE target.

        Given:
            A DISJOIN call whose target is an in-query CTE and whose explicit
            reference is a distinct registered table.
        When:
            Transpiling to SQL.
        Then:
            It should emit the coverage EXISTS subquery against the reference
            CTE.
        """
        # Arrange & act
        sql = transpile(
            "WITH bins AS (SELECT 1) SELECT * FROM DISJOIN(bins, reference := refs)",
            tables=["refs"],
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
        assert sql.count("EXISTS (") == 1

    def test_giqldisjoin_sql_should_succeed_when_outer_select_projects_cte_target_columns(
        self,
    ):
        """Test that an outer SELECT can project columns from a CTE target.

        Given:
            A DISJOIN whose target is a CTE and whose outer SELECT projects
            qualified columns from that CTE.
        When:
            Transpiling to SQL.
        Then:
            It should transpile cleanly and reference the canonical
            ``"chrom"`` / ``"start"`` / ``"end"`` columns without any
            coordinate-system arithmetic.
        """
        # Arrange & act
        sql = transpile(
            "WITH bins AS (SELECT 1) "
            'SELECT bins."chrom", bins."start", bins."end" FROM DISJOIN(bins)',
            tables=[],
        )

        # Assert
        assert '__giql_dj_tgt AS (SELECT * FROM "bins")' in sql
        assert 'bins."start" - 1' not in sql

    def test_giqldisjoin_sql_should_emit_canonical_output_columns_for_cte_target(self):
        """Test that DISJOIN over a CTE emits canonical 0-based half-open output.

        Given:
            A self-mode DISJOIN whose target is an in-query CTE with no
            registered table.
        When:
            Transpiling to SQL.
        Then:
            It should emit ``disjoin_start`` and ``disjoin_end`` projections
            free of any ``+ 1`` / ``- 1`` decanonicalization arithmetic.
        """
        # Arrange & act
        sql = transpile(
            "WITH bins AS (SELECT 1) SELECT * FROM DISJOIN(bins)",
            tables=[],
        )

        # Assert
        assert "s.seg_start AS disjoin_start" in sql
        assert "s.seg_end AS disjoin_end" in sql
        assert "seg_start + 1" not in sql
        assert "seg_start - 1" not in sql

    def test_giqldisjoin_sql_should_emit_canonical_output_when_cte_shadows_one_based_target(
        self,
    ):
        """Test that a CTE-shadowed 1-based target still emits canonical output.

        Given:
            A self-mode DISJOIN whose target is a CTE shadowing a registered
            1-based-closed table of the same name.
        When:
            Transpiling to SQL.
        Then:
            It should emit ``disjoin_start`` / ``disjoin_end`` without
            ``seg_start + 1`` decanonicalization — the registered table's
            config must not leak through the shadowed name.
        """
        # Arrange & act
        sql = transpile(
            "WITH features AS (SELECT 1) SELECT * FROM DISJOIN(features)",
            tables=[
                Table(
                    "features",
                    coordinate_system="1based",
                    interval_type="closed",
                ),
            ],
        )

        # Assert
        assert "seg_start + 1" not in sql

    def test_giqldisjoin_sql_should_accept_none_tables_argument_for_cte_target(self):
        """Test that ``tables=None`` is supported for a CTE-target DISJOIN.

        Given:
            A CTE-target DISJOIN transpiled with ``tables=None``.
        When:
            Transpiling to SQL.
        Then:
            It should succeed and emit the target CTE selection with canonical
            default columns.
        """
        # Arrange & act
        sql = transpile(
            "WITH bins AS (SELECT 1) SELECT * FROM DISJOIN(bins)",
            tables=None,
        )

        # Assert
        assert '__giql_dj_tgt AS (SELECT * FROM "bins")' in sql

    def test_giqldisjoin_sql_should_inline_subquery_reference_against_cte_target(self):
        """Test that a subquery reference composes with a CTE target.

        Given:
            A DISJOIN whose target is an in-query CTE and whose reference is
            a subquery against a registered table.
        When:
            Transpiling to SQL.
        Then:
            It should emit ``__giql_dj_tgt AS (SELECT * FROM \"bins\")``, alias
            the subquery as ``__giql_dj_rs``, and retain the coverage EXISTS
            clause.
        """
        # Arrange & act
        sql = transpile(
            "WITH bins AS (SELECT 1) "
            "SELECT * FROM DISJOIN(bins, reference := (SELECT * FROM refs))",
            tables=["refs"],
        )

        # Assert
        assert '__giql_dj_tgt AS (SELECT * FROM "bins")' in sql
        assert "AS __giql_dj_rs" in sql
        assert "EXISTS (" in sql

    def test_giqldisjoin_sql_should_resolve_cte_target_independent_of_unrelated_tables(
        self,
    ):
        """Test that an unrelated registered table does not interfere with CTE resolution.

        Given:
            A CTE target named ``bins`` and a registered, unrelated table
            ``refs`` not named anywhere in the DISJOIN call.
        When:
            Transpiling to SQL.
        Then:
            It should resolve the target to the CTE regardless of the other
            registered table's presence.
        """
        # Arrange & act
        sql = transpile(
            "WITH bins AS (SELECT * FROM refs) SELECT * FROM DISJOIN(bins)",
            tables=["refs"],
        )

        # Assert
        assert '__giql_dj_tgt AS (SELECT * FROM "bins")' in sql

    def test_giqldisjoin_sql_should_resolve_both_sides_to_cte_when_target_and_reference_share_shadowed_name(
        self,
    ):
        """Test that target and reference both resolve to a same-named CTE.

        Given:
            A CTE ``features`` and a 1-based-closed registered table
            ``features``; the DISJOIN names ``features`` on both sides.
        When:
            Transpiling to SQL.
        Then:
            It should select both target and reference from the CTE and emit
            no ``- 1`` coordinate-system arithmetic anywhere.
        """
        # Arrange & act
        sql = transpile(
            "WITH features AS (SELECT 1) "
            "SELECT * FROM DISJOIN(features, reference := features)",
            tables=[
                Table(
                    "features",
                    coordinate_system="1based",
                    interval_type="closed",
                ),
            ],
        )

        # Assert
        assert '__giql_dj_tgt AS (SELECT * FROM "features")' in sql
        assert '__giql_dj_ref AS (SELECT * FROM "features")' in sql
        assert '("start" - 1)' not in sql
        assert '(t."start" - 1)' not in sql

    def test_giqldisjoin_sql_should_skip_coordinate_shift_on_breakpoint_when_reference_is_shadowed_cte(
        self,
    ):
        """Test that a CTE-shadowed reference uses canonical breakpoint endpoints.

        Given:
            A 1-based-closed registered ``features`` table shadowed by a
            same-named CTE; DISJOIN names ``features`` on both sides.
        When:
            Transpiling to SQL.
        Then:
            The breakpoint CTE endpoints should contain no ``("start" - 1)``
            shift — the CTE is assumed canonical.
        """
        # Arrange & act
        sql = transpile(
            "WITH features AS (SELECT 1) "
            "SELECT * FROM DISJOIN(features, reference := features)",
            tables=[
                Table(
                    "features",
                    coordinate_system="1based",
                    interval_type="closed",
                ),
            ],
        )

        # Assert
        assert '("start" - 1)' not in sql

    def test_giqldisjoin_sql_should_skip_target_shift_when_target_is_shadowed_cte_and_reference_omitted(
        self,
    ):
        """Test that a CTE-shadowed target keeps canonical coordinates in self-mode.

        Given:
            A 1-based-closed registered ``features`` table shadowed by a
            same-named CTE; DISJOIN omits ``reference``.
        When:
            Transpiling to SQL.
        Then:
            The target endpoints should contain no ``(t."start" - 1)`` shift.
        """
        # Arrange & act
        sql = transpile(
            "WITH features AS (SELECT 1) SELECT * FROM DISJOIN(features)",
            tables=[
                Table(
                    "features",
                    coordinate_system="1based",
                    interval_type="closed",
                ),
            ],
        )

        # Assert
        assert '(t."start" - 1)' not in sql

    def test_giqldisjoin_sql_should_ignore_custom_columns_of_shadowed_registered_target(
        self,
    ):
        """Test that custom registered columns do not leak through a CTE target.

        Given:
            A registered ``features`` table with custom column names
            (``seqid`` / ``lo`` / ``hi``) shadowed by a same-named CTE.
        When:
            Transpiling a self-mode DISJOIN.
        Then:
            The generated SQL should reference canonical
            ``t."chrom"`` / ``t."start"`` / ``t."end"`` and contain no
            ``seqid`` / ``lo`` / ``hi`` references.
        """
        # Arrange & act
        sql = transpile(
            "WITH features AS (SELECT 1) SELECT * FROM DISJOIN(features)",
            tables=[
                Table("features", chrom_col="seqid", start_col="lo", end_col="hi"),
            ],
        )

        # Assert
        assert 't."chrom"' in sql
        assert 't."start"' in sql
        assert 't."end"' in sql
        assert "seqid" not in sql
        assert '"lo"' not in sql
        assert '"hi"' not in sql

    def test_giqldisjoin_sql_should_emit_canonical_output_when_cte_shadows_one_based_target_self_mode(
        self,
    ):
        """Test that a CTE-shadowed 1-based target emits canonical output columns.

        Given:
            A registered 1-based-closed ``features`` table shadowed by a
            same-named CTE in self-mode.
        When:
            Transpiling to SQL.
        Then:
            The output column aliases for ``disjoin_start`` / ``disjoin_end``
            should carry no ``+ 1`` arithmetic.
        """
        # Arrange & act
        sql = transpile(
            "WITH features AS (SELECT 1) SELECT * FROM DISJOIN(features)",
            tables=[
                Table(
                    "features",
                    coordinate_system="1based",
                    interval_type="closed",
                ),
            ],
        )

        # Assert
        assert "seg_start + 1" not in sql
        assert "seg_end + 1" not in sql

    def test_giqldisjoin_sql_should_prevent_all_registered_leakage_when_cte_shadows_target(
        self,
    ):
        """Test that a CTE shadow blocks every override of the registered target.

        Given:
            A registered ``features`` table combining custom column names and
            a 1-based-closed encoding, shadowed by a same-named CTE.
        When:
            Transpiling a self-mode DISJOIN.
        Then:
            The generated SQL should reference canonical columns, contain no
            ``- 1`` / ``+ 1`` shifts, and contain no ``seqid`` / ``lo`` /
            ``hi`` qualifiers anywhere.
        """
        # Arrange & act
        sql = transpile(
            "WITH features AS (SELECT 1) SELECT * FROM DISJOIN(features)",
            tables=[
                Table(
                    "features",
                    chrom_col="seqid",
                    start_col="lo",
                    end_col="hi",
                    coordinate_system="1based",
                    interval_type="closed",
                ),
            ],
        )

        # Assert
        assert '__giql_dj_tgt AS (SELECT * FROM "features")' in sql
        assert "seqid" not in sql
        assert '"lo"' not in sql
        assert '"hi"' not in sql
        assert "- 1" not in sql
        assert "+ 1" not in sql

    def test_giqldisjoin_sql_should_emit_one_exists_clause_when_explicit_self_reference_is_shadowed_cte(
        self,
    ):
        """Test that an explicit self-reference through a CTE emits exactly one EXISTS.

        Given:
            A registered ``features`` table shadowed by a same-named CTE; the
            DISJOIN names ``features`` for both target and reference.
        When:
            Transpiling to SQL.
        Then:
            Exactly one ``EXISTS (`` should appear — the CTE-vs-CTE path is
            conservatively treated as distinct, and target-side resolution
            does not regress.
        """
        # Arrange & act
        sql = transpile(
            "WITH features AS (SELECT 1) "
            "SELECT * FROM DISJOIN(features, reference := features)",
            tables=["features"],
        )

        # Assert
        assert sql.count("EXISTS (") == 1

    def test_giqldisjoin_sql_should_preserve_target_shift_when_reference_only_is_shadowed_cte(
        self,
    ):
        """Test asymmetric shadow on the reference side only.

        Given:
            Registered 1-based ``features`` and ``refs`` tables; a CTE
            ``refs`` shadows only the reference side.
        When:
            Transpiling a DISJOIN with ``reference := refs``.
        Then:
            The target side should keep its ``(t."start" - 1)`` shift while
            the breakpoint side (CTE) should contain no shift.
        """
        # Arrange & act
        sql = transpile(
            "WITH refs AS (SELECT 1) SELECT * FROM DISJOIN(features, reference := refs)",
            tables=[
                Table(
                    "features",
                    coordinate_system="1based",
                    interval_type="closed",
                ),
                Table(
                    "refs",
                    coordinate_system="1based",
                    interval_type="closed",
                ),
            ],
        )

        # Assert
        assert '(t."start" - 1)' in sql
        assert '__giql_dj_ref AS (SELECT * FROM "refs")' in sql
        # The breakpoint CTE pulls "start"/"end" without alias and would shift
        # them if the reference were the registered 1-based refs. The form
        # `("start" - 1) AS pos` (no `t.` qualifier) appears only when the
        # breakpoint side is non-canonical.
        assert '("start" - 1) AS pos' not in sql

    def test_giqldisjoin_sql_should_preserve_reference_shift_when_target_only_is_shadowed_cte(
        self,
    ):
        """Test asymmetric shadow on the target side only.

        Given:
            Registered 1-based ``features`` and ``refs`` tables; a CTE
            ``features`` shadows only the target side.
        When:
            Transpiling a DISJOIN with ``reference := refs``.
        Then:
            The target side should contain no shift, the breakpoint side
            should preserve its ``("start" - 1)`` shift, and EXISTS should
            be emitted.
        """
        # Arrange & act
        sql = transpile(
            "WITH features AS (SELECT 1) "
            "SELECT * FROM DISJOIN(features, reference := refs)",
            tables=[
                Table(
                    "features",
                    coordinate_system="1based",
                    interval_type="closed",
                ),
                Table(
                    "refs",
                    coordinate_system="1based",
                    interval_type="closed",
                ),
            ],
        )

        # Assert
        assert '(t."start" - 1)' not in sql
        assert '("start" - 1) AS pos' in sql
        assert "EXISTS (" in sql

    def test_giqldisjoin_sql_should_raise_with_disjoin_specific_message_for_unknown_target(
        self,
    ):
        """Test that the DISJOIN unknown-target error message is pinned.

        Given:
            A bare DISJOIN target name matching neither a CTE nor a registered
            table.
        When:
            Transpiling to SQL.
        Then:
            It should raise ``ValueError`` whose message contains
            ``DISJOIN target``, ``neither a registered table``, ``CTE
            defined in an enclosing WITH clause``, and a pointer to the
            public ``tables=`` parameter.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError) as excinfo:
            transpile("SELECT * FROM DISJOIN(missing)", tables=[])
        message = str(excinfo.value)
        assert "DISJOIN target" in message
        assert "neither a registered table" in message
        assert "CTE defined in an enclosing WITH clause" in message
        assert "tables=" in message

    def test_giqldisjoin_sql_should_reject_cte_target_with_reserved_prefix(self):
        """Test that a CTE target using the reserved prefix is rejected.

        Given:
            A CTE named with the reserved ``__giql_dj_`` prefix used as a
            DISJOIN target.
        When:
            Transpiling to SQL.
        Then:
            It should raise ``ValueError`` mentioning ``reserved`` and the
            reserved prefix string.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError) as excinfo:
            transpile(
                "WITH __giql_dj_tgt AS (SELECT 1) SELECT * FROM DISJOIN(__giql_dj_tgt)",
                tables=[],
            )
        message = str(excinfo.value)
        assert "reserved" in message
        assert "__giql_dj_" in message

    def test_giqldisjoin_sql_should_reject_reserved_prefix_even_with_registered_table(
        self,
    ):
        """Test that the reserved-prefix check fires even when a registered table shadows the CTE.

        Given:
            A CTE named ``__giql_dj_tgt`` that also has a same-named
            registered table.
        When:
            Transpiling to SQL.
        Then:
            It should raise ``ValueError`` matching ``reserved`` — the prefix
            guard runs post-resolution regardless of which branch resolved
            the name.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError, match="reserved"):
            transpile(
                "WITH __giql_dj_tgt AS (SELECT 1) SELECT * FROM DISJOIN(__giql_dj_tgt)",
                tables=["__giql_dj_tgt"],
            )

    def test_giqldisjoin_sql_should_reject_subquery_target(self):
        """Test that a subquery in target position is rejected.

        Given:
            A DISJOIN call whose target is a ``(SELECT ...)`` subquery.
        When:
            Transpiling to SQL.
        Then:
            It should raise ``ValueError`` matching ``must be a bare table
            or CTE name`` — subquery targets remain unsupported.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError, match="must be a bare table or CTE name"):
            transpile(
                "SELECT * FROM DISJOIN((SELECT * FROM features))",
                tables=["features"],
            )

    def test_giqldisjoin_sql_should_reject_nested_operator_target(self):
        """Test that a nested DISJOIN in target position is rejected.

        Given:
            A DISJOIN call whose target is another DISJOIN call.
        When:
            Transpiling to SQL.
        Then:
            It should raise ``ValueError`` matching ``must be a bare table
            or CTE name`` — nested-operator targets remain unsupported.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError, match="must be a bare table or CTE name"):
            transpile(
                "SELECT * FROM DISJOIN(DISJOIN(features))",
                tables=["features"],
            )

    def test_giqldisjoin_sql_should_resolve_named_cte_among_multiple_siblings(self):
        """Test that DISJOIN picks the named sibling CTE in a multi-CTE WITH.

        Given:
            A WITH clause defining two sibling CTEs ``a`` and ``b``; DISJOIN
            targets ``a``.
        When:
            Transpiling to SQL.
        Then:
            The target CTE selects from ``a`` only — ``b`` does not leak.
        """
        # Arrange & act
        sql = transpile(
            "WITH a AS (SELECT 1), b AS (SELECT 2) SELECT * FROM DISJOIN(a)",
            tables=[],
        )

        # Assert
        assert '__giql_dj_tgt AS (SELECT * FROM "a")' in sql
        assert '__giql_dj_tgt AS (SELECT * FROM "b")' not in sql

    def test_giqldisjoin_sql_should_fall_back_to_registered_table_when_cte_name_does_not_match(
        self,
    ):
        """Test that a non-matching CTE name does not block table fallback.

        Given:
            A WITH clause defining an unrelated CTE ``unrelated`` and a
            registered ``features`` table; DISJOIN targets ``features``.
        When:
            Transpiling to SQL.
        Then:
            The target CTE selects from the registered ``features``, not from
            ``unrelated``.
        """
        # Arrange & act
        sql = transpile(
            "WITH unrelated AS (SELECT 1) SELECT * FROM DISJOIN(features)",
            tables=["features"],
        )

        # Assert
        assert '__giql_dj_tgt AS (SELECT * FROM "features")' in sql
        assert '__giql_dj_tgt AS (SELECT * FROM "unrelated")' not in sql

    def test_giqldisjoin_sql_should_resolve_cte_target_from_derived_subquery(self):
        """Test that ancestor walk reaches a WITH past a derived subquery.

        Given:
            A query where the outer WITH defines ``x`` and DISJOIN sits inside
            a derived subquery.
        When:
            Transpiling to SQL.
        Then:
            The target CTE selects from the outer WITH's ``x``.
        """
        # Arrange & act
        sql = transpile(
            "WITH x AS (SELECT 1) SELECT * FROM (SELECT * FROM DISJOIN(x)) AS sub",
            tables=[],
        )

        # Assert
        assert '__giql_dj_tgt AS (SELECT * FROM "x")' in sql

    def test_giqldisjoin_sql_should_resolve_cte_target_from_join_right_side(self):
        """Test that ancestor walk reaches a WITH past a JOIN.

        Given:
            A WITH clause defining ``x`` and a DISJOIN sitting in the RHS of
            a JOIN.
        When:
            Transpiling to SQL.
        Then:
            The target CTE selects from the outer WITH's ``x``.
        """
        # Arrange & act
        sql = transpile(
            "WITH x AS (SELECT 1) "
            "SELECT * FROM features "
            "JOIN DISJOIN(x) AS d ON features.chrom = d.disjoin_chrom",
            tables=["features"],
        )

        # Assert
        assert '__giql_dj_tgt AS (SELECT * FROM "x")' in sql

    def test_giqldisjoin_sql_should_resolve_cte_target_when_with_attaches_to_union(self):
        """Test that ancestor walk picks up CTEs on a ``UNION``'s WITH clause.

        Given:
            A top-level ``WITH x AS (...) SELECT * FROM DISJOIN(x) UNION ALL
            SELECT * FROM features`` where the WITH attaches to the UNION
            node, not directly to the first SELECT.
        When:
            Transpiling to SQL.
        Then:
            DISJOIN resolves ``x`` to the outer-WITH CTE — regression guard
            for ``_enclosing_cte_names`` recognising non-Select ancestors.
        """
        # Arrange & act
        sql = transpile(
            "WITH x AS (SELECT 1) "
            "SELECT disjoin_chrom AS c FROM DISJOIN(x) "
            "UNION ALL SELECT chrom AS c FROM features",
            tables=["features"],
        )

        # Assert
        assert '__giql_dj_tgt AS (SELECT * FROM "x")' in sql

    def test_giqldisjoin_sql_should_reject_target_naming_inner_only_cte_sibling(self):
        """Test that an inner-subquery CTE is not visible to an outer DISJOIN.

        Given:
            A query where an inner subquery defines a CTE ``x`` and the outer
            DISJOIN names ``x`` — no registered table provides ``x``.
        When:
            Transpiling to SQL.
        Then:
            It should raise ``ValueError`` matching ``neither a registered
            table`` — negative scoping: ancestor walk does not descend into
            sibling subtrees.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError, match="neither a registered table"):
            transpile(
                "SELECT * FROM DISJOIN(x) "
                "WHERE EXISTS (WITH x AS (SELECT 1) SELECT 1 FROM x)",
                tables=[],
            )

    def test_giqldisjoin_sql_should_succeed_when_inner_with_redeclares_outer_cte_name(
        self,
    ):
        """Test that an inner WITH redeclaring an outer CTE name still resolves.

        Given:
            An outer ``WITH x AS (...)`` and an inner subquery that re-declares
            ``x`` inside which DISJOIN runs.
        When:
            Transpiling to SQL.
        Then:
            It should transpile successfully and emit ``__giql_dj_tgt AS
            (SELECT * FROM x)``.
        """
        # Arrange & act
        sql = transpile(
            "WITH x AS (SELECT 1) "
            "SELECT * FROM (WITH x AS (SELECT 2) SELECT * FROM DISJOIN(x)) AS sub",
            tables=[],
        )

        # Assert
        assert '__giql_dj_tgt AS (SELECT * FROM "x")' in sql

    def test_giqldisjoin_sql_should_match_quoted_cte_alias_with_bare_target_name(self):
        """Test that a quoted CTE alias matches a bare DISJOIN target.

        Given:
            A WITH clause defining ``"X"`` and a bare DISJOIN target ``X``.
        When:
            Transpiling to SQL.
        Then:
            It should transpile successfully — the quoted alias is collected
            without quotes, matching the bare target name.
        """
        # Arrange & act
        sql = transpile(
            'WITH "X" AS (SELECT 1) SELECT * FROM DISJOIN(X)',
            tables=[],
        )

        # Assert
        assert '__giql_dj_tgt AS (SELECT * FROM "X")' in sql

    def test_giqldisjoin_sql_should_match_quoted_cte_alias_with_quoted_target_name(
        self,
    ):
        """Test that a quoted DISJOIN target matches a quoted CTE alias.

        Given:
            A WITH clause defining ``"X"`` and a quoted DISJOIN target
            ``"X"``.
        When:
            Transpiling to SQL.
        Then:
            It should transpile successfully — target-name extraction
            normalises the quoted identifier to match the unquoted alias
            (regression guard for the quoted-target-name normalisation).
        """
        # Arrange & act
        sql = transpile(
            'WITH "X" AS (SELECT 1) SELECT * FROM DISJOIN("X")',
            tables=[],
        )

        # Assert
        assert '__giql_dj_tgt AS (SELECT * FROM "X")' in sql

    def test_giqldisjoin_sql_should_resolve_distinct_cte_pair_as_target_and_reference(
        self,
    ):
        """Test that two distinct in-query CTEs serve as target and reference.

        Given:
            A WITH clause defining two CTEs ``x`` and ``y``; DISJOIN targets
            ``x`` and references ``y``.
        When:
            Transpiling to SQL.
        Then:
            Both ``__giql_dj_tgt AS (SELECT * FROM \"x\")`` and ``__giql_dj_ref
            AS (SELECT * FROM y)`` should appear; exactly one ``EXISTS (``
            is emitted; the breakpoint CTE references canonical columns.
        """
        # Arrange & act
        sql = transpile(
            "WITH x AS (SELECT 1), y AS (SELECT 1) "
            "SELECT * FROM DISJOIN(x, reference := y)",
            tables=[],
        )

        # Assert
        assert '__giql_dj_tgt AS (SELECT * FROM "x")' in sql
        assert '__giql_dj_ref AS (SELECT * FROM "y")' in sql
        assert sql.count("EXISTS (") == 1
        assert '"chrom"' in sql
        assert '"start"' in sql
        assert '"end"' in sql

    def test_giqldisjoin_sql_should_emit_exists_when_same_cte_used_as_target_and_reference(
        self,
    ):
        """Test that a CTE-to-itself reference is conservatively distinct.

        Given:
            A CTE ``x`` used as both target and reference (with no registered
            table).
        When:
            Transpiling to SQL.
        Then:
            It should emit exactly one ``EXISTS (`` — the CTE branch
            conservatively treats the reference as distinct
            (``is_self_reference=False``).
        """
        # Arrange & act
        sql = transpile(
            "WITH x AS (SELECT 1) SELECT * FROM DISJOIN(x, reference := x)",
            tables=[],
        )

        # Assert
        assert '__giql_dj_tgt AS (SELECT * FROM "x")' in sql
        assert '__giql_dj_ref AS (SELECT * FROM "x")' in sql
        assert sql.count("EXISTS (") == 1

    def test_giqldisjoin_sql_should_use_custom_columns_of_registered_reference_with_cte_target(
        self,
    ):
        """Test that a CTE target + custom-column registered reference resolves cleanly.

        Given:
            A CTE target ``bins`` and a registered ``refs`` with custom
            column names ``seqid`` / ``lo`` / ``hi``.
        When:
            Transpiling to SQL.
        Then:
            The target side should use canonical columns; the reference side
            should use ``r."seqid"`` / ``r."lo"`` / ``r."hi"``.
        """
        # Arrange & act
        sql = transpile(
            "WITH bins AS (SELECT 1) SELECT * FROM DISJOIN(bins, reference := refs)",
            tables=[
                Table("refs", chrom_col="seqid", start_col="lo", end_col="hi"),
            ],
        )

        # Assert
        assert 't."chrom"' in sql
        assert 'r."seqid"' in sql
        assert 'r."lo"' in sql
        assert 'r."hi"' in sql
        assert "EXISTS (" in sql

    def test_giqldisjoin_sql_should_canonicalize_only_reference_when_target_is_cte(self):
        """Test that a 1-based registered reference is canonicalized; the CTE target is not.

        Given:
            A CTE target ``bins`` and a registered 1-based-closed ``refs``.
        When:
            Transpiling to SQL.
        Then:
            The breakpoint CTE should emit ``("start" - 1)`` for the
            reference; the target side should carry no shift; EXISTS is
            emitted.
        """
        # Arrange & act
        sql = transpile(
            "WITH bins AS (SELECT 1) SELECT * FROM DISJOIN(bins, reference := refs)",
            tables=[
                Table(
                    "refs",
                    coordinate_system="1based",
                    interval_type="closed",
                ),
            ],
        )

        # Assert
        assert '(t."start" - 1)' not in sql
        assert '("start" - 1) AS pos' in sql
        assert "EXISTS (" in sql

    def test_giqldisjoin_sql_should_emit_exists_when_omitted_reference_under_cte_shadow(
        self,
    ):
        """Test that omitted-reference under CTE shadowing preserves EXISTS.

        Given:
            A CTE ``features`` shadowing a registered ``features`` table;
            DISJOIN omits ``reference``.
        When:
            Transpiling to SQL.
        Then:
            EXISTS should be emitted: the CTE shadows the registered table,
            so target/reference are not provably the same registered
            relation and the optimisation cannot fire safely.
        """
        # Arrange & act
        sql = transpile(
            "WITH features AS (SELECT 1) SELECT * FROM DISJOIN(features)",
            tables=["features"],
        )

        # Assert
        assert "EXISTS (" in sql

    def test_giqldisjoin_sql_should_reject_literal_range_reference_with_cte_target(
        self,
    ):
        """Test that literal-range references remain rejected for CTE targets.

        Given:
            A CTE target and a literal-range ``reference``.
        When:
            Transpiling to SQL.
        Then:
            It should raise ``ValueError`` matching ``DISJOIN reference must
            be``.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError, match="DISJOIN reference must be"):
            transpile(
                "WITH x AS (SELECT 1) SELECT * FROM DISJOIN(x, reference := 'chr1:1-9')",
                tables=[],
            )

    def test_giqldisjoin_sql_should_reject_reserved_prefix_reference_with_cte_target(
        self,
    ):
        """Test that reserved-prefix references remain rejected for CTE targets.

        Given:
            A CTE target and a reference whose name uses the reserved
            ``__giql_dj_`` prefix.
        When:
            Transpiling to SQL.
        Then:
            It should raise ``ValueError`` matching ``reserved``.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError, match="reserved"):
            transpile(
                "WITH x AS (SELECT 1) "
                "SELECT * FROM DISJOIN(x, reference := __giql_dj_ref)",
                tables=[],
            )

    def test_giqldisjoin_sql_should_reject_qualified_reference(self):
        """Test that a qualified reference identifier is rejected.

        Given:
            A DISJOIN call whose reference is a qualified name like
            ``other.t`` rather than a bare identifier.
        When:
            Transpiling to SQL.
        Then:
            It should raise ``ValueError`` matching ``qualified`` — the bare
            identifier required by the resolver was not provided, and silently
            dropping the qualifier prefix masks the user's intent.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError, match="qualified"):
            transpile(
                "SELECT * FROM DISJOIN(features, reference := other.t)",
                tables=["features"],
            )

    def test_giqldisjoin_sql_should_reject_enclosing_cte_with_reserved_prefix(self):
        """Test that an enclosing CTE name colliding with the internal prefix is rejected.

        Given:
            A query whose enclosing ``WITH`` defines a CTE whose name starts
            with ``__giql_dj_``, alongside an ordinary DISJOIN-target CTE.
        When:
            Transpiling to SQL.
        Then:
            It should raise ``ValueError`` whose message names the offending
            CTE — pre-empting the opaque duplicate-CTE error the engine would
            otherwise produce once DISJOIN's internal ``__giql_dj_*`` CTEs
            are emitted alongside it.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError, match="__giql_dj_tgt"):
            transpile(
                "WITH __giql_dj_tgt AS (SELECT 1), x AS (SELECT 1) "
                "SELECT * FROM DISJOIN(x)",
                tables=[],
            )
