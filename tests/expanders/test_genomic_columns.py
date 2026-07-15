"""Genomic-column resolution for CLUSTER/MERGE over a derived-table FROM (#164).

CLUSTER and MERGE derive their partition/order/aggregation columns from the
enclosing FROM relation through the shared CLUSTER/MERGE resolution toolkit.
These tests drive the public ``transpile`` API and pin the resolution the #164
fix newly performs: tracing the genomic columns *through* a derived table or CTE
to the underlying registered mapping, honoring an explicit projection that
re-exposes the canonical columns, and raising a clear ``ValueError`` when the
source is opaque or combines multiple relations through a join — rather than
silently defaulting
to the canonical ``chrom`` / ``start`` / ``end`` names and emitting SQL that
references columns that do not exist.

The behavior is operator-neutral (both expanders call the same resolver), so the
resolver branches are exercised through CLUSTER, with dedicated MERGE tests
pinning that the shared resolver reaches MERGE and composes with its GROUP BY /
aggregation layer.
"""

import pytest

from giql.table import Table
from giql.transpile import transpile


def _query(op, source):
    """Build a CLUSTER or MERGE query over *source* for operator-parametrized tests."""
    if op == "CLUSTER":
        return f"SELECT *, CLUSTER(interval) AS cid FROM {source}"
    return f"SELECT MERGE(interval) FROM {source}"


@pytest.fixture
def custom():
    """A ``regions`` table declaring non-canonical genomic column names."""
    return Table("regions", chrom_col="ch", start_col="s", end_col="e")


class TestGenomicColumnResolution:
    """Resolution of genomic columns for CLUSTER/MERGE over a derived FROM (#164)."""

    def test_transpile_should_emit_custom_window_when_cluster_over_derived_table(
        self, custom
    ):
        """Test a CLUSTER over a ``SELECT *`` derived table uses the custom columns.

        Given:
            A CLUSTER whose FROM is a ``SELECT *`` subquery over a table declaring
            custom column names.
        When:
            Transpiling the query.
        Then:
            The window should partition/order/LAG on the underlying custom columns,
            not the canonical defaults — resolution traces through the subquery
            instead of silently defaulting (the #164 fix).
        """
        # Arrange
        query = "SELECT *, CLUSTER(interval) AS cid FROM (SELECT * FROM regions) AS sub"

        # Act
        sql = transpile(query, tables=[custom])

        # Assert
        assert 'PARTITION BY "ch" ORDER BY "s" NULLS LAST' in sql
        assert 'MAX("e")' in sql
        assert '"chrom"' not in sql and '"start"' not in sql and '"end"' not in sql

    def test_transpile_should_emit_custom_aggregation_when_merge_over_derived_table(
        self, custom
    ):
        """Test a MERGE over a ``SELECT *`` derived table uses the custom columns.

        Given:
            A MERGE whose FROM is a ``SELECT *`` subquery over a table declaring
            custom column names.
        When:
            Transpiling the query.
        Then:
            The clustered aggregation should group and MIN/MAX on the underlying
            custom columns — proving the shared resolver reaches MERGE and composes
            with its GROUP BY / aggregation layer.
        """
        # Arrange
        query = "SELECT MERGE(interval) FROM (SELECT * FROM regions) AS sub"

        # Act
        sql = transpile(query, tables=[custom])

        # Assert
        assert 'MIN("s") AS s' in sql and 'MAX("e") AS e' in sql
        assert 'GROUP BY "ch"' in sql
        assert '"chrom"' not in sql and '"start"' not in sql and '"end"' not in sql

    def test_transpile_should_resolve_custom_columns_when_from_is_qualified_star(
        self, custom
    ):
        """Test resolution through a derived table projecting a qualified star.

        Given:
            A CLUSTER over a derived table whose projection is ``alias.*`` rather
            than a bare ``*``.
        When:
            Transpiling the query.
        Then:
            The qualified star should pass its columns through, resolving to the
            underlying custom mapping.
        """
        # Arrange
        query = (
            "SELECT *, CLUSTER(interval) AS cid "
            "FROM (SELECT sub.* FROM regions sub) AS d"
        )

        # Act
        sql = transpile(query, tables=[custom])

        # Assert
        assert 'PARTITION BY "ch" ORDER BY "s" NULLS LAST' in sql
        assert '"chrom"' not in sql

    def test_transpile_should_resolve_custom_columns_when_star_mixed_with_projection(
        self, custom
    ):
        """Test a star mixed with an extra projection still passes columns through.

        Given:
            A CLUSTER over a derived table projecting ``alias.*`` alongside an extra
            computed column.
        When:
            Transpiling the query.
        Then:
            The presence of the star makes the source a passthrough, so it resolves
            to the underlying custom mapping.
        """
        # Arrange
        query = (
            "SELECT *, CLUSTER(interval) AS cid "
            "FROM (SELECT r.*, 1 AS extra FROM regions r) AS sub"
        )

        # Act
        sql = transpile(query, tables=[custom])

        # Assert
        assert 'PARTITION BY "ch" ORDER BY "s" NULLS LAST' in sql
        assert '"chrom"' not in sql

    def test_transpile_should_resolve_custom_columns_when_from_is_cte(self, custom):
        """Test resolution through a pass-through CTE to the underlying table.

        Given:
            A CLUSTER whose FROM references a ``SELECT *`` CTE over a custom-mapped
            table.
        When:
            Transpiling the query.
        Then:
            It should resolve through the CTE to the custom columns, and the
            enclosing WITH clause is preserved so the rewritten SQL stays
            executable (#174 — WITH preservation is pinned in the CLUSTER/MERGE
            expander suites; asserted here to confirm resolution and preservation
            compose).
        """
        # Arrange
        query = (
            "WITH sub AS (SELECT * FROM regions) "
            "SELECT *, CLUSTER(interval) AS cid FROM sub"
        )

        # Act
        sql = transpile(query, tables=[custom])

        # Assert
        assert 'PARTITION BY "ch" ORDER BY "s" NULLS LAST' in sql
        assert '"chrom"' not in sql
        assert sql.startswith("WITH sub AS (SELECT * FROM regions)")

    def test_transpile_should_resolve_custom_columns_when_from_is_nested_derived_table(
        self, custom
    ):
        """Test resolution recurses through nested derived tables.

        Given:
            A CLUSTER whose FROM is a ``SELECT *`` derived table wrapping another
            ``SELECT *`` derived table over a custom-mapped source.
        When:
            Transpiling the query.
        Then:
            Resolution should recurse through both levels to the custom columns.
        """
        # Arrange
        query = (
            "SELECT *, CLUSTER(interval) AS cid "
            "FROM (SELECT * FROM (SELECT * FROM regions) a) b"
        )

        # Act
        sql = transpile(query, tables=[custom])

        # Assert
        assert 'PARTITION BY "ch" ORDER BY "s" NULLS LAST' in sql
        assert '"chrom"' not in sql

    @pytest.mark.parametrize("setop", ["UNION ALL", "INTERSECT", "EXCEPT"])
    def test_transpile_should_resolve_left_arm_when_from_is_set_operation_derived_table(
        self, setop, custom
    ):
        """Test a set-operation-bodied derived table resolves via its left arm.

        Given:
            A CLUSTER over a derived table whose body combines a custom-mapped
            ``regions`` (left arm) with a differently-mapped ``regions2`` (right arm)
            through UNION ALL / INTERSECT / EXCEPT.
        When:
            Transpiling the query.
        Then:
            It should not raise and should resolve to the LEFT arm's mapping — the
            left table's physical columns appear in the emitted SQL and the right
            arm's do not, proving the left arm alone determines the resolution across
            every set-operation keyword (INTERSECT / EXCEPT no longer over-raise).
        """
        # Arrange
        regions2 = Table("regions2", chrom_col="seqid", start_col="lo", end_col="hi")
        query = (
            "SELECT *, CLUSTER(interval) AS cid "
            f"FROM (SELECT * FROM regions {setop} SELECT * FROM regions2) AS sub"
        )

        # Act
        sql = transpile(query, tables=[custom, regions2])

        # Assert
        assert 'PARTITION BY "ch" ORDER BY "s" NULLS LAST' in sql
        assert '"seqid"' not in sql and '"lo"' not in sql and '"hi"' not in sql

    def test_transpile_should_resolve_custom_columns_when_cte_defined_on_ancestor(
        self, custom
    ):
        """Test a CTE defined on an ancestor SELECT is found by the resolver.

        Given:
            A MERGE nested in a subquery whose FROM references a CTE defined on the
            enclosing root SELECT.
        When:
            Transpiling the query.
        Then:
            The resolver should walk up to the ancestor's WITH clause and resolve
            the custom columns.
        """
        # Arrange
        query = (
            "WITH c AS (SELECT * FROM regions) "
            "SELECT * FROM (SELECT MERGE(interval) FROM c) z"
        )

        # Act
        sql = transpile(query, tables=[custom])

        # Assert
        assert 'GROUP BY "ch"' in sql
        assert '"chrom"' not in sql

    def test_transpile_should_prefer_cte_body_over_registered_table_of_same_name(self):
        """Test a CTE shadowing a registered table resolves through the CTE body.

        Given:
            A CLUSTER over a CTE named identically to a registered table, whose body
            selects from a *different* registered table with its own mapping.
        When:
            Transpiling the query.
        Then:
            The CTE body should win (resolution checks CTEs before registered
            tables), so the columns come from the CTE's underlying source, not the
            shadowed same-named table.
        """
        # Arrange
        regions = Table("regions", chrom_col="ch", start_col="s", end_col="e")
        realsrc = Table("realsrc", chrom_col="rc", start_col="rs", end_col="re_")
        query = (
            "WITH regions AS (SELECT * FROM realsrc) "
            "SELECT *, CLUSTER(interval) AS cid FROM regions"
        )

        # Act
        sql = transpile(query, tables=[regions, realsrc])

        # Assert
        assert 'PARTITION BY "rc" ORDER BY "rs" NULLS LAST' in sql
        assert 'MAX("re_")' in sql
        assert '"ch"' not in sql

    def test_transpile_should_resolve_partial_custom_mapping_through_derived_table(self):
        """Test a partial custom mapping resolves field-by-field through a subquery.

        Given:
            A CLUSTER over a ``SELECT *`` derived table whose source customizes only
            the chromosome column and leaves start/end at their defaults.
        When:
            Transpiling the query.
        Then:
            The window should mix the custom chromosome with the default start/end,
            exactly as the mapping declares.
        """
        # Arrange
        regions = Table("regions", chrom_col="ch")
        query = "SELECT *, CLUSTER(interval) AS cid FROM (SELECT * FROM regions) AS sub"

        # Act
        sql = transpile(query, tables=[regions])

        # Assert
        assert 'PARTITION BY "ch" ORDER BY "start" NULLS LAST' in sql
        assert 'MAX("end")' in sql

    def test_transpile_should_partition_by_custom_strand_when_stranded_over_derived(
        self,
    ):
        """Test a stranded CLUSTER resolves the custom strand column through a subquery.

        Given:
            A stranded CLUSTER over a ``SELECT *`` derived table whose source
            declares a custom strand column.
        When:
            Transpiling the query.
        Then:
            The window partition should include the custom strand column.
        """
        # Arrange
        regions = Table(
            "regions", chrom_col="ch", start_col="s", end_col="e", strand_col="st"
        )
        query = (
            "SELECT *, CLUSTER(interval, stranded := true) AS cid "
            "FROM (SELECT * FROM regions) AS sub"
        )

        # Act
        sql = transpile(query, tables=[regions])

        # Assert
        assert 'PARTITION BY "ch", "st" ORDER BY "s" NULLS LAST' in sql

    def test_transpile_should_use_default_strand_when_table_declares_none_over_derived(
        self,
    ):
        """Test a stranded CLUSTER falls back to the default strand column.

        Given:
            A stranded CLUSTER over a ``SELECT *`` derived table whose source
            declares no strand column.
        When:
            Transpiling the query.
        Then:
            The window partition should fall back to the canonical ``strand``
            column alongside the custom chromosome.
        """
        # Arrange
        regions = Table(
            "regions", chrom_col="ch", start_col="s", end_col="e", strand_col=None
        )
        query = (
            "SELECT *, CLUSTER(interval, stranded := true) AS cid "
            "FROM (SELECT * FROM regions) AS sub"
        )

        # Act
        sql = transpile(query, tables=[regions])

        # Assert
        assert 'PARTITION BY "ch", "strand" ORDER BY "s" NULLS LAST' in sql

    @pytest.mark.parametrize("op", ["CLUSTER", "MERGE"])
    def test_transpile_should_resolve_canonical_when_derived_reexposes_canonical_names(
        self, op, custom
    ):
        """Test a derived table re-aliasing to canonical names resolves to canonical.

        Given:
            A CLUSTER/MERGE over a derived table whose explicit projection aliases
            the custom columns back to the canonical ``chrom`` / ``start`` / ``end``
            names.
        When:
            Transpiling the query.
        Then:
            It should not raise and should resolve to the canonical columns the
            subquery re-exposes — emitting the full canonical window/group shape
            while the underlying custom columns never leak.
        """
        # Arrange
        source = "(SELECT ch AS chrom, s AS start, e AS end FROM regions) AS sub"
        query = _query(op, source)

        # Act
        sql = transpile(query, tables=[custom])

        # Assert
        if op == "CLUSTER":
            assert 'PARTITION BY "chrom" ORDER BY "start"' in sql
        else:
            assert 'GROUP BY "chrom"' in sql
        assert '"ch"' not in sql

    def test_transpile_should_fall_back_to_canonical_when_derived_table_unregistered(
        self,
    ):
        """Test a derived table over an unregistered source defaults to canonical.

        Given:
            A CLUSTER over a ``SELECT *`` derived table whose source table is not
            registered.
        When:
            Transpiling the query.
        Then:
            It should not raise and should fall back to the canonical columns, as an
            unregistered relation is assumed canonical.
        """
        # Arrange
        query = "SELECT *, CLUSTER(interval) AS cid FROM (SELECT * FROM unknown) AS sub"

        # Act
        sql = transpile(query, tables=[])

        # Assert
        assert 'PARTITION BY "chrom" ORDER BY "start" NULLS LAST' in sql

    def test_transpile_should_fall_back_to_canonical_when_no_from_clause(self):
        """Test a CLUSTER with no FROM clause defaults to the canonical columns.

        Given:
            A CLUSTER whose enclosing SELECT has no FROM clause.
        When:
            Transpiling the query.
        Then:
            It should not raise and should fall back to the canonical columns, as
            there is no source relation to resolve.
        """
        # Arrange
        query = "SELECT CLUSTER(interval) AS cid"

        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        assert 'PARTITION BY "chrom" ORDER BY "start" NULLS LAST' in sql


class TestGenomicColumnResolutionErrors:
    """CLUSTER/MERGE raise when a derived FROM cannot be resolved (#164)."""

    @pytest.mark.parametrize("op", ["CLUSTER", "MERGE"])
    def test_transpile_should_raise_when_derived_table_hides_genomic_columns(
        self, op, custom
    ):
        """Test an opaque derived-table projection raises instead of defaulting.

        Given:
            A CLUSTER/MERGE over a derived table whose explicit projection exposes
            only non-canonical column names.
        When:
            Transpiling the query.
        Then:
            It should raise a ValueError explaining the source neither exposes the
            canonical columns nor resolves to a registered table — rather than
            silently emitting SQL that references non-existent columns.
        """
        # Arrange
        source = "(SELECT ch, s, e FROM regions) AS sub"
        query = _query(op, source)

        # Act & assert
        with pytest.raises(ValueError, match="neither exposes the canonical"):
            transpile(query, tables=[custom])

    def test_transpile_should_raise_when_from_is_a_values_relation(self):
        """Test a VALUES-relation FROM raises as an untraceable source.

        Given:
            A CLUSTER whose FROM is a derived table wrapping a VALUES clause.
        When:
            Transpiling the query.
        Then:
            It should raise a ValueError, as a VALUES body exposes no traceable
            genomic columns.
        """
        # Arrange
        query = (
            "SELECT *, CLUSTER(interval) AS cid FROM (VALUES (1, 2, 3)) AS v(a, b, c)"
        )

        # Act & assert
        with pytest.raises(ValueError, match="neither exposes the canonical"):
            transpile(query, tables=[])

    def test_transpile_should_raise_when_from_is_a_table_valued_function(self):
        """Test a table-valued-function FROM raises as an untraceable source.

        Given:
            A CLUSTER whose FROM is a table-valued function (``UNNEST``).
        When:
            Transpiling the query.
        Then:
            It should raise a ValueError, as a function source has no genomic column
            mapping.
        """
        # Arrange
        query = "SELECT *, CLUSTER(interval) AS cid FROM UNNEST([1, 2, 3])"

        # Act & assert
        with pytest.raises(ValueError, match="neither exposes the canonical"):
            transpile(query, tables=[])

    def test_transpile_should_raise_when_cte_is_recursive(self, custom):
        """Test a mutually-recursive CTE raises instead of looping forever.

        Given:
            A CLUSTER over a CTE whose body transitively selects from itself.
        When:
            Transpiling the query.
        Then:
            The cycle guard should terminate resolution and raise a ValueError
            rather than recursing without bound.
        """
        # Arrange
        query = (
            "WITH a AS (SELECT * FROM b), b AS (SELECT * FROM a) "
            "SELECT *, CLUSTER(interval) AS cid FROM a"
        )

        # Act & assert
        with pytest.raises(ValueError, match="neither exposes the canonical"):
            transpile(query, tables=[custom])

    def test_transpile_should_raise_when_from_clause_joins_two_relations(self, custom):
        """Test a join of two identically-mapped relations is rejected.

        Given:
            A CLUSTER whose FROM joins two tables that declare the same genomic
            column mapping.
        When:
            Transpiling the query.
        Then:
            It should raise a ValueError rejecting the join outright — the
            whole-query rewrite cannot preserve a join, so even agreeing mappings
            are not resolved.
        """
        # Arrange
        regions2 = Table("regions2", chrom_col="ch", start_col="s", end_col="e")
        query = (
            "SELECT *, CLUSTER(interval) AS cid FROM regions a JOIN regions2 b ON TRUE"
        )

        # Act & assert
        with pytest.raises(ValueError, match="over a join is not supported"):
            transpile(query, tables=[custom, regions2])

    def test_transpile_should_raise_when_join_mixes_concrete_and_opaque_sources(
        self, custom
    ):
        """Test a join of a concrete table and an opaque subquery is rejected.

        Given:
            A CLUSTER whose FROM joins a registered custom-mapped table with an
            opaque derived sibling.
        When:
            Transpiling the query.
        Then:
            It should raise the join rejection rather than silently preferring the
            concrete mapping and dropping the join partner.
        """
        # Arrange
        query = (
            "SELECT *, CLUSTER(interval) AS cid "
            "FROM regions JOIN (SELECT foo FROM regions) q ON TRUE"
        )

        # Act & assert
        with pytest.raises(ValueError, match="over a join is not supported"):
            transpile(query, tables=[custom])

    @pytest.mark.parametrize("op", ["CLUSTER", "MERGE"])
    def test_transpile_should_raise_when_from_clause_joins_conflicting_relations(
        self, op, custom
    ):
        """Test a join of two conflicting-mapped relations is rejected as a join.

        Given:
            A CLUSTER/MERGE whose FROM joins two tables declaring different genomic
            column mappings.
        When:
            Transpiling the query.
        Then:
            It should raise the join rejection — the join is refused before any
            conflicting-mapping comparison is even reached.
        """
        # Arrange
        other = Table("other", chrom_col="seqid", start_col="lo", end_col="hi")
        source = "regions JOIN other ON TRUE"
        query = _query(op, source)

        # Act & assert
        with pytest.raises(ValueError, match="over a join is not supported"):
            transpile(query, tables=[custom, other])

    def test_transpile_should_raise_when_qualified_star_projects_over_a_join(
        self, custom
    ):
        """Test a qualified ``alias.*`` projection does not sneak past a join.

        Given:
            A CLUSTER selecting a qualified ``a.*`` over a FROM that joins two
            relations.
        When:
            Transpiling the query.
        Then:
            It should still raise the join rejection — qualifying the star to one
            side does not exempt the FROM from being a join.
        """
        # Arrange
        regions2 = Table("regions2", chrom_col="ch", start_col="s", end_col="e")
        query = (
            "SELECT a.*, CLUSTER(interval) AS cid FROM regions a JOIN regions2 b ON TRUE"
        )

        # Act & assert
        with pytest.raises(ValueError, match="over a join is not supported"):
            transpile(query, tables=[custom, regions2])
