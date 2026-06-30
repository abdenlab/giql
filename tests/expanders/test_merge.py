"""Behavioral tests for the MERGE operator expander (#144).

MERGE migrated from the pre-pass ``MergeTransformer`` to a registered expander
(``giql.expanders.merge``) that rewrites the enclosing SELECT in place into the
clustered-aggregation form. These tests drive the public ``transpile`` API and
pin the behaviors the migration newly exercises (the registry pass, the
recursion-removal, and the error/limitation guards) that the legacy
transpilation suites did not already cover. The predicate byte-shape is pinned
separately in ``tests/test_cluster_predicate_transpilation.py``.
"""

import pytest

from giql.table import Table
from giql.transpile import transpile


class TestMergeExpander:
    """Expansion of MERGE through the operator-expander registry (#144)."""

    def test_transpile_should_raise_when_multiple_merge_in_one_select(self):
        """Test that two MERGE expressions in one SELECT are rejected.

        Given:
            A SELECT projecting two MERGE expressions.
        When:
            Transpiling the query.
        Then:
            It should raise ValueError naming the unsupported multiple-MERGE case.
        """
        # Arrange
        query = "SELECT MERGE(interval), MERGE(interval, 100) FROM peaks"

        # Act & assert
        with pytest.raises(
            ValueError, match="Multiple MERGE expressions not yet supported"
        ):
            transpile(query, tables=["peaks"])

    @pytest.mark.parametrize(
        "query",
        [
            "SELECT MERGE(interval), CLUSTER(interval) AS cid FROM peaks",
            "SELECT CLUSTER(interval) AS cid, MERGE(interval) FROM peaks",
        ],
    )
    def test_transpile_should_raise_when_cluster_and_merge_share_select(self, query):
        """Test that combining CLUSTER and MERGE in one SELECT is rejected.

        Given:
            A SELECT projecting both a MERGE and a CLUSTER (either order).
        When:
            Transpiling the query.
        Then:
            It should raise ValueError naming the unsupported combination, rather
            than silently emitting SQL with a leaked, unexpanded operator.
        """
        # Act & assert
        with pytest.raises(ValueError, match="CLUSTER and MERGE cannot be combined"):
            transpile(query, tables=["peaks"])

    def test_transpile_should_group_by_strand_when_merge_stranded(self):
        """Test that a stranded MERGE aggregates within strand.

        Given:
            A stranded MERGE query.
        When:
            Transpiling the query.
        Then:
            It should aggregate MIN(start)/MAX(end), group by chrom, strand, and
            the synthesized cluster id, project strand, and partition the inner
            windows by chrom and strand.
        """
        # Arrange
        query = "SELECT MERGE(interval, stranded := true) FROM peaks"

        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        assert 'MIN("start") AS start' in sql
        assert 'MAX("end") AS end' in sql
        assert 'GROUP BY "chrom", "strand", __giql_cluster_id' in sql
        assert 'PARTITION BY "chrom", "strand"' in sql
        assert "G_I_Q_L" not in sql

    def test_transpile_should_use_custom_columns_when_table_declares_them(self):
        """Test that a stranded MERGE honors a custom column mapping.

        Given:
            A stranded MERGE over a Table declaring custom chrom/start/end/strand
            column names.
        When:
            Transpiling the query.
        Then:
            The aggregation, GROUP BY, and window partitions should use the custom
            column names, never the canonical defaults.
        """
        # Arrange
        regions = Table(
            "regions", chrom_col="ch", start_col="s", end_col="e", strand_col="st"
        )
        query = "SELECT MERGE(interval, stranded := true) FROM regions"

        # Act
        sql = transpile(query, tables=[regions])

        # Assert
        assert 'MIN("s") AS s' in sql
        assert 'MAX("e") AS e' in sql
        assert 'GROUP BY "ch", "st", __giql_cluster_id' in sql
        assert 'PARTITION BY "ch", "st"' in sql
        assert '"chrom"' not in sql and '"start"' not in sql

    @pytest.mark.parametrize(
        "query",
        [
            "SELECT * FROM (SELECT MERGE(interval) FROM peaks) x",
            "WITH c AS (SELECT MERGE(interval) FROM peaks) SELECT * FROM c",
        ],
    )
    def test_transpile_should_expand_merge_nested_in_subquery_or_cte(self, query):
        """Test that a MERGE nested below the root SELECT still expands.

        Given:
            A MERGE inside a FROM-subquery, and a MERGE inside a WITH CTE.
        When:
            Transpiling the query.
        Then:
            The nested MERGE should expand into the clustered-aggregation form
            with no leaked operator — the pass walk replaces the manual recursion
            the legacy transformer performed.
        """
        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        assert "G_I_Q_L" not in sql
        assert "AS clustered" in sql
        assert "__giql_cluster_id" in sql

    def test_transpile_should_carry_where_into_clustered_subquery(self):
        """Test that a WHERE clause is pushed into MERGE's clustered subquery.

        Given:
            A MERGE query with a WHERE clause.
        When:
            Transpiling the query.
        Then:
            The WHERE predicate should appear inside the inner clustered subquery
            that feeds the aggregation.
        """
        # Arrange
        query = "SELECT MERGE(interval) AS m FROM peaks WHERE chrom = 'chr1'"

        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        assert "WHERE chrom = 'chr1'" in sql
        assert "AS clustered" in sql

    @pytest.mark.parametrize(
        "predicate, message",
        [
            ("depth = PREV(depth, score)", "exactly one column argument"),
            ("depth = PREV(PREV(depth))", "cannot be nested"),
        ],
    )
    def test_transpile_should_validate_prev_in_merge_predicate(self, predicate, message):
        """Test that MERGE reuses CLUSTER's PREV validation.

        Given:
            A MERGE predicate calling PREV with the wrong arity or a nested PREV.
        When:
            Transpiling the query.
        Then:
            It should raise the matching ValueError, proving the MERGE path
            reaches the shared predecessor-reference validation.
        """
        # Arrange
        query = f"SELECT MERGE(interval, predicate := {predicate}) FROM peaks"

        # Act & assert
        with pytest.raises(ValueError, match=message):
            transpile(query, tables=["peaks"])

    def test_transpile_should_expand_merge_in_projection_scalar_subquery(self):
        """Test that a MERGE in a projection scalar-subquery is expanded.

        Given:
            A MERGE inside a scalar subquery in the SELECT list (a shape the legacy
            transformer's manual recursion did not descend into, leaking an
            unexpanded operator).
        When:
            Transpiling the query.
        Then:
            The MERGE should expand into the clustered-aggregation form with no
            leaked operator — the pass walk reaches projection subqueries.
        """
        # Arrange
        query = "SELECT (SELECT MERGE(interval) FROM peaks) AS m FROM other"

        # Act
        sql = transpile(query, tables=["peaks", "other"])

        # Assert
        assert "G_I_Q_L" not in sql
        assert "AS clustered" in sql
        assert "__giql_cluster_id" in sql

    @pytest.mark.parametrize("predicate_op", ["INTERSECTS", "CONTAINS", "WITHIN"])
    @pytest.mark.parametrize(
        "projection",
        ["*, MERGE(interval) AS m", "MERGE(interval)"],
        ids=["aliased", "bare"],
    )
    def test_transpile_should_expand_spatial_predicate_copied_into_clustered(
        self, projection, predicate_op
    ):
        """Test that a spatial WHERE predicate survives the MERGE rewrite.

        Given:
            A MERGE query (aliased or bare) whose WHERE filters on a spatial
            predicate, which the rewrite copies into the inner clustered subquery.
        When:
            Transpiling the query.
        Then:
            The copied predicate should itself be expanded — no leaked, unexpanded
            operator — for both projection depths, pinning the #144 B1 regression
            where the aliased MERGE expanded before the predicate and stranded a
            live, unexpanded copy in the subquery.
        """
        # Arrange
        query = (
            f"SELECT {projection} FROM peaks a "
            f"WHERE a.interval {predicate_op} 'chr1:1-100'"
        )

        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        assert "G_I_Q_L" not in sql
        assert "AS clustered" in sql

    @pytest.mark.parametrize(
        "query",
        [
            "SELECT ABS(MERGE(interval)) FROM peaks",
            "SELECT MERGE(interval) + 1 AS m FROM peaks",
        ],
    )
    def test_transpile_should_raise_when_merge_nested_in_projection_expression(
        self, query
    ):
        """Test that a MERGE buried in a projection expression is rejected.

        Given:
            A MERGE nested inside a larger projection expression (a function call or
            arithmetic), which has no coherent whole-query rewrite.
        When:
            Transpiling the query.
        Then:
            It should raise ValueError requiring a top-level projection item,
            rather than leaking an unexpanded operator to the generator.
        """
        # Act & assert
        with pytest.raises(ValueError, match="must be a top-level projection item"):
            transpile(query, tables=["peaks"])

    def test_transpile_should_quote_group_by_chrom_when_chrom_is_reserved_word(self):
        """Test that the MERGE GROUP BY chrom term is quoted.

        Given:
            A MERGE over a Table whose chrom column is a SQL reserved word.
        When:
            Transpiling the query.
        Then:
            The GROUP BY chrom term should be quoted (like every other chrom
            reference), so the reserved-word column emits valid SQL (#144 A13).
        """
        # Arrange
        regions = Table("regions", chrom_col="order", start_col="s", end_col="e")
        query = "SELECT MERGE(interval) FROM regions"

        # Act
        sql = transpile(query, tables=[regions])

        # Assert
        assert 'GROUP BY "order"' in sql
        assert "GROUP BY order" not in sql
