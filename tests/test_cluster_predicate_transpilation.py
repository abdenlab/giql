"""Transpilation tests for the CLUSTER/MERGE predicate argument.

Tests verify that an optional ``predicate :=`` argument ANDs into the
cluster-boundary CASE, that ``PREV(col)`` calls are rewritten to quoted LAG
windows over the operator's partition/order, and that omitting the predicate
leaves the emitted SQL byte-identical to the pre-predicate behavior.
"""

import pytest

from giql import transpile


class TestClusterPredicateTranspilation:
    """Tests for CLUSTER predicate transpilation to SQL."""

    def test_transpile_without_predicate_is_unchanged(self):
        """Test that a predicate-free CLUSTER transpiles to the legacy shape.

        Given:
            A CLUSTER query with no predicate argument.
        When:
            Transpiling to SQL.
        Then:
            It should emit the bare adjacency CASE with no AND clause.
        """
        # Act
        sql = transpile(
            "SELECT *, CLUSTER(interval) AS cid FROM peaks", tables=["peaks"]
        )

        # Assert
        assert (
            'CASE WHEN MAX("end") OVER (PARTITION BY "chrom" ORDER BY "start" '
            "NULLS LAST ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) "
            '>= "start" THEN 0 ELSE 1 END' in sql
        )

    def test_transpile_with_predicate_ands_into_case(self):
        """Test that a predicate ANDs into the cluster-boundary CASE.

        Given:
            A CLUSTER query with predicate := depth = PREV(depth).
        When:
            Transpiling to SQL.
        Then:
            It should AND the predicate onto the adjacency test inside the CASE.
        """
        # Act
        sql = transpile(
            "SELECT *, CLUSTER(interval, predicate := depth = PREV(depth)) AS cid "
            "FROM peaks",
            tables=["peaks"],
        )

        # Assert
        assert ' >= "start" AND (' in sql

    def test_transpile_rewrites_prev_call_to_lag_window(self):
        """Test that a PREV() call becomes a quoted LAG over the cluster window.

        Given:
            A CLUSTER query with predicate := depth = PREV(depth).
        When:
            Transpiling to SQL.
        Then:
            It should rewrite PREV(depth) to a quoted LAG("depth") over the same
            partition/order as the adjacency window.
        """
        # Act
        sql = transpile(
            "SELECT *, CLUSTER(interval, predicate := depth = PREV(depth)) AS cid "
            "FROM peaks",
            tables=["peaks"],
        )

        # Assert
        assert (
            '"depth" = LAG("depth") OVER (PARTITION BY "chrom" ORDER BY "start" '
            "NULLS LAST)" in sql
        )

    def test_transpile_predicate_uses_stranded_partition(self):
        """Test that PREV() LAG windows honor the stranded partition.

        Given:
            A stranded CLUSTER query with a predicate referencing PREV(name).
        When:
            Transpiling to SQL.
        Then:
            It should partition the predicate's LAG window by chrom and strand.
        """
        # Act
        sql = transpile(
            "SELECT *, CLUSTER(interval, stranded := true, "
            "predicate := name = PREV(name)) AS cid FROM peaks",
            tables=["peaks"],
        )

        # Assert
        assert (
            'LAG("name") OVER (PARTITION BY "chrom", "strand" ORDER BY "start" '
            "NULLS LAST)" in sql
        )

    def test_transpile_quotes_reserved_word_predicate_columns(self):
        """Test that reserved-word predicate columns are quoted on both sides.

        Given:
            A CLUSTER query whose predicate references the reserved-word
            genomic columns end and start (end = PREV(start)).
        When:
            Transpiling to SQL.
        Then:
            It should quote both the current-row column and the LAG argument so
            the emitted SQL is valid.
        """
        # Act
        sql = transpile(
            "SELECT *, CLUSTER(interval, predicate := end = PREV(start)) AS cid "
            "FROM peaks",
            tables=["peaks"],
        )

        # Assert
        assert '"end" = LAG("start") OVER (' in sql

    def test_transpile_prev_with_wrong_arity_raises(self):
        """Test that a PREV() call with the wrong arity is rejected.

        Given:
            A CLUSTER predicate calling PREV() with two arguments.
        When:
            Transpiling to SQL.
        Then:
            It should raise a ValueError naming the one-argument requirement.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError, match="exactly one column argument"):
            transpile(
                "SELECT *, CLUSTER(interval, predicate := depth = PREV(depth, score)) "
                "AS cid FROM peaks",
                tables=["peaks"],
            )

    def test_transpile_nested_prev_raises(self):
        """Test that a nested PREV() call is rejected.

        Given:
            A CLUSTER predicate nesting PREV() inside another PREV().
        When:
            Transpiling to SQL.
        Then:
            It should raise a ValueError naming the no-nesting restriction.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError, match="cannot be nested"):
            transpile(
                "SELECT *, CLUSTER(interval, predicate := depth = PREV(PREV(depth))) "
                "AS cid FROM peaks",
                tables=["peaks"],
            )


class TestMergePredicateTranspilation:
    """Tests for MERGE predicate transpilation to SQL."""

    def test_transpile_without_predicate_is_unchanged(self):
        """Test that a predicate-free MERGE transpiles to the legacy shape.

        Given:
            A MERGE query with no predicate argument.
        When:
            Transpiling to SQL.
        Then:
            It should emit the bare adjacency CASE with no AND clause.
        """
        # Act
        sql = transpile("SELECT MERGE(interval) FROM peaks", tables=["peaks"])

        # Assert
        assert (
            'CASE WHEN MAX("end") OVER (PARTITION BY "chrom" ORDER BY "start" '
            "NULLS LAST ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) "
            '>= "start" THEN 0 ELSE 1 END' in sql
        )

    def test_transpile_predicate_inherited_through_cluster(self):
        """Test that MERGE inherits predicate-aware cluster boundaries.

        Given:
            A MERGE query with predicate := depth = PREV(depth).
        When:
            Transpiling to SQL.
        Then:
            It should AND the rewritten predicate into the underlying CLUSTER
            CASE that drives the GROUP BY.
        """
        # Act
        sql = transpile(
            "SELECT MERGE(interval, predicate := depth = PREV(depth)) FROM peaks",
            tables=["peaks"],
        )

        # Assert
        assert ' >= "start" AND (' in sql
        assert (
            '"depth" = LAG("depth") OVER (PARTITION BY "chrom" ORDER BY "start" '
            "NULLS LAST)" in sql
        )
