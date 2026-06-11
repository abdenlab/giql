"""Transpilation tests for the CLUSTER/MERGE predicate argument.

Tests verify that an optional ``predicate :=`` argument ANDs into the
cluster-boundary CASE, that ``prev.col`` references are rewritten to LAG
windows over the operator's partition/order, and that omitting the predicate
leaves the emitted SQL byte-identical to the pre-predicate behavior.
"""

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
            'CASE WHEN LAG("end") OVER (PARTITION BY "chrom" ORDER BY "start" '
            'NULLS LAST) >= "start" THEN 0 ELSE 1 END' in sql
        )

    def test_transpile_with_predicate_ands_into_case(self):
        """Test that a predicate ANDs into the cluster-boundary CASE.

        Given:
            A CLUSTER query with predicate := depth = prev.depth.
        When:
            Transpiling to SQL.
        Then:
            It should AND the predicate onto the adjacency test inside the CASE.
        """
        # Act
        sql = transpile(
            "SELECT *, CLUSTER(interval, predicate := depth = prev.depth) AS cid "
            "FROM peaks",
            tables=["peaks"],
        )

        # Assert
        assert ' >= "start" AND (' in sql

    def test_transpile_rewrites_prev_ref_to_lag_window(self):
        """Test that a prev. reference becomes a LAG over the cluster window.

        Given:
            A CLUSTER query with predicate := depth = prev.depth.
        When:
            Transpiling to SQL.
        Then:
            It should rewrite prev.depth to LAG(depth) over the same
            partition/order as the adjacency window.
        """
        # Act
        sql = transpile(
            "SELECT *, CLUSTER(interval, predicate := depth = prev.depth) AS cid "
            "FROM peaks",
            tables=["peaks"],
        )

        # Assert
        assert (
            'depth = LAG(depth) OVER (PARTITION BY "chrom" ORDER BY "start" '
            "NULLS LAST)" in sql
        )

    def test_transpile_predicate_uses_stranded_partition(self):
        """Test that prev. LAG windows honor the stranded partition.

        Given:
            A stranded CLUSTER query with a predicate referencing prev.name.
        When:
            Transpiling to SQL.
        Then:
            It should partition the predicate's LAG window by chrom and strand.
        """
        # Act
        sql = transpile(
            "SELECT *, CLUSTER(interval, stranded := true, "
            "predicate := name = prev.name) AS cid FROM peaks",
            tables=["peaks"],
        )

        # Assert
        assert (
            'LAG(name) OVER (PARTITION BY "chrom", "strand" ORDER BY "start" '
            "NULLS LAST)" in sql
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
            'CASE WHEN LAG("end") OVER (PARTITION BY "chrom" ORDER BY "start" '
            'NULLS LAST) >= "start" THEN 0 ELSE 1 END' in sql
        )

    def test_transpile_predicate_inherited_through_cluster(self):
        """Test that MERGE inherits predicate-aware cluster boundaries.

        Given:
            A MERGE query with predicate := depth = prev.depth.
        When:
            Transpiling to SQL.
        Then:
            It should AND the rewritten predicate into the underlying CLUSTER
            CASE that drives the GROUP BY.
        """
        # Act
        sql = transpile(
            "SELECT MERGE(interval, predicate := depth = prev.depth) FROM peaks",
            tables=["peaks"],
        )

        # Assert
        assert ' >= "start" AND (' in sql
        assert (
            'depth = LAG(depth) OVER (PARTITION BY "chrom" ORDER BY "start" '
            "NULLS LAST)" in sql
        )
