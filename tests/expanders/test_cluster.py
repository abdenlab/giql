"""Behavioral tests for the CLUSTER operator expander (#144).

CLUSTER migrated from the pre-pass ``ClusterTransformer`` to a registered
expander (``giql.expanders.cluster``) that rewrites the enclosing SELECT in place
into the two-level ``lag_calc`` form. These tests drive the public ``transpile``
API and pin the behaviors the migration newly exercises — the pass walk
replacing the transformer's manual recursion (nested placements), the
FROM-table column derivation, the distance/clause branches, and projection
shapes — that the legacy transpilation suites did not already cover. The
predicate byte-shape is pinned in ``tests/test_cluster_predicate_transpilation.py``.
"""

import pytest

from giql.table import Table
from giql.transpile import transpile


class TestClusterExpander:
    """Expansion of CLUSTER through the operator-expander registry (#144)."""

    @pytest.mark.parametrize(
        "query",
        [
            "SELECT * FROM (SELECT *, CLUSTER(interval) AS cid FROM peaks) x",
            "WITH c AS (SELECT *, CLUSTER(interval) AS cid FROM peaks) SELECT * FROM c",
        ],
    )
    def test_transpile_should_expand_cluster_nested_in_subquery_or_cte(self, query):
        """Test that a CLUSTER nested below the root SELECT still expands.

        Given:
            A CLUSTER inside a FROM-subquery, and a CLUSTER inside a WITH CTE.
        When:
            Transpiling the query.
        Then:
            The nested CLUSTER should expand into the two-level lag_calc form with
            no leaked operator — the pass walk replaces the manual recursion the
            legacy transformer performed.
        """
        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        assert "G_I_Q_L" not in sql
        assert "AS __giql_lag_calc" in sql
        assert "__giql_is_new_cluster" in sql

    def test_transpile_should_raise_when_multiple_cluster_in_one_select(self):
        """Test that two CLUSTER expressions in one SELECT are rejected.

        Given:
            A SELECT projecting two CLUSTER expressions.
        When:
            Transpiling the query.
        Then:
            It should raise ValueError naming the unsupported multiple-CLUSTER
            case, rather than chaining the rewrite into non-executable SQL (a
            duplicate __giql_lag_calc alias / __giql_is_new_cluster binder error).
        """
        # Arrange
        query = (
            "SELECT *, CLUSTER(interval) AS a, CLUSTER(interval, 100) AS b FROM peaks"
        )

        # Act & assert
        with pytest.raises(
            ValueError, match="Multiple CLUSTER expressions not yet supported"
        ):
            transpile(query, tables=["peaks"])

    def test_transpile_should_use_custom_columns_when_table_declares_them(self):
        """Test that a stranded CLUSTER honors a custom column mapping.

        Given:
            A stranded CLUSTER over a Table declaring custom column names.
        When:
            Transpiling the query.
        Then:
            The window partition and order should use the custom column names, not
            the canonical defaults.
        """
        # Arrange
        regions = Table(
            "regions", chrom_col="ch", start_col="s", end_col="e", strand_col="st"
        )
        query = "SELECT *, CLUSTER(interval, stranded := true) AS cid FROM regions"

        # Act
        sql = transpile(query, tables=[regions])

        # Assert
        assert 'PARTITION BY "ch", "st" ORDER BY "s"' in sql
        assert 'LAG("e")' in sql
        assert '"chrom"' not in sql and '"start"' not in sql

    def test_transpile_should_add_distance_offset_to_lag_when_distance_positive(self):
        """Test that a positive CLUSTER distance offsets the adjacency LAG.

        Given:
            A CLUSTER with a positive distance.
        When:
            Transpiling the query.
        Then:
            The adjacency should add the distance to the LAG before comparing to
            start.
        """
        # Arrange
        query = "SELECT *, CLUSTER(interval, 100) AS cid FROM peaks"

        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        window = 'OVER (PARTITION BY "chrom" ORDER BY "start" NULLS LAST)'
        assert f'LAG("end") {window} + 100 >= "start"' in sql

    def test_transpile_should_not_offset_lag_when_no_distance(self):
        """Test that a CLUSTER without distance uses a bare adjacency.

        Given:
            A CLUSTER with no distance argument.
        When:
            Transpiling the query.
        Then:
            The adjacency should compare the bare LAG to start with no offset.
        """
        # Arrange
        query = "SELECT *, CLUSTER(interval) AS cid FROM peaks"

        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        window = 'OVER (PARTITION BY "chrom" ORDER BY "start" NULLS LAST)'
        assert f'LAG("end") {window} >= "start"' in sql
        assert f"{window} + " not in sql

    def test_transpile_should_split_clauses_between_lag_calc_and_outer_query(self):
        """Test that CLUSTER places clauses at the correct query level.

        Given:
            A CLUSTER query carrying WHERE, GROUP BY, HAVING, and ORDER BY.
        When:
            Transpiling the query.
        Then:
            WHERE/GROUP BY/HAVING should land inside the inner lag_calc subquery
            and ORDER BY should attach to the outer query.
        """
        # Arrange
        query = (
            "SELECT chrom, CLUSTER(interval) AS cid FROM peaks "
            "WHERE chrom = 'chr1' GROUP BY chrom HAVING COUNT(*) > 1 ORDER BY chrom"
        )

        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        assert (
            "WHERE chrom = 'chr1' GROUP BY chrom HAVING COUNT(*) > 1) AS __giql_lag_calc"
            in sql
        )
        assert ") AS __giql_lag_calc ORDER BY chrom" in sql

    def test_transpile_should_expand_bare_cluster_without_alias(self):
        """Test that an un-aliased CLUSTER expands to a bare SUM window.

        Given:
            A CLUSTER projected without an AS alias.
        When:
            Transpiling the query.
        Then:
            The bare CLUSTER should be replaced by an un-aliased SUM window with no
            leaked operator.
        """
        # Arrange
        query = "SELECT *, CLUSTER(interval) FROM peaks"

        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        assert "G_I_Q_L" not in sql
        assert "SUM(__giql_is_new_cluster) OVER" in sql

    def test_transpile_should_keep_explicit_projection_columns_in_lag_calc(self):
        """Test that explicit (non-star) projection columns flow into lag_calc.

        Given:
            A CLUSTER query with an explicit column projection (not SELECT *).
        When:
            Transpiling the query.
        Then:
            The explicit columns should be projected by the inner lag_calc
            subquery feeding the cluster window.
        """
        # Arrange
        query = "SELECT chrom, start, CLUSTER(interval) AS c FROM peaks"

        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        assert "SELECT chrom, start," in sql
        assert "AS __giql_lag_calc" in sql
        assert "G_I_Q_L" not in sql

    # Note: CLUSTER combined with an INTERSECTS *join* in the same SELECT is not
    # tested here. The clause copy into lag_calc carries only FROM/WHERE/GROUP/
    # HAVING (never JOINs), so a join is dropped — pre-existing legacy behavior
    # that #144 preserves byte-for-byte, not a migration concern.

    def test_transpile_should_compose_distance_stranded_and_predicate(self):
        """Test that CLUSTER composes distance, stranded, and predicate together.

        Given:
            A CLUSTER with a distance, stranded mode, and a PREV-based predicate.
        When:
            Transpiling the query.
        Then:
            The output should offset the LAG by the distance, partition by chrom
            and strand, and rewrite PREV into a LAG window inside the adjacency.
        """
        # Arrange
        query = (
            "SELECT *, CLUSTER(interval, 1000, stranded := true, "
            "predicate := name = PREV(name)) AS cid FROM peaks"
        )

        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        assert "+ 1000 >= " in sql
        assert 'PARTITION BY "chrom", "strand"' in sql
        assert 'LAG("name")' in sql
        assert "G_I_Q_L" not in sql

    def test_transpile_should_expand_cluster_in_union_branch(self):
        """Test that a CLUSTER inside a UNION branch is expanded.

        Given:
            A CLUSTER in each branch of a UNION (a shape the legacy transformer's
            manual recursion did not descend into, leaking an unexpanded operator).
        When:
            Transpiling the query.
        Then:
            Both branches should expand to the lag_calc form with no leaked
            operator — the pass walk reaches UNION branches.
        """
        # Arrange
        query = (
            "SELECT *, CLUSTER(interval) AS c FROM peaks "
            "UNION ALL SELECT *, CLUSTER(interval) AS c FROM peaks"
        )

        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        assert "G_I_Q_L" not in sql
        assert sql.count("AS __giql_lag_calc") == 2

    @pytest.mark.parametrize("predicate_op", ["INTERSECTS", "CONTAINS", "WITHIN"])
    @pytest.mark.parametrize(
        "projection",
        ["*, CLUSTER(interval) AS cid", "CLUSTER(interval)"],
        ids=["aliased", "bare"],
    )
    def test_transpile_should_expand_spatial_predicate_copied_into_lag_calc(
        self, projection, predicate_op
    ):
        """Test that a spatial WHERE predicate survives the CLUSTER rewrite.

        Given:
            A CLUSTER query (aliased or bare) whose WHERE filters on a spatial
            predicate, which the rewrite copies into the inner lag_calc subquery.
        When:
            Transpiling the query.
        Then:
            The copied predicate should itself be expanded — no leaked, unexpanded
            operator — for both projection depths, pinning the #144 B1 regression
            where the aliased CLUSTER expanded before the predicate and stranded a
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
        assert "AS __giql_lag_calc" in sql

    @pytest.mark.parametrize(
        "query",
        [
            "SELECT ABS(CLUSTER(interval)) FROM peaks",
            "SELECT CLUSTER(interval) + 1 AS c FROM peaks",
        ],
    )
    def test_transpile_should_raise_when_cluster_nested_in_projection_expression(
        self, query
    ):
        """Test that a CLUSTER buried in a projection expression is rejected.

        Given:
            A CLUSTER nested inside a larger projection expression (a function call
            or arithmetic), which has no coherent whole-query rewrite.
        When:
            Transpiling the query.
        Then:
            It should raise ValueError requiring a top-level projection item,
            rather than leaking an unexpanded operator to the generator.
        """
        # Act & assert
        with pytest.raises(ValueError, match="must be a top-level projection item"):
            transpile(query, tables=["peaks"])

    def test_transpile_should_inject_lag_calc_columns_deterministically(self):
        """Test that the injected lag_calc column order is hash-seed independent.

        Given:
            An explicit-projection CLUSTER whose predicate forces several residual
            columns to be injected into lag_calc, transpiled in two child
            interpreters under differing PYTHONHASHSEED values.
        When:
            Comparing the two emitted strings.
        Then:
            They should be byte-identical, proving the injected-column order is
            sorted rather than set-iteration (PYTHONHASHSEED) dependent (#144 A2).
        """
        # Arrange
        import os
        import subprocess
        import sys

        code = (
            "from giql.transpile import transpile; "
            "print(transpile("
            '"SELECT chrom, CLUSTER(interval, stranded := true, '
            "predicate := name = PREV(score)) AS c FROM peaks\", tables=['peaks']))"
        )
        base_env = {
            k: v
            for k, v in os.environ.items()
            if not k.startswith("COV_CORE") and k != "COVERAGE_PROCESS_START"
        }

        def _run(seed: str) -> str:
            result = subprocess.run(
                [sys.executable, "-c", code],
                capture_output=True,
                text=True,
                env={**base_env, "PYTHONHASHSEED": seed},
            )
            assert result.returncode == 0, result.stderr
            return result.stdout

        # Act
        out_a = _run("0")
        out_b = _run("1")

        # Assert
        assert out_a == out_b
        assert "G_I_Q_L" not in out_a

    def test_transpile_should_namespace_synthesized_cluster_identifiers(self):
        """Test that CLUSTER's synthesized identifiers carry the reserved prefix.

        Given:
            A plain CLUSTER query.
        When:
            Transpiling the query.
        Then:
            The emitted SQL should alias the synthesized helpers as the reserved
            __giql_lag_calc and __giql_is_new_cluster identifiers and never emit the
            bare forms, guarding #161 against a partial revert that would re-collide
            with a user column named is_new_cluster.
        """
        # Arrange
        query = "SELECT *, CLUSTER(interval) AS cid FROM peaks"

        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        assert "AS __giql_lag_calc" in sql
        assert "AS __giql_is_new_cluster" in sql
        assert "AS lag_calc" not in sql
        assert "AS is_new_cluster" not in sql

    def test_transpile_should_disambiguate_cluster_when_source_has_is_new_cluster_column(
        self,
    ):
        """Test that CLUSTER runs correctly beside a user is_new_cluster column.

        Given:
            A table carrying a user column literally named is_new_cluster, seeded
            with a poison value that would corrupt the cumulative SUM if the outer
            window bound to it instead of CLUSTER's synthesized flag (the pre-#161
            collision, which executed silently and returned wrong ids).
        When:
            Transpiling a star-projected CLUSTER and executing it on DuckDB.
        Then:
            The query should execute and return the correct per-partition cluster
            ids, proving the synthesized flag now lives in the reserved __giql_
            namespace and no longer collides with the user column, which itself
            survives unchanged in the output.
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            "CREATE TABLE intervals ("
            'chrom VARCHAR, "start" INTEGER, "end" INTEGER, is_new_cluster INTEGER)'
        )
        conn.executemany(
            "INSERT INTO intervals VALUES (?, ?, ?, ?)",
            [
                ("chr1", 100, 200, 7),  # i1
                ("chr1", 150, 250, 7),  # overlaps i1 -> same cluster
                ("chr1", 400, 500, 7),  # separate -> new cluster
            ],
        )

        # Act
        sql = transpile(
            "SELECT *, CLUSTER(interval) AS cluster_id FROM intervals",
            tables=["intervals"],
        )
        cursor = conn.execute(sql)
        columns = [desc[0] for desc in cursor.description]
        rows = [dict(zip(columns, row)) for row in cursor.fetchall()]

        # Assert
        assert [r["cluster_id"] for r in rows] == [1, 1, 2]
        assert [r["is_new_cluster"] for r in rows] == [7, 7, 7]

    @pytest.mark.parametrize(
        "query, expected_with",
        [
            (
                "WITH sub AS (SELECT * FROM peaks) "
                "SELECT *, CLUSTER(interval) AS cid FROM sub",
                "WITH sub AS (SELECT * FROM peaks)",
            ),
            (
                "WITH a AS (SELECT * FROM peaks), sub AS (SELECT * FROM a) "
                "SELECT *, CLUSTER(interval) AS cid FROM sub",
                "WITH a AS (SELECT * FROM peaks), sub AS (SELECT * FROM a)",
            ),
        ],
        ids=["single_cte", "chained_ctes"],
    )
    def test_transpile_should_preserve_with_clause_when_cluster_over_cte(
        self, query, expected_with
    ):
        """Test that CLUSTER over a CTE FROM keeps the enclosing WITH clause.

        Given:
            A CLUSTER whose FROM references a CTE — a single CTE, and a chain of
            two CTEs where the driving relation is defined in terms of the first.
        When:
            Transpiling the query.
        Then:
            The emitted SQL should retain the enclosing WITH clause verbatim so the
            rewritten __giql_lag_calc subquery's ``FROM sub`` reference still
            resolves, rather than dropping the WITH and dangling an undefined CTE
            reference — the pre-#174 transplant defect.
        """
        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        assert sql.startswith(expected_with)
        assert "AS __giql_lag_calc" in sql
        assert "FROM sub) AS __giql_lag_calc" in sql

    def test_transpile_should_execute_when_cluster_over_cte(self):
        """Test that CLUSTER over a CTE FROM emits executable SQL.

        Given:
            A CLUSTER whose FROM references a pass-through CTE over a seeded table.
        When:
            Transpiling the query and executing it on DuckDB.
        Then:
            The query should execute (the preserved WITH keeps the CTE reference
            resolvable) and return the correct per-partition cluster ids, proving
            the #174 fix emits runnable SQL rather than referencing a dropped CTE.
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?)",
            [
                ("chr1", 100, 200),  # i1
                ("chr1", 150, 250),  # overlaps i1 -> same cluster
                ("chr1", 400, 500),  # separate -> new cluster
            ],
        )
        query = (
            "WITH sub AS (SELECT * FROM peaks) "
            "SELECT *, CLUSTER(interval) AS cid FROM sub"
        )

        # Act
        sql = transpile(query, tables=["peaks"])
        cursor = conn.execute(sql)
        columns = [desc[0] for desc in cursor.description]
        rows = [dict(zip(columns, row)) for row in cursor.fetchall()]

        # Assert
        # Sorted: the emitted CLUSTER query carries no outer ORDER BY, so assert the
        # cluster-id multiset rather than an insertion-order-dependent sequence.
        assert sorted(r["cid"] for r in rows) == [1, 1, 2]

    def test_transpile_should_not_hoist_outer_with_when_cluster_nested_in_cte_body(
        self,
    ):
        """Test that a CLUSTER inside a CTE body keeps the outer WITH at the root.

        Given:
            A CLUSTER nested inside the body of one CTE whose driving relation is a
            sibling CTE, under a root SELECT that carries the enclosing WITH.
        When:
            Transpiling the query.
        Then:
            The enclosing WITH should stay at the root and the rewritten
            __giql_lag_calc subquery should still reference the sibling CTE ``a`` —
            the inner transplant reads only its own SELECT's (absent) WITH, so it
            neither hoists the ancestor WITH nor drops the sibling CTE reference.
        """
        # Arrange
        query = (
            "WITH a AS (SELECT * FROM peaks), "
            "sub AS (SELECT *, CLUSTER(interval) AS cid FROM a) "
            "SELECT * FROM sub"
        )

        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        assert sql.startswith("WITH a AS (SELECT * FROM peaks), sub AS (")
        assert "FROM a) AS __giql_lag_calc" in sql
        assert sql.rstrip().endswith("SELECT * FROM sub")

    def test_transpile_should_not_emit_with_when_cluster_over_bare_table(self):
        """Test that a CLUSTER over a bare table grows no spurious WITH.

        Given:
            A CLUSTER over a bare registered table with no enclosing WITH.
        When:
            Transpiling the query.
        Then:
            The emitted SQL should not begin with a WITH clause — the transplant
            WITH-preservation branch is a no-op on the common, CTE-free path (#174).
        """
        # Act
        sql = transpile(
            "SELECT *, CLUSTER(interval) AS cid FROM peaks", tables=["peaks"]
        )

        # Assert
        assert not sql.lstrip().upper().startswith("WITH")
