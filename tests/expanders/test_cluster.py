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
        # dialect="duckdb" so the flag-hiding star emits DuckDB's EXCLUDE spelling
        # (#184); the generic EXCEPT form is not DuckDB-executable.
        sql = transpile(
            "SELECT *, CLUSTER(interval) AS cluster_id FROM intervals",
            tables=["intervals"],
            dialect="duckdb",
        )
        cursor = conn.execute(sql)
        columns = [desc[0] for desc in cursor.description]
        rows = [dict(zip(columns, row)) for row in cursor.fetchall()]

        # Assert
        assert [r["cluster_id"] for r in rows] == [1, 1, 2]
        assert [r["is_new_cluster"] for r in rows] == [7, 7, 7]
        # The synthesized flag is hidden; only the user's is_new_cluster surfaces.
        assert "__giql_is_new_cluster" not in columns

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
        # dialect="duckdb" so the flag-hiding star emits DuckDB's EXCLUDE spelling
        # (#184); the generic EXCEPT form is not DuckDB-executable.
        sql = transpile(query, tables=["peaks"], dialect="duckdb")
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

    @pytest.mark.parametrize(
        "clause",
        ["LIMIT 2", "LIMIT 2 OFFSET 1", "QUALIFY cid > 1"],
        ids=["limit", "offset", "qualify"],
    )
    def test_transpile_should_preserve_result_clause_when_cluster(self, clause):
        """Test that CLUSTER preserves an outer result-shaping clause.

        Given:
            A CLUSTER query carrying an outer LIMIT / OFFSET / QUALIFY clause.
        When:
            Transpiling the query.
        Then:
            The full clause should survive verbatim on the rewritten outer query
            rather than being silently dropped by the whole-query transplant (#181),
            so the emitted SQL returns the rows the original query asked for.
        """
        # Act
        sql = transpile(
            f"SELECT *, CLUSTER(interval) AS cid FROM peaks {clause}", tables=["peaks"]
        )

        # Assert
        assert clause in sql

    def test_transpile_should_preserve_distinct_when_cluster(self):
        """Test that CLUSTER preserves an outer DISTINCT.

        Given:
            A CLUSTER query with a DISTINCT projection.
        When:
            Transpiling the query.
        Then:
            The rewritten outer query should keep SELECT DISTINCT rather than
            dropping it in the transplant (#181).
        """
        # Act
        sql = transpile(
            "SELECT DISTINCT chrom, CLUSTER(interval) AS cid FROM peaks",
            tables=["peaks"],
        )

        # Assert
        assert "SELECT DISTINCT" in sql

    def test_transpile_should_dedup_when_cluster_with_distinct(self):
        """Test that a preserved DISTINCT dedups the CLUSTER output rows.

        Given:
            A CLUSTER over rows that collapse to one distinct (chrom, cluster-id)
            pair, projected with DISTINCT on an explicit column list.
        When:
            Transpiling the query and executing it on DuckDB.
        Then:
            The preserved DISTINCT should land on the outer query and dedup the final
            rows to the single distinct pair — proving DISTINCT is re-attached to the
            outer query, not an inner subquery (#181). An explicit projection is used
            so the dedup is not confounded by the pre-existing star-column leak (#184).
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?)",
            [("chr1", 10, 20), ("chr1", 12, 22), ("chr1", 14, 24)],
        )
        query = "SELECT DISTINCT chrom, CLUSTER(interval) AS cid FROM peaks"

        # Act
        sql = transpile(query, tables=["peaks"])
        rows = conn.execute(sql).fetchall()

        # Assert
        assert rows == [("chr1", 1)]

    def test_transpile_should_hide_flag_when_cluster_projects_star(self):
        """Test that a star-projected CLUSTER EXCEPTs the synthesized flag.

        Given:
            A top-level CLUSTER with a SELECT * projection.
        When:
            Transpiling the query.
        Then:
            The outer star should EXCEPT the synthesized __giql_is_new_cluster flag so
            it never surfaces in the query output (#184), and the bare flag name
            should not appear as a projected column.
        """
        # Act
        sql = transpile(
            "SELECT *, CLUSTER(interval) AS cid FROM peaks", tables=["peaks"]
        )

        # Assert
        assert 'EXCEPT ("__giql_is_new_cluster")' in sql

    def test_transpile_should_execute_star_without_flag_when_cluster_projects_star(
        self,
    ):
        """Test that a star-projected CLUSTER omits the flag from its output columns.

        Given:
            A star-projected CLUSTER over a seeded table.
        When:
            Transpiling for DuckDB and executing it.
        Then:
            The output columns should be exactly the user columns plus the cluster id,
            with the synthesized __giql_is_new_cluster flag hidden (#184), and the
            cluster ids should be correct.
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?)",
            [("chr1", 10, 20), ("chr1", 15, 30), ("chr1", 100, 110)],
        )

        # Act
        sql = transpile(
            "SELECT *, CLUSTER(interval) AS cid FROM peaks",
            tables=["peaks"],
            dialect="duckdb",
        )
        cursor = conn.execute(sql)
        columns = [desc[0] for desc in cursor.description]
        rows = [dict(zip(columns, row)) for row in cursor.fetchall()]

        # Assert
        assert columns == ["chrom", "start", "end", "cid"]
        assert sorted(r["cid"] for r in rows) == [1, 1, 2]

    def test_transpile_should_apply_user_except_inside_and_hide_flag_outside(self):
        """Test that a user's own EXCEPT survives while the flag is still hidden.

        Given:
            A star-projected CLUSTER whose projection carries a user EXCEPT on a column.
        When:
            Transpiling for DuckDB and executing it.
        Then:
            The user-excepted column should be dropped (applied inside the lag_calc
            subquery) and the synthesized flag should also be absent, leaving the
            remaining user columns plus the cluster id (#184).
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            "CREATE TABLE peaks "
            '(chrom VARCHAR, "start" INTEGER, "end" INTEGER, score INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?, ?)",
            [("chr1", 10, 20, 5), ("chr1", 15, 30, 5), ("chr1", 100, 110, 9)],
        )

        # Act
        sql = transpile(
            "SELECT * EXCEPT (score), CLUSTER(interval) AS cid FROM peaks",
            tables=["peaks"],
            dialect="duckdb",
        )
        columns = [desc[0] for desc in conn.execute(sql).description]

        # Assert
        assert columns == ["chrom", "start", "end", "cid"]

    def test_transpile_should_dedup_star_without_flag_when_cluster_distinct(self):
        """Test that DISTINCT over a star-projected CLUSTER dedups without the flag.

        Given:
            A CLUSTER over rows that collapse to one distinct interval, projected with
            DISTINCT over a bare star.
        When:
            Transpiling for DuckDB and executing it.
        Then:
            The output should dedup to the single distinct row rather than splitting on
            the now-hidden __giql_is_new_cluster flag — the star-leak fix (#184) also
            resolves the DISTINCT-star dedup confound noted in #181.
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?)",
            [("chr1", 10, 20), ("chr1", 10, 20), ("chr1", 10, 20)],
        )

        # Act
        sql = transpile(
            "SELECT DISTINCT *, CLUSTER(interval) AS cid FROM peaks",
            tables=["peaks"],
            dialect="duckdb",
        )
        rows = conn.execute(sql).fetchall()

        # Assert
        assert rows == [("chr1", 10, 20, 1)]

    def test_transpile_should_hide_flag_when_cluster_nested_below_root(self):
        """Test that a nested CLUSTER's flag does not leak through an enclosing star.

        Given:
            A star-projected CLUSTER nested inside a FROM-subquery, with the outer
            query projecting a bare star over it.
        When:
            Transpiling for DuckDB and executing it.
        Then:
            The nested CLUSTER should EXCEPT its synthesized flag so the enclosing
            star does not re-surface it — hiding runs for every CLUSTER expansion,
            not only a root one (#184).
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?)",
            [("chr1", 10, 20), ("chr1", 15, 30), ("chr1", 100, 110)],
        )
        query = "SELECT * FROM (SELECT *, CLUSTER(interval) AS cid FROM peaks) x"

        # Act
        sql = transpile(query, tables=["peaks"], dialect="duckdb")
        columns = [desc[0] for desc in conn.execute(sql).description]

        # Assert
        assert "__giql_is_new_cluster" not in columns
        assert columns == ["chrom", "start", "end", "cid"]

    def test_transpile_should_execute_when_cluster_projects_qualified_star(self):
        """Test that a CLUSTER over a qualified rel.* projection executes.

        Given:
            A CLUSTER whose projection is a qualified ``t.*`` over an aliased FROM.
        When:
            Transpiling for DuckDB and executing it.
        Then:
            The query should run rather than raising a binder error for the missing
            ``t`` alias — the outer FROM is the renamed __giql_lag_calc subquery, so the
            outer ``t.*`` is rewritten to a bare star that resolves against it — and the
            cluster ids should be correct (#185).
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?)",
            [("chr1", 10, 20), ("chr1", 15, 30), ("chr1", 100, 110)],
        )

        # Act
        sql = transpile(
            "SELECT t.*, CLUSTER(interval) AS cid FROM peaks t",
            tables=["peaks"],
            dialect="duckdb",
        )
        rows = conn.execute(sql).fetchall()

        # Assert
        assert rows == [("chr1", 10, 20, 1), ("chr1", 15, 30, 1), ("chr1", 100, 110, 2)]

    def test_transpile_should_not_duplicate_columns_when_qualified_star(self):
        """Test that a qualified rel.* CLUSTER does not duplicate required columns.

        Given:
            A CLUSTER over a qualified ``t.*`` projection, whose required chrom/start/end
            columns the flag-building subquery would re-inject if it failed to see the
            qualified star as a full-column projection.
        When:
            Transpiling for DuckDB and executing it.
        Then:
            The output columns should be exactly the table's columns plus the cluster id,
            with no duplicated chrom_1/start_1/end_1 columns and no synthesized flag
            (#185).
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?)",
            [("chr1", 10, 20), ("chr1", 15, 30), ("chr1", 100, 110)],
        )

        # Act
        sql = transpile(
            "SELECT t.*, CLUSTER(interval) AS cid FROM peaks t",
            tables=["peaks"],
            dialect="duckdb",
        )
        columns = [desc[0] for desc in conn.execute(sql).description]

        # Assert
        assert columns == ["chrom", "start", "end", "cid"]

    def test_transpile_should_hide_flag_when_cluster_projects_qualified_star(self):
        """Test that a qualified rel.* CLUSTER EXCEPTs the synthesized flag.

        Given:
            A top-level CLUSTER with a qualified ``t.*`` projection.
        When:
            Transpiling the query.
        Then:
            The outer star should be rewritten to a bare ``* EXCEPT`` that drops the
            synthesized __giql_is_new_cluster flag, extending the bare-star flag-hiding
            (#184) to the qualified form (#185).
        """
        # Act
        sql = transpile(
            "SELECT t.*, CLUSTER(interval) AS cid FROM peaks t", tables=["peaks"]
        )

        # Assert
        assert 'EXCEPT ("__giql_is_new_cluster")' in sql

    def test_transpile_should_not_reinject_columns_in_lag_calc_when_qualified_star(
        self,
    ):
        """Test that a qualified star does not re-inject required columns into lag_calc.

        Given:
            A CLUSTER over a qualified ``t.*`` projection, whose lag_calc subquery would
            re-inject the required chrom/start/end columns if it failed to recognize the
            qualified star as covering all columns.
        When:
            Transpiling the query.
        Then:
            The inner __giql_lag_calc projection should carry the qualified star directly
            followed by the synthesized flag, with no re-injected chrom/start/end columns
            between them — isolating the duplication defect in the emitted SQL,
            independent of any engine's tolerance for duplicate column names (#185).
        """
        # Act
        sql = transpile(
            "SELECT t.*, CLUSTER(interval) AS cid FROM peaks t", tables=["peaks"]
        )

        # Assert
        assert "t.*, CASE" in sql
        assert 't.*, "chrom"' not in sql

    def test_transpile_should_drop_qualifier_when_multipart_qualified_star(self):
        """Test that a multi-part db.t.* projection drops its qualifier at the outer.

        Given:
            A CLUSTER over a multi-part qualified ``db.t.*`` projection.
        When:
            Transpiling the query.
        Then:
            The outer projection should be a bare flag-hiding star with the dangling
            ``db.t`` qualifier dropped (CLUSTER is single-relation), while the inner
            lag_calc retains the original ``db.t.*`` against the aliased base table
            (#185).
        """
        # Act
        sql = transpile(
            "SELECT db.t.*, CLUSTER(interval) AS cid FROM peaks t", tables=["peaks"]
        )

        # Assert
        assert 'SELECT * EXCEPT ("__giql_is_new_cluster"), SUM(' in sql
        assert "db.t.*, CASE" in sql

    def test_transpile_should_apply_user_except_when_qualified_star(self):
        """Test that a user's EXCEPT on a qualified star survives while the flag hides.

        Given:
            A qualified ``t.* EXCEPT (score)`` CLUSTER projection over an aliased FROM.
        When:
            Transpiling for DuckDB and executing it.
        Then:
            The user-excepted column should be dropped (applied inside the lag_calc
            subquery, which retains the qualified projection) and the synthesized flag
            should also be absent, leaving the remaining user columns plus the cluster id
            (#185).
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            "CREATE TABLE peaks "
            '(chrom VARCHAR, "start" INTEGER, "end" INTEGER, score INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?, ?)",
            [("chr1", 10, 20, 5), ("chr1", 15, 30, 6), ("chr1", 100, 110, 9)],
        )

        # Act
        sql = transpile(
            "SELECT t.* EXCEPT (score), CLUSTER(interval) AS cid FROM peaks t",
            tables=["peaks"],
            dialect="duckdb",
        )
        columns = [desc[0] for desc in conn.execute(sql).description]

        # Assert
        assert columns == ["chrom", "start", "end", "cid"]

    def test_transpile_should_execute_when_cluster_with_limit_and_offset(self):
        """Test that CLUSTER with LIMIT/OFFSET returns the correct row window.

        Given:
            A CLUSTER query with an outer ORDER BY, LIMIT, and OFFSET over a seeded
            table.
        When:
            Transpiling the query and executing it on DuckDB.
        Then:
            The preserved LIMIT/OFFSET should return exactly the requested ordered
            window — proving the fix emits the right rows rather than the silently
            unbounded result the dropped clauses produced (#181).
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?)",
            [("chr1", 10, 20), ("chr1", 15, 30), ("chr1", 100, 110), ("chr1", 105, 120)],
        )
        query = (
            'SELECT chrom, start, "end", CLUSTER(interval) AS cid '
            "FROM peaks ORDER BY start LIMIT 2 OFFSET 1"
        )

        # Act
        sql = transpile(query, tables=["peaks"])
        rows = conn.execute(sql).fetchall()

        # Assert
        assert rows == [("chr1", 15, 30, 1), ("chr1", 100, 110, 2)]

    def test_transpile_should_execute_when_cluster_with_qualify(self):
        """Test that CLUSTER with QUALIFY filters on the cluster id.

        Given:
            A CLUSTER query with a QUALIFY predicate over the cluster-id alias.
        When:
            Transpiling the query and executing it on DuckDB.
        Then:
            The preserved QUALIFY should filter the rewritten outer query on the SUM
            window result, returning only the rows in clusters past the first (#181).
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?)",
            [("chr1", 10, 20), ("chr1", 15, 30), ("chr1", 100, 110), ("chr1", 105, 120)],
        )
        query = (
            'SELECT chrom, start, "end", CLUSTER(interval) AS cid '
            "FROM peaks QUALIFY cid > 1"
        )

        # Act
        sql = transpile(query, tables=["peaks"])
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert rows == [("chr1", 100, 110, 2), ("chr1", 105, 120, 2)]

    def test_transpile_should_preserve_with_and_limit_when_cluster_over_cte(self):
        """Test that CLUSTER over a CTE preserves both the WITH and an outer LIMIT.

        Given:
            A CLUSTER over a CTE FROM that also carries an outer LIMIT — the #174
            (WITH) and #181 (LIMIT) preservation cases composed in one query.
        When:
            Transpiling the query and executing it on DuckDB.
        Then:
            Both the enclosing WITH and the LIMIT should survive the transplant, so
            the query executes and returns exactly the requested ordered window.
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?)",
            [("chr1", 10, 20), ("chr1", 15, 30), ("chr1", 100, 110)],
        )
        query = (
            "WITH sub AS (SELECT * FROM peaks) "
            'SELECT chrom, start, "end", CLUSTER(interval) AS cid '
            "FROM sub ORDER BY start LIMIT 2"
        )

        # Act
        sql = transpile(query, tables=["peaks"])
        rows = conn.execute(sql).fetchall()

        # Assert
        assert sql.startswith("WITH sub AS (SELECT * FROM peaks)")
        assert rows == [("chr1", 10, 20, 1), ("chr1", 15, 30, 1)]

    @pytest.mark.parametrize(
        "clause",
        [
            "INTO newtbl FROM peaks",
            "FROM peaks FOR UPDATE",
            "FROM peaks SORT BY start",
            "FROM peaks WINDOW w AS (PARTITION BY chrom)",
        ],
        ids=["into", "for_update", "sort_by", "named_window"],
    )
    def test_transpile_should_raise_when_cluster_carries_unsupported_root_clause(
        self, clause
    ):
        """Test that an unhandled outer clause fails loud rather than dropping.

        Given:
            A CLUSTER query carrying a top-level clause the whole-query rewrite
            neither rebuilds nor preserves — SELECT INTO, FOR UPDATE, SORT BY, or a
            named WINDOW.
        When:
            Transpiling the query.
        Then:
            It should raise ValueError naming the unsupported clause rather than
            silently dropping it (the #174/#181 failure mode) — a loud, diagnosable
            error for a clause the transplant cannot carry.
        """
        # Act & assert
        with pytest.raises(ValueError, match="cannot carry these top-level clause"):
            transpile(
                f"SELECT chrom, CLUSTER(interval) AS cid {clause}", tables=["peaks"]
            )
