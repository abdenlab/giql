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
        assert "AS __giql_clustered" in sql
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
        assert "AS __giql_clustered" in sql

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
        assert "AS __giql_clustered" in sql
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
        assert "AS __giql_clustered" in sql

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

    def test_transpile_should_namespace_synthesized_merge_identifiers(self):
        """Test that MERGE's synthesized identifiers carry the reserved prefix.

        Given:
            A plain MERGE query.
        When:
            Transpiling the query.
        Then:
            The emitted SQL should alias the synthesized helpers as the reserved
            __giql_clustered, __giql_lag_calc, and __giql_is_new_cluster identifiers
            (the latter two composed from CLUSTER) and never emit the bare forms,
            guarding #161 against a partial revert that would re-collide with user
            columns named clustered / lag_calc / is_new_cluster.
        """
        # Arrange
        query = "SELECT MERGE(interval) FROM peaks"

        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        assert "AS __giql_clustered" in sql
        assert "AS __giql_lag_calc" in sql
        assert "AS __giql_is_new_cluster" in sql
        assert "AS clustered" not in sql
        assert "AS lag_calc" not in sql
        assert "AS is_new_cluster" not in sql

    def test_transpile_should_disambiguate_merge_when_source_has_is_new_cluster_column(
        self,
    ):
        """Test that MERGE merges correctly beside a user is_new_cluster column.

        Given:
            A table carrying a user column literally named is_new_cluster, seeded
            with a poison value. MERGE always composes CLUSTER over ``SELECT *``, so
            pre-#161 the outer SUM bound to the user column and produced wrong
            cluster ids, silently splitting rows that should merge.
        When:
            Transpiling the MERGE and executing it on DuckDB.
        Then:
            Overlapping intervals should collapse into one merged row and the
            separate interval into another, proving the synthesized flag now lives
            in the reserved __giql_ namespace and no longer collides.
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
                ("chr1", 1, 5, 7),  # overlaps the next -> one cluster
                ("chr1", 3, 8, 7),
                ("chr1", 20, 25, 7),  # separate -> its own cluster
            ],
        )

        # Act
        sql = transpile("SELECT MERGE(interval) FROM intervals", tables=["intervals"])
        cursor = conn.execute(sql)
        rows = cursor.fetchall()

        # Assert
        assert rows == [("chr1", 1, 8), ("chr1", 20, 25)]

    @pytest.mark.parametrize(
        "query, expected_with",
        [
            (
                "WITH sub AS (SELECT * FROM peaks) SELECT MERGE(interval) FROM sub",
                "WITH sub AS (SELECT * FROM peaks)",
            ),
            (
                "WITH a AS (SELECT * FROM peaks), sub AS (SELECT * FROM a) "
                "SELECT MERGE(interval) FROM sub",
                "WITH a AS (SELECT * FROM peaks), sub AS (SELECT * FROM a)",
            ),
        ],
        ids=["single_cte", "chained_ctes"],
    )
    def test_transpile_should_preserve_with_clause_when_merge_over_cte(
        self, query, expected_with
    ):
        """Test that MERGE over a CTE FROM keeps the enclosing WITH clause.

        Given:
            A MERGE whose FROM references a CTE — a single CTE, and a chain of two
            CTEs where the driving relation is defined in terms of the first.
        When:
            Transpiling the query.
        Then:
            The emitted SQL should retain the enclosing WITH clause verbatim so the
            rewritten __giql_clustered subquery's ``FROM sub`` reference still
            resolves, rather than dropping the WITH and dangling an undefined CTE
            reference — the pre-#174 transplant defect.
        """
        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        assert sql.startswith(expected_with)
        assert "AS __giql_clustered" in sql
        assert "FROM sub) AS __giql_lag_calc" in sql

    def test_transpile_should_execute_when_merge_over_cte(self):
        """Test that MERGE over a CTE FROM emits executable SQL.

        Given:
            A MERGE whose FROM references a pass-through CTE over a seeded table.
        When:
            Transpiling the query and executing it on DuckDB.
        Then:
            The query should execute (the preserved WITH keeps the CTE reference
            resolvable) and collapse overlapping intervals into merged rows,
            proving the #174 fix emits runnable SQL rather than referencing a
            dropped CTE.
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
                ("chr1", 1, 5),  # overlaps the next -> one merged row
                ("chr1", 3, 8),
                ("chr1", 20, 25),  # separate -> its own row
            ],
        )
        query = "WITH sub AS (SELECT * FROM peaks) SELECT MERGE(interval) FROM sub"

        # Act
        sql = transpile(query, tables=["peaks"])
        rows = conn.execute(sql).fetchall()

        # Assert
        assert rows == [("chr1", 1, 8), ("chr1", 20, 25)]

    def test_transpile_should_not_hoist_outer_with_when_merge_nested_in_cte_body(self):
        """Test that a MERGE inside a CTE body keeps the outer WITH at the root.

        Given:
            A MERGE nested inside a CTE body, under a root SELECT that carries the
            enclosing WITH.
        When:
            Transpiling the query and executing it on DuckDB.
        Then:
            The enclosing WITH should stay at the root, the MERGE should expand into
            the __giql_clustered form inside the CTE body, and the query should
            execute — the inner transplant reads only its own SELECT's (absent) WITH,
            so it does not hoist the ancestor WITH.
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?)",
            [("chr1", 1, 5), ("chr1", 3, 8), ("chr1", 20, 25)],
        )
        query = "WITH sub AS (SELECT MERGE(interval) FROM peaks) SELECT * FROM sub"

        # Act
        sql = transpile(query, tables=["peaks"])
        rows = conn.execute(sql).fetchall()

        # Assert
        assert sql.startswith("WITH sub AS (")
        assert "AS __giql_clustered" in sql
        assert sql.rstrip().endswith("SELECT * FROM sub")
        assert rows == [("chr1", 1, 8), ("chr1", 20, 25)]

    def test_transpile_should_not_emit_with_when_merge_over_bare_table(self):
        """Test that a MERGE over a bare table grows no spurious WITH.

        Given:
            A MERGE over a bare registered table with no enclosing WITH.
        When:
            Transpiling the query.
        Then:
            The emitted SQL should not begin with a WITH clause — the transplant
            WITH-preservation branch is a no-op on the common, CTE-free path (#174).
        """
        # Act
        sql = transpile("SELECT MERGE(interval) FROM peaks", tables=["peaks"])

        # Assert
        assert not sql.lstrip().upper().startswith("WITH")

    @pytest.mark.parametrize(
        "clause, expected",
        [
            ("LIMIT 1", "LIMIT 1"),
            ("LIMIT 5 OFFSET 1", "OFFSET 1"),
        ],
        ids=["limit", "offset"],
    )
    def test_transpile_should_preserve_result_clause_when_merge(self, clause, expected):
        """Test that MERGE preserves an outer result-shaping clause.

        Given:
            A MERGE query carrying an outer LIMIT / OFFSET clause.
        When:
            Transpiling the query.
        Then:
            The clause should survive on the rewritten outer aggregation query rather
            than being silently dropped by the whole-query transplant (#181).
        """
        # Act
        sql = transpile(f"SELECT MERGE(interval) FROM peaks {clause}", tables=["peaks"])

        # Assert
        assert expected in sql

    def test_transpile_should_execute_when_merge_with_limit_and_offset(self):
        """Test that MERGE with LIMIT/OFFSET returns the correct merged-row window.

        Given:
            A MERGE query with an outer LIMIT and OFFSET over a seeded table. MERGE
            already emits an ORDER BY chrom, start, so the window is deterministic.
        When:
            Transpiling the query and executing it on DuckDB.
        Then:
            The preserved LIMIT/OFFSET should return exactly the requested slice of
            the ordered merged rows rather than the silently unbounded result the
            dropped clauses produced (#181).
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
        query = "SELECT MERGE(interval) FROM peaks LIMIT 5 OFFSET 1"

        # Act
        sql = transpile(query, tables=["peaks"])
        rows = conn.execute(sql).fetchall()

        # Assert
        assert rows == [("chr1", 100, 120)]

    def test_transpile_should_default_order_by_chrom_start_when_merge_unordered(self):
        """Test that MERGE keeps its chrom, start default ORDER BY without a user one.

        Given:
            A MERGE query with no user ORDER BY.
        When:
            Transpiling the query.
        Then:
            The emitted SQL should carry MERGE's default ORDER BY on chromosome then
            start, so unordered MERGE output stays deterministic.
        """
        # Act
        sql = transpile("SELECT MERGE(interval) FROM peaks", tables=["peaks"])

        # Assert
        assert 'ORDER BY "chrom" NULLS LAST, "start" NULLS LAST' in sql

    def test_transpile_should_honor_user_order_by_when_merge_with_limit(self):
        """Test that MERGE honors the user's ORDER BY under an outer LIMIT.

        Given:
            A MERGE query with a user ORDER BY on the merged end descending and an
            outer LIMIT of one row, over a seeded table.
        When:
            Transpiling the query and executing it on DuckDB.
        Then:
            MERGE should order the merged rows by the user's ORDER BY (not its fixed
            chrom, start default) so the preserved LIMIT returns the largest-end
            merged row — otherwise the preserved LIMIT would silently return the
            wrong row over MERGE's forced ordering (#181).
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?)",
            [("chr1", 10, 20), ("chr1", 100, 110), ("chr1", 5000, 6000)],
        )
        query = 'SELECT MERGE(interval) FROM peaks ORDER BY "end" DESC LIMIT 1'

        # Act
        sql = transpile(query, tables=["peaks"])
        rows = conn.execute(sql).fetchall()

        # Assert
        assert rows == [("chr1", 5000, 6000)]

    # DISTINCT/QUALIFY over MERGE: transplant re-attaches every preserved clause
    # operator-agnostically, so CLUSTER's DISTINCT/QUALIFY tests already exercise the
    # shared re-attach path. MERGE's outer query is a GROUP BY aggregation, so a
    # DISTINCT there is redundant-but-valid (asserted below); a non-windowed QUALIFY
    # is ill-formed over an aggregation and is faithfully passed through to fail at the
    # engine rather than silently drop — hence no MERGE+QUALIFY executability test.
    def test_transpile_should_preserve_distinct_when_merge(self):
        """Test that MERGE preserves an outer DISTINCT.

        Given:
            A MERGE query with a DISTINCT projection.
        When:
            Transpiling the query.
        Then:
            The rewritten outer aggregation query should keep SELECT DISTINCT rather
            than dropping it in the transplant (#181).
        """
        # Act
        sql = transpile("SELECT DISTINCT MERGE(interval) FROM peaks", tables=["peaks"])

        # Assert
        assert "SELECT DISTINCT" in sql
