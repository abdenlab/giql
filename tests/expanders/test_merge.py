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

    @pytest.mark.parametrize(
        "query",
        [
            "SELECT *, MERGE(interval) FROM peaks",
            "SELECT MERGE(interval), * FROM peaks",
            "SELECT t.*, MERGE(interval) FROM peaks t",
            "SELECT t.*, MERGE(interval) AS m FROM peaks t",
            "SELECT db.t.*, MERGE(interval) FROM peaks t",
            'SELECT * EXCEPT ("start"), MERGE(interval) FROM peaks',
            "SELECT MERGE(interval, stranded := true), * FROM peaks",
        ],
        ids=[
            "bare-star-after",
            "bare-star-before",
            "qualified-star",
            "qualified-star-aliased-merge",
            "multipart-qualified-star",
            "star-except",
            "star-with-parameterized-merge",
        ],
    )
    def test_transpile_should_raise_when_merge_shares_select_with_star(self, query):
        """Test that a star projected alongside MERGE is rejected.

        Given:
            A SELECT projecting both a MERGE and a star — a bare ``*`` (either
            order), a qualified ``t.*``, a multi-part ``db.t.*``, a ``* EXCEPT``
            star (still an ``exp.Star``, carrying exclusion args), or a star beside
            a parameterized MERGE — optionally with the MERGE aliased.
        When:
            Transpiling the query.
        Then:
            It should raise ValueError naming the unsupported star-with-MERGE
            combination, rather than emitting non-executable SQL — a bare ``*``
            re-surfaces non-grouped columns under the synthesized GROUP BY and a
            qualified ``rel.*`` dangles an alias the aggregation no longer exposes
            (#189).
        """
        # Act & assert
        with pytest.raises(
            ValueError, match="MERGE cannot be combined with a star projection"
        ):
            transpile(query, tables=["peaks"])

    def test_transpile_should_not_reject_when_merge_shares_select_with_count_star(self):
        """Test that a COUNT(*) sibling does not trip the star guard.

        Given:
            A MERGE projected beside ``COUNT(*)`` — whose inner ``*`` is wrapped in
            an aggregate, not a bare or qualified star projection.
        When:
            Transpiling the query.
        Then:
            It should not raise the star-with-MERGE error and should expand into the
            clustered-aggregation form, proving the guard rejects only a bare or
            qualified star and not an aggregate-wrapped ``*`` (#189).
        """
        # Act
        sql = transpile(
            "SELECT MERGE(interval), COUNT(*) AS n FROM peaks", tables=["peaks"]
        )

        # Assert
        assert "AS __giql_clustered" in sql
        assert "COUNT(*)" in sql
        assert "G_I_Q_L" not in sql

    def test_transpile_should_reject_star_in_merge_own_scalar_subquery(self):
        """Test that the star guard keys on the MERGE's own SELECT scope.

        Given:
            A MERGE inside a scalar subquery whose own projection carries a star
            (``SELECT (SELECT *, MERGE(interval) FROM peaks) AS m FROM other``).
        When:
            Transpiling the query.
        Then:
            It should raise the star-with-MERGE error — the mirror of the
            FROM-subquery case — proving the guard inspects the MERGE's own
            enclosing SELECT, not just the outermost query (#189).
        """
        # Act & assert
        with pytest.raises(
            ValueError, match="MERGE cannot be combined with a star projection"
        ):
            transpile(
                "SELECT (SELECT *, MERGE(interval) FROM peaks) AS m FROM other",
                tables=["peaks", "other"],
            )

    def test_transpile_should_not_reject_star_inside_merge_from_subquery(self):
        """Test that the star guard fires only on the MERGE's own projection.

        Given:
            A MERGE whose FROM is a pass-through ``SELECT *`` subquery — a star that
            belongs to the source relation, not the MERGE's projection.
        When:
            Transpiling the query.
        Then:
            It should not raise the star-with-MERGE error and should expand into the
            clustered-aggregation form, proving the guard inspects the MERGE's own
            SELECT list rather than any star anywhere in the tree (#189).
        """
        # Arrange
        query = "SELECT MERGE(interval) FROM (SELECT * FROM peaks) x"

        # Act
        sql = transpile(query, tables=["peaks"])

        # Assert
        assert "AS __giql_clustered" in sql
        assert "G_I_Q_L" not in sql

    def test_transpile_should_dedup_grouping_key_when_merge_by_chromosome(self):
        """Test that an explicit chrom grouping key is not duplicated by MERGE.

        Given:
            The documented "merge by chromosome" shape — a MERGE projected beside an
            explicit ``chrom`` grouping key and a ``COUNT(*)`` aggregate, grouped by
            chrom — over a seeded table.
        When:
            Transpiling for DuckDB and executing it.
        Then:
            The output should carry chrom exactly once (MERGE already projects it) with
            the merged bounds and the per-region feature count, rather than surfacing
            chrom twice as the verbatim copy did (#192).
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?)",
            [("chr1", 10, 20), ("chr1", 15, 30), ("chr1", 100, 110), ("chr2", 5, 8)],
        )

        # Act
        sql = transpile(
            "SELECT chrom, MERGE(interval), COUNT(*) AS feature_count "
            "FROM peaks GROUP BY chrom ORDER BY chrom",
            tables=["peaks"],
            dialect="duckdb",
        )
        cursor = conn.execute(sql)
        columns = [desc[0] for desc in cursor.description]
        rows = cursor.fetchall()

        # Assert
        assert columns == ["chrom", "start", "end", "feature_count"]
        assert sorted(rows) == [
            ("chr1", 10, 30, 2),
            ("chr1", 100, 110, 1),
            ("chr2", 5, 8, 1),
        ]

    def test_transpile_should_dedup_strand_key_when_merge_stranded(self):
        """Test that a stranded MERGE does not duplicate an explicit strand key.

        Given:
            A stranded MERGE projected beside explicit ``chrom`` and ``strand`` grouping
            keys and a ``COUNT(*)`` aggregate, grouped by chrom and strand.
        When:
            Transpiling for DuckDB and executing it.
        Then:
            Both grouping keys should surface exactly once (the stranded MERGE already
            projects chrom and strand), with the merged bounds and per-region count
            (#192).
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            "CREATE TABLE peaks "
            '(chrom VARCHAR, "start" INTEGER, "end" INTEGER, strand VARCHAR)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?, ?)",
            [("chr1", 10, 20, "+"), ("chr1", 15, 30, "+"), ("chr1", 12, 40, "-")],
        )

        # Act
        sql = transpile(
            "SELECT chrom, strand, MERGE(interval, stranded := true), COUNT(*) AS n "
            "FROM peaks GROUP BY chrom, strand",
            tables=["peaks"],
            dialect="duckdb",
        )
        cursor = conn.execute(sql)
        columns = [desc[0] for desc in cursor.description]
        rows = cursor.fetchall()

        # Assert
        assert columns == ["chrom", "strand", "start", "end", "n"]
        assert sorted(rows) == [("chr1", "+", 10, 30, 2), ("chr1", "-", 12, 40, 1)]

    def test_transpile_should_keep_expression_over_grouping_key_when_merge(self):
        """Test that a non-aggregate expression over only grouping keys is kept.

        Given:
            A MERGE projected beside an expression computed purely from the grouping key
            (``UPPER(chrom) AS label``), grouped by chrom.
        When:
            Transpiling for DuckDB and executing it.
        Then:
            The expression should surface (it is functionally determined by the
            grouping, so it is well-formed under the GROUP BY) rather than being
            rejected as a raw column (#192).
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?)",
            [("chr1", 10, 20), ("chr1", 15, 30)],
        )

        # Act
        sql = transpile(
            "SELECT UPPER(chrom) AS label, MERGE(interval) FROM peaks GROUP BY chrom",
            tables=["peaks"],
            dialect="duckdb",
        )
        cursor = conn.execute(sql)
        columns = [d[0] for d in cursor.description]
        rows = [dict(zip(columns, r)) for r in cursor.fetchall()]

        # Assert
        assert rows == [{"chrom": "chr1", "start": 10, "end": 30, "label": "CHR1"}]

    @pytest.mark.parametrize(
        "query",
        [
            pytest.param(
                "SELECT chrom, score, MERGE(interval) FROM peaks GROUP BY chrom",
                id="raw-column-score",
            ),
            pytest.param(
                'SELECT "start", MERGE(interval) FROM peaks',
                id="raw-bounds-column",
            ),
            pytest.param(
                "SELECT MERGE(interval), COUNT(*) + score AS bad "
                "FROM peaks GROUP BY chrom",
                id="mixed-aggregate-and-raw",
            ),
        ],
    )
    def test_transpile_should_raise_when_merge_projects_raw_non_key_column(self, query):
        """Test that a raw non-aggregate, non-grouping-key projection is rejected.

        Given:
            A MERGE projected beside a column referencing, outside any aggregate, a
            column that is neither a grouping key nor aggregated — a bare ``score``, a
            raw interval-bound column, or a mixed ``COUNT(*) + score`` whose ``score``
            is added outside the aggregate.
        When:
            Transpiling the query.
        Then:
            It should raise a ValueError naming the non-aggregated column, rather than
            copying it verbatim into the aggregation where it is neither grouped nor
            aggregated (a silent duplicate or a raw engine binder error) (#192).
        """
        # Act & assert
        with pytest.raises(
            ValueError, match="MERGE cannot project the non-aggregated column"
        ):
            transpile(query, tables=["peaks"])

    @pytest.mark.parametrize(
        "query",
        [
            pytest.param(
                "SELECT chrom AS start, MERGE(interval) FROM peaks",
                id="alias-onto-start",
            ),
            pytest.param(
                "SELECT chrom AS end, MERGE(interval) FROM peaks",
                id="alias-onto-end",
            ),
            pytest.param(
                "SELECT UPPER(chrom) AS chrom, MERGE(interval) "
                "FROM peaks GROUP BY chrom",
                id="expr-onto-chrom",
            ),
        ],
    )
    def test_transpile_should_raise_when_projection_collides_with_synthesized_name(
        self, query
    ):
        """Test that an item aliased onto a synthesized column name is rejected.

        Given:
            A MERGE projected beside an item whose output name collides with a column
            MERGE synthesizes — ``chrom AS start`` / ``chrom AS end`` (onto the
            aggregated bounds) or ``UPPER(chrom) AS chrom`` (a different value onto the
            grouping key).
        When:
            Transpiling the query.
        Then:
            It should raise a ValueError naming the collision, rather than emitting two
            columns of that name — which silently duplicates one and, for ``start`` /
            ``end``, hijacks the default ``ORDER BY "chrom", "start"`` and returns rows
            in the wrong order (#192).
        """
        # Act & assert
        with pytest.raises(
            ValueError, match="collides with the chrom/start/end columns"
        ):
            transpile(query, tables=["peaks"])

    def test_transpile_should_dedup_when_grouping_key_self_aliased(self):
        """Test that a grouping key aliased to its own name is not duplicated.

        Given:
            A MERGE projected beside a no-op self-aliased grouping key
            (``chrom AS chrom``).
        When:
            Transpiling for DuckDB and executing it.
        Then:
            chrom should surface exactly once (MERGE already projects it), proving the
            dedup catches the aliased form and not only a bare column (#192).
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?)",
            [("chr1", 10, 20), ("chr1", 15, 30)],
        )

        # Act
        sql = transpile(
            "SELECT chrom AS chrom, MERGE(interval) FROM peaks",
            tables=["peaks"],
            dialect="duckdb",
        )
        cursor = conn.execute(sql)
        columns = [desc[0] for desc in cursor.description]
        rows = cursor.fetchall()

        # Assert
        assert columns == ["chrom", "start", "end"]
        assert rows == [("chr1", 10, 30)]

    def test_transpile_should_raise_when_group_by_conflicts_with_merge_grouping(self):
        """Test that a GROUP BY over a non-grouping-key column is rejected.

        Given:
            A MERGE whose enclosing query groups by a column MERGE does not group by
            (``GROUP BY score``).
        When:
            Transpiling the query.
        Then:
            It should raise a ValueError naming the offending GROUP BY column, rather
            than silently discarding the user's GROUP BY for MERGE's synthesized one
            (#192).
        """
        # Act & assert
        with pytest.raises(ValueError, match="MERGE cannot honor GROUP BY score"):
            transpile(
                "SELECT MERGE(interval), COUNT(*) FROM peaks GROUP BY score",
                tables=["peaks"],
            )

    @pytest.mark.parametrize(
        "query",
        [
            pytest.param(
                "SELECT MERGE(interval), -MAX(start) AS Start FROM peaks",
                id="aggregate-alias-Start",
            ),
            pytest.param(
                "SELECT 1 AS START, MERGE(interval) FROM peaks",
                id="literal-alias-START",
            ),
            pytest.param(
                "SELECT chrom AS End, MERGE(interval) FROM peaks",
                id="grouping-key-alias-End",
            ),
        ],
    )
    def test_transpile_should_raise_when_case_variant_alias_collides_with_synthesized(
        self, query
    ):
        """Test that a case-variant alias colliding with a synthesized name is rejected.

        Given:
            A MERGE projected beside an item whose alias differs only in case from a
            column MERGE synthesizes — ``AS Start`` / ``AS START`` / ``AS End`` against
            the aggregated ``start`` / ``end`` bounds.
        When:
            Transpiling the query.
        Then:
            It should raise a ValueError naming the collision: SQL binds unquoted
            identifiers case-insensitively, so the alias folds onto the synthesized
            column and would otherwise duplicate it and silently hijack the default
            ``ORDER BY "chrom", "start"``, returning rows in the wrong order (#192).
        """
        # Act & assert
        with pytest.raises(
            ValueError, match="collides with the chrom/start/end columns"
        ):
            transpile(query, tables=["peaks"])

    @pytest.mark.parametrize(
        "query",
        [
            pytest.param(
                "SELECT UPPER(chrom) AS c, MERGE(interval) FROM peaks "
                "GROUP BY UPPER(chrom)",
                id="group-by-expression",
            ),
            pytest.param(
                "SELECT chrom, MERGE(interval) FROM peaks GROUP BY 1",
                id="group-by-ordinal",
            ),
        ],
    )
    def test_transpile_should_raise_when_group_by_is_expression_or_ordinal(self, query):
        """Test that a GROUP BY expression or ordinal is rejected, not silently dropped.

        Given:
            A MERGE whose enclosing query groups by an expression (``GROUP BY
            UPPER(chrom)``) or an ordinal (``GROUP BY 1``) rather than a bare
            grouping-key column.
        When:
            Transpiling the query.
        Then:
            It should raise a ValueError: MERGE cannot verify a non-column GROUP BY names
            a grouping key, so honoring it would silently substitute MERGE's synthesized
            grouping and change the result (a case-folding ``GROUP BY UPPER(chrom)``
            would no longer merge case-insensitively) (#192).
        """
        # Act & assert
        with pytest.raises(ValueError, match="MERGE cannot honor GROUP BY"):
            transpile(query, tables=["peaks"])

    def test_transpile_should_accept_case_variant_grouping_key(self):
        """Test that a grouping key referenced in a different case is honored.

        Given:
            A MERGE grouped by ``CHROM`` — the ``chrom`` grouping key in a different
            case — over a seeded table.
        When:
            Transpiling for DuckDB and executing it.
        Then:
            The GROUP BY should be accepted (SQL folds the unquoted identifier onto the
            grouping key) and the merged regions returned, rather than falsely rejected
            as a non-grouping column (#192).
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
            "SELECT chrom, MERGE(interval) FROM peaks GROUP BY CHROM",
            tables=["peaks"],
            dialect="duckdb",
        )
        cursor = conn.execute(sql)
        columns = [desc[0] for desc in cursor.description]
        rows = cursor.fetchall()

        # Assert
        assert columns == ["chrom", "start", "end"]
        assert sorted(rows) == [("chr1", 10, 30), ("chr1", 100, 110)]

    def test_transpile_should_keep_unaliased_cast_over_grouping_key(self):
        """Test that an unaliased cast over the grouping key is kept, not rejected.

        Given:
            A MERGE projected beside an unaliased ``CAST(chrom AS VARCHAR)`` over the
            grouping key, grouped by chrom.
        When:
            Transpiling for DuckDB and executing it.
        Then:
            The cast should surface (it is functionally determined by the grouping key,
            and its emitted column name is the cast text, not ``chrom``) rather than
            falsely rejected as colliding with the synthesized ``chrom`` on the strength
            of sqlglot's output-name heuristic (#192).
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?)",
            [("chr1", 10, 20), ("chr1", 15, 30)],
        )

        # Act
        sql = transpile(
            "SELECT CAST(chrom AS VARCHAR), MERGE(interval) FROM peaks GROUP BY chrom",
            tables=["peaks"],
            dialect="duckdb",
        )
        rows = conn.execute(sql).fetchall()

        # Assert
        assert rows == [("chr1", 10, 30, "chr1")]

    def test_transpile_should_raise_when_windowed_aggregate_over_raw_column(self):
        """Test that a window aggregate over a raw column is rejected loudly.

        Given:
            A MERGE projected beside a window aggregate whose argument is a raw
            non-grouping column (``SUM(score) OVER (PARTITION BY chrom)``).
        When:
            Transpiling the query.
        Then:
            It should raise a ValueError naming ``score``: a window aggregate runs
            after the GROUP BY and does not collapse the group, so its raw argument is
            still ungrouped — MERGE rejects it with a clear diagnostic rather than
            emitting SQL the engine rejects with a raw binder error (#192).
        """
        # Act & assert
        with pytest.raises(
            ValueError, match="MERGE cannot project the non-aggregated column"
        ):
            transpile(
                "SELECT SUM(score) OVER (PARTITION BY chrom) AS s, MERGE(interval) "
                "FROM peaks GROUP BY chrom",
                tables=["peaks"],
            )

    def test_transpile_should_honor_having_over_merged_regions(self):
        """Test that a HAVING filters merged regions instead of being dropped.

        Given:
            A MERGE with ``HAVING COUNT(*) > 1`` over a seeded table where one merged
            region is built from two input intervals and another from a single interval.
        When:
            Transpiling for DuckDB and executing it.
        Then:
            Only the multi-interval region should be returned — the HAVING is honored
            against the per-region grouping rather than silently dropped by the rewrite,
            which had returned the unfiltered regions (#192).
        """
        # Arrange
        duckdb = pytest.importorskip("duckdb")
        conn = duckdb.connect(":memory:")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.executemany(
            "INSERT INTO peaks VALUES (?, ?, ?)",
            [("chr1", 0, 10), ("chr1", 5, 15), ("chr1", 100, 110)],
        )

        # Act
        sql = transpile(
            "SELECT chrom, COUNT(*) AS n, MERGE(interval) FROM peaks "
            "HAVING COUNT(*) > 1",
            tables=["peaks"],
            dialect="duckdb",
        )
        cursor = conn.execute(sql)
        columns = [desc[0] for desc in cursor.description]
        rows = cursor.fetchall()

        # Assert
        assert columns == ["chrom", "start", "end", "n"]
        assert rows == [("chr1", 0, 15, 2)]

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
        ["MERGE(interval) AS m, COUNT(*) AS n", "MERGE(interval)"],
        ids=["aliased", "bare"],
    )
    def test_transpile_should_expand_spatial_predicate_copied_into_clustered(
        self, projection, predicate_op
    ):
        """Test that a spatial WHERE predicate survives the MERGE rewrite.

        Given:
            A MERGE query (aliased alongside an extra aggregate, or bare) whose
            WHERE filters on a spatial predicate, which the rewrite copies into the
            inner clustered subquery.
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

    def test_transpile_should_not_hide_flag_when_merge(self):
        """Test that MERGE's composed CLUSTER does not add a flag-hiding star.

        Given:
            A plain MERGE query.
        When:
            Transpiling the query.
        Then:
            The emitted SQL should not carry a flag-hiding ``* EXCEPT`` / ``EXCLUDE``
            wrapper: MERGE composes CLUSTER with hide_reserved disabled because its
            explicit outer projection never surfaces the __giql_is_new_cluster flag,
            so no needless exclusion is added to the intermediate subquery (#184).
        """
        # Act
        sql = transpile("SELECT MERGE(interval) FROM peaks", tables=["peaks"])

        # Assert
        assert 'EXCEPT ("__giql_is_new_cluster")' not in sql
