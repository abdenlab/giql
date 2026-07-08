"""Tests for the DuckDB IEJoin dialect path on column-to-column INTERSECTS joins."""

import re

import pytest

from giql import Table
from giql import transpile

duckdb = pytest.importorskip("duckdb")
hypothesis = pytest.importorskip("hypothesis")

from hypothesis import HealthCheck  # noqa: E402
from hypothesis import given  # noqa: E402
from hypothesis import settings  # noqa: E402
from hypothesis import strategies as st  # noqa: E402


def _make_table(conn, name: str, rows: list[tuple]) -> None:
    conn.execute(f"""
        CREATE TABLE {name} (
            chrom VARCHAR,
            "start" INTEGER,
            "end" INTEGER,
            name VARCHAR,
            score INTEGER,
            strand VARCHAR
        )
    """)
    if rows:
        conn.executemany(f"INSERT INTO {name} VALUES (?, ?, ?, ?, ?, ?)", rows)


def _python_overlap(peaks: list[tuple], genes: list[tuple]) -> list[tuple]:
    """Reference half-open overlap on (chrom, start, end) triples.

    Returns a *sorted multiset* (list) so the caller can assert exact row
    multiplicity — a `set` would collapse duplicate input rows and hide
    multiplicity bugs in the dialect's rewrite.
    """
    out: list[tuple] = []
    for pc, ps, pe in peaks:
        for gc, gs, ge in genes:
            if pc == gc and pe > gs and ge > ps:
                out.append((pc, ps, pe, gc, gs, ge))
    return sorted(out)


def _python_semi_overlap(peaks: list[tuple], genes: list[tuple]) -> list[tuple]:
    """Reference SEMI-join overlap: distinct left rows with any match."""
    matched: set[tuple] = set()
    for pc, ps, pe in peaks:
        for gc, gs, ge in genes:
            if pc == gc and pe > gs and ge > ps:
                matched.add((pc, ps, pe))
                break
    return sorted(matched)


def _python_anti_overlap(peaks: list[tuple], genes: list[tuple]) -> list[tuple]:
    """Reference ANTI-join overlap: distinct left rows with no match."""
    matched = set(_python_semi_overlap(peaks, genes))
    return sorted({(pc, ps, pe) for (pc, ps, pe) in peaks} - matched)


def _python_count_overlaps(peaks: list[tuple], genes: list[tuple]) -> list[tuple]:
    """Reference count_overlaps as ``(chrom, start, end, count)`` per distinct left key.

    Mirrors ``COUNT(b.col)`` under ``GROUP BY`` on the left keys: each distinct
    left interval's count is the number of overlapping right intervals summed
    across the duplicate left rows sharing that key (so a key appearing twice
    with one overlap counts 2, matching SQL group-by-over-a-LEFT-join). Left keys
    with no overlap are present with count 0.
    """
    out: list[tuple] = []
    for key in set(peaks):
        pc, ps, pe = key
        duplicates = sum(1 for p in peaks if p == key)
        b_overlaps = sum(1 for (gc, gs, ge) in genes if pc == gc and pe > gs and ge > ps)
        out.append((pc, ps, pe, duplicates * b_overlaps))
    return sorted(out)


def _explain_dynamic_sql(conn, sql: str) -> str:
    """Return the DuckDB EXPLAIN text of the per-chromosome dynamic SQL for *sql*.

    Runs the leading ``SET VARIABLE`` statement, reads the aggregated
    per-chromosome ``UNION ALL`` string back out of the session variable, and
    returns the plain EXPLAIN for it. The outer ``query(getvariable(...))``
    wrapper is opaque to EXPLAIN, so the dynamic SQL must be planned directly.
    """
    statements = [s for s in sql.split(";\n") if s.strip()]
    for statement in statements[:-1]:
        conn.execute(statement)
    var_name = re.search(r"__giql_iejoin_[0-9a-f]+", sql).group(0)
    dynamic_sql = conn.execute(f"SELECT getvariable('{var_name}')").fetchone()[0]
    return conn.execute("EXPLAIN " + dynamic_sql).fetchone()[1]


@pytest.fixture
def conn():
    c = duckdb.connect(":memory:")
    yield c
    c.close()


@pytest.fixture
def peaks_genes(conn):
    _make_table(
        conn,
        "peaks",
        [
            ("chr1", 100, 200, "p1", 10, "+"),
            ("chr1", 300, 400, "p2", 20, "+"),
            ("chr1", 500, 600, "p3", 25, "+"),
            ("chr2", 100, 200, "p4", 30, "-"),
            ("chr2", 800, 900, "p5", 35, "-"),
        ],
    )
    _make_table(
        conn,
        "genes",
        [
            ("chr1", 150, 250, "g1", 1, "+"),
            ("chr1", 500, 600, "g2", 2, "-"),
            ("chr1", 700, 800, "g3", 3, "+"),
            ("chr2", 50, 150, "g4", 4, "-"),
            ("chr2", 250, 350, "g5", 5, "+"),
        ],
    )
    return conn


# Truth table for peaks_genes overlapping pairs (half-open):
# - chr1: peaks[100,200) vs genes[150,250) -> overlap
# - chr1: peaks[500,600) vs genes[500,600) -> overlap
# - chr2: peaks[100,200) vs genes[50,150)  -> overlap
_EXPECTED_OVERLAPS_PEAKS_GENES = [
    ("chr1", 100, 200, "chr1", 150, 250),
    ("chr1", 500, 600, "chr1", 500, 600),
    ("chr2", 100, 200, "chr2", 50, 150),
]


class TestTranspileDuckDBIEJoinSQLStructure:
    """SQL-shape and routing assertions for the DuckDB IEJoin dialect path."""

    def test_transpile_should_keep_naive_predicate_path_when_dialect_is_none(self):
        """Test that omitting ``dialect`` keeps the generic naive-predicate plan.

        Given:
            A column-to-column INTERSECTS JOIN.
        When:
            ``transpile`` is called without specifying ``dialect``.
        Then:
            It should emit SQL that contains no ``SET VARIABLE`` /
            ``getvariable`` scaffolding (the IEJoin path is opt-in).
        """
        # Arrange
        query = """
            SELECT a.*, b.*
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """

        # Act
        sql = transpile(query, tables=["peaks", "genes"])

        # Assert
        assert "SET VARIABLE" not in sql.upper()
        assert "getvariable" not in sql

    def test_transpile_should_route_outer_join_intersects_to_naive_predicate_plan_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that LEFT JOIN INTERSECTS falls back to the naive-predicate plan.

        Given:
            A LEFT JOIN INTERSECTS query and a peak with no overlapping
            gene on the join chromosome.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            It should preserve LEFT JOIN semantics by returning the
            unmatched peak with NULL gene columns (which the IEJoin path
            cannot do, so the naive-predicate fallback fires).
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 0, "+"),
                ("chr1", 5000, 6000, "p_lonely", 0, "+"),
            ],
        )
        _make_table(conn, "genes", [("chr1", 150, 250, "g1", 0, "+")])
        sql = transpile(
            """
            SELECT a.start AS a_start, b.start AS b_start
            FROM peaks a
            LEFT JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall(), key=lambda r: r[0])

        # Assert
        assert rows == [(100, 150), (5000, None)]

    def test_query_should_honor_extra_join_on_predicate_alongside_where_intersects_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that an extra non-INTERSECTS join predicate is honored.

        Given:
            A query whose JOIN ON carries a non-INTERSECTS predicate
            (``a.score > b.score``) alongside a WHERE INTERSECTS, and
            input rows where one pair overlaps but violates the score
            predicate.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            The dialect path should inline the extra predicate into each
            per-chromosome subquery's join ON, excluding the pair that
            violates the predicate and including the pair that satisfies
            it. (Note: the naive-predicate plan still drops this predicate per
            bug #94; the dialect now sidesteps that bug entirely by
            handling extras directly.)
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p_high", 100, "+"),
                ("chr1", 300, 400, "p_low", 1, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g_low", 5, "+"),
                ("chr1", 350, 450, "g_high", 50, "+"),
            ],
        )
        sql = transpile(
            """
            SELECT a.start AS a_start, b.start AS b_start
            FROM peaks a JOIN genes b ON a.score > b.score
            WHERE a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        # Only ``p_high`` (score 100) overlaps ``g_low`` (score 5) AND beats
        # its score; ``p_low`` (score 1) overlaps ``g_high`` (score 50) but
        # fails the score predicate.
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert rows == [(100, 150)]

    def test_transpile_should_inline_extra_join_on_predicate_when_dialect_is_duckdb(
        self,
    ):
        """Test that an extra ON predicate appears inside each per-chromosome subquery.

        Given:
            A column-to-column INTERSECTS INNER join whose ON also
            carries ``a.score > 20`` (a single-side predicate).
        When:
            The query is transpiled with ``dialect='duckdb'``.
        Then:
            The emitted SQL should engage the dialect path (no fallback)
            and inline the extra predicate into the per-chromosome
            subquery's ON, rewritten against the inner ``a``/``b`` scope.
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.start FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval AND a.score > 20",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        # The score predicate must show up inside the dynamic SQL literal so
        # DuckDB evaluates it inside each per-chromosome IEJoin candidate set.
        assert "a.score > 20" in sql

    def test_query_should_filter_by_extra_where_predicate_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that a WHERE predicate alongside an ON-INTERSECTS filters rows.

        Given:
            Two overlap rows on the same chromosome and an extra
            ``WHERE b.score < 30`` predicate that should exclude one
            of them.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            Only the row whose ``b.score`` satisfies the predicate is
            returned.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),
                ("chr1", 300, 400, "p2", 2, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g_low", 5, "+"),
                ("chr1", 350, 450, "g_high", 50, "+"),
            ],
        )
        sql = transpile(
            "SELECT a.start AS s FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "WHERE b.score < 30",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert rows == [(100,)]

    def test_query_should_combine_multiple_extra_predicates_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that two extra predicates AND together correctly under the dialect.

        Given:
            Three peaks/genes with varying scores, and ``ON ... AND a.score > 20``
            plus ``WHERE b.score < 40`` — only one pair satisfies both.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            Only the pair satisfying both extras is returned.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 10, "+"),  # score too low for ON
                ("chr1", 300, 400, "p2", 25, "+"),  # score passes ON
                ("chr1", 500, 600, "p3", 30, "+"),  # score passes ON
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 5, "+"),
                ("chr1", 350, 450, "g2", 5, "+"),  # b.score passes WHERE
                ("chr1", 550, 650, "g3", 50, "+"),  # b.score fails WHERE
            ],
        )
        sql = transpile(
            "SELECT a.start AS s FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval AND a.score > 20 "
            "WHERE b.score < 40",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert rows == [(300,)]

    def test_transpile_should_apply_anti_join_where_residual_as_outer_wrapper_filter(
        self,
    ):
        """Test that an ANTI-join WHERE residual filters outside the per-chrom join.

        Given:
            An ``ANTI JOIN ... ON a.interval INTERSECTS b.interval`` whose
            top-level ``WHERE a.start >= 550`` is a post-join filter.
        When:
            The query is transpiled with ``dialect='duckdb'``.
        Then:
            The residual is applied as a ``WHERE`` on the outer wrapper
            relation, never inlined into the per-chromosome ``NOT EXISTS``
            subquery (which would invert the anti-join for rows failing the
            residual, #200).
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.start FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval "
            "WHERE a.start >= 550",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        per_chrom_template, _, outer_select = sql.partition("query(getvariable")
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert "NOT EXISTS" in per_chrom_template
        # The residual predicate sits in the outer wrapper filter, never in the
        # per-chromosome NOT EXISTS template that DuckDB expands via string_agg.
        # Anchor on ``>= 550`` (the rendered comparison) rather than the bare
        # ``550`` literal: the per-chromosome template embeds the random
        # ``uuid4`` var-name token, whose hex could coincidentally contain the
        # substring ``550`` — but never the operator-and-space form.
        assert "AS __giql_iejoin_wrapper WHERE" in outer_select
        assert ">= 550" in outer_select
        assert ">= 550" not in per_chrom_template

    def test_transpile_should_raise_when_anti_join_where_residual_references_right_side(
        self,
    ):
        """Test that an ANTI-join WHERE residual on the right side raises.

        Given:
            An ``ANTI JOIN`` whose ``WHERE`` references ``b.score`` — a right-
            side column the left-only outer wrapper filter cannot resolve.
        When:
            The query is transpiled with ``dialect='duckdb'``.
        Then:
            A ``ValueError`` is raised rather than the residual being silently
            inlined or mis-bound.
        """
        # Arrange & act & assert
        with pytest.raises(ValueError, match="only projects left-side columns"):
            transpile(
                "SELECT a.start FROM peaks a "
                "ANTI JOIN genes b ON a.interval INTERSECTS b.interval "
                "WHERE b.score > 5",
                tables=["peaks", "genes"],
                dialect="duckdb",
            )

    def test_transpile_should_raise_when_semi_join_where_residual_references_right_side(
        self,
    ):
        """Test that a SEMI-join WHERE residual on the right side raises.

        Given:
            A ``SEMI JOIN`` whose ``WHERE`` references ``b.score`` — a right-
            side column the left-only outer wrapper filter cannot resolve.
        When:
            The query is transpiled with ``dialect='duckdb'``.
        Then:
            A ``ValueError`` is raised, mirroring the ANTI case (both are
            left-only outer projections).
        """
        # Arrange & act & assert
        with pytest.raises(ValueError, match="only projects left-side columns"):
            transpile(
                "SELECT a.start FROM peaks a "
                "SEMI JOIN genes b ON a.interval INTERSECTS b.interval "
                "WHERE b.score > 5",
                tables=["peaks", "genes"],
                dialect="duckdb",
            )

    def test_transpile_should_route_to_naive_predicate_plan_when_extra_predicate_uses_or(
        self,
    ):
        """Test that an OR-wrapped predicate falls back to the naive-predicate plan.

        Given:
            A column-to-column INTERSECTS INNER join whose JOIN ON wraps
            the INTERSECTS and a sibling predicate in an OR.
        When:
            The query is transpiled with ``dialect='duckdb'``.
        Then:
            The dialect cannot peel the INTERSECTS out of the OR tree
            (only AND-conjunctions are decomposable safely) and routes
            to the naive-predicate plan.
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.start AS s FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval OR a.score > 0",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql

    def test_transpile_should_raise_when_extra_predicate_is_unqualified(
        self,
    ):
        """Test that an unqualified-column extra predicate raises a ValueError.

        Given:
            A column-to-column INTERSECTS INNER join carrying
            ``WHERE strand = '+'`` (the column reference has no table
            qualifier).
        When:
            The query is transpiled with ``dialect='duckdb'``.
        Then:
            The dialect cannot safely scope the predicate to either
            ``a`` or ``b`` inside the inner subquery, and the naive-predicate
            plan can't handle it correctly either (issue #94), so the
            dialect raises with an actionable message naming the join's
            valid aliases instead of silently routing to a broken plan.
        """
        # Act & assert
        with pytest.raises(ValueError, match="qualified with 'a' or 'b'"):
            transpile(
                "SELECT a.start AS s FROM peaks a "
                "JOIN genes b ON a.interval INTERSECTS b.interval "
                "WHERE strand = '+'",
                tables=["peaks", "genes"],
                dialect="duckdb",
            )

    def test_transpile_should_rewrite_residual_aliases_to_inner_a_b_when_user_uses_other_aliases(
        self,
    ):
        """Test that residuals are rewritten when the user picks aliases other than a/b.

        Given:
            A column-to-column INTERSECTS INNER join whose user aliases
            are ``p`` and ``g`` (instead of the dialect's hardcoded inner
            ``a`` / ``b``) and an extra ON predicate ``p.score > g.score``.
        When:
            The query is transpiled with ``dialect='duckdb'``.
        Then:
            The emitted inner SQL should reference ``a`` and ``b`` rather
            than ``p`` and ``g`` in the extra predicate (otherwise the
            inner subquery would fail to resolve the column refs).
        """
        # Arrange & act
        sql = transpile(
            "SELECT p.start AS s FROM peaks p "
            "JOIN genes g ON p.interval INTERSECTS g.interval AND p.score > g.score",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert "a.score > b.score" in sql
        assert "p.score" not in sql
        assert "g.score" not in sql

    def test_transpile_should_handle_custom_chrom_col_when_dialect_is_duckdb(self, conn):
        """Test that custom column names round-trip through the IEJoin path.

        Given:
            Two DuckDB tables defined with ``chromosome`` / ``start_pos`` /
            ``end_pos`` instead of the default column names, and a
            qualified-projection query referring to those columns.
        When:
            The query is transpiled with the matching ``Table`` configs and
            ``dialect='duckdb'`` and executed against the connection.
        Then:
            It should return the overlap rows correctly using the custom
            column names on both sides of the join.
        """
        # Arrange
        conn.execute(
            "CREATE TABLE peaks (chromosome VARCHAR, start_pos INTEGER, end_pos INTEGER)"
        )
        conn.execute("INSERT INTO peaks VALUES ('chr1', 100, 200), ('chr1', 1000, 1100)")
        conn.execute(
            "CREATE TABLE genes (chromosome VARCHAR, start_pos INTEGER, end_pos INTEGER)"
        )
        conn.execute("INSERT INTO genes VALUES ('chr1', 150, 250), ('chr1', 5000, 6000)")
        sql = transpile(
            """
            SELECT a.chromosome AS a_chrom, a.start_pos AS a_start,
                   b.start_pos AS b_start
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table(
                    "peaks",
                    chrom_col="chromosome",
                    start_col="start_pos",
                    end_col="end_pos",
                ),
                Table(
                    "genes",
                    chrom_col="chromosome",
                    start_col="start_pos",
                    end_col="end_pos",
                ),
            ],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert rows == [("chr1", 100, 150)]

    def test_transpile_should_passthrough_literal_intersects_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that literal-range INTERSECTS bypasses the IEJoin scaffolding.

        Given:
            A WHERE literal-range INTERSECTS query (not column-to-column).
        When:
            The query is transpiled with ``dialect=None`` and with
            ``dialect='duckdb'`` and both are executed.
        Then:
            It should produce the same SQL and the same result rows under
            both dialects (the IEJoin path is a no-op for literal ranges).
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 0, "+"),
                ("chr1", 1500, 1700, "p2", 0, "+"),
                ("chr1", 5000, 6000, "p3", 0, "+"),
            ],
        )
        query = "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'"

        # Act
        sql_default = transpile(query, tables=["peaks"])
        sql_duckdb = transpile(query, tables=["peaks"], dialect="duckdb")
        rows_default = sorted(conn.execute(sql_default).fetchall())
        rows_duckdb = sorted(conn.execute(sql_duckdb).fetchall())

        # Assert
        assert "SET VARIABLE" not in sql_duckdb.upper()
        assert "getvariable" not in sql_duckdb
        assert rows_default == rows_duckdb
        assert rows_duckdb == [("chr1", 1500, 1700, "p2", 0, "+")]

    def test_transpile_should_raise_when_dialect_is_unknown(self):
        """Test that an unknown dialect string raises ``ValueError``.

        Given:
            A dialect name other than ``'duckdb'`` or ``None``.
        When:
            ``transpile`` is called.
        Then:
            It should raise ``ValueError`` with a message that names the
            unknown dialect.
        """
        # Arrange
        query = "SELECT * FROM peaks"

        # Act & assert — message must name both the failure mode and the
        # offending value so the user can recover without guessing.
        with pytest.raises(ValueError, match=r"[Uu]nknown dialect.*'postgres'"):
            transpile(query, tables=["peaks"], dialect="postgres")

    def test_transpile_should_decline_bare_star_to_naive_plan_when_dialect_is_duckdb(
        self,
    ):
        """Test that bare ``SELECT *`` declines to the naive-predicate plan.

        Given:
            A column-to-column INTERSECTS JOIN whose projection is a bare
            ``*``.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should decline the IEJoin (no ``SET VARIABLE`` scaffolding)
            so DuckDB expands the star against the live schema (#202),
            rather than raising or silently narrowing the projection.
        """
        # Arrange
        query = """
            SELECT *
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE" not in sql.upper()
        assert "getvariable" not in sql

    # --- DI-001..DI-009: new SQL-structure / contract tests --------------

    def test_transpile_should_be_round_trip_executable_when_dialect_is_duckdb(
        self, peaks_genes
    ):
        """Test end-to-end execution of an INNER JOIN INTERSECTS under DuckDB.

        Given:
            A column-to-column INTERSECTS INNER JOIN with default column
            names against a multi-chromosome fixture with a known truth
            table of overlapping pairs.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed
            against an in-memory DuckDB connection.
        Then:
            It should return exactly the expected set of overlapping
            (peak, gene) tuples.
        """
        # Arrange
        sql = transpile(
            """
            SELECT a.chrom AS a_chrom, a.start AS a_start, a.end AS a_end,
                   b.chrom AS b_chrom, b.start AS b_start, b.end AS b_end
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(peaks_genes.execute(sql).fetchall())

        # Assert
        assert rows == sorted(_EXPECTED_OVERLAPS_PEAKS_GENES)

    def test_transpile_should_route_to_naive_predicate_plan_when_query_has_two_intersects_predicates(
        self, conn
    ):
        """Test that two INTERSECTS predicates fall back to the naive-predicate plan.

        Given:
            A three-table chained INTERSECTS query (peaks JOIN genes
            JOIN genes2) with ON-INTERSECTS predicates on each join.
        When:
            The query is transpiled both with ``dialect='duckdb'`` and
            with ``dialect=None`` and both are executed.
        Then:
            The DuckDB-dialect path should fall back to the naive-predicate plan
            and return the same multiset of triply-overlapping rows as the
            default ``dialect=None`` path.
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p", 0, "+")])
        _make_table(conn, "genes", [("chr1", 150, 250, "g", 0, "+")])
        conn.execute(
            "CREATE TABLE genes2 ("
            'chrom VARCHAR, "start" INTEGER, "end" INTEGER, '
            "name VARCHAR, score INTEGER, strand VARCHAR)"
        )
        conn.execute("INSERT INTO genes2 VALUES ('chr1', 175, 350, 'g2', 0, '+')")
        query = """
            SELECT a.start AS a_s, b.start AS b_s, c.start AS c_s
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            JOIN genes2 c ON b.interval INTERSECTS c.interval
            """
        sql_dd = transpile(
            query,
            tables=["peaks", "genes", "genes2"],
            dialect="duckdb",
        )
        sql_default = transpile(query, tables=["peaks", "genes", "genes2"])

        # Act
        rows_dd = sorted(conn.execute(sql_dd).fetchall())
        rows_default = sorted(conn.execute(sql_default).fetchall())

        # Assert
        assert rows_dd == rows_default
        assert rows_dd == [(100, 150, 175)]

    def test_query_should_route_to_naive_predicate_plan_and_execute_correctly_when_sibling_contains_present(
        self, conn
    ):
        """Test that a sibling CONTAINS defers an ANTI IEJoin to the naive plan.

        Given:
            An ANTI JOIN INTERSECTS query whose WHERE carries a sibling CONTAINS
            predicate, and inputs where inlining that residual into the
            per-chromosome ANTI ON would invert the anti-join semantics.
        When:
            The query is transpiled with ``dialect='duckdb'`` and with
            ``dialect=None`` and both are executed.
        Then:
            The DuckDB path should defer to the naive-predicate plan (no IEJoin
            SET VARIABLE) and return the same rows as the default plan.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 50, 250, "p1", 0, "+"),
                ("chr1", 120, 140, "p2", 0, "+"),
                ("chr1", 350, 500, "p3", 0, "+"),
                ("chr1", 10, 500, "p4", 0, "+"),
            ],
        )
        _make_table(conn, "genes", [("chr1", 300, 400, "g1", 0, "+")])
        query = (
            "SELECT a.start FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval "
            "WHERE a.interval CONTAINS 'chr1:100-200'"
        )
        sql_dd = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_default = transpile(query, tables=["peaks", "genes"])

        # Act
        rows_dd = sorted(conn.execute(sql_dd).fetchall())
        rows_default = sorted(conn.execute(sql_default).fetchall())

        # Assert
        assert "SET VARIABLE" not in sql_dd.upper()
        assert rows_dd == rows_default
        assert rows_dd == [(50,)]

    def test_transpile_should_route_to_naive_predicate_plan_when_residual_intersects_is_literal(
        self,
    ):
        """Test that a sibling literal INTERSECTS defers the join to the naive plan.

        Given:
            A column-to-column INTERSECTS join whose WHERE carries a sibling
            literal-range INTERSECTS predicate.
        When:
            The query is transpiled with ``dialect='duckdb'``.
        Then:
            It should defer to the naive-predicate plan rather than engaging the
            IEJoin (which cannot see the residual once pass 3 expands it).
        """
        # Arrange
        query = (
            "SELECT a.start FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "WHERE a.interval INTERSECTS 'chr1:100-200'"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE" not in sql.upper()

    def test_transpile_should_route_self_join_to_fallback_plan_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that a self-join INTERSECTS falls back and still executes.

        Given:
            A query joining ``peaks`` against itself on INTERSECTS.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            It should fall back to the naive-predicate plan and return the
            self-overlap pairs without raising.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 0, "+"),
                ("chr1", 150, 300, "p2", 0, "+"),
                ("chr1", 1000, 1100, "p3", 0, "+"),
            ],
        )
        sql = transpile(
            """
            SELECT a.start AS a_s, b.start AS b_s
            FROM peaks a JOIN peaks b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        # Each peak overlaps itself; p1[100,200) and p2[150,300) overlap each
        # other; p3 only overlaps itself.
        assert rows == [
            (100, 100),
            (100, 150),
            (150, 100),
            (150, 150),
            (1000, 1000),
        ]

    def test_transpile_should_emit_order_by_and_limit_on_outer_select_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that ORDER BY and LIMIT survive on the dialect's outer SELECT.

        Given:
            A column-to-column INTERSECTS INNER join carrying a top-level
            ``ORDER BY a.start LIMIT 1``.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            It should engage the dialect path (``SET VARIABLE
            __giql_iejoin_`` is emitted), append the user's ``ORDER BY``
            / ``LIMIT`` to the outer SELECT, and return exactly the
            first row by ascending ``a.start``.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),
                ("chr1", 300, 400, "p2", 2, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 1, "+"),
                ("chr1", 350, 450, "g2", 2, "+"),
            ],
        )
        sql = transpile(
            "SELECT a.start AS s FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "ORDER BY a.start LIMIT 1",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        # Structural check: ORDER BY must sit on the outer wrapper (after
        # ``query(getvariable(...)``) and before ``LIMIT`` — pins the clause
        # ordering contract of ``_build_sql``.
        wrapper_idx = sql.index("query(getvariable(")
        order_idx = sql.index("ORDER BY", wrapper_idx)
        limit_idx = sql.index("LIMIT 1", wrapper_idx)
        assert wrapper_idx < order_idx < limit_idx
        assert rows == [(100,)]

    def test_transpile_should_emit_distinct_on_outer_select_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that SELECT DISTINCT survives on the dialect's outer SELECT.

        Given:
            A column-to-column INTERSECTS INNER join projecting
            ``DISTINCT a.chrom`` over fixture data with two overlapping
            peaks sharing the same chromosome.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            It should engage the dialect path (``SET VARIABLE
            __giql_iejoin_`` is emitted), prepend ``DISTINCT`` to the
            outer SELECT, and return exactly one distinct chromosome row.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),
                ("chr1", 300, 400, "p2", 2, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 1, "+"),
                ("chr1", 350, 450, "g2", 2, "+"),
            ],
        )
        sql = transpile(
            "SELECT DISTINCT a.chrom AS c FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert "SELECT DISTINCT" in sql
        assert rows == [("chr1",)]

    def test_query_should_return_rows_in_order_by_order_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that ORDER BY DESC over multiple chromosomes orders rows.

        Given:
            A column-to-column INTERSECTS INNER join with two overlap
            rows on different chromosomes and a top-level
            ``ORDER BY a.start DESC``.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            The returned rows should be sorted by ``a.start`` descending.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),
                ("chr2", 500, 600, "p2", 2, "-"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 1, "+"),
                ("chr2", 550, 650, "g2", 2, "-"),
            ],
        )
        sql = transpile(
            "SELECT a.start AS s, a.chrom AS c FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "ORDER BY a.start DESC",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert rows == [(500, "chr2"), (100, "chr1")]

    def test_query_should_respect_offset_when_dialect_is_duckdb(self, conn):
        """Test that LIMIT N OFFSET M skips the first M rows.

        Given:
            Two overlap rows after sorting by ``a.start``.
        When:
            The query is transpiled with ``dialect='duckdb'``, a
            ``LIMIT 1 OFFSET 1`` is included, and the SQL is executed.
        Then:
            The second-by-order row is returned.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),
                ("chr1", 300, 400, "p2", 2, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 1, "+"),
                ("chr1", 350, 450, "g2", 2, "+"),
            ],
        )
        sql = transpile(
            "SELECT a.start AS s FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "ORDER BY a.start LIMIT 1 OFFSET 1",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert rows == [(300,)]

    def test_query_should_resolve_order_by_user_alias_when_dialect_is_duckdb(self, conn):
        """Test that ORDER BY can reference a SELECT-list user alias directly.

        Given:
            A column-to-column INTERSECTS INNER join with
            ``SELECT a.start AS s ... ORDER BY s``.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            The order-by reference to the bare user alias should resolve
            against the outer wrapper and return rows in ascending order
            by ``a.start``.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 300, 400, "p2", 2, "+"),
                ("chr1", 100, 200, "p1", 1, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 1, "+"),
                ("chr1", 350, 450, "g2", 2, "+"),
            ],
        )
        sql = transpile(
            "SELECT a.start AS s FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "ORDER BY s",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert rows == [(100,), (300,)]

    def test_query_should_auto_project_order_by_column_when_not_in_select_under_dialect(
        self, conn
    ):
        """Test that an ORDER BY column not in the SELECT is auto-projected internally.

        Given:
            Two overlap rows sortable by ``b.score`` (which is not in
            the user's SELECT list).
        When:
            The query is transpiled with ``dialect='duckdb'`` and
            ``ORDER BY b.score``, and the SQL is executed.
        Then:
            The dialect should silently project ``b.score`` into the
            inner subquery so the outer ORDER BY can resolve it; the
            user's SELECT list remains unchanged, and rows return in
            ``b.score`` order.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),
                ("chr1", 300, 400, "p2", 2, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g_high", 99, "+"),
                ("chr1", 350, 450, "g_low", 1, "+"),
            ],
        )
        sql = transpile(
            "SELECT a.start AS s FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "ORDER BY b.score",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert: ordered by g_low (b.score=1) first then g_high (b.score=99)
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert [r[0] for r in rows] == [300, 100]
        # The user's projection is still single-column ``s``.
        assert all(len(r) == 1 for r in rows)

    def test_query_should_allow_order_by_when_outer_name_is_shared_under_dialect(
        self, conn
    ):
        """Test that ORDER BY on a shared outer name resolves unambiguously.

        Given:
            A column-to-column INTERSECTS INNER join projecting both
            ``a.start`` and ``b.start`` (both would otherwise surface
            under the outer name ``"start"``) and ``ORDER BY a.start``.
        When:
            The query is transpiled with ``dialect='duckdb'`` and
            executed.
        Then:
            The rewriter resolves the reference to the inner alias for
            ``a.start`` rather than to a colliding outer name, so the
            query succeeds and orders by ``a.start``.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 300, 400, "p2", 2, "+"),
                ("chr1", 100, 200, "p1", 1, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 1, "+"),
                ("chr1", 350, 450, "g2", 2, "+"),
            ],
        )
        sql = transpile(
            "SELECT a.start, b.start FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "ORDER BY a.start",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert rows == [(100, 150), (300, 350)]

    def test_query_should_group_overlap_pairs_by_chrom_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that an outer GROUP BY with COUNT(*) under the dialect groups rows.

        Given:
            A column-to-column INTERSECTS INNER join with a top-level
            ``GROUP BY a.chrom`` and a ``COUNT(*)`` aggregate over two
            chromosomes each having a single overlap.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            The dialect path should engage and return one row per
            chromosome with the matching count.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),
                ("chr2", 100, 200, "p2", 2, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 1, "+"),
                ("chr2", 150, 250, "g2", 2, "+"),
            ],
        )
        sql = transpile(
            "SELECT a.chrom AS c, COUNT(*) AS n FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "GROUP BY a.chrom",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert "GROUP BY" in sql
        assert rows == [("chr1", 1), ("chr2", 1)]

    def test_query_should_filter_groups_via_having_when_dialect_is_duckdb(self, conn):
        """Test that an outer HAVING filter survives under the dialect.

        Given:
            A column-to-column INTERSECTS INNER join with a
            ``GROUP BY ... HAVING COUNT(*) > 1`` over data where exactly
            one chromosome has more than one overlap.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            The dialect path should engage and return only the
            chromosome whose overlap count exceeds the threshold.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),
                ("chr1", 110, 210, "p2", 2, "+"),
                ("chr2", 100, 200, "p3", 3, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 1, "+"),
                ("chr2", 150, 250, "g2", 2, "+"),
            ],
        )
        sql = transpile(
            "SELECT a.chrom AS c, COUNT(*) AS n FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "GROUP BY a.chrom HAVING COUNT(*) > 1",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert "HAVING" in sql
        assert rows == [("chr1", 2)]

    def test_query_should_compute_sum_aggregate_with_group_by_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that SUM(a.score) over GROUP BY a.chrom works under the dialect.

        Given:
            A column-to-column INTERSECTS INNER join over two chromosomes
            with multiple overlaps per side, plus a ``GROUP BY a.chrom``
            and ``SUM(a.score)`` aggregate.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            The dialect path should engage and return one row per
            chromosome with the per-chromosome sum.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 10, "+"),
                ("chr1", 300, 400, "p2", 5, "+"),
                ("chr2", 100, 200, "p3", 7, "-"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 0, "+"),
                ("chr1", 350, 450, "g2", 0, "+"),
                ("chr2", 150, 250, "g3", 0, "-"),
            ],
        )
        sql = transpile(
            "SELECT a.chrom AS c, SUM(a.score) AS s FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "GROUP BY a.chrom",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert rows == [("chr1", 15), ("chr2", 7)]

    def test_query_should_handle_count_distinct_aggregate_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that COUNT(DISTINCT a.x) preserves DISTINCT through the rewrite.

        Given:
            Overlap data where the same ``a.name`` appears twice on the
            same chromosome and a ``COUNT(DISTINCT a.name)`` aggregate.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            The dialect path should engage and the count should reflect
            the number of distinct names (not the total row count).
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p_dup", 1, "+"),
                ("chr1", 110, 210, "p_dup", 2, "+"),
                ("chr1", 300, 400, "p_other", 3, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 0, "+"),
                ("chr1", 350, 450, "g2", 0, "+"),
            ],
        )
        sql = transpile(
            "SELECT a.chrom AS c, COUNT(DISTINCT a.name) AS n FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "GROUP BY a.chrom",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        # The two p_dup peaks both overlap g1 — under DISTINCT they
        # collapse to 1; p_other overlaps g2; total distinct = 2.
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert "DISTINCT" in sql.upper()
        assert rows == [("chr1", 2)]

    def test_query_should_match_naive_predicate_plan_for_avg_aggregate(self, conn):
        """Test that AVG aggregate matches the naive-predicate plan executing-side.

        Given:
            Two chromosomes with multiple overlap pairs and varied
            ``a.score`` values that produce a non-integer per-chrom
            average.
        When:
            A ``SELECT a.chrom, AVG(a.score) ... GROUP BY a.chrom`` is
            transpiled under both ``dialect=None`` (naive-predicate plan) and
            ``dialect='duckdb'`` and executed.
        Then:
            Both plans return the same per-chrom averages (compared via
            ``pytest.approx`` since AVG yields floats). This is the
            execution-side AVG-equivalence check that the aggregate PBT
            intentionally skips (PBT excludes AVG to keep its equality
            assertion strict).
        """
        # Arrange — scores chosen so chr1's mean is 31.0 (integer) and
        # chr2's mean is 12.5 (genuinely non-integer), forcing
        # ``pytest.approx`` to actually do float comparison.
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 10, "+"),
                ("chr1", 110, 210, "p2", 32, "+"),
                ("chr1", 120, 220, "p3", 51, "+"),
                ("chr2", 100, 200, "p4", 7, "-"),
                ("chr2", 110, 210, "p5", 18, "-"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 0, "+"),
                ("chr2", 150, 250, "g2", 0, "-"),
            ],
        )
        query = (
            "SELECT a.chrom AS c, AVG(a.score) AS v FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "GROUP BY a.chrom ORDER BY a.chrom"
        )
        sql_default = transpile(query, tables=["peaks", "genes"])
        sql_duckdb = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Act
        rows_default = conn.execute(sql_default).fetchall()
        rows_duckdb = conn.execute(sql_duckdb).fetchall()

        # Assert — value-equivalence with float tolerance.
        assert len(rows_default) == len(rows_duckdb)
        for (c_d, v_d), (c_q, v_q) in zip(rows_default, rows_duckdb):
            assert c_d == c_q
            assert v_d == pytest.approx(v_q)

    def test_query_should_combine_count_star_and_sum_when_dialect_is_duckdb(self, conn):
        """Test that multiple aggregates coexist correctly under the dialect.

        Given:
            Two chromosomes with multiple overlaps and a SELECT list that
            mixes ``COUNT(*)`` and ``SUM(b.score)``.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            The dialect should engage, project the underlying column for
            ``SUM(b.score)`` into the inner subquery, and emit both
            aggregates against the wrapper relation.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),
                ("chr1", 110, 210, "p2", 2, "+"),
                ("chr2", 100, 200, "p3", 3, "-"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 10, "+"),
                ("chr2", 150, 250, "g2", 20, "-"),
            ],
        )
        sql = transpile(
            "SELECT a.chrom AS c, COUNT(*) AS n, SUM(b.score) AS bs "
            "FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval "
            "GROUP BY a.chrom",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        # chr1: two overlaps, each with g1 (score 10) -> sum 20; chr2: 1 with g2 (20)
        assert rows == [("chr1", 2, 20), ("chr2", 1, 20)]

    def test_transpile_should_raise_when_aggregate_argument_is_unqualified(self):
        """Test that an aggregate with an unqualified-column argument raises.

        Given:
            A SELECT-list aggregate ``SUM(score)`` whose argument is not
            table-qualified.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            A ``ValueError`` should be raised that names the offending
            projection and asks for a qualified reference.
        """
        # Act & assert
        with pytest.raises(ValueError, match=r"qualified"):
            transpile(
                "SELECT a.chrom AS c, SUM(score) AS s FROM peaks a "
                "JOIN genes b ON a.interval INTERSECTS b.interval "
                "GROUP BY a.chrom",
                tables=["peaks", "genes"],
                dialect="duckdb",
            )

    def test_transpile_should_route_to_naive_predicate_plan_when_outer_join_intersects_lives_in_where(
        self,
    ):
        """Test that LEFT JOIN ON TRUE WHERE INTERSECTS falls back.

        Given:
            A ``LEFT JOIN ... ON TRUE WHERE a.interval INTERSECTS
            b.interval`` query — the INTERSECTS lives in WHERE while
            the join keeps its LEFT side modifier.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The dialect must NOT engage (rewriting as an inner IEJoin
            would silently drop LEFT JOIN's unmatched-row guarantee).
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.start AS s FROM peaks a "
            "LEFT JOIN genes b ON TRUE "
            "WHERE a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql

    def test_transpile_should_route_to_naive_predicate_plan_when_full_outer_join_intersects_lives_in_where(
        self,
    ):
        """Test that FULL OUTER JOIN ... WHERE INTERSECTS falls back.

        Given:
            A ``FULL OUTER JOIN ... ON TRUE WHERE a.interval INTERSECTS
            b.interval`` — same shape as LEFT, but FULL.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The dialect must not engage; the naive-predicate plan handles the
            outer-join semantics.
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.start AS s FROM peaks a "
            "FULL OUTER JOIN genes b ON TRUE "
            "WHERE a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql

    def test_transpile_should_route_to_naive_predicate_plan_when_extra_predicate_contains_window_function(
        self,
    ):
        """Test that window functions inside extras force a fallback.

        Given:
            An INTERSECTS join carrying ``WHERE a.score > (ROW_NUMBER()
            OVER (PARTITION BY a.chrom))`` — a window function in an
            extra predicate. DuckDB forbids window functions in JOIN ON
            clauses, so the dialect cannot inline this residual.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The dialect must reject the residual and route to the
            naive-predicate plan (rather than emit SQL that DuckDB rejects at
            runtime).
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.start AS s FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "WHERE a.score > (ROW_NUMBER() OVER (PARTITION BY a.chrom))",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql

    def test_transpile_should_decline_window_aggregate_to_naive_plan(self):
        """Test that a window-aggregate projection declines to the naive plan.

        Given:
            A SELECT list containing ``SUM(a.score) OVER (PARTITION BY
            a.chrom)`` — a window aggregate over qualified columns the IEJoin
            projection rebuild cannot express.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should decline the IEJoin (no ``SET VARIABLE`` scaffolding) so
            the naive-predicate plan evaluates the window aggregate (#205).
        """
        # Act
        sql = transpile(
            "SELECT a.chrom AS c, "
            "SUM(a.score) OVER (PARTITION BY a.chrom) AS s "
            "FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE" not in sql.upper()
        assert "getvariable" not in sql

    def test_transpile_should_decline_paren_wrapped_window_aggregate_to_naive_plan(
        self,
    ):
        """Test that ``(SUM(a.score) OVER (...))`` declines to the naive plan.

        Given:
            A paren-wrapped window aggregate over qualified columns in the
            SELECT list.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The dispatcher peels the Paren and declines the IEJoin (no ``SET
            VARIABLE`` scaffolding) so the naive-predicate plan handles the
            window aggregate (#205).
        """
        # Act
        sql = transpile(
            "SELECT a.chrom AS c, "
            "(SUM(a.score) OVER (PARTITION BY a.chrom)) AS s "
            "FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE" not in sql.upper()
        assert "getvariable" not in sql

    def test_transpile_should_decline_aggregate_in_expression_to_naive_plan(self):
        """Test that arithmetic over an aggregate declines to the naive plan.

        Given:
            A SELECT list containing ``COUNT(*) * 2`` — an aggregate nested in
            an arithmetic expression the IEJoin projection rebuild cannot
            express.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should decline the IEJoin (no ``SET VARIABLE`` scaffolding) so
            the naive-predicate plan evaluates the expression (#205).
        """
        # Act
        sql = transpile(
            "SELECT a.chrom AS c, COUNT(*) * 2 AS doubled "
            "FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval "
            "GROUP BY a.chrom",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE" not in sql.upper()
        assert "getvariable" not in sql

    def test_transpile_should_route_to_naive_predicate_plan_when_modifier_has_subquery(
        self,
    ):
        """Test that a subquery inside ORDER BY routes to the naive-predicate plan.

        Given:
            A query whose ``ORDER BY`` clause embeds a scalar subquery
            (``(SELECT MAX(score) FROM peaks)``). The modifier rewriter
            is not scope-aware, so the dialect should fall back to the
            naive-predicate plan rather than rewrite column refs inside the
            nested scope.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The dialect must not engage; the naive-predicate plan handles the
            subquery natively.
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.start AS s FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "ORDER BY (SELECT MAX(score) FROM peaks)",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql

    def test_transpile_should_route_to_naive_predicate_plan_when_modifier_has_exists(
        self,
    ):
        """Test that EXISTS inside ORDER BY routes to the naive-predicate plan.

        Given:
            A query whose ``ORDER BY`` clause embeds an ``EXISTS
            (SELECT ...)`` clause. sqlglot parses EXISTS as
            ``Exists(Select(...))`` without an ``exp.Subquery``
            wrapper, so the gate must also check for ``exp.Select`` —
            otherwise EXISTS slips past and the non-scope-aware
            modifier rewriter corrupts the EXISTS scope.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The dialect must not engage.
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.start AS s FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "ORDER BY EXISTS (SELECT 1 FROM peaks p WHERE p.score > 0)",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql

    def test_transpile_should_route_to_naive_predicate_plan_when_having_has_subquery(
        self,
    ):
        """Test that a subquery inside HAVING routes to the naive-predicate plan.

        Given:
            A query whose HAVING compares against a scalar subquery
            (``HAVING SUM(a.score) > (SELECT AVG(score) FROM peaks)``).
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The dialect must not engage.
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.chrom AS c, SUM(a.score) AS s FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "GROUP BY a.chrom "
            "HAVING SUM(a.score) > (SELECT AVG(score) FROM peaks)",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql

    def test_transpile_should_accept_paren_wrapped_aggregate_in_select(
        self,
    ):
        """Test that ``(SUM(a.score)) AS s`` is accepted, not rejected.

        Given:
            A SELECT-list projection that wraps an aggregate in an
            otherwise-redundant ``exp.Paren`` (``(SUM(a.score)) AS s``).
            Previously this hit the arithmetic-over-aggregate raise
            with misleading guidance, because ``exp.Paren`` is not an
            ``exp.AggFunc``.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The dialect should peel the paren and route the projection
            through the aggregate branch, emitting valid SQL with the
            inner-alias rewrite.
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.chrom AS c, (SUM(a.score)) AS s "
            "FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval "
            "GROUP BY a.chrom",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        # The aggregate's argument should be rewritten to an inner alias
        # rather than left as ``a.score`` on the outer SELECT.
        outer_select = sql.split(";\n", 1)[-1]
        assert "a.score" not in outer_select
        assert re.search(r"SUM\(__giql_p\d+\)", outer_select) is not None

    def test_transpile_should_decline_scalar_subquery_to_naive_plan(self):
        """Test that a scalar subquery in the SELECT list declines to the naive plan.

        Given:
            A SELECT list containing a scalar subquery
            (``(SELECT SUM(score) FROM peaks)``) — the subquery's columns
            resolve against its own scope, and the IEJoin projection rebuild
            cannot express it.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should decline the IEJoin (no ``SET VARIABLE`` scaffolding) so
            the naive-predicate plan evaluates the subquery (#205).
        """
        # Act
        sql = transpile(
            "SELECT a.chrom AS c, "
            "(SELECT SUM(score) FROM peaks) AS s "
            "FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE" not in sql.upper()
        assert "getvariable" not in sql

    def test_transpile_should_route_to_naive_predicate_plan_when_query_has_distinct_on(
        self,
    ):
        """Test that DISTINCT ON (...) falls back to the naive-predicate plan.

        Given:
            A column-to-column INTERSECTS INNER join with a
            ``SELECT DISTINCT ON (a.chrom) ...`` projection.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The dialect should not engage; the SQL should not contain
            the ``SET VARIABLE __giql_iejoin_`` marker.
        """
        # Arrange & act
        sql = transpile(
            "SELECT DISTINCT ON (a.chrom) a.chrom, a.start FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql

    def test_transpile_should_route_to_naive_predicate_plan_and_execute_correctly_when_query_has_top_level_with_clause(
        self, conn
    ):
        """Test that a top-level WITH clause falls back to the naive-predicate plan.

        Given:
            A query whose INTERSECTS join is wrapped in a top-level
            ``WITH big AS (SELECT * FROM genes)`` whose CTE is then used
            as one of the joined operands.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            It should fall back to the naive-predicate plan, preserve the CTE
            definition, and execute against the materialized CTE. The
            dialect path would otherwise emit a `FROM big` referencing
            a missing table.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [("chr1", 100, 200, "p1", 1, "+")],
        )
        _make_table(
            conn,
            "genes",
            [("chr1", 150, 250, "g1", 1, "+")],
        )
        sql = transpile(
            "WITH big AS (SELECT * FROM genes) "
            "SELECT a.start AS s FROM peaks a "
            "JOIN big b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql
        assert rows == [(100,)]

    def test_transpile_should_route_to_naive_predicate_plan_when_three_table_cross_join(
        self, conn
    ):
        """Test that a 3-table cross-join falls back to the naive-predicate plan.

        Given:
            A 3-table comma-style FROM clause where the INTERSECTS
            predicate only connects two of the three tables. The third
            table contributes a real cross-product factor to the result.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            It should fall back to the naive-predicate plan and return a row
            count consistent with the full 3-way join — the dialect
            path would silently drop the third table.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [("chr1", 100, 200, "p1", 1, "+")],
        )
        _make_table(
            conn,
            "genes",
            [("chr1", 150, 250, "g1", 1, "+")],
        )
        _make_table(
            conn,
            "extra",
            [
                ("chr1", 0, 10, "e1", 0, "+"),
                ("chr1", 0, 10, "e2", 0, "+"),
            ],
        )
        sql = transpile(
            "SELECT a.start AS s FROM peaks a, genes b, extra c "
            "WHERE a.interval INTERSECTS b.interval",
            tables=["peaks", "genes", "extra"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql
        # One (peaks, genes) overlap pair cross-joined with 2 extra rows
        # = 2 result rows. The dialect path would emit only 1 (extra dropped).
        assert len(rows) == 2

    def test_transpile_should_route_to_naive_predicate_plan_when_target_is_giql_operator(
        self, conn
    ):
        """Test that a GIQL table-function in JOIN position falls back.

        Given:
            A query whose INTERSECTS join uses ``DISJOIN(genes)`` as the
            right-hand operand. ``DISJOIN(...)`` parses as an
            ``exp.Table`` with an empty ``name``; the pre-fix dialect
            path would interpolate that empty name into broken SQL.
        When:
            The query is transpiled with ``dialect='duckdb'``.
        Then:
            The dialect path should return ``None`` (fall back) and the
            emitted SQL should contain no ``SET VARIABLE`` block and no
            empty ``FROM `` interpolation. (Execution is the naive-predicate
            path's domain and is covered by DISJOIN's own tests.)
        """
        # Arrange & Act
        sql = transpile(
            "SELECT a.start AS s FROM peaks a "
            "JOIN DISJOIN(genes) b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql
        # No empty-name table interpolation leaked into the output.
        assert "FROM  " not in sql
        assert "FROM (SELECT * FROM  " not in sql

    def test_transpile_should_double_quote_table_name_and_execute_correctly_when_table_name_is_a_reserved_word(
        self, conn
    ):
        """Test that the dialect double-quotes a reserved-word table name.

        Given:
            A column-to-column INTERSECTS INNER join whose left table is
            named ``select`` — a SQL reserved word that must be quoted
            wherever the dialect interpolates it.
        When:
            The query is transpiled with ``dialect='duckdb'`` and the
            generated SQL is executed against a real DuckDB table named
            ``select``.
        Then:
            The dialect should engage (this *is* the supported shape),
            the emitted SQL should reference ``"select"`` (quoted) and
            never bare ``FROM select``, and the query should execute to
            the expected overlap row.
        """
        # Arrange
        conn.execute(
            'CREATE TABLE "select" '
            '(chrom VARCHAR, "start" BIGINT, "end" BIGINT, name VARCHAR, '
            "score INTEGER, strand VARCHAR)"
        )
        conn.execute(
            'INSERT INTO "select" VALUES (?, ?, ?, ?, ?, ?)',
            ("chr1", 100, 200, "p", 1, "+"),
        )
        conn.execute(
            "CREATE TABLE genes "
            '(chrom VARCHAR, "start" BIGINT, "end" BIGINT, name VARCHAR, '
            "score INTEGER, strand VARCHAR)"
        )
        conn.execute(
            "INSERT INTO genes VALUES (?, ?, ?, ?, ?, ?)",
            ("chr1", 150, 250, "g", 1, "+"),
        )
        sql = transpile(
            'SELECT a.start AS s FROM "select" a '
            "JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=[Table("select"), "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert '"select"' in sql
        # Bare unquoted identifier in a FROM position would have broken DuckDB.
        assert "FROM select " not in sql
        assert "FROM select," not in sql
        assert rows == [(100,)]

    def test_transpile_should_emit_iejoin_block_and_return_overlap_when_implicit_cross_join_has_swapped_alias_order(
        self, conn
    ):
        """Test that an implicit cross-join with swapped FROM order works.

        Given:
            An implicit cross-join where the FROM clause lists ``genes``
            first and ``peaks`` second, while the INTERSECTS predicate
            references ``a.interval`` (peaks) before ``b.interval``
            (genes).
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            It should emit the dialect's ``SET VARIABLE`` IEJoin block
            and return the correct overlap pair (the IEJoin path must
            not depend on FROM order matching INTERSECTS argument order).
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p", 0, "+")])
        _make_table(conn, "genes", [("chr1", 150, 250, "g", 0, "+")])
        sql = transpile(
            """
            SELECT a.start AS a_s, b.start AS b_s
            FROM genes b, peaks a
            WHERE a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert "SET VARIABLE" in sql.upper()
        assert rows == [(100, 150)]

    def test_transpile_should_raise_when_select_list_has_only_unqualified_columns(
        self,
    ):
        """Test that unqualified column projections are rejected under DuckDB.

        Given:
            A query whose SELECT list contains bare column names (no
            table qualifier).
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should raise ``ValueError`` whose message mentions
            "qualified".
        """
        # Arrange
        query = """
            SELECT chrom, start
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """

        # Act & assert — message should name the offending column and the
        # ``a.col``/``b.col`` form so the user can fix the query directly.
        with pytest.raises(ValueError) as excinfo:
            transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        message = str(excinfo.value)
        assert "qualified" in message
        assert "chrom" in message
        assert "a.col" in message or "b.col" in message

    def test_transpile_should_raise_when_projection_references_unknown_table_alias(
        self,
    ):
        """Test that a projection referring to an out-of-scope alias is rejected.

        Given:
            A SELECT list that qualifies a column with a table alias
            (``c``) not present in the FROM/JOIN clauses.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should raise ``ValueError`` whose message names the
            unknown qualifier.
        """
        # Arrange
        query = """
            SELECT c.col
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """

        # Act & assert
        with pytest.raises(ValueError, match="[Uu]nknown.*qualifier|'c'"):
            transpile(query, tables=["peaks", "genes"], dialect="duckdb")

    def test_transpile_should_decline_expression_form_projection_to_naive_plan(self):
        """Test that an arithmetic projection expression declines to the naive plan.

        Given:
            A SELECT list containing an expression over a qualified column
            (``a.start + 1``) rather than a bare qualified column.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should decline the IEJoin (no ``SET VARIABLE`` scaffolding) so
            the naive-predicate plan evaluates the expression (#205).
        """
        # Arrange
        query = """
            SELECT a.start + 1
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE" not in sql.upper()
        assert "getvariable" not in sql

    def test_transpile_should_decline_star_with_unknown_qualifier_to_naive_plan(
        self,
    ):
        """Test that ``c.*`` against an out-of-scope alias declines to naive.

        Given:
            A SELECT list with a qualified-star projection ``c.*`` whose
            qualifier ``c`` is not present in the FROM/JOIN clauses.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should decline the IEJoin (no ``SET VARIABLE`` scaffolding)
            like every star (#202) and defer the out-of-scope qualifier to
            DuckDB's binder, which rejects it at execution — matching the
            naive-predicate plan rather than raising at transpile time.
        """
        # Arrange
        query = """
            SELECT c.*
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE" not in sql.upper()
        assert "getvariable" not in sql

    def test_transpile_should_use_unique_variable_names_across_interleaved_calls(
        self, conn
    ):
        """Test that two interleaved transpile outputs do not collide.

        Given:
            Two ``transpile(..., dialect='duckdb')`` outputs whose
            ``SET VARIABLE`` and ``SELECT`` statements are interleaved
            (both ``SET``, then both ``SELECT``) in a single DuckDB
            session, with the two queries projecting different columns
            so that a fixed-name collision would surface one query's
            result in the other's slot.
        When:
            The interleaved statements are executed.
        Then:
            Each ``SELECT`` should return the result of its own query —
            proving the dialect emits a per-call unique variable name.
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p", 0, "+")])
        _make_table(conn, "genes", [("chr1", 150, 250, "g", 0, "+")])
        sql_a = transpile(
            "SELECT a.chrom AS v FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )
        sql_b = transpile(
            "SELECT a.start AS v FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )
        # The dialect emits exactly one `;\n` between the SET and the
        # SELECT, so split each output on that to interleave the four
        # statements deliberately. Under a single fixed variable name
        # the second `SET` would overwrite the first before either
        # `SELECT` ran, surfacing query B's row in query A's slot.
        set_a, sep_a, select_a = sql_a.partition(";\n")
        set_b, sep_b, select_b = sql_b.partition(";\n")
        assert sep_a and sep_b
        assert set_a.startswith("SET VARIABLE")
        assert set_b.startswith("SET VARIABLE")
        assert select_a.startswith("SELECT")
        assert select_b.startswith("SELECT")

        # Act
        conn.execute(set_a)
        conn.execute(set_b)
        row_a = conn.execute(select_a).fetchall()
        row_b = conn.execute(select_b).fetchall()

        # Assert
        assert row_a == [("chr1",)]
        assert row_b == [(100,)]

    def test_transpile_should_preserve_order_by_modifiers_when_dialect_is_duckdb(
        self,
    ):
        """Test that DESC and NULLS FIRST modifiers survive the rewrite.

        Given:
            A column-to-column INTERSECTS join with
            ``ORDER BY a.start DESC NULLS FIRST``.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The emitted outer SELECT should preserve both ``DESC`` and
            ``NULLS FIRST`` in the rewritten ORDER BY clause. (sqlglot
            strips its own defaults like ``ASC NULLS FIRST``, so the
            test uses the non-default-for-DESC variant to assert that
            both modifiers round-trip through the rewriter.)
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.start AS s FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "ORDER BY a.start DESC NULLS FIRST",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        outer_select = sql.split(";\n", 1)[-1]
        order_idx = outer_select.index("ORDER BY")
        assert "DESC" in outer_select[order_idx:]
        assert "NULLS FIRST" in outer_select[order_idx:]

    def test_transpile_should_preserve_multiple_order_by_expressions_when_dialect_is_duckdb(
        self,
    ):
        """Test that multi-column ORDER BY survives, with per-term modifiers.

        Given:
            A column-to-column INTERSECTS join with
            ``ORDER BY a.start, b.end DESC``.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The outer ORDER BY clause should carry two comma-separated
            terms, ``DESC`` should attach to the second term only, and
            neither the user's ``a.``/``b.`` qualifiers should remain on
            the outer ORDER BY (they get rewritten to inner aliases).
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.start AS s, b.end AS e FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "ORDER BY a.start, b.end DESC",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        wrapper_idx = sql.index("query(getvariable(")
        order_clause = sql[sql.index("ORDER BY", wrapper_idx) :]
        # Pin the DESC attachment positionally so a regression that
        # accidentally moved DESC to the first term (or attached it to
        # both) would fail. Two inner-alias placeholders, comma-separated,
        # with DESC on the second.
        assert (
            re.search(
                r"ORDER BY\s+__giql_p\d+(?:\s+NULLS\s+(?:FIRST|LAST))?\s*,"
                r"\s*__giql_p\d+\s+DESC",
                order_clause,
            )
            is not None
        )
        # The rewriter binds user-qualified refs to inner aliases on the
        # outer SELECT, so ``a.start`` / ``b.end`` should not appear in
        # the outer ORDER BY.
        assert "a.start" not in order_clause
        assert "b.end" not in order_clause

    def test_transpile_should_inline_cross_side_on_predicate_when_dialect_is_duckdb(
        self,
    ):
        """Test that a cross-side ON predicate is inlined into each subquery.

        Given:
            A join ``ON a.interval INTERSECTS b.interval AND a.score > b.score``
            referencing both sides of the inner subquery.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The cross-side predicate ``a.score > b.score`` should appear
            in the dynamic SQL literal (inside the per-chromosome
            subquery template), proving both ``a`` and ``b`` are in
            scope inside the inlined extra predicate.
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.start FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "AND a.score > b.score",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert "a.score > b.score" in sql

    def test_transpile_should_inline_between_predicate_when_dialect_is_duckdb(
        self,
    ):
        """Test that a ``BETWEEN`` extra predicate is inlined.

        Given:
            A join carrying an extra ``BETWEEN`` predicate on a
            qualified column.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The ``BETWEEN 10 AND 50`` substring should appear in the
            dynamic SQL literal.
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.start FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "AND a.score BETWEEN 10 AND 50",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert — the predicate must land inside the per-chromosome
        # template (SET VARIABLE half), not on the outer wrapper SELECT.
        assert "SET VARIABLE __giql_iejoin_" in sql
        set_part, outer_part = sql.split(";\n", 1)
        assert "BETWEEN 10 AND 50" in set_part
        assert "BETWEEN" not in outer_part

    def test_transpile_should_inline_is_null_predicate_when_dialect_is_duckdb(
        self,
    ):
        """Test that an ``IS NULL`` extra predicate is inlined.

        Given:
            A join carrying ``AND a.strand IS NULL`` as an extra
            predicate.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The ``IS NULL`` substring should appear in the dynamic SQL
            literal.
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.start FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "AND a.strand IS NULL",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert — the predicate must land inside the per-chromosome
        # template (SET VARIABLE half), not on the outer wrapper SELECT.
        assert "SET VARIABLE __giql_iejoin_" in sql
        set_part, outer_part = sql.split(";\n", 1)
        assert "IS NULL" in set_part
        assert "IS NULL" not in outer_part

    def test_transpile_should_route_to_naive_predicate_plan_when_paren_wraps_intersects(
        self,
    ):
        """Test that a Paren-wrapped INTERSECTS routes to the naive-predicate plan.

        Given:
            A join ``ON (a.interval INTERSECTS b.interval) AND a.score > 0``
            where the INTERSECTS is wrapped in explicit parentheses.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The dialect only decomposes AND-conjunctions when stripping
            the INTERSECTS out of a join clause; a paren wrapper leaves
            the INTERSECTS embedded in the residual and the dialect
            must fall back to the naive-predicate plan.
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.start FROM peaks a "
            "JOIN genes b ON (a.interval INTERSECTS b.interval) "
            "AND a.score > 0",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql

    def test_transpile_should_route_to_naive_predicate_plan_when_residual_has_subquery(
        self,
    ):
        """Test that a subquery in the residual predicate forces a fallback.

        Given:
            A join INTERSECTS plus
            ``WHERE a.score > (SELECT MAX(score) FROM peaks)``.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The dialect's residual-safety check rejects predicates
            containing a subquery, so the query falls back to the
            naive-predicate plan.
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.start FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "WHERE a.score > (SELECT MAX(score) FROM peaks)",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql

    def test_transpile_should_rewrite_group_by_multiple_columns_to_inner_aliases_when_dialect_is_duckdb(
        self,
    ):
        """Test that multi-column GROUP BY rewrites every key to an inner alias.

        Given:
            A query ``SELECT a.chrom, a.score, COUNT(*) ... GROUP BY a.chrom, a.score``.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The outer GROUP BY clause should reference two distinct
            ``__giql_p<n>`` placeholders (one per key) and the user's
            ``a.chrom``/``a.score`` should not appear in the outer
            GROUP BY.
        """

        # Arrange & act
        sql = transpile(
            "SELECT a.chrom AS c, a.score AS s, COUNT(*) AS n FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "GROUP BY a.chrom, a.score",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        wrapper_idx = sql.index("query(getvariable(")
        outer = sql[wrapper_idx:]
        group_clause = outer[outer.index("GROUP BY") :]
        # Two distinct inner-alias placeholders should appear in GROUP BY.
        placeholders = set(
            re.findall(r"__giql_p\d+", group_clause.split("HAVING")[0].split("ORDER")[0])
        )
        assert len(placeholders) >= 2
        assert "a.chrom" not in group_clause.split("ORDER")[0].split("HAVING")[0]
        assert "a.score" not in group_clause.split("ORDER")[0].split("HAVING")[0]

    def test_transpile_should_emit_sum_min_max_avg_aggregates_with_inner_alias_when_dialect_is_duckdb(
        self,
    ):
        """Test that four aggregates over the same column reuse one inner alias.

        Given:
            A query ``SELECT SUM(a.score), MIN(a.score), MAX(a.score),
            AVG(a.score) ... GROUP BY a.chrom``.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The inner SELECT should project ``a."score"`` exactly once,
            and all four outer aggregates should reference the same
            ``__giql_p<n>`` placeholder.
        """

        # Arrange & act
        sql = transpile(
            "SELECT a.chrom AS c, "
            "SUM(a.score) AS s, MIN(a.score) AS mn, "
            "MAX(a.score) AS mx, AVG(a.score) AS av FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "GROUP BY a.chrom",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        # Behavior-level dedup check: ``a."score"`` should be allocated to
        # exactly one inner alias across the whole SQL string, regardless
        # of which ``__giql_p<n>`` index that alias happens to get
        # (allocation order depends on sqlglot's traversal of the SELECT
        # list and is not part of the dialect's public contract).
        score_inner_allocs = re.findall(r'a\."score"\s+AS\s+(__giql_p\d+)', sql)
        assert len(set(score_inner_allocs)) == 1
        outer_select = sql.split(";\n", 1)[-1]
        # Each of the four aggregates appears once with an inner-alias arg,
        # and all four reference the same alias (dedup propagated to the
        # outer SELECT).
        agg_refs = re.findall(r"(SUM|MIN|MAX|AVG)\((__giql_p\d+)\)", outer_select)
        assert len(agg_refs) == 4
        assert len({alias for _, alias in agg_refs}) == 1
        # User aliases preserved on the outer SELECT.
        for label in ('"s"', '"mn"', '"mx"', '"av"'):
            assert label in outer_select

    def test_transpile_should_auto_project_group_by_column_into_inner_select_when_not_in_outer_select_when_dialect_is_duckdb(
        self,
    ):
        """Test that a GROUP BY column absent from SELECT is auto-projected.

        Given:
            A query ``SELECT COUNT(*) AS n FROM peaks a JOIN genes b ON
            ... GROUP BY a.chrom`` where ``a.chrom`` only appears in
            GROUP BY.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The inner subquery should project ``a."chrom"`` (so the
            wrapper relation has it available); the outer SELECT list
            should be ``COUNT(*)`` only; and GROUP BY should reference
            the inner alias rather than ``a.chrom``.
        """
        # Arrange & act
        sql = transpile(
            "SELECT COUNT(*) AS n FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "GROUP BY a.chrom",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert 'a."chrom"' in sql  # inner projection auto-added
        outer_select = sql.split(";\n", 1)[-1]
        select_clause = outer_select[: outer_select.index("FROM query(")]
        # Only COUNT(*) — no qualified column projection.
        assert "COUNT(*)" in select_clause
        assert "a.chrom" not in select_clause
        # GROUP BY references inner alias.
        group_clause = outer_select[outer_select.index("GROUP BY") :]
        assert "a.chrom" not in group_clause
        assert "__giql_p" in group_clause

    def test_transpile_should_emit_aggregate_of_expression_with_rewritten_inner_aliases_when_dialect_is_duckdb(
        self,
    ):
        """Test that ``SUM(a.start + a.end)`` rewrites both inner columns.

        Given:
            A query whose aggregate argument is an expression over two
            qualified columns: ``SUM(a.start + a.end)``.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The inner subquery should project both ``a."start"`` and
            ``a."end"`` under distinct inner aliases, the outer
            aggregate should reference both inner aliases inside the
            arithmetic expression, and the original ``a.start`` /
            ``a.end`` should not survive in the outer SELECT.
        """

        # Arrange & act
        sql = transpile(
            "SELECT a.chrom AS c, SUM(a.start + a.end) AS s FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "GROUP BY a.chrom",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        # Both genomic columns projected into the inner subquery.
        assert 'a."start"' in sql
        assert 'a."end"' in sql
        outer_select = sql.split(";\n", 1)[-1]
        # Outer SUM references two inner aliases under arithmetic.
        sum_match = re.search(r"SUM\(__giql_p\d+\s*\+\s*__giql_p\d+\)", outer_select)
        assert sum_match is not None
        # The two operand aliases must be distinct.
        operand_aliases = set(re.findall(r"__giql_p\d+", sum_match.group(0)))
        assert len(operand_aliases) == 2
        # No raw a.start / a.end in the outer SELECT.
        assert "a.start" not in outer_select
        assert "a.end" not in outer_select

    def test_transpile_should_route_right_outer_join_intersects_to_naive_predicate_plan_when_dialect_is_duckdb(
        self,
    ):
        """Test that ``RIGHT JOIN INTERSECTS`` falls back to the naive-predicate plan.

        Given:
            A query with a ``RIGHT JOIN`` whose ON contains a
            column-to-column INTERSECTS.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            ``_has_outer_join_intersects`` should detect the
            outer-join side and route to the naive-predicate plan.
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.start FROM peaks a "
            "RIGHT JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql

    def test_transpile_should_route_to_naive_predicate_plan_when_join_target_is_subquery(
        self,
    ):
        """Test that a subquery JOIN target routes to the naive-predicate plan.

        Given:
            A query whose JOIN target is a parenthesised subquery,
            ``JOIN (SELECT * FROM genes) b ON a.interval INTERSECTS
            b.interval``.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The dialect's "must be a base table" gate should fire and
            route to the naive-predicate plan.
        """
        # Arrange & act
        sql = transpile(
            "SELECT a.start FROM peaks a "
            "JOIN (SELECT * FROM genes) b "
            "ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql

    def test_transpile_should_decline_literal_projection_to_naive_plan(self):
        """Test that a literal-only SELECT-list entry declines to the naive plan.

        Given:
            A query whose SELECT list contains a bare integer literal
            (``SELECT 100 FROM peaks a JOIN ...``) — a constant projection
            with no column the IEJoin projection rebuild handles.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should decline the IEJoin (no ``SET VARIABLE`` scaffolding) so
            the naive-predicate plan projects the literal (#205).
        """
        # Act
        sql = transpile(
            "SELECT 100 FROM peaks a JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE" not in sql.upper()
        assert "getvariable" not in sql

    def test_transpile_should_decline_aliased_star_to_naive_plan(
        self,
    ):
        """Test that ``SELECT a.* AS x`` declines to the naive-predicate plan.

        Given:
            A query with a star projection that also has a user alias
            (``SELECT a.* AS x``).
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should decline the IEJoin (no ``SET VARIABLE`` scaffolding)
            like every star (#202) and defer the aliased-star shape to
            DuckDB — matching the naive-predicate plan rather than raising
            at transpile time.
        """
        # Act
        sql = transpile(
            "SELECT a.* AS x FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE" not in sql.upper()
        assert "getvariable" not in sql

    def test_transpile_should_surface_value_error_for_unqualified_projection_under_dialect(
        self,
    ):
        """Test that the dialect surfaces a plain ``ValueError`` on an unqualified column.

        Given:
            A query that triggers the dialect's unqualified-projection
            rejection (an unqualified column ``SELECT score``). Bare
            ``SELECT *`` no longer raises — stars decline to the naive
            plan (#202) — so an unqualified column is the shape that still
            exercises this contract.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The exception surfaced to the caller should be a
            ``ValueError`` and the docs-promised behavior holds: the
            ``"qualified"`` keyword appears in the message, and the
            class is reachable via ``isinstance(..., ValueError)`` (no
            internal subclass leaks across that contract).
        """
        # Act & assert
        with pytest.raises(ValueError) as excinfo:
            transpile(
                "SELECT score FROM peaks a "
                "JOIN genes b ON a.interval INTERSECTS b.interval",
                tables=["peaks", "genes"],
                dialect="duckdb",
            )
        assert isinstance(excinfo.value, ValueError)
        assert "qualified" in str(excinfo.value)

    def test_transpile_should_fall_back_to_naive_predicate_plan_when_join_is_natural(
        self,
    ):
        """Test that NATURAL JOIN falls back to the naive-predicate plan.

        Given:
            A query joining two tables with ``NATURAL JOIN`` and an
            INTERSECTS predicate in WHERE.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should emit SQL with no ``SET VARIABLE`` scaffolding
            (the dialect cannot enumerate NATURAL's implicit shared
            columns at transpile time, so it falls back).
        """
        # Arrange
        query = (
            "SELECT a.start FROM peaks a "
            "NATURAL JOIN genes b "
            "WHERE a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql

    def test_transpile_should_engage_dialect_when_using_join_targets_chrom_only(self):
        """Test that ``USING(chrom)`` engages the dialect path.

        Given:
            A query joining two tables with ``USING(chrom)`` and an
            INTERSECTS predicate in WHERE — the per-chromosome partition
            IS the equi-join that USING(chrom) requests.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should emit a ``SET VARIABLE __giql_iejoin_`` block.
        """
        # Arrange
        query = (
            "SELECT a.start FROM peaks a JOIN genes b USING (chrom) "
            "WHERE a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql

    def test_transpile_should_fall_back_when_using_join_targets_non_chrom_column(self):
        """Test that ``USING(<non-chrom>)`` falls back to the naive-predicate plan.

        Given:
            A query with ``USING(strand)`` plus INTERSECTS — the partition
            cannot satisfy a strand equi-join.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should fall back to the naive-predicate plan (no
            ``SET VARIABLE __giql_iejoin_`` block).
        """
        # Arrange
        query = (
            "SELECT a.start FROM peaks a JOIN genes b USING (strand) "
            "WHERE a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql

    def test_transpile_should_fall_back_when_using_join_targets_multiple_columns(self):
        """Test that multi-column ``USING(chrom, strand)`` falls back.

        Given:
            A query with ``USING(chrom, strand)`` plus INTERSECTS.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should fall back to the naive-predicate plan (multi-column USING
            inline support is a documented follow-up).
        """
        # Arrange
        query = (
            "SELECT a.start FROM peaks a "
            "JOIN genes b USING (chrom, strand) "
            "WHERE a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql

    def test_transpile_should_engage_dialect_when_join_is_implicit_cross_join(self):
        """Test that CROSS JOIN + WHERE INTERSECTS continues to engage the dialect.

        Given:
            A query joining two tables with ``CROSS JOIN`` and an
            INTERSECTS predicate in WHERE — equivalent to the implicit
            cross-join shape the dialect was designed for.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should emit a ``SET VARIABLE __giql_iejoin_`` block
            (regression guard so the new gate doesn't reject the
            legitimate cross-join shape).
        """
        # Arrange
        query = (
            "SELECT a.start FROM peaks a CROSS JOIN genes b "
            "WHERE a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql

    def test_transpile_should_fall_back_when_aggregate_argument_contains_subquery(
        self,
    ):
        """Test that a subquery inside an aggregate argument falls back.

        Given:
            A query whose SELECT list contains ``SUM((SELECT a.start
            FROM peaks))`` — the wrapper-relation rewriter would walk
            into the subquery and substitute wrapper aliases for refs
            that only resolve in the subquery's own scope.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should fall back to the naive-predicate plan.
        """
        # Arrange
        query = (
            "SELECT SUM((SELECT a.start FROM peaks)) "
            "FROM peaks a JOIN genes b ON a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql

    def test_transpile_should_fall_back_when_select_aggregate_contains_correlated_subquery(
        self,
    ):
        """Test that a correlated subquery inside an aggregate argument falls back.

        Given:
            ``SUM((SELECT a.start FROM genes WHERE genes.chrom = a.chrom))``
            in the SELECT list — same scope hazard as the uncorrelated
            case.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should fall back to the naive-predicate plan.
        """
        # Arrange
        query = (
            "SELECT SUM((SELECT a.start FROM genes "
            "WHERE genes.chrom = a.chrom)) "
            "FROM peaks a JOIN genes b ON a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE __giql_iejoin_" not in sql

    def test_transpile_should_resolve_projections_when_aliases_differ_in_case(
        self,
    ):
        """Test that mixed-case aliases are normalized to match references.

        Given:
            A query with ``peaks A`` / ``genes B`` (uppercase aliases)
            and references via ``A.start`` / ``B.interval``.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should emit a ``SET VARIABLE __giql_iejoin_`` block
            (case-folded alias matching aligns with DuckDB's identifier
            semantics, not Python's strict equality).
        """
        # Arrange
        query = (
            "SELECT A.start FROM peaks A "
            "JOIN genes B ON A.interval INTERSECTS B.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql

    def test_transpile_should_reject_extra_predicate_with_catalog_or_schema_qualified_column_when_dialect_is_duckdb(  # noqa: E501
        self,
    ):
        """Test that a catalog/schema-qualified column in an extra predicate raises.

        Given:
            A query whose extra predicate uses
            ``mycat.myschema.a.score > 0`` — the alias-only rewriter
            would leave the catalog/schema qualifier intact in the
            inner subquery where it addresses a different relation
            than intended.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should raise ``ValueError`` with a message naming the
            catalog/schema qualifier.
        """
        # Arrange
        query = (
            "SELECT a.start FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "AND mycat.myschema.a.score > 0"
        )

        # Act & assert
        with pytest.raises(ValueError) as excinfo:
            transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        assert "catalog/schema qualifier" in str(excinfo.value)

    def test_transpile_should_engage_dialect_for_semi_join_with_intersects(self):
        """Test that SEMI JOIN with INTERSECTS engages the dialect path.

        Given:
            A query joining two tables with ``SEMI JOIN`` and an
            INTERSECTS predicate.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should emit a ``SET VARIABLE __giql_iejoin_`` block.
        """
        # Arrange
        query = (
            "SELECT a.start FROM peaks a "
            "SEMI JOIN genes b ON a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql

    def test_transpile_should_emit_exists_when_join_is_semi(self):
        """Test that SEMI JOIN's per-chromosome template uses a correlated EXISTS.

        Given:
            A SEMI JOIN + INTERSECTS query.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The emitted SQL should express the semi-join as a correlated
            ``WHERE EXISTS`` subquery (which DuckDB plans through IE_JOIN),
            not a bare ``SEMI JOIN`` keyword (which DuckDB plans as a
            nested loop, #208).
        """
        # Arrange
        query = (
            "SELECT a.start FROM peaks a "
            "SEMI JOIN genes b ON a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "EXISTS" in sql
        assert "SEMI JOIN" not in sql

    def test_transpile_should_emit_not_exists_when_join_is_anti(self):
        """Test that ANTI JOIN's per-chromosome template uses a correlated NOT EXISTS.

        Given:
            An ANTI JOIN + INTERSECTS query.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The emitted SQL should express the anti-join as a correlated
            ``WHERE NOT EXISTS`` subquery (which DuckDB plans through
            IE_JOIN), not a bare ``ANTI JOIN`` keyword (which DuckDB plans
            as a nested loop, #208).
        """
        # Arrange
        query = (
            "SELECT a.start FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "NOT EXISTS" in sql
        assert "ANTI JOIN" not in sql

    def test_transpile_should_reject_right_side_projection_when_join_is_semi(self):
        """Test that SEMI JOIN raises on ``b.col`` in the outer SELECT.

        Given:
            A SEMI JOIN + INTERSECTS query whose outer SELECT projects
            ``b.start`` — but SEMI returns left-side rows only.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should raise ``ValueError`` naming the left-only
            constraint.
        """
        # Arrange
        query = (
            "SELECT b.start FROM peaks a "
            "SEMI JOIN genes b ON a.interval INTERSECTS b.interval"
        )

        # Act & assert
        with pytest.raises(ValueError) as excinfo:
            transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        assert "left-side" in str(excinfo.value)

    def test_transpile_should_engage_dialect_for_anti_join_with_intersects(self):
        """Test that ANTI JOIN with INTERSECTS engages the dialect path.

        Given:
            A query joining two tables with ``ANTI JOIN`` and an
            INTERSECTS predicate.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should emit a ``SET VARIABLE __giql_iejoin_`` block.
        """
        # Arrange
        query = (
            "SELECT a.start FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql

    def test_transpile_should_use_left_only_chrom_partition_when_join_is_anti(self):
        """Test that ANTI JOIN's partition source is left-distinct only.

        Given:
            An ANTI JOIN + INTERSECTS query.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            The emitted SET VARIABLE block should NOT contain
            ``INTERSECT`` (the partition source is left-distinct only,
            not the chromosome-INTERSECT used by INNER / SEMI).
        """
        # Arrange
        query = (
            "SELECT a.start FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        set_var_stmt = sql.split(";\n", 1)[0]

        # Assert
        assert "INTERSECT" not in set_var_stmt

    def test_transpile_should_reject_right_side_projection_when_join_is_anti(self):
        """Test that ANTI JOIN raises on ``b.col`` in the outer SELECT.

        Given:
            An ANTI JOIN + INTERSECTS query whose outer SELECT projects
            ``b.start``.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should raise ``ValueError`` naming the left-only
            constraint.
        """
        # Arrange
        query = (
            "SELECT b.start FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval"
        )

        # Act & assert
        with pytest.raises(ValueError) as excinfo:
            transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        assert "left-side" in str(excinfo.value)


class TestTranspileDuckDBIEJoinCountOverlaps:
    """The count_overlaps fast path: LEFT JOIN + COUNT(b.col) + GROUP BY (#209)."""

    def test_transpile_should_emit_zero_fill_wrapper_when_left_join_count(self):
        """Test that a LEFT-join COUNT over INTERSECTS emits the zero-fill fast path.

        Given:
            A ``SELECT a.cols, COUNT(b.col) ... FROM a LEFT JOIN b ON
            a.interval INTERSECTS b.interval GROUP BY a.cols`` query.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should emit the INNER IEJoin count wrapped in a zero-fill
            LEFT join (a ``SET VARIABLE`` block, a ``__giql_counts`` CTE, and
            a ``COALESCE(..., 0)``), never the naive ``a.chrom = b.chrom``
            hash-join predicate.
        """
        # Arrange
        query = (
            'SELECT a.chrom, a.start, a."end", COUNT(b.chrom) AS n '
            "FROM peaks a LEFT JOIN genes b ON a.interval INTERSECTS b.interval "
            'GROUP BY a.chrom, a.start, a."end"'
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert "__giql_counts" in sql
        assert "COALESCE" in sql
        assert 'a."chrom" = b."chrom"' not in sql

    def test_transpile_should_decline_count_overlaps_when_count_star(self):
        """Test that COUNT(*) over a LEFT join declines to the naive plan.

        Given:
            A LEFT-join query aggregating ``COUNT(*)`` rather than a
            right-side column.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should decline the fast path (no ``__giql_counts`` wrapper),
            because COUNT(*) over a LEFT join counts the null-padded row for
            unmatched left rows, which the zero-fill rewrite cannot express.
        """
        # Arrange
        query = (
            "SELECT a.chrom, COUNT(*) AS n "
            "FROM a LEFT JOIN b ON a.interval INTERSECTS b.interval GROUP BY a.chrom"
        )

        # Act
        sql = transpile(query, tables=["a", "b"], dialect="duckdb")

        # Assert
        assert "__giql_counts" not in sql

    def test_transpile_should_decline_count_overlaps_when_non_count_aggregate(self):
        """Test that a non-COUNT aggregate over a LEFT join declines to the naive plan.

        Given:
            A LEFT-join query aggregating ``SUM(b.score)``.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should decline the fast path (no ``__giql_counts`` wrapper);
            the fast path supports only ``COUNT`` in v1.
        """
        # Arrange
        query = (
            "SELECT a.chrom, SUM(b.score) AS s "
            "FROM a LEFT JOIN b ON a.interval INTERSECTS b.interval GROUP BY a.chrom"
        )

        # Act
        sql = transpile(query, tables=["a", "b"], dialect="duckdb")

        # Assert
        assert "__giql_counts" not in sql

    def test_transpile_should_decline_count_overlaps_when_no_group_by(self):
        """Test that a LEFT-join COUNT without GROUP BY declines to the naive plan.

        Given:
            A LEFT-join COUNT query with no GROUP BY clause.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should decline the fast path (no ``__giql_counts`` wrapper).
        """
        # Arrange
        query = (
            "SELECT COUNT(b.chrom) AS n "
            "FROM a LEFT JOIN b ON a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["a", "b"], dialect="duckdb")

        # Assert
        assert "__giql_counts" not in sql

    def test_transpile_should_decline_count_overlaps_when_order_by(self):
        """Test that a LEFT-join COUNT with ORDER BY declines to the naive plan.

        Given:
            A LEFT-join COUNT query carrying a top-level ORDER BY.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should decline the fast path (no ``__giql_counts`` wrapper);
            result-ordering clauses are out of scope for v1.
        """
        # Arrange
        query = (
            "SELECT a.chrom, COUNT(b.chrom) AS n "
            "FROM a LEFT JOIN b ON a.interval INTERSECTS b.interval "
            "GROUP BY a.chrom ORDER BY a.chrom"
        )

        # Act
        sql = transpile(query, tables=["a", "b"], dialect="duckdb")

        # Assert
        assert "__giql_counts" not in sql

    def test_transpile_should_decline_count_overlaps_when_self_join(self):
        """Test that a self-join LEFT COUNT declines to the naive plan.

        Given:
            A LEFT-join COUNT whose two sides are the same table (a self
            overlap count).
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should decline the fast path (no ``__giql_counts`` wrapper),
            because the reused INNER path declines the self-join and the
            count builder falls back to the naive plan.
        """
        # Arrange
        query = (
            "SELECT a.chrom, COUNT(b.chrom) AS n "
            "FROM t a LEFT JOIN t b ON a.interval INTERSECTS b.interval GROUP BY a.chrom"
        )

        # Act
        sql = transpile(query, tables=["t"], dialect="duckdb")

        # Assert
        assert "__giql_counts" not in sql

    @settings(
        max_examples=25,
        deadline=None,
        suppress_health_check=[HealthCheck.function_scoped_fixture],
    )
    @given(
        peaks=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2", "chr3"]),
                st.integers(min_value=0, max_value=100),
                st.integers(min_value=1, max_value=30),
            ),
            min_size=0,
            max_size=8,
        ),
        genes=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2", "chr3"]),
                st.integers(min_value=0, max_value=100),
                st.integers(min_value=1, max_value=30),
            ),
            min_size=0,
            max_size=8,
        ),
    )
    def test_count_overlaps_should_match_python_reference_for_random_inputs(
        self, conn, peaks, genes
    ):
        """Test the count_overlaps fast path against a Python-native reference.

        Given:
            A Hypothesis-generated pair of small interval lists.
        When:
            A ``COUNT(b.chrom)`` LEFT-join over INTERSECTS is transpiled with
            ``dialect='duckdb'`` and executed.
        Then:
            The returned ``(chrom, start, end, count)`` rows should equal the
            Python count_overlaps reference, including left rows with no
            overlap (count 0), chromosomes present on only one side, an empty
            right table, and duplicate left rows.
        """
        # Arrange
        peak_rows = [(c, s, s + length) for (c, s, length) in peaks]
        gene_rows = [(c, s, s + length) for (c, s, length) in genes]
        conn.execute("DROP TABLE IF EXISTS peaks")
        conn.execute("DROP TABLE IF EXISTS genes")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        if peak_rows:
            conn.executemany("INSERT INTO peaks VALUES (?, ?, ?)", peak_rows)
        conn.execute(
            'CREATE TABLE genes (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        if gene_rows:
            conn.executemany("INSERT INTO genes VALUES (?, ?, ?)", gene_rows)
        sql = transpile(
            'SELECT a.chrom, a.start, a."end", COUNT(b.chrom) AS n '
            "FROM peaks a LEFT JOIN genes b ON a.interval INTERSECTS b.interval "
            'GROUP BY a.chrom, a.start, a."end"',
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())
        expected = _python_count_overlaps(peak_rows, gene_rows)

        # Assert
        assert rows == expected

    def test_count_overlaps_plan_should_avoid_nested_loop_when_executed_on_duckdb(
        self, conn
    ):
        """Test that the count_overlaps fast path is not planned as a nested loop.

        Given:
            A dense single-chromosome peaks/genes pair large enough that a
            naive LEFT-join count would be planned as a quadratic
            ``BLOCKWISE_NL_JOIN``.
        When:
            The ``dialect='duckdb'`` count query is generated and its inner
            per-chromosome dynamic SQL is planned with ``EXPLAIN``.
        Then:
            It should avoid ``BLOCKWISE_NL_JOIN`` — the INNER IEJoin count
            routes through DuckDB's range-join family (#209).
        """
        # Arrange
        _make_table(
            conn, "peaks", [("chr1", i * 7, i * 7 + 5, "p", 1, "+") for i in range(4000)]
        )
        _make_table(
            conn,
            "genes",
            [("chr1", i * 7 + 2, i * 7 + 9, "g", 1, "+") for i in range(4000)],
        )
        sql = transpile(
            'SELECT a.chrom, a.start, a."end", COUNT(b.chrom) AS n '
            "FROM peaks a LEFT JOIN genes b ON a.interval INTERSECTS b.interval "
            'GROUP BY a.chrom, a.start, a."end"',
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        plan = _explain_dynamic_sql(conn, sql)

        # Assert
        assert "BLOCKWISE_NL_JOIN" not in plan


class TestTranspileDuckDBIEJoinExecution:
    """End-to-end execution tests against in-memory DuckDB."""

    def test_semi_join_plan_should_avoid_nested_loop_when_executed_on_duckdb(self, conn):
        """Test that the SEMI-join EXISTS rewrite is not planned as a nested loop.

        Given:
            A dense single-chromosome peaks/genes pair large enough that a
            bare ``SEMI JOIN`` inequality overlap would be planned as a
            quadratic ``BLOCKWISE_NL_JOIN``.
        When:
            The ``dialect='duckdb'`` SEMI-join SQL is generated and its
            per-chromosome dynamic SQL is planned with ``EXPLAIN``.
        Then:
            It should avoid ``BLOCKWISE_NL_JOIN`` — the correlated
            ``WHERE EXISTS`` form routes through DuckDB's range-join family
            instead of a nested loop (#208).
        """
        # Arrange
        _make_table(
            conn, "peaks", [("chr1", i * 7, i * 7 + 5, "p", 1, "+") for i in range(4000)]
        )
        _make_table(
            conn,
            "genes",
            [("chr1", i * 7 + 2, i * 7 + 9, "g", 1, "+") for i in range(4000)],
        )
        sql = transpile(
            "SELECT a.start FROM peaks a "
            "SEMI JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        plan = _explain_dynamic_sql(conn, sql)

        # Assert
        assert "BLOCKWISE_NL_JOIN" not in plan

    def test_anti_join_plan_should_avoid_nested_loop_when_executed_on_duckdb(self, conn):
        """Test that the ANTI-join NOT EXISTS rewrite is not planned as a nested loop.

        Given:
            A dense single-chromosome peaks/genes pair large enough that a
            bare ``ANTI JOIN`` inequality overlap would be planned as a
            quadratic ``BLOCKWISE_NL_JOIN``.
        When:
            The ``dialect='duckdb'`` ANTI-join SQL is generated and its
            per-chromosome dynamic SQL is planned with ``EXPLAIN``.
        Then:
            It should avoid ``BLOCKWISE_NL_JOIN`` — the correlated
            ``WHERE NOT EXISTS`` form routes through DuckDB's range-join
            family instead of a nested loop (#208).
        """
        # Arrange
        _make_table(
            conn, "peaks", [("chr1", i * 7, i * 7 + 5, "p", 1, "+") for i in range(4000)]
        )
        _make_table(
            conn,
            "genes",
            [("chr1", i * 7 + 2, i * 7 + 9, "g", 1, "+") for i in range(4000)],
        )
        sql = transpile(
            "SELECT a.start FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        plan = _explain_dynamic_sql(conn, sql)

        # Assert
        assert "BLOCKWISE_NL_JOIN" not in plan

    def test_query_should_return_overlapping_pairs_when_executed_on_duckdb(
        self, peaks_genes
    ):
        """Test the multi-chromosome inner-join overlap truth table.

        Given:
            A peaks/genes fixture with a non-trivial truth table of
            overlapping pairs across multiple chromosomes.
        When:
            The DuckDB-dialect SQL is executed against the connection.
        Then:
            It should return exactly the expected sorted set of
            overlapping pairs.
        """
        # Arrange
        sql = transpile(
            """
            SELECT a.chrom AS a_chrom, a.start AS a_start, a.end AS a_end,
                   b.chrom AS b_chrom, b.start AS b_start, b.end AS b_end
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(peaks_genes.execute(sql).fetchall())

        # Assert
        assert rows == sorted(_EXPECTED_OVERLAPS_PEAKS_GENES)

    def test_query_should_return_empty_set_when_no_chromosomes_intersect(self, conn):
        """Test the empty-schema fallback when no chromosomes intersect.

        Given:
            Two tables with disjoint chromosome sets.
        When:
            The DuckDB-dialect SQL is executed.
        Then:
            It should return an empty result without raising (the
            COALESCE empty-schema fallback fires).
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p1", 0, "+")])
        _make_table(conn, "genes", [("chr2", 100, 200, "g1", 0, "+")])
        sql = transpile(
            """
            SELECT a.chrom AS a_chrom, a.start AS a_start, a.end AS a_end,
                   b.chrom AS b_chrom, b.start AS b_start, b.end AS b_end
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert rows == []

    def test_query_should_match_naive_predicate_plan_results_when_executed_on_duckdb(
        self, peaks_genes
    ):
        """Test that the naive-predicate plan and the IEJoin path agree on results.

        Given:
            The ``peaks_genes`` fixture and a column-to-column
            INTERSECTS INNER JOIN.
        When:
            The query is transpiled both with ``dialect=None`` (naive-predicate
            plan) and with ``dialect='duckdb'`` (IEJoin path) and both
            are executed.
        Then:
            Both should produce the same sorted set of overlapping
            (peak, gene) rows.
        """
        # Arrange
        query = """
            SELECT a.chrom AS a_chrom, a.start AS a_start, a.end AS a_end,
                   b.chrom AS b_chrom, b.start AS b_start, b.end AS b_end
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """
        generic_sql = transpile(query, tables=["peaks", "genes"])
        duckdb_sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Act
        generic_rows = sorted(peaks_genes.execute(generic_sql).fetchall())
        duckdb_rows = sorted(peaks_genes.execute(duckdb_sql).fetchall())

        # Assert
        assert generic_rows == duckdb_rows
        assert duckdb_rows == sorted(_EXPECTED_OVERLAPS_PEAKS_GENES)

    def test_query_should_handle_chrom_with_single_quote_in_name(self, conn):
        """Test that single quotes in chrom names are escaped correctly.

        Given:
            Intervals on a chromosome whose name contains a single
            quote.
        When:
            The DuckDB-dialect SQL is executed.
        Then:
            It should escape the quote inside the dynamic SQL builder
            and return the matching pair byte-for-byte.
        """
        # Arrange
        _make_table(conn, "peaks", [("chr'1", 100, 200, "p1", 0, "+")])
        _make_table(conn, "genes", [("chr'1", 150, 250, "g1", 0, "+")])
        sql = transpile(
            """
            SELECT a.chrom AS a_chrom, a.start AS a_start, a.end AS a_end,
                   b.chrom AS b_chrom, b.start AS b_start, b.end AS b_end
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert rows == [("chr'1", 100, 200, "chr'1", 150, 250)]

    def test_query_should_handle_one_based_closed_intervals(self, conn):
        """Test that 1-based closed intervals report touching endpoints as overlapping.

        Given:
            Tables declared as 1-based closed-closed intervals where
            ``peaks[100, 200]`` and ``genes[200, 300]`` share endpoint
            200.
        When:
            The DuckDB-dialect SQL is executed.
        Then:
            It should report the pair as overlapping (closed-closed
            shares endpoint 200).
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p1", 0, "+")])
        _make_table(conn, "genes", [("chr1", 200, 300, "g1", 0, "+")])
        sql = transpile(
            """
            SELECT a.chrom AS a_chrom, a.start AS a_start, a.end AS a_end,
                   b.chrom AS b_chrom, b.start AS b_start, b.end AS b_end
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", coordinate_system="1based", interval_type="closed"),
                Table("genes", coordinate_system="1based", interval_type="closed"),
            ],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert ("chr1", 100, 200, "chr1", 200, 300) in rows

    def test_query_should_not_report_non_touching_one_based_closed_as_overlapping(
        self, conn
    ):
        """Test that strictly disjoint 1-based closed intervals do not match.

        Given:
            Tables declared as 1-based closed-closed with non-touching
            intervals (``peaks[100, 199]`` and ``genes[200, 300]``).
        When:
            The DuckDB-dialect SQL is executed.
        Then:
            It should report no overlap.
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 199, "p1", 0, "+")])
        _make_table(conn, "genes", [("chr1", 200, 300, "g1", 0, "+")])
        sql = transpile(
            """
            SELECT a.chrom AS a_chrom, a.start AS a_start, a.end AS a_end,
                   b.chrom AS b_chrom, b.start AS b_start, b.end AS b_end
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", coordinate_system="1based", interval_type="closed"),
                Table("genes", coordinate_system="1based", interval_type="closed"),
            ],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert rows == []

    def test_query_should_handle_zero_based_closed_intervals(self, conn):
        """Test that 0-based closed intervals match abutting endpoints.

        Given:
            Tables declared as 0-based closed-closed intervals where
            ``peaks[100, 200]`` abuts ``genes[200, 300]``.
        When:
            The DuckDB-dialect SQL is executed.
        Then:
            It should report the pair as overlapping (closed-closed
            shares endpoint 200).
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p1", 0, "+")])
        _make_table(conn, "genes", [("chr1", 200, 300, "g1", 0, "+")])
        sql = transpile(
            """
            SELECT a.chrom AS a_chrom, a.start AS a_start, a.end AS a_end,
                   b.chrom AS b_chrom, b.start AS b_start, b.end AS b_end
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", coordinate_system="0based", interval_type="closed"),
                Table("genes", coordinate_system="0based", interval_type="closed"),
            ],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert ("chr1", 100, 200, "chr1", 200, 300) in rows

    def test_query_should_handle_implicit_cross_join_when_dialect_is_duckdb(
        self, peaks_genes
    ):
        """Test that an implicit cross-join INTERSECTS produces the truth table.

        Given:
            An implicit cross-join (``FROM peaks a, genes b WHERE
            a.interval INTERSECTS b.interval``) over the multi-chromosome
            ``peaks_genes`` fixture.
        When:
            The DuckDB-dialect SQL is executed.
        Then:
            It should return exactly the expected sorted truth table of
            overlapping pairs.
        """
        # Arrange
        sql = transpile(
            """
            SELECT a.chrom AS a_chrom, a.start AS a_start, a.end AS a_end,
                   b.chrom AS b_chrom, b.start AS b_start, b.end AS b_end
            FROM peaks a, genes b
            WHERE a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(peaks_genes.execute(sql).fetchall())

        # Assert
        assert rows == sorted(_EXPECTED_OVERLAPS_PEAKS_GENES)

    # --- DX-001..DX-011: behavioral truth tables -------------------------

    def test_query_should_preserve_all_columns_for_qualified_star_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that ``a.*, b.*`` preserves every base-table column, matching the naive plan.

        Given:
            A ``SELECT a.*, b.*`` projection over an INTERSECTS join, where
            each table carries non-genomic columns (``name`` / ``score``)
            beyond the configured chrom / start / end / strand.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed,
            and separately with ``dialect=None``.
        Then:
            The DuckDB path should decline the IEJoin (#202) so the star
            expands against the live schema — the result carries every
            base-table column from both sides (not just the genomic subset),
            identical in columns and rows to the naive-predicate plan.
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p1", 7, "+")])
        _make_table(conn, "genes", [("chr1", 150, 250, "g1", 9, "-")])
        query = """
            SELECT a.*, b.*
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """
        sql_dd = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act
        result_dd = conn.execute(sql_dd)
        cols_dd = [d[0] for d in result_dd.description]
        rows_dd = result_dd.fetchall()
        result_naive = conn.execute(sql_naive)
        cols_naive = [d[0] for d in result_naive.description]
        rows_naive = result_naive.fetchall()

        # Assert
        assert "SET VARIABLE" not in sql_dd.upper()
        assert cols_dd == ["chrom", "start", "end", "name", "score", "strand"] * 2
        assert rows_dd == [
            ("chr1", 100, 200, "p1", 7, "+", "chr1", 150, 250, "g1", 9, "-")
        ]
        assert cols_dd == cols_naive
        assert rows_dd == rows_naive

    def test_query_should_propagate_user_alias_on_qualified_column(self, conn):
        """Test that a ``a.start AS s`` alias survives the IEJoin rewrite.

        Given:
            A ``SELECT a.start AS s`` projection over an INTERSECTS
            join.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            The result column should be named ``s`` and carry the value
            of ``a.start`` for the matching pair.
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p", 0, "+")])
        _make_table(conn, "genes", [("chr1", 150, 250, "g", 0, "+")])
        sql = transpile(
            """
            SELECT a.start AS s
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        result = conn.execute(sql)
        column_names = [d[0] for d in result.description]
        rows = result.fetchall()

        # Assert
        assert column_names == ["s"]
        assert rows == [(100,)]

    def test_query_should_treat_touching_half_open_intervals_as_non_overlapping(
        self, conn
    ):
        """Test that half-open intervals that share an endpoint do not overlap.

        Given:
            Two half-open tables where ``peaks[100, 200)`` ends exactly
            where ``genes[200, 300)`` begins.
        When:
            The DuckDB-dialect SQL is executed.
        Then:
            It should return an empty result (touching half-open
            intervals do not overlap).
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p", 0, "+")])
        _make_table(conn, "genes", [("chr1", 200, 300, "g", 0, "+")])
        sql = transpile(
            """
            SELECT a.start AS a_s, b.start AS b_s
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", coordinate_system="0based", interval_type="half_open"),
                Table("genes", coordinate_system="0based", interval_type="half_open"),
            ],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert rows == []

    def test_query_should_handle_mixed_coordinate_systems(self, conn):
        """Test semantic equivalence across mixed coordinate systems.

        Given:
            ``peaks`` declared as 1-based closed (``[100, 200]``) and
            ``genes`` declared as 0-based half-open (``[99, 200)``)
            describing the same genomic span.
        When:
            The DuckDB-dialect SQL is executed.
        Then:
            It should report the pair as overlapping (both canonicalize
            to the same 0-based half-open span).
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p", 0, "+")])
        _make_table(conn, "genes", [("chr1", 99, 200, "g", 0, "+")])
        sql = transpile(
            """
            SELECT a.start AS a_s, b.start AS b_s
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table("peaks", coordinate_system="1based", interval_type="closed"),
                Table("genes", coordinate_system="0based", interval_type="half_open"),
            ],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert rows == [(100, 99)]

    def test_query_should_apply_one_based_offset_when_only_left_table_is_one_based(
        self, conn
    ):
        """Test that the 1-based offset is applied per side, not globally.

        Given:
            ``peaks`` declared as 1-based half-open and ``genes`` declared
            as 0-based half-open, with values that overlap iff the offset
            is applied to the peaks side
            (``peaks[100, 101)`` and ``genes[99, 100)``).
        When:
            The DuckDB-dialect SQL is executed.
        Then:
            It should report the overlap; if the offset were skipped no
            overlap would be reported.
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 101, "p", 0, "+")])
        _make_table(conn, "genes", [("chr1", 99, 100, "g", 0, "+")])
        sql_offset = transpile(
            """
            SELECT a.start AS a_s, b.start AS b_s
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table(
                    "peaks",
                    coordinate_system="1based",
                    interval_type="half_open",
                ),
                Table(
                    "genes",
                    coordinate_system="0based",
                    interval_type="half_open",
                ),
            ],
            dialect="duckdb",
        )
        sql_no_offset = transpile(
            """
            SELECT a.start AS a_s, b.start AS b_s
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table(
                    "peaks",
                    coordinate_system="0based",
                    interval_type="half_open",
                ),
                Table(
                    "genes",
                    coordinate_system="0based",
                    interval_type="half_open",
                ),
            ],
            dialect="duckdb",
        )

        # Act
        rows_offset = conn.execute(sql_offset).fetchall()
        rows_no_offset = conn.execute(sql_no_offset).fetchall()

        # Assert
        assert rows_offset == [(100, 99)]
        assert rows_no_offset == []

    def test_query_should_return_empty_when_left_table_is_empty(self, conn):
        """Test the empty-left input case for the IEJoin path.

        Given:
            ``peaks`` empty and ``genes`` populated with one row.
        When:
            The DuckDB-dialect SQL is executed.
        Then:
            It should return an empty result without raising.
        """
        # Arrange
        _make_table(conn, "peaks", [])
        _make_table(conn, "genes", [("chr1", 150, 250, "g1", 0, "+")])
        sql = transpile(
            """
            SELECT a.start AS a_s, b.start AS b_s
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert rows == []

    def test_query_should_return_empty_when_right_table_is_empty(self, conn):
        """Test the empty-right input case for the IEJoin path.

        Given:
            ``peaks`` populated with one row and ``genes`` empty.
        When:
            The DuckDB-dialect SQL is executed.
        Then:
            It should return an empty result without raising.
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p1", 0, "+")])
        _make_table(conn, "genes", [])
        sql = transpile(
            """
            SELECT a.start AS a_s, b.start AS b_s
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert rows == []

    def test_query_should_return_typed_empty_result_when_chromosome_intersection_is_empty(
        self, conn
    ):
        """Test that the empty-schema fallback exposes correctly typed columns.

        Given:
            Two tables whose chromosome sets are disjoint so the
            chromosome-intersection ``string_agg`` is NULL.
        When:
            The DuckDB-dialect SQL is executed.
        Then:
            It should fetch an empty result whose column descriptions
            match the source-table schema declared by ``_make_table``
            (chrom as VARCHAR, start/end as INTEGER).
        """
        # Arrange
        _make_table(conn, "peaks", [("chrA", 100, 200, "p", 0, "+")])
        _make_table(conn, "genes", [("chrB", 100, 200, "g", 0, "+")])
        sql = transpile(
            """
            SELECT a.chrom AS a_chrom, a.start AS a_start, a.end AS a_end,
                   b.chrom AS b_chrom, b.start AS b_start, b.end AS b_end
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        result = conn.execute(sql)
        rows = result.fetchall()
        type_by_name = {d[0]: str(d[1]) for d in result.description}

        # Assert
        assert rows == []
        assert type_by_name["a_chrom"] == "VARCHAR"
        assert type_by_name["b_chrom"] == "VARCHAR"
        assert type_by_name["a_start"] == "INTEGER"
        assert type_by_name["a_end"] == "INTEGER"
        assert type_by_name["b_start"] == "INTEGER"
        assert type_by_name["b_end"] == "INTEGER"

    def test_query_should_type_non_genomic_columns_in_empty_schema_fallback(self, conn):
        """Test that non-genomic projected columns retain their real types when empty.

        Given:
            Two tables whose chromosome sets are disjoint, projecting a
            non-genomic column whose declared type is INTEGER (not VARCHAR
            or BIGINT).
        When:
            The DuckDB-dialect SQL is executed.
        Then:
            It should fetch an empty result whose ``score`` column is typed
            as INTEGER, matching the source table schema.
        """
        # Arrange
        conn.execute(
            "CREATE TABLE peaks (chrom VARCHAR, start BIGINT, "
            '"end" BIGINT, score INTEGER)'
        )
        conn.execute(
            "CREATE TABLE genes (chrom VARCHAR, start BIGINT, "
            '"end" BIGINT, score INTEGER)'
        )
        conn.execute("INSERT INTO peaks VALUES ('chrA', 100, 200, 7)")
        conn.execute("INSERT INTO genes VALUES ('chrB', 100, 200, 9)")
        sql = transpile(
            """
            SELECT a.score AS a_score, b.score AS b_score
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        result = conn.execute(sql)
        rows = result.fetchall()
        type_by_name = {d[0]: str(d[1]) for d in result.description}

        # Assert
        assert rows == []
        assert type_by_name["a_score"] == "INTEGER"
        assert type_by_name["b_score"] == "INTEGER"

    def test_query_should_preserve_outer_join_semantics_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that LEFT JOIN INTERSECTS preserves unmatched left rows.

        Given:
            A LEFT JOIN of ``peaks`` against ``genes`` on INTERSECTS,
            with one peak that has no matching gene.
        When:
            The query is transpiled with ``dialect='duckdb'`` and
            executed.
        Then:
            The unmatched peak should appear in the result with NULL
            gene columns (the IEJoin path falls back to the naive-predicate plan
            for outer joins).
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p_match", 0, "+"),
                ("chr1", 5000, 6000, "p_lonely", 0, "+"),
            ],
        )
        _make_table(conn, "genes", [("chr1", 150, 250, "g1", 0, "+")])
        sql = transpile(
            """
            SELECT a.start AS a_start, b.start AS b_start
            FROM peaks a
            LEFT JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall(), key=lambda r: r[0])

        # Assert
        assert rows == [(100, 150), (5000, None)]

    def test_query_should_preserve_extra_join_predicate_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that an extra non-INTERSECTS predicate is honored end-to-end.

        Given:
            A query whose JOIN ON carries
            ``a.interval INTERSECTS b.interval AND a.score > b.score``,
            with input rows where one pair overlaps and satisfies the
            score predicate and another overlaps but violates it.
        When:
            The query is transpiled with ``dialect='duckdb'`` and
            executed.
        Then:
            Only the pair that overlaps and satisfies the extra
            predicate should be returned.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p_high", 100, "+"),
                ("chr1", 300, 400, "p_low", 1, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g_low", 5, "+"),
                ("chr1", 350, 450, "g_high", 50, "+"),
            ],
        )
        sql = transpile(
            """
            SELECT a.start AS a_start, b.start AS b_start
            FROM peaks a
            JOIN genes b
              ON a.interval INTERSECTS b.interval
             AND a.score > b.score
            """,
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert rows == [(100, 150)]

    def test_query_should_return_correct_rows_for_literal_intersects_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that literal-range INTERSECTS produces the same rows as default.

        Given:
            A WHERE literal-range INTERSECTS query over ``peaks``.
        When:
            The query is transpiled with ``dialect=None`` and with
            ``dialect='duckdb'`` and both are executed.
        Then:
            Both should return the same sorted rows.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 0, "+"),
                ("chr1", 1500, 1700, "p2", 0, "+"),
                ("chr1", 5000, 6000, "p3", 0, "+"),
            ],
        )
        query = "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'"

        # Act
        rows_default = sorted(
            conn.execute(transpile(query, tables=["peaks"])).fetchall()
        )
        rows_duckdb = sorted(
            conn.execute(transpile(query, tables=["peaks"], dialect="duckdb")).fetchall()
        )

        # Assert
        assert rows_default == rows_duckdb
        assert rows_duckdb == [("chr1", 1500, 1700, "p2", 0, "+")]

    def test_query_should_partition_on_custom_chrom_col_when_table_uses_custom_columns(
        self, conn
    ):
        """Test that the IEJoin partition uses the configured chrom column.

        Given:
            DuckDB tables with ``chromosome``/``start_pos``/``end_pos``
            custom column names and rows on multiple chromosomes (only
            one of which has overlapping intervals).
        When:
            The query is transpiled with the matching ``Table`` configs
            and ``dialect='duckdb'`` and executed.
        Then:
            It should return only the rows where the custom-named
            chromosome column matches and the ranges overlap.
        """
        # Arrange
        conn.execute(
            "CREATE TABLE peaks (chromosome VARCHAR, start_pos INTEGER, end_pos INTEGER)"
        )
        conn.execute(
            "INSERT INTO peaks VALUES "
            "('chr1', 100, 200), "
            "('chr2', 100, 200), "
            "('chr3', 100, 200)"
        )
        conn.execute(
            "CREATE TABLE genes (chromosome VARCHAR, start_pos INTEGER, end_pos INTEGER)"
        )
        conn.execute(
            "INSERT INTO genes VALUES "
            "('chr1', 150, 250), "
            "('chr2', 1000, 2000), "
            "('chr3', 5000, 6000)"
        )
        sql = transpile(
            """
            SELECT a.chromosome AS a_chrom, a.start_pos AS a_start,
                   b.start_pos AS b_start
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=[
                Table(
                    "peaks",
                    chrom_col="chromosome",
                    start_col="start_pos",
                    end_col="end_pos",
                ),
                Table(
                    "genes",
                    chrom_col="chromosome",
                    start_col="start_pos",
                    end_col="end_pos",
                ),
            ],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert rows == [("chr1", 100, 150)]

    def test_query_should_preserve_row_multiplicity_for_duplicate_input_rows(self, conn):
        """Test that duplicate input rows surface as duplicate output rows.

        Given:
            Two identical ``peaks`` rows that each overlap the same single
            ``genes`` row, projecting genomic columns from both sides.
        When:
            The DuckDB-dialect SQL is executed.
        Then:
            The result should contain *two* identical overlap rows, not
            one — a regression detector for a rewrite that silently
            collapses duplicates via an unintended ``DISTINCT``. (The
            property-based oracle below sorts a multiset, so this is the
            canonical-example pin for the same invariant.)
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p_dup", 1, "+"),
                ("chr1", 100, 200, "p_dup", 1, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [("chr1", 150, 250, "g1", 10, "-")],
        )
        sql = transpile(
            "SELECT a.chrom AS c, a.start AS s, a.end AS e "
            "FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert rows == [("chr1", 100, 200), ("chr1", 100, 200)]

    def test_query_should_order_lexicographically_when_order_by_has_two_columns(
        self, conn
    ):
        """Test that ORDER BY on two keys orders lexicographically.

        Given:
            Multiple overlap rows where the first ORDER BY key has ties
            and the second key discriminates among them.
        When:
            The query is transpiled with ``dialect='duckdb'`` and
            executed.
        Then:
            Rows should come back in lexicographic order on the two
            keys.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),
                ("chr1", 100, 200, "p2", 2, "+"),
                ("chr1", 300, 400, "p3", 3, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g_a", 10, "+"),
                ("chr1", 350, 450, "g_b", 20, "+"),
            ],
        )
        sql = transpile(
            "SELECT a.start AS s, b.score AS gs FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "ORDER BY a.start ASC, b.score DESC",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        # Two chr1 peaks at start=100 overlap g_a (score 10); chr1 peak at
        # start=300 overlaps g_b (score 20). Within the start=100 ties,
        # b.score DESC is a no-op (only one matching gene), so the keyed
        # order is: (100, 10), (100, 10), then (300, 20).
        assert rows == [(100, 10), (100, 10), (300, 20)]

    def test_query_should_skip_first_m_and_return_remainder_when_offset_only(self, conn):
        """Test that OFFSET past the first ordered row returns the rest.

        Given:
            Three ordered overlap rows on chr1.
        When:
            The query is transpiled with ``LIMIT 5 OFFSET 1`` and
            executed.
        Then:
            The trailing two rows (after skipping the first) are
            returned, in order.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),
                ("chr1", 300, 400, "p2", 2, "+"),
                ("chr1", 500, 600, "p3", 3, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 1, "+"),
                ("chr1", 350, 450, "g2", 2, "+"),
                ("chr1", 550, 650, "g3", 3, "+"),
            ],
        )
        sql = transpile(
            "SELECT a.start AS s FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "ORDER BY a.start LIMIT 5 OFFSET 1",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert rows == [(300,), (500,)]

    def test_query_should_dedup_on_tuple_when_select_distinct_has_multiple_columns(
        self, conn
    ):
        """Test that ``SELECT DISTINCT`` over a multi-column list dedups on tuples.

        Given:
            Multiple overlap rows that produce duplicate
            ``(a.chrom, b.chrom)`` tuples but distinct ``a.start``
            values.
        When:
            The query is transpiled with ``SELECT DISTINCT a.chrom AS
            ca, b.chrom AS cb`` and executed.
        Then:
            The result should collapse to the set of distinct tuples,
            not deduplicate on either column alone.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),
                ("chr1", 300, 400, "p2", 2, "+"),
                ("chr2", 100, 200, "p3", 3, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 1, "+"),
                ("chr1", 350, 450, "g2", 2, "+"),
                ("chr2", 150, 250, "g3", 3, "+"),
            ],
        )
        sql = transpile(
            "SELECT DISTINCT a.chrom AS ca, b.chrom AS cb FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = set(conn.execute(sql).fetchall())

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        # Two chr1 join rows collapse to one ('chr1','chr1') tuple; chr2
        # contributes one ('chr2','chr2') tuple. Total: 2 distinct tuples.
        assert rows == {("chr1", "chr1"), ("chr2", "chr2")}

    def test_query_should_filter_cross_side_predicate_in_where_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that a cross-side predicate in WHERE filters rows correctly.

        Given:
            Three overlapping pairs on chr1 where only those whose
            peak's ``start`` is strictly less than the gene's ``start``
            satisfy the cross-side WHERE predicate. (Variant of the
            ON-side cross-predicate test — exercises the WHERE-extras
            extraction path rather than the JOIN-ON extras path.)
        When:
            The query is transpiled with ``dialect='duckdb'`` and
            executed.
        Then:
            Only the pair satisfying the cross-side WHERE predicate is
            returned.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                (
                    "chr1",
                    100,
                    200,
                    "p_lt",
                    1,
                    "+",
                ),  # a.start (100) < b.start (150) → keep
                (
                    "chr1",
                    300,
                    400,
                    "p_eq",
                    2,
                    "+",
                ),  # a.start (300) < b.start (350) → keep
                (
                    "chr1",
                    500,
                    600,
                    "p_gt",
                    3,
                    "+",
                ),  # a.start (500) > b.start (490) → drop
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 0, "+"),
                ("chr1", 350, 450, "g2", 0, "+"),
                ("chr1", 490, 590, "g3", 0, "+"),
            ],
        )
        sql = transpile(
            "SELECT a.start AS a_s, b.start AS b_s FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "WHERE a.start < b.start",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert rows == [(100, 150), (300, 350)]

    def test_query_should_filter_by_between_predicate_when_dialect_is_duckdb(self, conn):
        """Test that ``BETWEEN`` filters overlap rows.

        Given:
            Three overlap pairs with varied ``a.score``; the predicate
            ``a.score BETWEEN 10 AND 50`` excludes the row at score=1
            and the row at score=99.
        When:
            The query is transpiled with ``dialect='duckdb'`` and
            executed.
        Then:
            Only the score=20 row should be returned.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p_low", 1, "+"),
                ("chr1", 300, 400, "p_mid", 20, "+"),
                ("chr1", 500, 600, "p_high", 99, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 0, "+"),
                ("chr1", 350, 450, "g2", 0, "+"),
                ("chr1", 550, 650, "g3", 0, "+"),
            ],
        )
        sql = transpile(
            "SELECT a.start AS s FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "AND a.score BETWEEN 10 AND 50",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert rows == [(300,)]

    def test_query_should_honor_is_null_predicate_when_dialect_is_duckdb(self, conn):
        """Test that ``WHERE a.score IS NULL`` filters overlap rows.

        Given:
            Two overlap pairs on chr1 — one with NULL ``a.score`` and
            one with a non-NULL score.
        When:
            The query is transpiled with ``WHERE a.score IS NULL`` and
            executed.
        Then:
            Only the NULL-score row should be returned.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p_null", None, "+"),
                ("chr1", 300, 400, "p_set", 99, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 0, "+"),
                ("chr1", 350, 450, "g2", 0, "+"),
            ],
        )
        sql = transpile(
            "SELECT a.start AS s FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "WHERE a.score IS NULL",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert rows == [(100,)]

    def test_query_should_rewrite_user_aliases_p_g_to_inner_a_b_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that user aliases other than ``a``/``b`` get rewritten at runtime.

        Given:
            A join written as ``FROM peaks p JOIN genes g`` (user
            aliases ``p`` and ``g``, not ``a`` and ``b``) with an extra
            cross-side ON predicate ``p.score > g.score``.
        When:
            The query is transpiled with ``dialect='duckdb'`` and
            executed.
        Then:
            DuckDB should not raise an alias-resolution error — the
            dialect rewrote ``p`` / ``g`` to the inner subquery's
            hardcoded ``a`` / ``b`` aliases — and the filter should
            return only the pair satisfying the cross-side predicate.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 100, "+"),
                ("chr1", 300, 400, "p2", 5, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g_low", 10, "+"),
                ("chr1", 350, 450, "g_high", 50, "+"),
            ],
        )
        sql = transpile(
            "SELECT p.start AS p_s, g.start AS g_s FROM peaks p "
            "JOIN genes g ON p.interval INTERSECTS g.interval "
            "AND p.score > g.score",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert rows == [(100, 150)]

    def test_query_should_count_non_null_column_values_per_group_when_dialect_is_duckdb(
        self, conn
    ):
        """Test that ``COUNT(a.name)`` counts only non-NULL names per group.

        Given:
            Three overlap rows on chr1 — two with non-NULL ``a.name``
            and one with NULL — plus one overlap on chr2 with non-NULL
            ``a.name``.
        When:
            The query is transpiled with ``COUNT(a.name) GROUP BY
            a.chrom`` and executed.
        Then:
            The chr1 count should be 2 (excluding the NULL), and the
            chr2 count should be 1.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),
                ("chr1", 110, 210, None, 2, "+"),
                ("chr1", 120, 220, "p3", 3, "+"),
                ("chr2", 100, 200, "p4", 4, "-"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 1, "+"),
                ("chr2", 150, 250, "g2", 2, "-"),
            ],
        )
        sql = transpile(
            "SELECT a.chrom AS c, COUNT(a.name) AS n FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "GROUP BY a.chrom",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert rows == [("chr1", 2), ("chr2", 1)]

    def test_query_should_group_by_multiple_columns_when_dialect_is_duckdb(self, conn):
        """Test that multi-column GROUP BY groups by tuple.

        Given:
            Overlap rows whose ``(a.chrom, a.strand)`` tuples are
            heterogeneous: two ``(chr1, +)`` overlaps, one ``(chr1, -)``
            overlap.
        When:
            The query is transpiled with ``GROUP BY a.chrom, a.strand``
            and executed.
        Then:
            One row per distinct ``(chrom, strand)`` tuple should be
            returned.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),
                ("chr1", 110, 210, "p2", 2, "+"),
                ("chr1", 300, 400, "p3", 3, "-"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 1, "+"),
                ("chr1", 350, 450, "g2", 2, "+"),
            ],
        )
        sql = transpile(
            "SELECT a.chrom AS c, a.strand AS s, COUNT(*) AS n FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "GROUP BY a.chrom, a.strand",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert rows == [("chr1", "+", 2), ("chr1", "-", 1)]

    def test_query_should_filter_groups_when_having_combines_two_conditions(self, conn):
        """Test that ``HAVING`` with two AND-ed aggregate conditions filters groups.

        Given:
            Three chromosomes; chr1 has many low-score overlaps, chr2
            has one high-score overlap, chr3 has many high-score
            overlaps. Only chr3 satisfies both
            ``COUNT(*) > 1 AND SUM(a.score) > 100``.
        When:
            The query is transpiled with ``dialect='duckdb'`` and
            executed.
        Then:
            Only chr3 should be returned.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                # chr1: 3 overlaps, all score 1 → COUNT > 1, SUM = 3
                ("chr1", 100, 200, "p1", 1, "+"),
                ("chr1", 110, 210, "p2", 1, "+"),
                ("chr1", 120, 220, "p3", 1, "+"),
                # chr2: 1 overlap, score 99 → COUNT == 1 (fails)
                ("chr2", 100, 200, "p4", 99, "+"),
                # chr3: 2 overlaps, scores 60 each → COUNT > 1, SUM = 120
                ("chr3", 100, 200, "p5", 60, "+"),
                ("chr3", 110, 210, "p6", 60, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 0, "+"),
                ("chr2", 150, 250, "g2", 0, "+"),
                ("chr3", 150, 250, "g3", 0, "+"),
            ],
        )
        sql = transpile(
            "SELECT a.chrom AS c, COUNT(*) AS n, SUM(a.score) AS s "
            "FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval "
            "GROUP BY a.chrom HAVING COUNT(*) > 1 AND SUM(a.score) > 100",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        assert rows == [("chr3", 2, 120)]

    def test_query_should_compose_all_tier1_features_in_single_pipeline_when_dialect_is_duckdb(
        self, conn
    ):
        """Test the kitchen-sink composition of every Tier 1 feature.

        Given:
            A query combining an extra ON predicate (``a.score > 20``),
            an extra WHERE predicate (``b.score < 50``), multi-column
            ``GROUP BY``, a non-trivial ``HAVING`` (``SUM(a.score) >
            40`` — actually filters one of the groups out), ``ORDER BY
            n DESC``, and ``LIMIT 10`` — over data crafted so the
            counts after both filters are 2/1/0 on chr1/chr2/chr3 but
            chr2's surviving sum is below the HAVING threshold.
        When:
            The query is transpiled with ``dialect='duckdb'`` and
            executed.
        Then:
            Only chr1 survives — its 2-overlap SUM(a.score) of 70 beats
            the threshold; chr2's single 30 fails it; chr3 was already
            dropped by the ON predicate. LIMIT 10 has no effect (one
            group survives). Aggregates SUM and COUNT(*) coexist and
            the ORDER BY references the user-aliased COUNT alias.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                # chr1: two peaks both pass a.score > 20; sum = 70
                ("chr1", 100, 200, "p1", 30, "+"),
                ("chr1", 110, 210, "p2", 40, "+"),
                # chr2: one peak passes a.score > 20; sum = 30 (fails HAVING)
                ("chr2", 100, 200, "p3", 30, "-"),
                # chr3: peak fails a.score > 20
                ("chr3", 100, 200, "p4", 10, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                # chr1 gene passes b.score < 50 (=10)
                ("chr1", 150, 250, "g1", 10, "+"),
                # chr2 gene passes b.score < 50 (=20)
                ("chr2", 150, 250, "g2", 20, "-"),
                # chr3 gene would pass but no overlap-survivor exists
                ("chr3", 150, 250, "g3", 5, "+"),
            ],
        )
        sql = transpile(
            "SELECT a.chrom AS c, COUNT(*) AS n, SUM(a.score) AS s "
            "FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval "
            "AND a.score > 20 "
            "WHERE b.score < 50 "
            "GROUP BY a.chrom HAVING SUM(a.score) > 40 "
            "ORDER BY n DESC LIMIT 10",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
        # chr1: 2 overlaps, sum(a.score)=70 → passes HAVING
        # chr2: 1 overlap, sum(a.score)=30 → fails HAVING
        # chr3: 0 overlaps (peak filtered by ON) → no group
        assert rows == [("chr1", 2, 70)]

    @pytest.mark.parametrize(
        ("coordinate_system", "interval_type"),
        [
            ("0based", "half_open"),
            ("0based", "closed"),
            ("1based", "half_open"),
            ("1based", "closed"),
        ],
    )
    def test_query_should_match_naive_predicate_plan_for_coordinate_system_combinations(
        self, conn, coordinate_system, interval_type
    ):
        """Test cross-plan equivalence across all 4 coord-system × interval-type combos.

        Given:
            A fixed 6-row truth table of peaks/genes with overlaps that
            depend on how endpoints are interpreted.
        When:
            The same logical query is transpiled twice — once with
            ``dialect=None`` (naive-predicate plan) and once with
            ``dialect="duckdb"`` — under each of the 4
            ``(coordinate_system, interval_type)`` combinations.
        Then:
            Both plans should return the same sorted multiset for every
            combination, proving canonicalization composes with the
            dialect rewrite consistently.
        """
        # Arrange — fixture uses strict (non-touching) overlaps so the
        # overlap set is identical across all 4 coord-system × interval
        # combos. The test then asserts that canonicalization composes
        # with the dialect rewrite without changing which pairs match.
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 0, "+"),
                ("chr1", 300, 400, "p2", 0, "+"),
                ("chr1", 700, 800, "p3", 0, "+"),
                ("chr2", 100, 200, "p4", 0, "-"),
                ("chr2", 500, 600, "p5", 0, "-"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 0, "+"),
                ("chr1", 350, 450, "g2", 0, "+"),
                ("chr2", 50, 180, "g3", 0, "-"),
                ("chr2", 550, 700, "g4", 0, "-"),
            ],
        )
        tables = [
            Table(
                "peaks",
                coordinate_system=coordinate_system,
                interval_type=interval_type,
            ),
            Table(
                "genes",
                coordinate_system=coordinate_system,
                interval_type=interval_type,
            ),
        ]
        query = (
            "SELECT a.chrom AS c, a.start AS s "
            "FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval"
        )
        sql_default = transpile(query, tables=tables)
        sql_duckdb = transpile(query, tables=tables, dialect="duckdb")
        # Python-native ground truth. The fixture uses strict (non-
        # touching) overlaps, so the matching set is identical across
        # all 4 coord-system × interval combos — that's the property
        # the test is pinning. Asserting against this expected list
        # (not just `rows_default == rows_duckdb`) catches a regression
        # in which both plans share the same bug.
        expected = [
            ("chr1", 100),  # peaks[100,200) ∩ genes[150,250)
            ("chr1", 300),  # peaks[300,400) ∩ genes[350,450)
            ("chr2", 100),  # peaks[100,200) ∩ genes[50,180)
            ("chr2", 500),  # peaks[500,600) ∩ genes[550,700)
        ]

        # Act
        rows_default = sorted(conn.execute(sql_default).fetchall())
        rows_duckdb = sorted(conn.execute(sql_duckdb).fetchall())

        # Assert
        assert rows_duckdb == expected
        assert rows_default == rows_duckdb

    @settings(
        max_examples=25,
        deadline=None,
        suppress_health_check=[HealthCheck.function_scoped_fixture],
    )
    @given(
        peaks=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2", "chr3"]),
                st.integers(min_value=0, max_value=200),
                st.integers(min_value=1, max_value=50),
            ),
            min_size=0,
            max_size=8,
        ),
        genes=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2", "chr3"]),
                st.integers(min_value=0, max_value=200),
                st.integers(min_value=1, max_value=50),
            ),
            min_size=0,
            max_size=8,
        ),
    )
    def test_query_should_match_python_native_overlap_for_random_inputs(
        self, conn, peaks, genes
    ):
        """Test the IEJoin path against a Python-native overlap reference.

        Given:
            A Hypothesis-generated pair of small interval lists drawn
            from a tiny chromosome alphabet and small int starts/lengths
            (half-open coordinates on both sides).
        When:
            The DuckDB-dialect SQL is executed and compared against a
            Python-native overlap reference
            (``a.end > b.start AND b.end > a.start AND a.chrom == b.chrom``).
        Then:
            The set of overlapping ``(chrom, start, end, chrom, start,
            end)`` tuples should be identical between the two methods.
        """
        # Arrange
        peak_rows = [(c, s, s + length) for (c, s, length) in peaks]
        gene_rows = [(c, s, s + length) for (c, s, length) in genes]
        conn.execute("DROP TABLE IF EXISTS peaks")
        conn.execute("DROP TABLE IF EXISTS genes")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        if peak_rows:
            conn.executemany("INSERT INTO peaks VALUES (?, ?, ?)", peak_rows)
        conn.execute(
            'CREATE TABLE genes (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        if gene_rows:
            conn.executemany("INSERT INTO genes VALUES (?, ?, ?)", gene_rows)
        sql = transpile(
            """
            SELECT a.chrom AS a_chrom, a.start AS a_start, a.end AS a_end,
                   b.chrom AS b_chrom, b.start AS b_start, b.end AS b_end
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        sql_rows = sorted(conn.execute(sql).fetchall())
        expected = _python_overlap(peak_rows, gene_rows)

        # Assert: sorted-multiset equality — a `set` here would collapse
        # duplicate input rows and silently hide multiplicity bugs.
        assert sql_rows == expected

    @settings(
        max_examples=25,
        deadline=None,
        suppress_health_check=[HealthCheck.function_scoped_fixture],
    )
    @given(
        peaks=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2"]),
                st.integers(min_value=0, max_value=100),
                st.integers(min_value=1, max_value=30),
            ),
            min_size=0,
            max_size=6,
        ),
        genes=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2"]),
                st.integers(min_value=0, max_value=100),
                st.integers(min_value=1, max_value=30),
            ),
            min_size=0,
            max_size=6,
        ),
    )
    def test_query_should_match_naive_predicate_plan_for_random_inputs(
        self, conn, peaks, genes
    ):
        """Test that the naive-predicate plan and the IEJoin path agree on random inputs.

        Given:
            A Hypothesis-generated pair of small interval lists.
        When:
            The same logical query is transpiled both with
            ``dialect=None`` and with ``dialect='duckdb'`` and executed.
        Then:
            Both plans should return the same multiset of overlapping
            rows.
        """
        # Arrange
        peak_rows = [(c, s, s + length) for (c, s, length) in peaks]
        gene_rows = [(c, s, s + length) for (c, s, length) in genes]
        conn.execute("DROP TABLE IF EXISTS peaks")
        conn.execute("DROP TABLE IF EXISTS genes")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        if peak_rows:
            conn.executemany("INSERT INTO peaks VALUES (?, ?, ?)", peak_rows)
        conn.execute(
            'CREATE TABLE genes (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        if gene_rows:
            conn.executemany("INSERT INTO genes VALUES (?, ?, ?)", gene_rows)
        query = """
            SELECT a.chrom AS a_chrom, a.start AS a_start, a.end AS a_end,
                   b.chrom AS b_chrom, b.start AS b_start, b.end AS b_end
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """
        sql_default = transpile(query, tables=["peaks", "genes"])
        sql_duckdb = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Act
        rows_default = sorted(conn.execute(sql_default).fetchall())
        rows_duckdb = sorted(conn.execute(sql_duckdb).fetchall())

        # Assert
        assert rows_default == rows_duckdb

    @settings(
        max_examples=30,
        deadline=None,
        suppress_health_check=[HealthCheck.function_scoped_fixture],
    )
    @given(
        peaks=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2"]),
                st.integers(min_value=0, max_value=100),
                st.integers(min_value=1, max_value=30),
                st.integers(min_value=0, max_value=50),
            ),
            min_size=1,
            max_size=6,
        ),
        genes=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2"]),
                st.integers(min_value=0, max_value=100),
                st.integers(min_value=1, max_value=30),
                st.integers(min_value=0, max_value=50),
            ),
            min_size=1,
            max_size=6,
        ),
        score_threshold=st.integers(min_value=0, max_value=50),
    )
    def test_query_should_match_naive_predicate_plan_when_where_group_order_combined(
        self, conn, peaks, genes, score_threshold
    ):
        """Test cross-plan equivalence with extras, GROUP BY, and ORDER BY.

        Given:
            Hypothesis-generated peak / gene interval lists carrying
            scores, and a random ``score_threshold``.
        When:
            The same logical query — combining a WHERE extra predicate,
            a GROUP BY, a COUNT(*) aggregate, and ORDER BY — is
            transpiled both with ``dialect=None`` and ``dialect='duckdb'``
            and executed.
        Then:
            Both plans should return the same multiset of grouped rows.
        """
        # Arrange
        peak_rows = [(c, s, s + length, score) for (c, s, length, score) in peaks]
        gene_rows = [(c, s, s + length, score) for (c, s, length, score) in genes]
        conn.execute("DROP TABLE IF EXISTS peaks")
        conn.execute("DROP TABLE IF EXISTS genes")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, '
            '"end" INTEGER, score INTEGER)'
        )
        conn.executemany("INSERT INTO peaks VALUES (?, ?, ?, ?)", peak_rows)
        conn.execute(
            'CREATE TABLE genes (chrom VARCHAR, "start" INTEGER, '
            '"end" INTEGER, score INTEGER)'
        )
        conn.executemany("INSERT INTO genes VALUES (?, ?, ?, ?)", gene_rows)
        query = (
            "SELECT a.chrom AS c, COUNT(*) AS n "
            "FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval "
            f"WHERE a.score >= {score_threshold} "
            "GROUP BY a.chrom "
            "ORDER BY a.chrom"
        )
        sql_default = transpile(query, tables=["peaks", "genes"])
        sql_duckdb = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Act
        rows_default = conn.execute(sql_default).fetchall()
        rows_duckdb = conn.execute(sql_duckdb).fetchall()

        # Assert
        assert rows_default == rows_duckdb

    @settings(
        max_examples=25,
        deadline=None,
        suppress_health_check=[HealthCheck.function_scoped_fixture],
    )
    @given(
        chrom_name=st.text(
            alphabet=st.characters(
                min_codepoint=0x20,
                max_codepoint=0x7E,
                blacklist_characters=["\x00", "%", "?"],
            ),
            min_size=1,
            max_size=10,
        ),
    )
    def test_query_should_round_trip_arbitrary_chrom_names_under_dialect(
        self, conn, chrom_name
    ):
        """Test that arbitrary printable chrom names round-trip the IEJoin path.

        Given:
            A Hypothesis-drawn printable-ASCII chromosome name (length
            1..10, excluding NUL, ``%``, ``?``) that may include single
            quotes, double quotes, or backslashes, seeded as the chrom
            value of a single overlapping pair.
        When:
            The DuckDB-dialect SQL is generated and executed.
        Then:
            The pair should be returned exactly once and the chrom value
            should round-trip byte-for-byte.
        """
        # Arrange
        conn.execute("DROP TABLE IF EXISTS peaks")
        conn.execute("DROP TABLE IF EXISTS genes")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.execute(
            'CREATE TABLE genes (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        conn.execute("INSERT INTO peaks VALUES (?, ?, ?)", (chrom_name, 100, 200))
        conn.execute("INSERT INTO genes VALUES (?, ?, ?)", (chrom_name, 150, 250))
        sql = transpile(
            """
            SELECT a.chrom AS a_chrom, a.start AS a_start
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """,
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = conn.execute(sql).fetchall()

        # Assert
        assert rows == [(chrom_name, 100)]

    @settings(
        max_examples=30,
        deadline=None,
        suppress_health_check=[HealthCheck.function_scoped_fixture],
    )
    @given(
        peaks=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2"]),
                st.integers(min_value=0, max_value=100),
                st.integers(min_value=1, max_value=30),
                # Nullable score so IS NULL / IS NOT NULL aren't
                # tautological against an always-populated column.
                st.one_of(st.none(), st.integers(min_value=0, max_value=50)),
            ),
            min_size=1,
            max_size=5,
        ),
        genes=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2"]),
                st.integers(min_value=0, max_value=100),
                st.integers(min_value=1, max_value=30),
                st.one_of(st.none(), st.integers(min_value=0, max_value=50)),
            ),
            min_size=1,
            max_size=5,
        ),
        predicate_kind=st.sampled_from(
            ["cmp", "between", "in_list", "cross_side", "is_null"]
        ),
        threshold=st.integers(min_value=0, max_value=50),
    )
    def test_query_should_match_naive_predicate_plan_for_random_extra_where_predicates(
        self, conn, peaks, genes, predicate_kind, threshold
    ):
        """Test cross-plan equivalence under random WHERE extras.

        Given:
            Hypothesis-generated peak / gene interval lists with
            (possibly NULL) score columns, plus a sampled predicate
            kind from ``{cmp, between, in_list, cross_side, is_null}``.
        When:
            The composed query is transpiled under both ``dialect=None``
            and ``dialect="duckdb"`` and executed.
        Then:
            Both plans should return the same sorted multiset.
        """
        # Arrange
        peak_rows = [(c, s, s + length, score) for (c, s, length, score) in peaks]
        gene_rows = [(c, s, s + length, score) for (c, s, length, score) in genes]
        conn.execute("DROP TABLE IF EXISTS peaks")
        conn.execute("DROP TABLE IF EXISTS genes")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, '
            '"end" INTEGER, score INTEGER)'
        )
        conn.executemany("INSERT INTO peaks VALUES (?, ?, ?, ?)", peak_rows)
        conn.execute(
            'CREATE TABLE genes (chrom VARCHAR, "start" INTEGER, '
            '"end" INTEGER, score INTEGER)'
        )
        conn.executemany("INSERT INTO genes VALUES (?, ?, ?, ?)", gene_rows)
        where_clauses = {
            "cmp": f"WHERE a.score >= {threshold}",
            "between": f"WHERE a.score BETWEEN {threshold} AND {threshold + 20}",
            "in_list": "WHERE a.score IN (0, 10, 20, 30, 40, 50)",
            "cross_side": "WHERE a.score >= b.score",
            "is_null": "WHERE a.score IS NOT NULL",
        }
        where = where_clauses[predicate_kind]
        query = (
            "SELECT a.chrom AS c, a.start AS s "
            "FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval "
            f"{where}"
        )
        sql_default = transpile(query, tables=["peaks", "genes"])
        sql_duckdb = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Act
        rows_default = sorted(conn.execute(sql_default).fetchall())
        rows_duckdb = sorted(conn.execute(sql_duckdb).fetchall())

        # Assert
        assert rows_default == rows_duckdb

    @settings(
        max_examples=30,
        deadline=None,
        suppress_health_check=[HealthCheck.function_scoped_fixture],
    )
    @given(
        peaks=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2"]),
                st.integers(min_value=0, max_value=100),
                st.integers(min_value=1, max_value=30),
                st.integers(min_value=0, max_value=50),
            ),
            min_size=1,
            max_size=5,
        ),
        genes=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2"]),
                st.integers(min_value=0, max_value=100),
                st.integers(min_value=1, max_value=30),
                st.integers(min_value=0, max_value=50),
            ),
            min_size=1,
            max_size=5,
        ),
        aggregate_kind=st.sampled_from(
            [
                "count_star",
                "sum_a_score",
                "min_b_start",
                "max_b_end",
                "count_distinct_b_chrom",
            ]
        ),
    )
    def test_query_should_match_naive_predicate_plan_for_random_group_by_aggregates(
        self, conn, peaks, genes, aggregate_kind
    ):
        """Test cross-plan equivalence for GROUP BY + random aggregate.

        Given:
            Hypothesis-generated peak / gene interval lists with score
            columns, plus a sampled aggregate from a small whitelist.
        When:
            A query with ``GROUP BY a.chrom`` and the chosen aggregate
            is transpiled under both ``dialect=None`` and
            ``dialect="duckdb"`` and executed.
        Then:
            Both plans should return the same row list (sorted by chrom
            for stable comparison; ``AVG`` excluded from this PBT to
            keep equality strict — covered by an example test instead).
        """
        # Arrange
        peak_rows = [(c, s, s + length, score) for (c, s, length, score) in peaks]
        gene_rows = [(c, s, s + length, score) for (c, s, length, score) in genes]
        conn.execute("DROP TABLE IF EXISTS peaks")
        conn.execute("DROP TABLE IF EXISTS genes")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, '
            '"end" INTEGER, score INTEGER)'
        )
        conn.executemany("INSERT INTO peaks VALUES (?, ?, ?, ?)", peak_rows)
        conn.execute(
            'CREATE TABLE genes (chrom VARCHAR, "start" INTEGER, '
            '"end" INTEGER, score INTEGER)'
        )
        conn.executemany("INSERT INTO genes VALUES (?, ?, ?, ?)", gene_rows)
        aggregates = {
            "count_star": "COUNT(*)",
            "sum_a_score": "SUM(a.score)",
            "min_b_start": "MIN(b.start)",
            "max_b_end": "MAX(b.end)",
            "count_distinct_b_chrom": "COUNT(DISTINCT b.chrom)",
        }
        agg = aggregates[aggregate_kind]
        query = (
            f"SELECT a.chrom AS c, {agg} AS v "
            "FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval "
            "GROUP BY a.chrom ORDER BY a.chrom"
        )
        sql_default = transpile(query, tables=["peaks", "genes"])
        sql_duckdb = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Act
        rows_default = conn.execute(sql_default).fetchall()
        rows_duckdb = conn.execute(sql_duckdb).fetchall()

        # Assert — ORDER BY makes positional comparison meaningful.
        assert rows_default == rows_duckdb

    @settings(
        # Full-composite has the largest input space (extras + GROUP BY
        # + HAVING + ORDER BY + LIMIT + OFFSET), so we run more examples.
        max_examples=60,
        deadline=None,
        suppress_health_check=[HealthCheck.function_scoped_fixture],
    )
    @given(
        # Allow ``min_size=0`` on peaks so empty-input cross-plan
        # equivalence is exercised in combination with the full Tier 1
        # stack (the simpler PBTs use min_size=1).
        peaks=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2"]),
                st.integers(min_value=0, max_value=100),
                st.integers(min_value=1, max_value=30),
                st.integers(min_value=0, max_value=50),
            ),
            min_size=0,
            max_size=6,
        ),
        genes=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2"]),
                st.integers(min_value=0, max_value=100),
                st.integers(min_value=1, max_value=30),
                st.integers(min_value=0, max_value=50),
            ),
            min_size=0,
            max_size=6,
        ),
        score_threshold=st.integers(min_value=0, max_value=50),
        having_min=st.integers(min_value=0, max_value=5),
        # Include LIMIT 0 so the trivial-truncation edge case is
        # exercised against both plans.
        limit=st.integers(min_value=0, max_value=10),
    )
    def test_query_should_match_naive_predicate_plan_for_full_composite_query(
        self, conn, peaks, genes, score_threshold, having_min, limit
    ):
        """Test cross-plan equivalence under the full Tier 1 composition.

        Given:
            Hypothesis-generated peaks/genes with random thresholds for
            an extra WHERE predicate, a HAVING minimum count, and a
            LIMIT.
        When:
            A query exercising extra ON + extra WHERE + GROUP BY +
            HAVING + ORDER BY + LIMIT is transpiled under both
            ``dialect=None`` and ``dialect="duckdb"`` and executed.
        Then:
            Both plans should return the same ordered row list (``==``,
            because ``ORDER BY`` makes order part of the contract).
        """
        # Arrange
        peak_rows = [(c, s, s + length, score) for (c, s, length, score) in peaks]
        gene_rows = [(c, s, s + length, score) for (c, s, length, score) in genes]
        conn.execute("DROP TABLE IF EXISTS peaks")
        conn.execute("DROP TABLE IF EXISTS genes")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, '
            '"end" INTEGER, score INTEGER)'
        )
        # DuckDB's executemany rejects empty parameter lists; skip the
        # INSERT step when the strategy generates no rows.
        if peak_rows:
            conn.executemany("INSERT INTO peaks VALUES (?, ?, ?, ?)", peak_rows)
        conn.execute(
            'CREATE TABLE genes (chrom VARCHAR, "start" INTEGER, '
            '"end" INTEGER, score INTEGER)'
        )
        if gene_rows:
            conn.executemany("INSERT INTO genes VALUES (?, ?, ?, ?)", gene_rows)
        query = (
            "SELECT a.chrom AS c, COUNT(*) AS n "
            "FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval AND a.score >= b.score "
            f"WHERE a.score >= {score_threshold} "
            "GROUP BY a.chrom "
            f"HAVING COUNT(*) >= {having_min} "
            "ORDER BY a.chrom "
            f"LIMIT {limit}"
        )
        sql_default = transpile(query, tables=["peaks", "genes"])
        sql_duckdb = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Act
        rows_default = conn.execute(sql_default).fetchall()
        rows_duckdb = conn.execute(sql_duckdb).fetchall()

        # Assert
        assert rows_default == rows_duckdb

    def test_query_should_match_inner_join_result_when_using_join_targets_chrom_only(
        self, peaks_genes
    ):
        """Test that USING(chrom) returns the same rows as plain INNER JOIN.

        Given:
            The shared peaks_genes fixture and two queries — one using
            ``JOIN ... ON a.interval INTERSECTS b.interval`` and another
            using ``JOIN ... USING (chrom) WHERE a.interval INTERSECTS
            b.interval`` — both with ``dialect='duckdb'``.
        When:
            Both are executed.
        Then:
            They should return identical sorted row sets, because
            USING(chrom) is redundant with the dialect's per-chromosome
            partition.
        """
        # Arrange
        plain_sql = transpile(
            "SELECT a.chrom, a.start, a.end "
            "FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )
        using_sql = transpile(
            "SELECT a.chrom, a.start, a.end "
            "FROM peaks a JOIN genes b USING (chrom) "
            "WHERE a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        plain_rows = sorted(peaks_genes.execute(plain_sql).fetchall())
        using_rows = sorted(peaks_genes.execute(using_sql).fetchall())

        # Assert
        assert plain_rows == using_rows

    def test_query_should_match_naive_predicate_plan_for_mixed_case_aliases(self, conn):
        """Test that mixed-case aliases under the dialect match the naive-predicate plan.

        Given:
            The shared peaks_genes fixture loaded onto ``conn`` and a
            query using uppercase aliases (``peaks A`` / ``genes B``).
        When:
            The same logical query is transpiled both with
            ``dialect=None`` and ``dialect='duckdb'`` and executed.
        Then:
            Both plans should return the same multiset of rows.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 10, "+"),
                ("chr1", 500, 600, "p2", 20, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 1, "+"),
                ("chr1", 500, 600, "g2", 2, "+"),
            ],
        )
        query = (
            "SELECT A.chrom, A.start, B.start "
            "FROM peaks A JOIN genes B "
            "ON A.interval INTERSECTS B.interval"
        )
        sql_default = transpile(query, tables=["peaks", "genes"])
        sql_duckdb = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Act
        rows_default = sorted(conn.execute(sql_default).fetchall())
        rows_duckdb = sorted(conn.execute(sql_duckdb).fetchall())

        # Assert
        assert rows_default == rows_duckdb

    def test_query_should_preserve_rows_from_left_only_chromosomes_when_join_is_anti(
        self, conn
    ):
        """Test that ANTI JOIN preserves left-only chromosomes (C1b regression).

        Given:
            Peaks on chr1 and chr3; genes on chr1 only. Under ANTI JOIN
            with INTERSECTS, chr3 peaks must be returned (they have no
            overlapping right-side row), and chr1 peaks that don't
            overlap any gene must also be returned.
        When:
            The query is transpiled with ``dialect='duckdb'`` and
            executed.
        Then:
            It should return the chr3 peak and any chr1 peaks with no
            overlap — the chrom-INTERSECT partition would have excluded
            chr3 entirely, so this is the critical regression for
            C1b's left-distinct partition swap.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 10, 20, "p1", 0, "+"),  # no overlap with chr1 gene
                ("chr1", 100, 200, "p2", 0, "+"),  # overlaps chr1 gene
                ("chr3", 1, 1000, "p3", 0, "+"),  # chr3 — no genes there
            ],
        )
        _make_table(
            conn,
            "genes",
            [("chr1", 50, 150, "g1", 0, "+")],
        )
        sql = transpile(
            "SELECT a.chrom, a.start, a.end "
            "FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        rows = sorted(conn.execute(sql).fetchall())

        # Assert
        assert rows == sorted([("chr1", 10, 20), ("chr3", 1, 1000)])

    @settings(
        max_examples=25,
        deadline=None,
        suppress_health_check=[HealthCheck.function_scoped_fixture],
    )
    @given(
        peaks=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2", "chr3"]),
                st.integers(min_value=0, max_value=200),
                st.integers(min_value=1, max_value=50),
            ),
            min_size=0,
            max_size=8,
        ),
        genes=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2", "chr3"]),
                st.integers(min_value=0, max_value=200),
                st.integers(min_value=1, max_value=50),
            ),
            min_size=0,
            max_size=8,
        ),
    )
    def test_query_should_match_python_native_semi_overlap_for_random_inputs(
        self, conn, peaks, genes
    ):
        """Test SEMI JOIN against a Python-native semi-overlap reference.

        Given:
            A Hypothesis-generated pair of small interval lists drawn
            from a tiny chromosome alphabet.
        When:
            A ``SEMI JOIN ... ON a.interval INTERSECTS b.interval`` is
            transpiled with ``dialect='duckdb'`` and executed.
        Then:
            The set of returned ``(chrom, start, end)`` tuples should
            equal the Python reference (distinct left rows that have
            at least one overlapping right row).
        """
        # Arrange
        peak_rows = [(c, s, s + length) for (c, s, length) in peaks]
        gene_rows = [(c, s, s + length) for (c, s, length) in genes]
        conn.execute("DROP TABLE IF EXISTS peaks")
        conn.execute("DROP TABLE IF EXISTS genes")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        if peak_rows:
            conn.executemany("INSERT INTO peaks VALUES (?, ?, ?)", peak_rows)
        conn.execute(
            'CREATE TABLE genes (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        if gene_rows:
            conn.executemany("INSERT INTO genes VALUES (?, ?, ?)", gene_rows)
        sql = transpile(
            "SELECT a.chrom, a.start, a.end "
            "FROM peaks a "
            "SEMI JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        sql_rows = sorted(set(conn.execute(sql).fetchall()))
        expected = _python_semi_overlap(peak_rows, gene_rows)

        # Assert
        assert sql_rows == expected

    @settings(
        max_examples=25,
        deadline=None,
        suppress_health_check=[HealthCheck.function_scoped_fixture],
    )
    @given(
        peaks=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2"]),
                st.integers(min_value=0, max_value=100),
                st.integers(min_value=1, max_value=30),
            ),
            min_size=0,
            max_size=6,
        ),
        genes=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2"]),
                st.integers(min_value=0, max_value=100),
                st.integers(min_value=1, max_value=30),
            ),
            min_size=0,
            max_size=6,
        ),
    )
    def test_query_should_match_python_native_semi_overlap_for_narrow_random_inputs(
        self, conn, peaks, genes
    ):
        """Test SEMI JOIN under the dialect against a Python-native reference.

        Given:
            A Hypothesis-generated pair of small interval lists.
        When:
            A ``SEMI JOIN ... ON a.interval INTERSECTS b.interval`` is
            transpiled with ``dialect='duckdb'`` and executed.
        Then:
            The set of returned rows should equal the Python reference of
            distinct left rows with at least one overlapping right row.
        """
        # Arrange
        peak_rows = [(c, s, s + length) for (c, s, length) in peaks]
        gene_rows = [(c, s, s + length) for (c, s, length) in genes]
        conn.execute("DROP TABLE IF EXISTS peaks")
        conn.execute("DROP TABLE IF EXISTS genes")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        if peak_rows:
            conn.executemany("INSERT INTO peaks VALUES (?, ?, ?)", peak_rows)
        conn.execute(
            'CREATE TABLE genes (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        if gene_rows:
            conn.executemany("INSERT INTO genes VALUES (?, ?, ?)", gene_rows)
        query = (
            "SELECT a.chrom, a.start, a.end "
            "FROM peaks a "
            "SEMI JOIN genes b ON a.interval INTERSECTS b.interval"
        )
        sql_duckdb = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Act
        rows_duckdb = sorted(set(conn.execute(sql_duckdb).fetchall()))
        expected = _python_semi_overlap(peak_rows, gene_rows)

        # Assert
        assert rows_duckdb == expected

    @settings(
        max_examples=25,
        deadline=None,
        suppress_health_check=[HealthCheck.function_scoped_fixture],
    )
    @given(
        peaks=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2", "chr3"]),
                st.integers(min_value=0, max_value=200),
                st.integers(min_value=1, max_value=50),
            ),
            min_size=0,
            max_size=8,
        ),
        genes=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2", "chr3"]),
                st.integers(min_value=0, max_value=200),
                st.integers(min_value=1, max_value=50),
            ),
            min_size=0,
            max_size=8,
        ),
    )
    def test_query_should_match_python_native_anti_overlap_for_random_inputs(
        self, conn, peaks, genes
    ):
        """Test ANTI JOIN against a Python-native anti-overlap reference.

        Given:
            A Hypothesis-generated pair of small interval lists.
        When:
            An ``ANTI JOIN ... ON a.interval INTERSECTS b.interval`` is
            transpiled with ``dialect='duckdb'`` and executed.
        Then:
            The set of returned ``(chrom, start, end)`` tuples should
            equal the Python reference (distinct left rows with no
            overlapping right row, including rows on chromosomes
            absent from the right table).
        """
        # Arrange
        peak_rows = [(c, s, s + length) for (c, s, length) in peaks]
        gene_rows = [(c, s, s + length) for (c, s, length) in genes]
        conn.execute("DROP TABLE IF EXISTS peaks")
        conn.execute("DROP TABLE IF EXISTS genes")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        if peak_rows:
            conn.executemany("INSERT INTO peaks VALUES (?, ?, ?)", peak_rows)
        conn.execute(
            'CREATE TABLE genes (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        if gene_rows:
            conn.executemany("INSERT INTO genes VALUES (?, ?, ?)", gene_rows)
        sql = transpile(
            "SELECT a.chrom, a.start, a.end "
            "FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Act
        sql_rows = sorted(set(conn.execute(sql).fetchall()))
        expected = _python_anti_overlap(peak_rows, gene_rows)

        # Assert
        assert sql_rows == expected

    @settings(
        max_examples=25,
        deadline=None,
        suppress_health_check=[HealthCheck.function_scoped_fixture],
    )
    @given(
        peaks=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2"]),
                st.integers(min_value=0, max_value=100),
                st.integers(min_value=1, max_value=30),
            ),
            min_size=0,
            max_size=6,
        ),
        genes=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2"]),
                st.integers(min_value=0, max_value=100),
                st.integers(min_value=1, max_value=30),
            ),
            min_size=0,
            max_size=6,
        ),
    )
    def test_query_should_match_python_native_anti_overlap_for_narrow_random_inputs(
        self, conn, peaks, genes
    ):
        """Test ANTI JOIN under the dialect against a Python-native reference.

        Given:
            A Hypothesis-generated pair of small interval lists.
        When:
            An ``ANTI JOIN ... ON a.interval INTERSECTS b.interval`` is
            transpiled with ``dialect='duckdb'`` and executed.
        Then:
            The set of returned rows should equal the Python reference of
            distinct left rows with no overlapping right row, including
            rows on chromosomes absent from the right table and every left
            row when the right table is empty.
        """
        # Arrange
        peak_rows = [(c, s, s + length) for (c, s, length) in peaks]
        gene_rows = [(c, s, s + length) for (c, s, length) in genes]
        conn.execute("DROP TABLE IF EXISTS peaks")
        conn.execute("DROP TABLE IF EXISTS genes")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        if peak_rows:
            conn.executemany("INSERT INTO peaks VALUES (?, ?, ?)", peak_rows)
        conn.execute(
            'CREATE TABLE genes (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        if gene_rows:
            conn.executemany("INSERT INTO genes VALUES (?, ?, ?)", gene_rows)
        query = (
            "SELECT a.chrom, a.start, a.end "
            "FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval"
        )
        sql_duckdb = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Act
        rows_duckdb = sorted(set(conn.execute(sql_duckdb).fetchall()))
        expected = _python_anti_overlap(peak_rows, gene_rows)

        # Assert
        assert rows_duckdb == expected

    def test_query_should_match_naive_plan_for_anti_join_with_where_residual(self, conn):
        """Test that an ANTI-join WHERE residual filters without inverting the join.

        Given:
            Three left rows — two overlap a right row (excluded by the
            anti-join), one does not — and a top-level ``WHERE a.start >= 550``
            post-join filter.
        When:
            The query is transpiled with ``dialect='duckdb'`` and
            ``dialect=None`` and executed.
        Then:
            Both plans return only the non-overlapping row that also passes
            the residual, proving the filter is applied after the anti-join
            rather than folded into its per-chromosome ON (#200); the buggy
            inlining would resurrect the two overlapping rows below 550.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),  # overlaps g1; anti-excluded
                ("chr1", 500, 600, "p2", 2, "+"),  # overlaps g2; anti-excluded
                ("chr1", 700, 800, "p3", 3, "+"),  # no overlap; passes residual
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 1, "+"),
                ("chr1", 500, 600, "g2", 2, "+"),
            ],
        )
        query = (
            "SELECT a.start FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval "
            "WHERE a.start >= 550"
        )
        sql_duckdb = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act
        rows_duckdb = sorted(conn.execute(sql_duckdb).fetchall())
        rows_naive = sorted(conn.execute(sql_naive).fetchall())

        # Assert
        assert rows_duckdb == [(700,)]
        assert rows_duckdb == rows_naive

    def test_query_should_filter_anti_join_by_where_residual_column_absent_from_select(
        self, conn
    ):
        """Test that an ANTI-join WHERE residual can filter on an unselected column.

        Given:
            An ``ANTI JOIN`` selecting only ``a.start`` but filtering on
            ``a.end`` — a column present only in the ``WHERE`` residual.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            The wrapper relation still exposes ``a.end`` so the outer filter
            resolves, and the result matches the naive-predicate plan.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),  # overlaps g1; anti-excluded
                ("chr1", 700, 800, "p3", 3, "+"),  # no overlap; a.end 800 > 750
            ],
        )
        _make_table(conn, "genes", [("chr1", 150, 250, "g1", 1, "+")])
        query = (
            "SELECT a.start FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval "
            "WHERE a.end > 750"
        )
        sql_duckdb = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act
        rows_duckdb = sorted(conn.execute(sql_duckdb).fetchall())
        rows_naive = sorted(conn.execute(sql_naive).fetchall())

        # Assert
        assert rows_duckdb == [(700,)]
        assert rows_duckdb == rows_naive

    def test_query_should_match_naive_plan_for_semi_join_with_where_residual(self, conn):
        """Test that a SEMI-join WHERE residual filters consistently with the naive plan.

        Given:
            Two overlapping left rows (kept by the semi-join) and a top-level
            ``WHERE a.start >= 300`` post-join filter that keeps one of them.
        When:
            The query is transpiled with ``dialect='duckdb'`` and
            ``dialect=None`` and executed.
        Then:
            Both plans return only the surviving overlapping row. SEMI's
            result is identical whether the residual is inlined or applied
            outside, so this locks the new SEMI outer-filter code path
            against regressions rather than proving a behavior change (#200).
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),  # overlaps g1; start < 300
                ("chr1", 500, 600, "p2", 2, "+"),  # overlaps g2; start >= 300
                ("chr1", 700, 800, "p3", 3, "+"),  # no overlap; semi-excluded
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 1, "+"),
                ("chr1", 500, 600, "g2", 2, "+"),
            ],
        )
        query = (
            "SELECT a.start FROM peaks a "
            "SEMI JOIN genes b ON a.interval INTERSECTS b.interval "
            "WHERE a.start >= 300"
        )
        sql_duckdb = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act
        rows_duckdb = sorted(conn.execute(sql_duckdb).fetchall())
        rows_naive = sorted(conn.execute(sql_naive).fetchall())

        # Assert
        assert rows_duckdb == [(500,)]
        assert rows_duckdb == rows_naive

    @settings(
        max_examples=25,
        deadline=None,
        suppress_health_check=[HealthCheck.function_scoped_fixture],
    )
    @given(
        peaks=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2"]),
                st.integers(min_value=0, max_value=100),
                st.integers(min_value=1, max_value=30),
            ),
            min_size=0,
            max_size=6,
        ),
        genes=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2"]),
                st.integers(min_value=0, max_value=100),
                st.integers(min_value=1, max_value=30),
            ),
            min_size=0,
            max_size=6,
        ),
    )
    def test_query_should_match_python_native_anti_overlap_with_where_residual_random_inputs(
        self, conn, peaks, genes
    ):
        """Test ANTI JOIN with a WHERE residual against a Python-native reference.

        Given:
            A Hypothesis-generated pair of small interval lists.
        When:
            An ``ANTI JOIN ... ON INTERSECTS ... WHERE a.start >= 50`` is
            transpiled with ``dialect='duckdb'`` and executed.
        Then:
            The set of returned rows should equal the Python anti-overlap
            reference filtered to ``start >= 50`` — the WHERE filter never
            resurrects an anti-excluded row on either side of the 50
            threshold (#200).
        """
        # Arrange
        peak_rows = [(c, s, s + length) for (c, s, length) in peaks]
        gene_rows = [(c, s, s + length) for (c, s, length) in genes]
        conn.execute("DROP TABLE IF EXISTS peaks")
        conn.execute("DROP TABLE IF EXISTS genes")
        conn.execute(
            'CREATE TABLE peaks (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        if peak_rows:
            conn.executemany("INSERT INTO peaks VALUES (?, ?, ?)", peak_rows)
        conn.execute(
            'CREATE TABLE genes (chrom VARCHAR, "start" INTEGER, "end" INTEGER)'
        )
        if gene_rows:
            conn.executemany("INSERT INTO genes VALUES (?, ?, ?)", gene_rows)
        query = (
            "SELECT a.chrom, a.start, a.end "
            "FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval "
            "WHERE a.start >= 50"
        )
        sql_duckdb = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Act
        rows_duckdb = sorted(set(conn.execute(sql_duckdb).fetchall()))
        expected = sorted(
            row for row in _python_anti_overlap(peak_rows, gene_rows) if row[1] >= 50
        )

        # Assert
        assert rows_duckdb == expected

    def test_query_should_inline_anti_join_on_residual_referencing_right_side(
        self, conn
    ):
        """Test that an ANTI-join ON residual on the right side stays inlined.

        Given:
            An ``ANTI JOIN ... ON a.interval INTERSECTS b.interval AND
            b.score > 100`` — a genuine join condition referencing the right
            side, which no gene satisfies.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            The residual stays inlined in the per-chromosome ON (it is a join
            condition, not a post-join filter), so every left row survives the
            anti-join and the result matches the naive plan. This guards the
            provenance split from over-correcting the ON path — only WHERE
            residuals move to the outer filter for the left-only shapes (#200).
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),
                ("chr1", 500, 600, "p2", 2, "+"),
                ("chr1", 700, 800, "p3", 3, "+"),
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 10, "+"),
                ("chr1", 500, 600, "g2", 50, "+"),
            ],
        )
        query = (
            "SELECT a.start FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval "
            "AND b.score > 100"
        )
        sql_duckdb = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act
        rows_duckdb = sorted(conn.execute(sql_duckdb).fetchall())
        rows_naive = sorted(conn.execute(sql_naive).fetchall())

        # Assert
        assert rows_duckdb == [(100,), (500,), (700,)]
        assert rows_duckdb == rows_naive

    def test_query_should_match_naive_plan_for_anti_join_with_compound_where_residual(
        self, conn
    ):
        """Test that a compound ANTI-join WHERE residual renders as one outer filter.

        Given:
            An ``ANTI JOIN`` with a two-clause post-join filter
            ``WHERE a.start >= 550 AND a.end < 900`` over two columns.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            Both clauses are applied on the outer wrapper relation and the
            result matches the naive plan (#200).
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),  # overlaps g1; anti-excluded
                ("chr1", 500, 600, "p2", 2, "+"),  # overlaps g2; anti-excluded
                ("chr1", 700, 800, "p3", 3, "+"),  # no overlap; passes both clauses
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 1, "+"),
                ("chr1", 500, 600, "g2", 2, "+"),
            ],
        )
        query = (
            "SELECT a.start FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval "
            "WHERE a.start >= 550 AND a.end < 900"
        )
        sql_duckdb = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act
        rows_duckdb = sorted(conn.execute(sql_duckdb).fetchall())
        rows_naive = sorted(conn.execute(sql_naive).fetchall())

        # Assert
        assert rows_duckdb == [(700,)]
        assert rows_duckdb == rows_naive

    def test_query_should_group_anti_join_with_where_residual(self, conn):
        """Test that an ANTI-join WHERE residual composes with GROUP BY on the wrapper.

        Given:
            An ``ANTI JOIN`` grouping by ``a.chrom`` with a post-join
            ``WHERE a.start >= 550`` filter — the outer filter and the
            aggregation both sit on the wrapper relation.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            The filter is applied before the aggregation and the grouped
            counts match the naive plan (#200).
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),  # overlaps g1; anti-excluded
                ("chr1", 500, 600, "p2", 2, "+"),  # overlaps g2; anti-excluded
                ("chr1", 700, 800, "p3", 3, "+"),  # no overlap; passes filter
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 1, "+"),
                ("chr1", 500, 600, "g2", 2, "+"),
            ],
        )
        query = (
            "SELECT a.chrom, COUNT(a.start) AS n FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval "
            "WHERE a.start >= 550 GROUP BY a.chrom"
        )
        sql_duckdb = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act
        rows_duckdb = sorted(conn.execute(sql_duckdb).fetchall())
        rows_naive = sorted(conn.execute(sql_naive).fetchall())

        # Assert
        assert rows_duckdb == [("chr1", 1)]
        assert rows_duckdb == rows_naive

    def test_query_should_order_anti_join_with_where_residual(self, conn):
        """Test that an ANTI-join WHERE residual composes with ORDER BY on the wrapper.

        Given:
            An ``ANTI JOIN`` with a post-join ``WHERE a.start >= 550`` filter and
            an order-sensitive ``ORDER BY a.start DESC`` — both the filter and the
            sort sit on the wrapper relation.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            The filtered rows come back in the requested descending order,
            matching the naive plan (rows compared without re-sorting so the
            order itself is asserted) (#200).
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [
                ("chr1", 100, 200, "p1", 1, "+"),  # overlaps g1; anti-excluded
                ("chr1", 500, 600, "p2", 2, "+"),  # overlaps g2; anti-excluded
                ("chr1", 700, 800, "p3", 3, "+"),  # no overlap; passes filter
                ("chr1", 900, 1000, "p4", 4, "+"),  # no overlap; passes filter
            ],
        )
        _make_table(
            conn,
            "genes",
            [
                ("chr1", 150, 250, "g1", 1, "+"),
                ("chr1", 500, 600, "g2", 2, "+"),
            ],
        )
        query = (
            "SELECT a.start FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval "
            "WHERE a.start >= 550 ORDER BY a.start DESC"
        )
        sql_duckdb = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act
        rows_duckdb = conn.execute(sql_duckdb).fetchall()
        rows_naive = conn.execute(sql_naive).fetchall()

        # Assert
        assert rows_duckdb == [(900,), (700,)]
        assert rows_duckdb == rows_naive


class TestTranspileDuckDBIEJoinKwargs:
    """Tests for the ``dialect`` kwarg surface.

    Covers the ``@overload`` set on :func:`giql.transpile.transpile`, the
    ``dialect=None``-equals-omitted equivalence, and the unknown-dialect error.
    """

    def test_transpile_should_match_default_when_dialect_none_is_explicit(self):
        """Test that ``dialect=None`` is equivalent to omitting the kwarg.

        Given:
            A column-to-column INTERSECTS JOIN.
        When:
            ``transpile`` is called once with ``dialect=None`` explicit
            and once with the kwarg omitted entirely.
        Then:
            Both calls should produce byte-identical SQL output.
        """
        # Arrange
        query = """
            SELECT a.chrom, a.start, b.start
            FROM peaks a
            JOIN genes b ON a.interval INTERSECTS b.interval
            """

        # Act
        sql_omitted = transpile(query, tables=["peaks", "genes"])
        sql_explicit_none = transpile(query, tables=["peaks", "genes"], dialect=None)

        # Assert
        assert sql_omitted == sql_explicit_none

    def test_transpile_should_echo_offending_value_in_unknown_dialect_error(self):
        """Test that the unknown-dialect error names the offending value.

        Given:
            An unknown dialect string ``'postgres'``.
        When:
            ``transpile`` is called.
        Then:
            It should raise ``ValueError`` whose message contains both
            ``Unknown dialect`` and the literal ``'postgres'``.
        """
        # Arrange
        query = "SELECT * FROM peaks"

        # Act & assert
        with pytest.raises(ValueError) as excinfo:
            transpile(query, tables=["peaks"], dialect="postgres")
        message = str(excinfo.value)
        assert "Unknown dialect" in message
        assert "'postgres'" in message


class TestTranspileDuckDBIEJoinLeftOnlyWhereIntersects:
    """Regression tests for #201.

    A ``SEMI`` / ``ANTI`` join whose column-to-column ``INTERSECTS`` lives in
    the top-level ``WHERE`` (out of scope for the right table after a left-only
    join) must decline the IEJoin rewrite so the naive-predicate plan surfaces
    the same binder error the reference plans raise — rather than relocating
    the predicate into the join and inventing anti/semi-overlap results.
    """

    @pytest.mark.parametrize("kind", ["ANTI", "SEMI"])
    @pytest.mark.parametrize(
        "on_clause", ["ON TRUE", "ON FALSE", "ON a.score > b.score"]
    )
    def test_transpile_should_decline_left_only_where_intersects_to_naive_plan(
        self, kind, on_clause
    ):
        """Test that a SEMI/ANTI join with a WHERE-INTERSECTS declines to naive.

        Given:
            A ``SEMI`` / ``ANTI`` join carrying an arbitrary ``ON`` clause and
            the column-to-column ``INTERSECTS`` in the top-level ``WHERE``.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should decline the IEJoin (no ``SET VARIABLE`` scaffolding) so
            the naive-predicate plan handles the out-of-scope WHERE reference,
            regardless of the ignored ``ON`` content.
        """
        # Arrange
        query = (
            f"SELECT a.name FROM peaks a {kind} JOIN genes b {on_clause} "
            "WHERE a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE" not in sql.upper()
        assert "getvariable" not in sql

    @pytest.mark.parametrize("kind", ["ANTI", "SEMI"])
    def test_query_should_raise_binder_error_for_where_intersects_like_naive_plan(
        self, conn, kind
    ):
        """Test that DuckDB rejects the WHERE-INTERSECTS shape as the naive plan does.

        Given:
            A ``SEMI`` / ``ANTI`` join with ``ON TRUE`` and the ``INTERSECTS``
            in the ``WHERE`` (right table out of scope), plus fixture data.
        When:
            The query is transpiled with ``dialect='duckdb'`` and with
            ``dialect=None`` and both are executed.
        Then:
            Both should raise a DuckDB binder error naming the out-of-scope
            right table — the dialect no longer diverges by silently inventing
            results.
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p1", 7, "+")])
        _make_table(conn, "genes", [("chr1", 150, 250, "g1", 9, "-")])
        query = (
            f"SELECT a.name FROM peaks a {kind} JOIN genes b ON TRUE "
            "WHERE a.interval INTERSECTS b.interval"
        )
        sql_dd = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act & assert
        with pytest.raises(duckdb.Error, match='(?i)referenced table "?b"? not found'):
            conn.execute(sql_dd)
        with pytest.raises(duckdb.Error, match='(?i)referenced table "?b"? not found'):
            conn.execute(sql_naive)

    @pytest.mark.parametrize("kind", ["ANTI", "SEMI"])
    def test_transpile_should_still_engage_iejoin_for_idiomatic_on_intersects(
        self, kind
    ):
        """Test that the idiomatic SEMI/ANTI ON-INTERSECTS form still uses the IEJoin.

        Given:
            The well-formed shape with the ``INTERSECTS`` in the join's own
            ``ON`` clause (the right table is legitimately in scope there).
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should still engage the IEJoin path (``SET VARIABLE
            __giql_iejoin_`` emitted) — the #201 decline must not over-reach
            and disable the supported idiomatic form.
        """
        # Arrange
        query = (
            f"SELECT a.name FROM peaks a {kind} JOIN genes b "
            "ON a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql

    def test_query_should_not_change_idiomatic_anti_result_after_decline_gate(
        self, peaks_genes
    ):
        """Test that the idiomatic ANTI ON-INTERSECTS still returns correct rows.

        Given:
            The idiomatic ``ANTI JOIN ... ON a.interval INTERSECTS b.interval``
            over the shared fixture (peaks with no overlapping gene survive).
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed.
        Then:
            It should return exactly the left rows with no overlap, matching
            the naive-predicate plan — confirming the new gate left the
            supported path untouched.
        """
        # Arrange
        conn = peaks_genes
        query = (
            "SELECT a.name FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval"
        )
        sql_dd = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act
        rows_dd = sorted(conn.execute(sql_dd).fetchall())
        rows_naive = sorted(conn.execute(sql_naive).fetchall())

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql_dd
        assert rows_dd == rows_naive
        assert rows_dd == [("p2",), ("p5",)]


class TestTranspileDuckDBIEJoinStarProjectionFallback:
    """Regression tests for #202.

    Star projections (bare ``*``, ``a.*``, ``b.*``, ``a.* AS x``) must decline
    the IEJoin rewrite so DuckDB expands the real star against the live schema,
    keeping every base-table column and staying identical to the naive plan —
    rather than silently narrowing the projection to the configured genomic
    columns.
    """

    @pytest.mark.parametrize(
        "projection",
        ["a.*", "b.*", "a.*, b.*", "*", "a.* AS x", "a.chrom, b.*"],
    )
    def test_transpile_should_decline_star_projection_to_naive_plan(self, projection):
        """Test that every star projection shape declines the IEJoin.

        Given:
            A column-to-column INTERSECTS INNER join whose SELECT list carries
            a star in some form (bare, qualified, aliased, or mixed).
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should decline the IEJoin (no ``SET VARIABLE`` scaffolding) so
            the naive-predicate plan expands the star against the live schema.
        """
        # Arrange
        query = (
            f"SELECT {projection} FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE" not in sql.upper()
        assert "getvariable" not in sql

    def test_query_should_match_naive_plan_for_a_star_under_anti(self, conn):
        """Test the #202 repro: ``a.*`` under ANTI preserves every column.

        Given:
            The exact issue reproduction — ``SELECT a.*`` over an ANTI join
            with one non-overlapping left row carrying ``name`` / ``score``
            beyond the configured genomic columns.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed,
            and separately with ``dialect=None``.
        Then:
            The DuckDB result should carry all six base-table columns
            (``name`` / ``score`` no longer dropped) and equal the naive plan
            in both columns and rows.
        """
        # Arrange
        _make_table(
            conn,
            "peaks",
            [("chr1", 100, 200, "p1", 7, "+"), ("chr1", 700, 800, "p3", 30, "+")],
        )
        _make_table(conn, "genes", [("chr1", 150, 250, "g1", 9, "-")])
        query = (
            "SELECT a.* FROM peaks a "
            "ANTI JOIN genes b ON a.interval INTERSECTS b.interval"
        )
        sql_dd = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act
        res_dd = conn.execute(sql_dd)
        cols_dd = [d[0] for d in res_dd.description]
        rows_dd = res_dd.fetchall()
        res_naive = conn.execute(sql_naive)
        cols_naive = [d[0] for d in res_naive.description]
        rows_naive = res_naive.fetchall()

        # Assert
        assert "SET VARIABLE" not in sql_dd.upper()
        assert cols_dd == ["chrom", "start", "end", "name", "score", "strand"]
        assert rows_dd == [("chr1", 700, 800, "p3", 30, "+")]
        assert cols_dd == cols_naive
        assert rows_dd == rows_naive

    def test_query_should_match_naive_plan_for_bare_star(self, conn):
        """Test that bare ``SELECT *`` returns both tables' full columns.

        Given:
            A ``SELECT *`` over an INNER INTERSECTS join, each table carrying
            non-genomic columns.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed,
            and separately with ``dialect=None``.
        Then:
            The DuckDB result should carry all twelve columns (both tables'
            full schema) and equal the naive plan in columns and rows.
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p1", 7, "+")])
        _make_table(conn, "genes", [("chr1", 150, 250, "g1", 9, "-")])
        query = "SELECT * FROM peaks a JOIN genes b ON a.interval INTERSECTS b.interval"
        sql_dd = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act
        res_dd = conn.execute(sql_dd)
        cols_dd = [d[0] for d in res_dd.description]
        rows_dd = res_dd.fetchall()
        res_naive = conn.execute(sql_naive)
        cols_naive = [d[0] for d in res_naive.description]
        rows_naive = res_naive.fetchall()

        # Assert
        assert "SET VARIABLE" not in sql_dd.upper()
        assert cols_dd == ["chrom", "start", "end", "name", "score", "strand"] * 2
        assert rows_dd == [
            ("chr1", 100, 200, "p1", 7, "+", "chr1", 150, 250, "g1", 9, "-")
        ]
        assert cols_dd == cols_naive
        assert rows_dd == rows_naive

    def test_transpile_should_still_engage_iejoin_for_explicit_columns(self):
        """Test that explicit (non-star) qualified projections still use the IEJoin.

        Given:
            A projection listing explicit qualified columns (no star).
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should still engage the IEJoin path (``SET VARIABLE
            __giql_iejoin_`` emitted) — the #202 star decline must not
            over-reach and disable non-star projections.
        """
        # Arrange
        query = (
            "SELECT a.chrom, a.start, a.name, b.score FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql

    def test_query_should_match_naive_plan_for_a_star_under_semi(self, conn):
        """Test that ``a.*`` under SEMI preserves every column, matching the naive plan.

        Given:
            A ``SELECT a.*`` over a SEMI join, where the one matching left
            row carries ``name`` / ``score`` beyond the genomic columns.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed,
            and separately with ``dialect=None``.
        Then:
            The DuckDB result should carry all six base-table columns and
            equal the naive plan in both columns and rows — the star decline
            preserves the full projection under SEMI just as under ANTI/INNER.
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p1", 7, "+")])
        _make_table(conn, "genes", [("chr1", 150, 250, "g1", 9, "-")])
        query = (
            "SELECT a.* FROM peaks a "
            "SEMI JOIN genes b ON a.interval INTERSECTS b.interval"
        )
        sql_dd = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act
        res_dd = conn.execute(sql_dd)
        cols_dd = [d[0] for d in res_dd.description]
        rows_dd = res_dd.fetchall()
        res_naive = conn.execute(sql_naive)
        cols_naive = [d[0] for d in res_naive.description]
        rows_naive = res_naive.fetchall()

        # Assert
        assert "SET VARIABLE" not in sql_dd.upper()
        assert cols_dd == ["chrom", "start", "end", "name", "score", "strand"]
        assert rows_dd == [("chr1", 100, 200, "p1", 7, "+")]
        assert cols_dd == cols_naive
        assert rows_dd == rows_naive

    def test_query_should_match_naive_plan_for_aliased_star(self, conn):
        """Test that ``a.* AS x`` executes to the same result as the naive plan.

        Given:
            A ``SELECT a.* AS x`` projection over an INNER INTERSECTS join.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed,
            and separately with ``dialect=None``.
        Then:
            The declined DuckDB plan should be genuinely executable and
            produce the same columns and rows as the naive plan (DuckDB
            applies the alias uniformly across the expanded star), confirming
            the fallback is not a broken shape.
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p1", 7, "+")])
        _make_table(conn, "genes", [("chr1", 150, 250, "g1", 9, "-")])
        query = (
            "SELECT a.* AS x FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval"
        )
        sql_dd = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act
        res_dd = conn.execute(sql_dd)
        cols_dd = [d[0] for d in res_dd.description]
        rows_dd = res_dd.fetchall()
        res_naive = conn.execute(sql_naive)
        cols_naive = [d[0] for d in res_naive.description]
        rows_naive = res_naive.fetchall()

        # Assert
        assert "SET VARIABLE" not in sql_dd.upper()
        assert cols_dd == cols_naive
        assert rows_dd == rows_naive
        assert rows_dd == [("chr1", 100, 200, "p1", 7, "+")]

    @pytest.mark.parametrize("kind", ["ANTI", "SEMI"])
    def test_query_should_raise_for_right_star_under_left_only_join_like_naive_plan(
        self, conn, kind
    ):
        """Test that a right-side ``b.*`` under SEMI/ANTI errors as the naive plan does.

        Given:
            A ``SELECT b.*`` projection over a SEMI / ANTI join — the right
            table is out of scope in a left-only join's output.
        When:
            The query is transpiled with ``dialect='duckdb'`` and with
            ``dialect=None`` and both are executed.
        Then:
            Both should raise the same DuckDB binder error naming the
            out-of-scope right table — the star decline defers the rejection
            to DuckDB rather than raising a bespoke transpile-time error.
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p1", 7, "+")])
        _make_table(conn, "genes", [("chr1", 150, 250, "g1", 9, "-")])
        query = (
            f"SELECT b.* FROM peaks a {kind} JOIN genes b "
            "ON a.interval INTERSECTS b.interval"
        )
        sql_dd = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act & assert
        assert "SET VARIABLE" not in sql_dd.upper()
        with pytest.raises(duckdb.Error, match='(?i)referenced table "?b"? not found'):
            conn.execute(sql_dd)
        with pytest.raises(duckdb.Error, match='(?i)referenced table "?b"? not found'):
            conn.execute(sql_naive)

    def test_query_should_raise_for_unknown_qualifier_star_like_naive_plan(self, conn):
        """Test that ``c.*`` errors at execution as the naive plan does.

        Given:
            A ``SELECT c.*`` projection whose qualifier ``c`` is not one of
            the join's two tables.
        When:
            The query is transpiled with ``dialect='duckdb'`` and with
            ``dialect=None`` and both are executed.
        Then:
            Both should raise the same DuckDB binder error naming the unknown
            table — the star decline defers the rejection to DuckDB rather
            than raising at transpile time.
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p1", 7, "+")])
        _make_table(conn, "genes", [("chr1", 150, 250, "g1", 9, "-")])
        query = (
            "SELECT c.* FROM peaks a JOIN genes b ON a.interval INTERSECTS b.interval"
        )
        sql_dd = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act & assert
        assert "SET VARIABLE" not in sql_dd.upper()
        with pytest.raises(duckdb.Error, match='(?i)referenced table "?c"? not found'):
            conn.execute(sql_dd)
        with pytest.raises(duckdb.Error, match='(?i)referenced table "?c"? not found'):
            conn.execute(sql_naive)


class TestTranspileDuckDBIEJoinUnsupportedProjectionFallback:
    """Regression tests for #204 and #205.

    Projections the naive-predicate plan compiles but the IEJoin projection
    rebuild cannot express — expressions (``a.start + 1``), window aggregates,
    ``FILTER`` clauses, scalar subqueries, aggregates nested in expressions
    (``COUNT(*) * 2``), and stars nested in an aggregate argument
    (``COUNT(a.*)`` / ``MIN(COLUMNS(*))``) — must decline the IEJoin and fall
    back to the naive plan, keeping ``dialect="duckdb"`` consistent with every
    other backend rather than hard-erroring (#205) or miscompiling (#204). A
    projection with an out-of-scope column reference stays a clean transpile
    error (the naive plan rejects it too).
    """

    @pytest.mark.parametrize(
        "projection",
        [
            "COUNT(a.*)",
            "COUNT(b.*)",
            "COUNT(a.* EXCLUDE (name))",
            "MIN(COLUMNS(*))",
            "a.start + 1",
            "SUM(a.score) OVER ()",
            "SUM(a.score) FILTER (WHERE a.score > 0)",
            "COUNT(*) * 2",
            "(SELECT 1) AS s",
            "100",
            "a.name, COUNT(a.*)",
        ],
    )
    def test_transpile_should_decline_unsupported_projection_to_naive_plan(
        self, projection
    ):
        """Test that every unsupported-but-naive-valid projection declines.

        Given:
            A column-to-column INTERSECTS INNER join whose SELECT list carries
            a projection the IEJoin cannot rebuild but the naive plan handles.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should decline the IEJoin (no ``SET VARIABLE`` scaffolding) so
            the naive-predicate plan compiles the projection.
        """
        # Arrange
        query = (
            f"SELECT {projection} FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE" not in sql.upper()
        assert "getvariable" not in sql

    def test_query_should_match_naive_plan_for_count_qualified_star(self, conn):
        """Test the #204 repro: ``COUNT(a.*)`` returns the naive count, not a crash.

        Given:
            A ``SELECT COUNT(a.*)`` over an INTERSECTS join with one
            overlapping pair.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed,
            and separately with ``dialect=None``.
        Then:
            The DuckDB result should equal the naive plan's count instead of
            crashing with a binder error on a synthesized ``a."*"`` column.
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p1", 7, "+")])
        _make_table(conn, "genes", [("chr1", 150, 250, "g1", 9, "-")])
        query = (
            "SELECT COUNT(a.*) FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval"
        )
        sql_dd = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act
        rows_dd = conn.execute(sql_dd).fetchall()
        rows_naive = conn.execute(sql_naive).fetchall()

        # Assert
        assert "SET VARIABLE" not in sql_dd.upper()
        assert rows_dd == rows_naive
        assert rows_dd == [(1,)]

    def test_query_should_match_naive_plan_for_min_columns_star(self, conn):
        """Test the #204 repro: ``MIN(COLUMNS(*))`` matches the naive per-column min.

        Given:
            A ``SELECT MIN(COLUMNS(*))`` over an INTERSECTS join.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed,
            and separately with ``dialect=None``.
        Then:
            The DuckDB result should equal the naive plan's per-column min
            across both tables' columns, not the previous bogus single scalar.
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p1", 7, "+")])
        _make_table(conn, "genes", [("chr1", 150, 250, "g1", 9, "-")])
        query = (
            "SELECT MIN(COLUMNS(*)) FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval"
        )
        sql_dd = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act
        res_dd = conn.execute(sql_dd)
        cols_dd = [d[0] for d in res_dd.description]
        rows_dd = res_dd.fetchall()
        res_naive = conn.execute(sql_naive)
        cols_naive = [d[0] for d in res_naive.description]
        rows_naive = res_naive.fetchall()

        # Assert
        assert "SET VARIABLE" not in sql_dd.upper()
        assert len(cols_dd) == 12
        assert cols_dd == cols_naive
        assert rows_dd == rows_naive

    @pytest.mark.parametrize(
        "projection, expected",
        [
            ("a.start + 1", [(101,)]),
            ("SUM(a.score) OVER ()", [(7,)]),
            ("SUM(a.score) FILTER (WHERE a.score > 0)", [(7,)]),
            ("COUNT(*) * 2", [(2,)]),
            ("(SELECT 1) AS s", [(1,)]),
            ("COUNT(a.*)", [(1,)]),
        ],
    )
    def test_query_should_match_naive_plan_for_unsupported_projection(
        self, conn, projection, expected
    ):
        """Test that a declined unsupported projection executes like the naive plan.

        Given:
            A projection the IEJoin cannot rebuild but the naive plan handles —
            an expression, a window aggregate, a FILTER clause, an
            arithmetic-over-aggregate, a scalar subquery, or an aggregate over a
            qualified star — over an INTERSECTS join with one overlapping pair.
        When:
            The query is transpiled with ``dialect='duckdb'`` and executed, and
            separately with ``dialect=None``.
        Then:
            The declined DuckDB plan should return the same rows as the naive
            plan (#204, #205), confirming the fallback is genuinely executable.
        """
        # Arrange
        _make_table(conn, "peaks", [("chr1", 100, 200, "p1", 7, "+")])
        _make_table(conn, "genes", [("chr1", 150, 250, "g1", 9, "-")])
        query = (
            f"SELECT {projection} FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval"
        )
        sql_dd = transpile(query, tables=["peaks", "genes"], dialect="duckdb")
        sql_naive = transpile(query, tables=["peaks", "genes"])

        # Act
        rows_dd = conn.execute(sql_dd).fetchall()
        rows_naive = conn.execute(sql_naive).fetchall()

        # Assert
        assert "SET VARIABLE" not in sql_dd.upper()
        assert rows_dd == rows_naive
        assert rows_dd == expected

    @pytest.mark.parametrize(
        "projection",
        [
            "score + 1",
            "SUM(score) OVER ()",
            "SUM(score) FILTER (WHERE score > 0)",
        ],
    )
    def test_transpile_should_raise_when_unsupported_wrapper_has_unqualified_column(
        self, projection
    ):
        """Test that an unqualified column inside an unsupported wrapper raises.

        Given:
            An expression / window aggregate / FILTER clause whose argument
            references an unqualified column the dialect cannot attribute to a
            join side.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should raise ``ValueError`` guiding the user to qualify the
            column (rather than declining) — every wrapper kind reaches the same
            diagnostic.
        """
        # Act & assert
        with pytest.raises(ValueError, match="qualified"):
            transpile(
                f"SELECT {projection} FROM peaks a "
                "JOIN genes b ON a.interval INTERSECTS b.interval",
                tables=["peaks", "genes"],
                dialect="duckdb",
            )

    @pytest.mark.parametrize("kind", ["SEMI", "ANTI"])
    def test_transpile_should_raise_left_side_only_for_right_ref_in_wrapper(self, kind):
        """Test that a right-side column in a wrapper names the left-only rule.

        Given:
            A ``SUM(b.score) OVER ()`` projection under a SEMI / ANTI join — a
            right-side reference the left-only output cannot resolve.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should raise the dedicated left-side-only ``ValueError`` (naming
            the right side) rather than the generic "qualify the column"
            message, which would steer the user toward another invalid form.
        """
        # Act & assert
        with pytest.raises(ValueError, match="left-side"):
            transpile(
                f"SELECT SUM(b.score) OVER () FROM peaks a {kind} JOIN genes b "
                "ON a.interval INTERSECTS b.interval",
                tables=["peaks", "genes"],
                dialect="duckdb",
            )

    @pytest.mark.parametrize("kind", ["SEMI", "ANTI"])
    def test_transpile_should_decline_unsupported_projection_under_left_only_join(
        self, kind
    ):
        """Test that an unsupported projection under SEMI/ANTI declines to naive.

        Given:
            An ``a.start + 1`` expression projection under a SEMI / ANTI join.
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should decline the IEJoin (no ``SET VARIABLE`` scaffolding) so
            the naive-predicate plan handles it — the decline applies under
            left-only joins as under INNER.
        """
        # Arrange
        query = (
            f"SELECT a.start + 1 FROM peaks a {kind} JOIN genes b "
            "ON a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="duckdb")

        # Assert
        assert "SET VARIABLE" not in sql.upper()
        assert "getvariable" not in sql

    def test_transpile_should_still_engage_iejoin_for_plain_aggregate(self):
        """Test that a plain qualified aggregate still engages the IEJoin.

        Given:
            A supported ``SUM(a.score)`` aggregate projection (no window,
            FILTER, expression wrapper, or star argument).
        When:
            ``transpile`` is called with ``dialect='duckdb'``.
        Then:
            It should still engage the IEJoin path (``SET VARIABLE
            __giql_iejoin_`` emitted) — the #204 / #205 declines must not
            over-reach and disable plain aggregates.
        """
        # Act
        sql = transpile(
            "SELECT SUM(a.score) FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval",
            tables=["peaks", "genes"],
            dialect="duckdb",
        )

        # Assert
        assert "SET VARIABLE __giql_iejoin_" in sql
