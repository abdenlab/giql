"""Tests for the CanonicalizeCoordinates normalization pass (pass 2).

These tests pin three guarantees of :mod:`giql.canonicalizer`:

* With every operator's ``GIQL_CANONICALIZE`` flag off (the state as of issue
  #121) the pass is a strict no-op and the emitted SQL is byte-identical for
  both canonical and non-canonical tables.
* With an operator opted in, the synthesized ``__giql_canon_*`` wrapper CTE
  projects each declared encoding to canonical 0-based half-open with the correct
  arithmetic, identity encodings pass through unwrapped, generated aliases avoid
  user-name collisions, identical wrapper bodies deduplicate to one CTE, and the
  CTEs are inserted in DAG order on the outermost ``WITH``.
"""

import itertools

import pytest
from sqlglot import exp
from sqlglot import parse_one

from giql.canonicalizer import CANON_PREFIX
from giql.canonicalizer import canonicalize_coordinates
from giql.dialect import GIQLDialect
from giql.expander import ExpandOperators
from giql.expressions import GIQLDisjoin
from giql.expressions import GIQLDistance
from giql.expressions import GIQLNearest
from giql.expressions import Intersects
from giql.resolver import META_KEY
from giql.resolver import resolve_operator_refs
from giql.table import Table
from giql.table import Tables
from giql.targets import DuckDBTarget
from giql.targets import GenericTarget
from giql.transpile import transpile

hypothesis = pytest.importorskip("hypothesis")
from hypothesis import given  # noqa: E402
from hypothesis import settings  # noqa: E402
from hypothesis import strategies as st  # noqa: E402

_COORD_SYSTEMS = ("0based", "1based")
_INTERVAL_TYPES = ("half_open", "closed")


def _tables(encoding=("0based", "half_open"), names=("variants",)) -> Tables:
    """Build a Tables container registering each name with the given encoding."""
    coordinate_system, interval_type = encoding
    container = Tables()
    for name in names:
        container.register(
            name,
            Table(
                name,
                coordinate_system=coordinate_system,
                interval_type=interval_type,
            ),
        )
    return container


def _prepare(query: str, tables: Tables) -> exp.Expression:
    """Parse *query* and run passes 1 and 2 over the resulting AST."""
    ast = parse_one(query, dialect=GIQLDialect)
    ast = resolve_operator_refs(ast, tables)
    return canonicalize_coordinates(ast)


def _resolve(query: str, tables: Tables) -> exp.Expression:
    """Parse *query* and run pass 1 only, leaving pass 2 for the caller."""
    return resolve_operator_refs(parse_one(query, dialect=GIQLDialect), tables)


def _canon_ctes(ast: exp.Expression) -> list[exp.CTE]:
    """Return every synthesized canonical CTE in *ast*, in declaration order."""
    with_ = ast.args.get("with_")
    if with_ is None:
        return []
    return [c for c in with_.expressions if c.alias.startswith(CANON_PREFIX)]


def _disjoin(ast: exp.Expression) -> GIQLDisjoin:
    """Return the single GIQLDisjoin node reachable from *ast*."""
    return next(n for n in ast.walk() if isinstance(n, GIQLDisjoin))


@pytest.fixture
def disjoin_opted_in(monkeypatch):
    """Opt GIQLDisjoin into canonicalization for the duration of a test."""
    monkeypatch.setattr(GIQLDisjoin, "GIQL_CANONICALIZE", True, raising=False)
    return GIQLDisjoin


class TestNoOpWhenFlagsOff:
    """The pass touches nothing while every GIQL_CANONICALIZE flag is off."""

    def test_canonical_table_sql_unchanged(self):
        """Test that a canonical table's SQL is byte-identical through the pass.

        Given:
            A DISJOIN over a canonical (0based/half_open) registered table and no
            operator opted into canonicalization.
        When:
            Transpiling the query.
        Then:
            The pass is inert and the SQL matches a run with the pass bypassed.
        """
        # Arrange
        query = "SELECT * FROM DISJOIN(variants)"
        tables = _tables(("0based", "half_open"))
        ast = resolve_operator_refs(parse_one(query, dialect=GIQLDialect), tables)
        # Bypass pass 2 (canonicalization) but keep pass 3 (expansion): DISJOIN
        # now goes through its registered expander, so the pass-bypassed reference
        # must too, isolating pass 2's contribution to nothing.
        ast = ExpandOperators(GenericTarget(), tables).transform(ast)
        expected = ast.sql()

        # Act
        actual = transpile(query, tables=[Table("variants")])

        # Assert
        assert actual == expected
        assert CANON_PREFIX not in actual

    def test_non_canonical_table_sql_unchanged(self, monkeypatch):
        """Test that a non-canonical table's SQL is byte-identical through the pass.

        Given:
            A DISJOIN over a non-canonical (1based/closed) registered table with
            its GIQL_CANONICALIZE flag explicitly toggled off (the gating an
            unmigrated operator relies on).
        When:
            Transpiling the query and comparing against a pass-bypassed run.
        Then:
            No wrapper CTE is synthesized and the SQL is unchanged.
        """
        # Arrange
        monkeypatch.setattr(GIQLDisjoin, "GIQL_CANONICALIZE", False, raising=False)
        query = "SELECT * FROM DISJOIN(variants)"
        variants = Table("variants", coordinate_system="1based", interval_type="closed")
        tables = Tables()
        tables.register("variants", variants)
        ast = resolve_operator_refs(parse_one(query, dialect=GIQLDialect), tables)
        # Bypass pass 2 (canonicalization) but keep pass 3 (expansion): DISJOIN
        # now expands through its registry entry, so the pass-bypassed reference
        # must run pass 3 too for the byte-identical comparison to isolate pass 2.
        ast = ExpandOperators(GenericTarget(), tables).transform(ast)
        expected = ast.sql()

        # Act
        actual = transpile(query, tables=[variants])

        # Assert
        assert actual == expected
        assert CANON_PREFIX not in actual

    def test_pass_returns_expression_with_no_with_added(self, monkeypatch):
        """Test that the no-op pass adds no WITH clause to the AST.

        Given:
            A non-canonical DISJOIN AST with its GIQL_CANONICALIZE flag explicitly
            toggled off.
        When:
            Running canonicalize_coordinates directly.
        Then:
            No canonical CTE is present on the returned tree.
        """
        # Arrange
        monkeypatch.setattr(GIQLDisjoin, "GIQL_CANONICALIZE", False, raising=False)
        tables = _tables(("1based", "closed"))

        # Act
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)

        # Assert
        assert _canon_ctes(ast) == []


class TestProjectionArithmetic:
    """The wrapper CTE replicates canonical_start/canonical_end per encoding."""

    def test_zero_based_closed_projection(self, disjoin_opted_in):
        """Test the canonical projection for a 0based/closed table.

        Given:
            A DISJOIN over a 0based/closed table with an operator opted in.
        When:
            Running pass 2.
        Then:
            The wrapper projects start unchanged and end + 1.
        """
        # Arrange
        tables = _tables(("0based", "closed"))

        # Act
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)

        # Assert
        body = _canon_ctes(ast)[0].this.sql()
        assert '"start" AS "start"' in body
        assert '("end" + 1) AS "end"' in body

    def test_one_based_half_open_projection(self, disjoin_opted_in):
        """Test the canonical projection for a 1based/half_open table.

        Given:
            A DISJOIN over a 1based/half_open table with an operator opted in.
        When:
            Running pass 2.
        Then:
            The wrapper projects start - 1 and end - 1.
        """
        # Arrange
        tables = _tables(("1based", "half_open"))

        # Act
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)

        # Assert
        body = _canon_ctes(ast)[0].this.sql()
        assert '("start" - 1) AS "start"' in body
        assert '("end" - 1) AS "end"' in body

    def test_one_based_closed_projection(self, disjoin_opted_in):
        """Test the canonical projection for a 1based/closed table.

        Given:
            A DISJOIN over a 1based/closed table with an operator opted in.
        When:
            Running pass 2.
        Then:
            The wrapper projects start - 1 and end unchanged.
        """
        # Arrange
        tables = _tables(("1based", "closed"))

        # Act
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)

        # Assert
        body = _canon_ctes(ast)[0].this.sql()
        assert '("start" - 1) AS "start"' in body
        assert '"end" AS "end"' in body

    def test_projection_preserves_custom_column_names(self, disjoin_opted_in):
        """Test the projection uses a table's physical column names.

        Given:
            A 1based/closed table with custom chrom/start/end column names.
        When:
            Running pass 2.
        Then:
            The wrapper exposes the canonical interval under those same names;
            chrom (never offset) flows through the star untouched.
        """
        # Arrange
        variants = Table(
            "variants",
            chrom_col="chr",
            start_col="lo",
            end_col="hi",
            coordinate_system="1based",
            interval_type="closed",
        )
        tables = Tables()
        tables.register("variants", variants)

        # Act
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)

        # Assert
        body = _canon_ctes(ast)[0].this.sql()
        assert body.startswith("SELECT *")
        assert '("lo" - 1) AS "lo"' in body
        assert '"hi" AS "hi"' in body


class TestPassThrough:
    """Identity-encoded and non-table references pass through unwrapped."""

    def test_identity_encoding_not_wrapped(self, disjoin_opted_in):
        """Test that a canonical table is not wrapped even when opted in.

        Given:
            A DISJOIN over a canonical (0based/half_open) table with the operator
            opted into canonicalization.
        When:
            Running pass 2.
        Then:
            No wrapper CTE is synthesized and the slot metadata is unchanged.
        """
        # Arrange
        tables = _tables(("0based", "half_open"))

        # Act
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)

        # Assert
        assert _canon_ctes(ast) == []
        ref = _disjoin(ast).meta[META_KEY].slot("this")
        assert ref.kind == "registered_table"

    def test_cte_reference_not_wrapped(self, disjoin_opted_in):
        """Test that a CTE reference passes through unwrapped.

        Given:
            A DISJOIN whose reference resolves to an enclosing CTE while its
            target is a non-canonical registered table.
        When:
            Running pass 2.
        Then:
            Only the target is wrapped; the CTE reference is left as a CTE ref.
        """
        # Arrange
        query = (
            "WITH foo AS (SELECT * FROM other) "
            "SELECT * FROM DISJOIN(variants, reference := foo)"
        )
        tables = _tables(("1based", "closed"))

        # Act
        ast = _prepare(query, tables)

        # Assert
        slots = _disjoin(ast).meta[META_KEY]
        assert slots.slot("this").kind == "cte"  # target rewritten to canon CTE
        assert slots.slot("this").name.startswith(CANON_PREFIX)
        assert slots.slot("reference").name == "foo"  # user CTE untouched

    def test_subquery_reference_not_wrapped(self, disjoin_opted_in):
        """Test that a subquery reference passes through unwrapped.

        Given:
            A DISJOIN whose reference is a subquery and whose target is a
            non-canonical registered table.
        When:
            Running pass 2.
        Then:
            The subquery reference stays a subquery ref carrying no Table config.
        """
        # Arrange
        query = "SELECT * FROM DISJOIN(variants, reference := (SELECT * FROM other))"
        tables = _tables(("1based", "closed"))

        # Act
        ast = _prepare(query, tables)

        # Assert
        ref = _disjoin(ast).meta[META_KEY].slot("reference")
        assert ref.kind == "subquery"
        assert ref.table is None


class TestSlotRewrite:
    """A wrapped slot is rewritten in both the AST and its metadata."""

    def test_metadata_rewritten_to_canonical_cte(self, disjoin_opted_in):
        """Test that a wrapped slot's ResolvedRef points at the canonical CTE.

        Given:
            A DISJOIN over a non-canonical table with the operator opted in.
        When:
            Running pass 2.
        Then:
            The slot's ResolvedRef becomes a CTE ref naming the canon CTE and
            carrying no Table config.
        """
        # Arrange
        tables = _tables(("1based", "closed"))

        # Act
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)

        # Assert
        ref = _disjoin(ast).meta[META_KEY].slot("this")
        assert ref.kind == "cte"
        assert ref.name.startswith(CANON_PREFIX)
        assert ref.table is None
        assert ref.cols == ("chrom", "start", "end")

    def test_ast_slot_points_at_canonical_cte(self, disjoin_opted_in):
        """Test that a wrapped slot's AST node names the canonical CTE.

        Given:
            A DISJOIN over a non-canonical table with the operator opted in.
        When:
            Running pass 2.
        Then:
            The operator's target AST node references the canon CTE by name.
        """
        # Arrange
        tables = _tables(("1based", "closed"))

        # Act
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)

        # Assert
        target = _disjoin(ast).this
        assert isinstance(target, exp.Table)
        assert target.name.startswith(CANON_PREFIX)


class TestDeduplication:
    """Identical wrapper bodies collapse to a single canonical CTE."""

    def test_same_table_two_slots_dedup_to_one_cte(self, disjoin_opted_in):
        """Test that two slots over the same table share one canonical CTE.

        Given:
            A self-referencing DISJOIN over a non-canonical table (both target
            and reference resolve to the same table) with the operator opted in.
        When:
            Running pass 2.
        Then:
            Exactly one canonical CTE is synthesized and both slots point at it.
        """
        # Arrange
        tables = _tables(("1based", "closed"))

        # Act
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)

        # Assert
        ctes = _canon_ctes(ast)
        assert len(ctes) == 1
        slots = _disjoin(ast).meta[META_KEY]
        assert slots.slot("this").name == slots.slot("reference").name == ctes[0].alias

    def test_distinct_tables_get_distinct_ctes(self, disjoin_opted_in):
        """Test that two distinct non-canonical tables get distinct CTEs.

        Given:
            A DISJOIN whose target and reference are two distinct non-canonical
            registered tables, with the operator opted in.
        When:
            Running pass 2.
        Then:
            Two distinct canonical CTEs are synthesized.
        """
        # Arrange
        tables = _tables(("1based", "closed"), names=("variants", "genes"))

        # Act
        ast = _prepare("SELECT * FROM DISJOIN(variants, reference := genes)", tables)

        # Assert
        ctes = _canon_ctes(ast)
        assert len(ctes) == 2
        slots = _disjoin(ast).meta[META_KEY]
        assert slots.slot("this").name != slots.slot("reference").name


class TestAliasCollision:
    """Generated aliases avoid colliding with user identifiers."""

    def test_alias_avoids_user_cte_collision(self, disjoin_opted_in):
        """Test that a generated alias dodges a same-named user CTE.

        Given:
            A query containing a user CTE literally named __giql_canon_0 and a
            non-canonical DISJOIN, with the operator opted in.
        When:
            Running pass 2.
        Then:
            The synthesized CTE takes a distinct, non-colliding alias.
        """
        # Arrange
        query = (
            f"WITH {CANON_PREFIX}0 AS (SELECT * FROM other) "
            "SELECT * FROM DISJOIN(variants)"
        )
        tables = _tables(("1based", "closed"))

        # Act
        ast = _prepare(query, tables)

        # Assert
        synthesized = _disjoin(ast).meta[META_KEY].slot("this").name
        assert synthesized.startswith(CANON_PREFIX)
        assert synthesized != f"{CANON_PREFIX}0"

    def test_alias_avoids_user_table_collision(self, disjoin_opted_in):
        """Test that a generated alias dodges a same-named user table.

        Given:
            A query joining a user table literally named __giql_canon_0 alongside
            a non-canonical DISJOIN, with the operator opted in.
        When:
            Running pass 2.
        Then:
            The synthesized CTE alias does not collide with that table name.
        """
        # Arrange
        query = (
            f"SELECT * FROM DISJOIN(variants) AS d JOIN {CANON_PREFIX}0 AS u ON 1 = 1"
        )
        tables = _tables(("1based", "closed"))

        # Act
        ast = _prepare(query, tables)

        # Assert
        synthesized = _disjoin(ast).meta[META_KEY].slot("this").name
        assert synthesized != f"{CANON_PREFIX}0"


class TestDagOrdering:
    """Canonical CTEs are inserted before the user CTEs that may reference them."""

    def test_canon_cte_prepended_before_existing_cte(self, disjoin_opted_in):
        """Test that a synthesized CTE precedes existing user CTEs in the WITH.

        Given:
            A query with a user CTE and a non-canonical DISJOIN, opted in.
        When:
            Running pass 2.
        Then:
            The canonical CTE is the first entry in the outermost WITH clause.
        """
        # Arrange
        query = "WITH userfoo AS (SELECT * FROM other) SELECT * FROM DISJOIN(variants)"
        tables = _tables(("1based", "closed"))

        # Act
        ast = _prepare(query, tables)

        # Assert
        with_exprs = ast.args["with_"].expressions
        assert with_exprs[0].alias.startswith(CANON_PREFIX)
        assert "userfoo" in {c.alias for c in with_exprs}

    def test_canon_cte_created_when_no_existing_with(self, disjoin_opted_in):
        """Test that a fresh WITH clause is created when none exists.

        Given:
            A non-canonical DISJOIN with no enclosing WITH, opted in.
        When:
            Running pass 2.
        Then:
            A WITH clause is attached carrying the canonical CTE.
        """
        # Arrange
        tables = _tables(("1based", "closed"))

        # Act
        ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)

        # Assert
        assert ast.args.get("with_") is not None
        assert len(_canon_ctes(ast)) == 1


class TestEncodingInvariants:
    """Property-based coverage over the full encoding domain."""

    @settings(max_examples=25, deadline=None)
    @given(
        coordinate_system=st.sampled_from(_COORD_SYSTEMS),
        interval_type=st.sampled_from(_INTERVAL_TYPES),
    )
    def test_wrap_iff_non_canonical(self, coordinate_system, interval_type):
        """Test that a table is wrapped exactly when its encoding is non-canonical.

        Given:
            A DISJOIN over a registered table with any encoding combination, with
            the operator opted into canonicalization.
        When:
            Running passes 1 and 2.
        Then:
            A canonical CTE exists iff the encoding is not 0based/half_open, and a
            wrapped slot always becomes a Table-free CTE ref.
        """
        # Arrange
        # DISJOIN is opted into canonicalization by default (issue #122); a plain
        # save/restore toggle pins it on, since @given re-runs the body many times
        # under one function-scoped fixture, making monkeypatch unsafe.
        previous = GIQLDisjoin.GIQL_CANONICALIZE
        GIQLDisjoin.GIQL_CANONICALIZE = True
        tables = _tables((coordinate_system, interval_type))
        is_canonical = coordinate_system == "0based" and interval_type == "half_open"

        # Act
        try:
            ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)
        finally:
            GIQLDisjoin.GIQL_CANONICALIZE = previous

        # Assert
        ctes = _canon_ctes(ast)
        ref = _disjoin(ast).meta[META_KEY].slot("this")
        if is_canonical:
            assert ctes == []
            assert ref.kind == "registered_table"
        else:
            assert len(ctes) == 1
            assert ref.kind == "cte"
            assert ref.table is None

    @settings(max_examples=25, deadline=None)
    @given(
        coordinate_system=st.sampled_from(_COORD_SYSTEMS),
        interval_type=st.sampled_from(_INTERVAL_TYPES),
    )
    def test_flags_off_is_always_noop(self, coordinate_system, interval_type):
        """Test that the pass never mutates the tree while flags are off.

        Given:
            A DISJOIN over a registered table with any encoding and its
            GIQL_CANONICALIZE flag explicitly toggled off (the gating an
            unmigrated operator relies on).
        When:
            Running passes 1 and 2.
        Then:
            No canonical CTE is ever synthesized.
        """
        # Arrange
        # Save/restore the class flag off: @given re-runs the body many times
        # under one function-scoped fixture, making monkeypatch unsafe.
        previous = GIQLDisjoin.GIQL_CANONICALIZE
        GIQLDisjoin.GIQL_CANONICALIZE = False
        tables = _tables((coordinate_system, interval_type))

        # Act
        try:
            ast = _prepare("SELECT * FROM DISJOIN(variants)", tables)
        finally:
            GIQLDisjoin.GIQL_CANONICALIZE = previous

        # Assert
        assert _canon_ctes(ast) == []


def test_all_encoding_pairs_covered():
    """Test that the sampled encoding domain spans every coordinate/type pair.

    Given:
        The coordinate-system and interval-type domains under test.
    When:
        Enumerating their product.
    Then:
        Exactly four encoding combinations exist, one of them canonical.
    """
    # Arrange
    pairs = list(itertools.product(_COORD_SYSTEMS, _INTERVAL_TYPES))

    # Act
    canonical = [p for p in pairs if p == ("0based", "half_open")]

    # Assert
    assert len(pairs) == 4
    assert len(canonical) == 1


def _two_tables(encoding) -> list[Table]:
    """Build two registered tables under the same (non-default) encoding."""
    coordinate_system, interval_type = encoding
    return [
        Table(name, coordinate_system=coordinate_system, interval_type=interval_type)
        for name in ("intervals_a", "intervals_b")
    ]


class TestColumnOperandCanonicalization:
    """Pass 2 canonicalizes column / interval operands in place (issue #123).

    A column operand references an alias bound in the enclosing query's FROM,
    shared with the user's own projection, so it cannot be wrapped in a canonical
    CTE without changing unrelated columns. Pass 2 therefore rewrites the operand's
    resolution metadata to carry the canonical arithmetic inline and the emitter
    consumes it verbatim — no in-emitter canonicalization, no wrapper CTE.
    """

    def test_distance_operand_canonicalized_inline_without_wrapper_cte(self):
        """Test a non-canonical DISTANCE operand canonicalizes inline, no wrapper.

        Given:
            Two 1-based closed tables and a DISTANCE between their columns.
        When:
            The query is transpiled.
        Then:
            The CASE arithmetic should wrap each side's start as (start - 1)
            inline and synthesize no __giql_canon_* wrapper CTE.
        """
        # Arrange
        sql = (
            "SELECT DISTANCE(a.interval, b.interval) AS dist "
            "FROM intervals_a a, intervals_b b"
        )

        # Act
        output = transpile(sql, tables=_two_tables(("1based", "closed")))

        # Assert
        assert CANON_PREFIX not in output
        assert '(a."start" - 1)' in output
        assert '(b."start" - 1)' in output

    def test_canonical_distance_operand_is_byte_identical(self):
        """Test a canonical DISTANCE operand emits byte-identical SQL.

        Given:
            Two default (0-based half-open) tables and a DISTANCE between their
            columns, transpiled with the operator opted in and with its
            GIQL_CANONICALIZE flag toggled off.
        When:
            The two transpilations are compared.
        Then:
            They should be byte-identical and carry no wrapper CTE — the pass is
            inert for an already-canonical operand.
        """
        # Arrange
        sql = (
            "SELECT DISTANCE(a.interval, b.interval) AS dist "
            "FROM intervals_a a, intervals_b b"
        )
        tables = [Table("intervals_a"), Table("intervals_b")]

        # Act
        opted_in = transpile(sql, tables=tables)
        previous = GIQLDistance.GIQL_CANONICALIZE
        GIQLDistance.GIQL_CANONICALIZE = False
        try:
            flag_off = transpile(sql, tables=tables)
        finally:
            GIQLDistance.GIQL_CANONICALIZE = previous

        # Assert
        assert opted_in == flag_off
        assert CANON_PREFIX not in opted_in

    def test_predicate_operand_canonicalized_inline_without_wrapper_cte(self):
        """Test a non-canonical INTERSECTS operand canonicalizes inline, no wrapper.

        Given:
            A 1-based closed table and a literal-range INTERSECTS predicate.
        When:
            The query is transpiled.
        Then:
            The predicate should wrap the table-side start as (start - 1) inline
            and synthesize no __giql_canon_* wrapper CTE.
        """
        # Arrange
        sql = "SELECT * FROM variants WHERE interval INTERSECTS 'chr1:100-200'"

        # Act
        output = transpile(
            sql,
            tables=[
                Table("variants", coordinate_system="1based", interval_type="closed")
            ],
        )

        # Assert
        assert CANON_PREFIX not in output
        assert '("start" - 1) < 200' in output

    def test_metadata_blanked_after_in_place_canonicalization(self):
        """Test a canonicalized column operand carries a blanked Table.

        Given:
            A 1-based closed INTERSECTS predicate annotated by pass 1.
        When:
            Pass 2 runs.
        Then:
            The operand's ResolvedColumn should carry the canonical arithmetic and
            its Table should be blanked so the emitter applies no further wrapping.
        """
        # Arrange
        tables = _tables(("1based", "closed"), names=("variants",))
        query = "SELECT * FROM variants WHERE interval INTERSECTS 'chr1:100-200'"
        ast = resolve_operator_refs(parse_one(query, dialect=GIQLDialect), tables)

        # Act
        ast = canonicalize_coordinates(ast)

        # Assert
        node = next(n for n in ast.walk() if isinstance(n, Intersects))
        column = node.meta[META_KEY].column("this")
        assert column.table is None
        assert column.start == '("start" - 1)'


class TestNearestTargetCanonicalization:
    """NEAREST's registered-table target is wrapped and its row round-trips."""

    def test_non_canonical_target_wrapped_and_row_round_tripped(self):
        """Test a non-canonical NEAREST target is wrapped and its row de-canonicalized.

        Given:
            A 1-based closed NEAREST target table and a literal reference.
        When:
            The query is transpiled.
        Then:
            The target should be wrapped in a __giql_canon_* CTE, the distance
            CASE should read the bare canonical columns, and the passed-through
            row should de-canonicalize the interval back to the declared encoding
            — via a star REPLACE on DuckDB and the portable star EXCEPT form on
            the generic / DataFusion family (#145).
        """
        # Arrange
        sql = "SELECT * FROM NEAREST(genes, reference := 'chr1:1000-2000', k := 1)"
        tables = [Table("genes", coordinate_system="1based", interval_type="closed")]

        # Act
        duckdb_output = transpile(sql, tables=tables, dialect="duckdb")
        generic_output = transpile(sql, tables=tables)

        # Assert
        assert f"{CANON_PREFIX}0 AS (SELECT * REPLACE" in duckdb_output
        assert f"{CANON_PREFIX}0 AS (SELECT * EXCEPT" in generic_output
        for output in (duckdb_output, generic_output):
            assert 'WHEN 1000 < __giql_canon_0."end"' in output
        assert (
            '__giql_canon_0.* REPLACE ((__giql_canon_0."start" + 1) AS "start"'
        ) in duckdb_output
        assert (
            '__giql_canon_0.* EXCEPT ("start", "end"), '
            '(__giql_canon_0."start" + 1) AS "start"'
        ) in generic_output

    def test_canonical_target_not_wrapped(self):
        """Test a canonical NEAREST target is left unwrapped.

        Given:
            A canonical 0-based half-open NEAREST target and a literal reference,
            transpiled with NEAREST opted in and with its flag toggled off.
        When:
            The two transpilations are compared.
        Then:
            They should be byte-identical with no wrapper CTE — the identity fast
            path.
        """
        # Arrange
        sql = "SELECT * FROM NEAREST(genes, reference := 'chr1:1000-2000', k := 1)"
        tables = [Table("genes")]

        # Act
        opted_in = transpile(sql, tables=tables)
        previous = GIQLNearest.GIQL_CANONICALIZE
        GIQLNearest.GIQL_CANONICALIZE = False
        try:
            flag_off = transpile(sql, tables=tables)
        finally:
            GIQLNearest.GIQL_CANONICALIZE = previous

        # Assert
        assert opted_in == flag_off
        assert CANON_PREFIX not in opted_in


class TestCanonicalProjectionCapabilities:
    """The wrapper-CTE emit strategy is chosen from target capabilities (#145)."""

    def test_canonicalize_coordinates_should_emit_except_without_star_replace(
        self, disjoin_opted_in
    ):
        """Test the wrapper uses the portable EXCEPT form without star-replace support.

        Given:
            A non-canonical DISJOIN AST and a Capabilities lacking
            supports_star_replace (the generic / DataFusion family).
        When:
            canonicalize_coordinates is called directly with those capabilities.
        Then:
            The synthesized wrapper CTE should project the portable
            SELECT * EXCEPT (...) form and never a star REPLACE.
        """
        # Arrange
        ast = _resolve("SELECT * FROM DISJOIN(variants)", _tables(("1based", "closed")))

        # Act
        canonicalize_coordinates(ast, GenericTarget().capabilities)

        # Assert
        body = _canon_ctes(ast)[0].this.sql()
        assert "* EXCEPT" in body
        assert "REPLACE" not in body
        assert '("start" - 1) AS "start"' in body

    def test_canonicalize_coordinates_should_emit_replace_with_star_replace(
        self, disjoin_opted_in
    ):
        """Test the wrapper uses star REPLACE on a REPLACE-capable target.

        Given:
            A non-canonical DISJOIN AST and a Capabilities with
            supports_star_replace (DuckDB).
        When:
            canonicalize_coordinates is called directly with those capabilities.
        Then:
            The synthesized wrapper CTE should substitute the interval columns in
            place via a star REPLACE and never the EXCEPT form.
        """
        # Arrange
        ast = _resolve("SELECT * FROM DISJOIN(variants)", _tables(("1based", "closed")))

        # Act
        canonicalize_coordinates(ast, DuckDBTarget().capabilities)

        # Assert
        body = _canon_ctes(ast)[0].this.sql()
        assert "* REPLACE" in body
        assert "EXCEPT" not in body
        assert '("start" - 1) AS "start"' in body

    def test_canonicalize_coordinates_should_emit_replace_when_capabilities_none(
        self, disjoin_opted_in
    ):
        """Test a direct caller with no capabilities keeps the historical REPLACE form.

        Given:
            A non-canonical DISJOIN AST.
        When:
            canonicalize_coordinates is called with no capabilities argument.
        Then:
            The wrapper should default to the star REPLACE form, preserving the
            historical behavior for direct callers.
        """
        # Arrange
        ast = _resolve("SELECT * FROM DISJOIN(variants)", _tables(("1based", "closed")))

        # Act
        canonicalize_coordinates(ast)

        # Assert
        body = _canon_ctes(ast)[0].this.sql()
        assert "* REPLACE" in body
        assert "EXCEPT" not in body

    def test_canonicalize_coordinates_should_not_re_wrap_on_a_second_pass(
        self, disjoin_opted_in
    ):
        """Test the pass is idempotent over the EXCEPT wrapper it already inserted.

        Given:
            A non-canonical DISJOIN AST canonicalized once with EXCEPT-form
            capabilities.
        When:
            canonicalize_coordinates is run a second time over the same AST.
        Then:
            It should not synthesize a second wrapper — the rewritten slot now
            points at the canonical CTE, so exactly one __giql_canon_ CTE remains.
        """
        # Arrange
        caps = GenericTarget().capabilities
        ast = _resolve("SELECT * FROM DISJOIN(variants)", _tables(("1based", "closed")))
        canonicalize_coordinates(ast, caps)

        # Act
        canonicalize_coordinates(ast, caps)

        # Assert
        assert len(_canon_ctes(ast)) == 1

    def test_transpile_should_emit_except_wrapper_for_datafusion_dialect(self):
        """Test the DataFusion dialect threads its capabilities into the wrapper.

        Given:
            A non-canonical DISJOIN over a 1-based closed table.
        When:
            Transpiling with dialect="datafusion".
        Then:
            The wrapper CTE and passthrough should use the portable EXCEPT form
            (DataFusion lacks star REPLACE) with no leaked operator.
        """
        # Arrange
        tables = [Table("variants", coordinate_system="1based", interval_type="closed")]

        # Act
        sql = transpile(
            "SELECT * FROM DISJOIN(variants)", tables=tables, dialect="datafusion"
        )

        # Assert
        assert "SELECT * EXCEPT" in sql
        assert "REPLACE" not in sql
        assert "G_I_Q_L" not in sql

    def test_transpile_should_emit_except_wrapper_when_no_dialect(self):
        """Test the generic (no-dialect) default emits the portable EXCEPT wrapper.

        Given:
            A non-canonical DISJOIN over a 1-based closed table.
        When:
            Transpiling with no dialect (the generic default target).
        Then:
            The wrapper CTE should use the portable EXCEPT form, not star REPLACE
            — pinning that the no-dialect default resolves GenericTarget
            (supports_star_replace=False) and emits the not-DuckDB-runnable
            portable form, guarding against a regression back to REPLACE.
        """
        # Arrange
        tables = [Table("variants", coordinate_system="1based", interval_type="closed")]

        # Act
        sql = transpile("SELECT * FROM DISJOIN(variants)", tables=tables)

        # Assert
        assert "SELECT * EXCEPT" in sql
        assert "REPLACE" not in sql
