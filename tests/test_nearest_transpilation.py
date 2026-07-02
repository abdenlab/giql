"""Transpilation tests for NEAREST operator SQL generation.

Tests verify that NEAREST() is correctly transpiled to SQL
(LATERAL joins for correlated queries, ORDER BY + LIMIT for standalone).
"""

import re

import pytest
from sqlglot import parse_one

import giql.expanders  # noqa: F401  (side-effect: registers the NEAREST expander)
from giql import Table
from giql.canonicalizer import canonicalize_coordinates
from giql.dialect import GIQLDialect
from giql.expander import ExpandOperators
from giql.resolver import resolve_operator_refs
from giql.table import Tables
from giql.targets import DataFusionTarget
from giql.targets import GenericTarget


def _generate_for_target(sql: str, tables: Tables, target) -> str:
    """Parse, run passes 1-3 against *target*, then generate SQL.

    Drives the expander for a specific :class:`~giql.targets.Target` so a
    capability-dependent shape (e.g. DataFusion's decorrelated window fallback,
    chosen because ``supports_lateral`` is False) can be asserted without an
    engine.
    """
    ast = parse_one(sql, dialect=GIQLDialect)
    ast = resolve_operator_refs(ast, tables)
    ast = canonicalize_coordinates(ast)
    ast = ExpandOperators(target, tables).transform(ast)
    return ast.sql()


def _generate(sql: str, tables: Tables) -> str:
    """Parse, run normalization passes 1-3, then generate SQL.

    Operator resolution, coordinate canonicalization, and operator expansion
    moved out of the emitter into the ResolveOperatorRefs / CanonicalizeCoordinates
    / ExpandOperators passes (epics #114, #137). NEAREST is now produced by its
    registered expander (issue #142) rather than a ``giqlnearest_sql`` emitter, so
    these tests must run pass 3 before generating, exactly as
    :func:`giql.transpile.transpile` does, rather than calling ``generate`` on a
    bare parsed AST.
    """
    ast = parse_one(sql, dialect=GIQLDialect)
    ast = resolve_operator_refs(ast, tables)
    ast = canonicalize_coordinates(ast)
    ast = ExpandOperators(GenericTarget(), tables).transform(ast)
    return ast.sql()


@pytest.fixture
def tables_with_peaks_and_genes():
    """Tables container with peaks and genes tables."""
    tables = Tables()
    tables.register("peaks", Table("peaks"))
    tables.register("genes", Table("genes"))
    return tables


class TestNearestTranspilation:
    """Tests for NEAREST transpilation to SQL."""

    def test_nearest_basic_k3(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k := 3)
        WHEN transpiling to SQL
        THEN should generate LATERAL join with DISTANCE and LIMIT 3
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference := peaks.interval, k := 3)
        """

        output = _generate(sql, tables_with_peaks_and_genes)

        # Expectations:
        # - LATERAL subquery
        # - DISTANCE(...) AS distance in SELECT
        # - WHERE peaks.chromosome = genes.chromosome (pre-filter)
        # - ORDER BY distance
        # - LIMIT 3
        assert "LATERAL" in output.upper()
        assert "CASE" in output or "DISTANCE" in output  # Distance calculation
        assert " AS distance" in output or " as distance" in output.lower()
        assert "LIMIT 3" in output
        assert "ORDER BY" in output

    def test_nearest_with_max_distance(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k := 5, max_distance := 100000)
        WHEN transpiling to SQL
        THEN should generate LATERAL join with distance filter
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference := peaks.interval, k := 5, max_distance := 100000)
        """

        output = _generate(sql, tables_with_peaks_and_genes)

        # Expectations:
        # - LATERAL subquery
        # - Distance filter: <= 100000
        # - LIMIT 5
        assert "LATERAL" in output.upper()
        assert "100000" in output
        assert "LIMIT 5" in output

    def test_nearest_standalone_literal(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQL query with literal reference NEAREST(genes, reference := 'chr1:1000-2000', k := 3)
        WHEN transpiling to SQL
        THEN should generate standalone query without LATERAL
        """
        sql = """
        SELECT *
        FROM NEAREST(genes, reference := 'chr1:1000-2000', k := 3)
        """

        output = _generate(sql, tables_with_peaks_and_genes)

        # Expectations:
        # - No LATERAL (standalone mode)
        # - Distance calculation with literal 'chr1', 1000, 2000
        # - ORDER BY distance
        # - LIMIT 3
        assert "LATERAL" not in output.upper()
        assert "chr1" in output.lower()
        assert "LIMIT 3" in output

    def test_nearest_with_stranded(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k := 3, stranded := true)
        WHEN transpiling to SQL
        THEN should generate SQL with strand filtering
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference := peaks.interval, k := 3, stranded := true)
        """

        output = _generate(sql, tables_with_peaks_and_genes)

        # Expectations:
        # - LATERAL subquery
        # - Strand filtering in WHERE clause
        # - LIMIT 3
        assert "LATERAL" in output.upper()
        assert "strand" in output.lower()
        assert "LIMIT 3" in output

    def test_nearest_with_signed(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQL query with NEAREST(genes, k := 3, signed := true)
        WHEN transpiling to SQL
        THEN should generate SQL with signed distance column
            (negative for upstream, positive for downstream)
        """
        sql = """
        SELECT *
        FROM peaks
        CROSS JOIN LATERAL NEAREST(genes, reference := peaks.interval, k := 3, signed := true)
        """

        output = _generate(sql, tables_with_peaks_and_genes)

        # Expectations:
        # - LATERAL subquery
        # - Signed distance calculation (includes negation for upstream)
        # - LIMIT 3
        assert "LATERAL" in output.upper()
        assert "LIMIT 3" in output
        # Check for signed distance: the ELSE branch should have a negation
        # for upstream features (B before A)
        assert "ELSE -(" in output, (
            f"Expected signed distance with negation for upstream, got:\n{output}"
        )


class TestNearestDataFusionFallbackShape:
    """Engine-free transpile-shape checks for the DataFusion window fallback (A8).

    A correlated NEAREST on the DataFusion target (``supports_lateral`` is False)
    expands to the decorrelated window-function form. These assert its structural
    invariants without running an engine: the window is present, the top-k filter
    is a ``<= k`` predicate, no correlated ``LATERAL`` survives, and the candidate
    cross-join and the window live at separate query levels.
    """

    def test_fallback_emits_window_with_topk_and_no_lateral(
        self, tables_with_peaks_and_genes
    ):
        """Test the DataFusion fallback emits a windowed top-k with no LATERAL.

        Given:
            A correlated NEAREST(genes, k := 1) on the DataFusion target.
        When:
            Transpiling.
        Then:
            It should emit a ROW_NUMBER() window, a `<= 1` top-k predicate, no
            surviving LATERAL, and the cross-join and window at separate query
            levels.
        """
        # Arrange
        sql = (
            "SELECT * FROM peaks "
            "CROSS JOIN LATERAL NEAREST(genes, reference := peaks.interval, k := 1) AS b"
        )

        # Act
        output = _generate_for_target(
            sql, tables_with_peaks_and_genes, DataFusionTarget()
        )

        # Assert
        assert "ROW_NUMBER(" in output.upper()
        assert "OVER (" in output.upper()
        assert "<= 1" in output
        assert "LATERAL" not in output.upper()
        # The candidate cross-join sits one level below the window: the window's
        # FROM is a parenthesized subquery, so a CROSS JOIN appears nested inside.
        assert "CROSS JOIN" in output.upper()

    def test_fallback_stranded_emits_window_and_strand_match(
        self, tables_with_peaks_and_genes
    ):
        """Test the stranded DataFusion fallback keeps a strand match, no LATERAL.

        Given:
            A stranded correlated NEAREST on the DataFusion target.
        When:
            Transpiling.
        Then:
            It should emit the window form, keep a strand equality in the
            candidate WHERE, and surface no LATERAL.
        """
        # Arrange
        sql = (
            "SELECT * FROM peaks CROSS JOIN LATERAL "
            "NEAREST(genes, reference := peaks.interval, k := 1, stranded := true) AS b"
        )

        # Act
        output = _generate_for_target(
            sql, tables_with_peaks_and_genes, DataFusionTarget()
        )

        # Assert
        assert "ROW_NUMBER(" in output.upper()
        assert "LATERAL" not in output.upper()
        assert 'peaks."strand"' in output
        assert 'genes."strand"' in output

    def test_fallback_k_greater_than_one_uses_k_in_topk_predicate(
        self, tables_with_peaks_and_genes
    ):
        """Test the DataFusion fallback carries the requested k in its top-k filter.

        Given:
            A correlated NEAREST(genes, k := 3) on the DataFusion target.
        When:
            Transpiling.
        Then:
            The top-k predicate should carry the requested k (`<= 3`) rather than a
            LIMIT, and no LATERAL should survive.
        """
        # Arrange
        sql = (
            "SELECT * FROM peaks "
            "CROSS JOIN LATERAL NEAREST(genes, reference := peaks.interval, k := 3) AS b"
        )

        # Act
        output = _generate_for_target(
            sql, tables_with_peaks_and_genes, DataFusionTarget()
        )

        # Assert
        assert "ROW_NUMBER(" in output.upper()
        assert "<= 3" in output
        assert "LATERAL" not in output.upper()

    def test_fallback_preserves_signed_and_max_distance(
        self, tables_with_peaks_and_genes
    ):
        """Test the DataFusion fallback keeps signed distance and the max filter.

        Given:
            A correlated NEAREST(signed := true, max_distance := 100000) on the
            DataFusion target (the decorrelated window path).
        When:
            Transpiling.
        Then:
            The window form should retain the signed distance CASE (`ELSE -(`)
            and the `<= 100000` max-distance filter — the `signed`/`max_distance`
            branches of the shared distance/filter builder on the fallback path.
        """
        # Arrange
        sql = (
            "SELECT * FROM peaks CROSS JOIN LATERAL NEAREST("
            "genes, reference := peaks.interval, signed := true, "
            "max_distance := 100000) AS b"
        )

        # Act
        output = _generate_for_target(
            sql, tables_with_peaks_and_genes, DataFusionTarget()
        )

        # Assert
        assert "ROW_NUMBER(" in output.upper()
        assert "LATERAL" not in output.upper()
        assert "ELSE -(" in output
        assert "<= 100000" in output


class TestNearestUnaliasedCorrelatedFallback:
    """The fallback synthesizes a LATERAL alias when the user omits one (B3)."""

    def test_expand_nearest_should_synthesize_alias_when_correlated_lateral_unaliased_on_nonlateral_target(  # noqa: E501
        self, tables_with_peaks_and_genes
    ):
        """Test an unaliased correlated NEAREST transpiles on a non-LATERAL target.

        Given:
            A correlated NEAREST whose surrounding CROSS JOIN LATERAL carries no
            table alias — legitimate GIQL that transpiles fine on lateral-capable
            engines — targeted at DataFusion (``supports_lateral`` is False), which
            takes the decorrelated window fallback.
        When:
            Transpiling for the DataFusion target.
        Then:
            It should not raise: the fallback synthesizes an alias via
            ``ctx.alias()`` instead of asserting one is present, and emits the
            decorrelated window form (a synthesized ``__giql_x_`` alias, no
            LATERAL).
        """
        # Arrange
        sql = (
            "SELECT a.start FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1)"
        )

        # Act
        output = _generate_for_target(
            sql, tables_with_peaks_and_genes, DataFusionTarget()
        )

        # Assert
        assert "ROW_NUMBER(" in output.upper()
        assert "LATERAL" not in output.upper()
        assert "__giql_x_" in output

    def test_expand_nearest_should_not_assert_on_unaliased_lateral_under_O(self):
        """Test the unaliased correlated fallback survives ``python -O``.

        Given:
            A fresh ``python -O`` interpreter (asserts stripped), in which an
            asserted alias precondition would degrade to a ``NoneType`` deref
            rather than a clear error.
        When:
            Transpiling an unaliased correlated NEAREST for the DataFusion dialect.
        Then:
            It should transpile without raising and emit the window fallback,
            proving the alias is synthesized (not asserted) so the optimized
            interpreter cannot strip the guard into an opaque crash.
        """
        # Arrange
        import subprocess
        import sys

        code = (
            "from giql import transpile, Table; "
            "sql = transpile("
            "'SELECT a.start FROM peaks a CROSS JOIN LATERAL "
            "NEAREST(genes, reference := a.interval, k := 1)', "
            "tables=[Table('peaks'), Table('genes')], dialect='datafusion'); "
            "assert 'ROW_NUMBER(' in sql.upper(), sql; "
            "assert 'LATERAL' not in sql.upper(), sql; "
            "print('ok')"
        )

        # Act
        result = subprocess.run(
            [sys.executable, "-O", "-c", code],
            capture_output=True,
            text=True,
        )

        # Assert
        assert result.returncode == 0, result.stderr
        assert result.stdout.strip() == "ok"


class TestNearestFallbackDetachContract:
    """The fallback detaches its NEAREST node so the pass's replace is a no-op (A10)."""

    def test_fallback_detaches_node_and_rewritten_join_survives_pass(
        self, tables_with_peaks_and_genes
    ):
        """Test the fallback detaches the NEAREST subtree and leaves the rewritten join.

        Given:
            A correlated NEAREST on the DataFusion target, captured before the
            ExpandOperators pass runs.
        When:
            Running the pass (which dispatches to the decorrelated fallback).
        Then:
            The original NEAREST node should be detached from the returned tree —
            its surrounding LATERAL is swapped out, so it no longer reaches the
            result root, making the pass's own ``node.replace`` a no-op — and the
            rewritten plain JOIN (no LATERAL, ROW_NUMBER window present) should
            survive in the returned tree.
        """
        # Arrange
        from giql.expressions import GIQLNearest

        sql = (
            "SELECT a.start FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) AS b"
        )
        ast = parse_one(sql, dialect=GIQLDialect)
        ast = resolve_operator_refs(ast, tables_with_peaks_and_genes)
        ast = canonicalize_coordinates(ast)
        nearest = ast.find(GIQLNearest)
        assert nearest is not None and nearest.root() is ast

        # Act
        result = ExpandOperators(
            DataFusionTarget(), tables_with_peaks_and_genes
        ).transform(ast)

        # Assert
        # The fallback swaps out the LATERAL holding the NEAREST, so the node's
        # subtree is detached: it no longer reaches the result root, which is what
        # makes the pass's ``node.replace`` a no-op.
        assert nearest.root() is not result
        assert not list(result.find_all(GIQLNearest))
        output = result.sql()
        assert "LATERAL" not in output.upper()
        assert "ROW_NUMBER(" in output.upper()


#: The reserved rank/key columns the decorrelated fallback exposes on its join and
#: that the statement finalizer must project away from a surfacing star (#160).
_RESERVED_FALLBACK_COLUMNS = (
    "__giql_x_rk_chrom",
    "__giql_x_rk_start",
    "__giql_x_rk_end",
    "__giql_x_rn",
)

#: The extra reserved key the fallback exposes in the stranded case (strand joins
#: the partition/reference key), on top of ``_RESERVED_FALLBACK_COLUMNS``.
_STRAND_RESERVED_COLUMN = "__giql_x_rk_strand"


def _wrapper_except_names(output: str) -> set[str] | None:
    """Return the excepted column names of a top-level ``SELECT * EXCEPT (...)``.

    Extracts the names from the leading ``SELECT * EXCEPT (...)`` clause so the wrap
    can be asserted by exact set rather than substring membership (which would pass
    with an extra, missing, or misquoted column). Returns ``None`` when *output* is
    not such a wrapper. The clause is parsed by string slice because sqlglot's
    default parser does not round-trip a ``* EXCEPT`` list back into a ``Star``
    except-set.
    """
    prefix = "SELECT * EXCEPT ("
    if not output.startswith(prefix):
        return None
    inner = output[len(prefix) :]
    names = inner[: inner.index(")")].split(",")
    return {name.strip().strip('"') for name in names}


class TestNearestFallbackReservedColumnProjection:
    """The statement finalizer that hides the fallback's reserved columns (#160).

    The DataFusion decorrelated fallback must expose reserved rank/key columns on
    the rewritten join, which a ``SELECT *`` / ``b.*`` would leak. A registered
    statement finalizer wraps the enclosing ``SELECT`` in ``SELECT * EXCEPT (...)``
    when — and only when — a star projection surfaces them. These assert the wrap
    shape without an engine; the cross-target result identity is covered by the
    oracle.
    """

    def test_expand_nearest_should_wrap_star_except_reserved_when_b_star_on_datafusion(
        self, tables_with_peaks_and_genes
    ):
        """Test a correlated ``SELECT b.*`` fallback is wrapped in ``* EXCEPT``.

        Given:
            A correlated NEAREST projected as ``SELECT b.*`` on the DataFusion
            target (no LATERAL plan, so the decorrelated fallback runs).
        When:
            Transpiling.
        Then:
            The output should be a top-level ``SELECT * EXCEPT (...)`` wrapper whose
            ``EXCEPT`` name-set equals exactly the reserved columns, hiding them
            from output.
        """
        # Arrange
        sql = (
            "SELECT b.* FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b"
        )

        # Act
        output = _generate_for_target(
            sql, tables_with_peaks_and_genes, DataFusionTarget()
        )

        # Assert
        assert _wrapper_except_names(output) == set(_RESERVED_FALLBACK_COLUMNS)

    def test_expand_nearest_should_wrap_star_except_reserved_when_unqualified_star_on_datafusion(  # noqa: E501
        self, tables_with_peaks_and_genes
    ):
        """Test a correlated unqualified ``SELECT *`` fallback is wrapped.

        Given:
            A correlated NEAREST projected as unqualified ``SELECT *`` on the
            DataFusion target — the star pulls the join relation's reserved
            columns.
        When:
            Transpiling.
        Then:
            The output should be a top-level ``SELECT * EXCEPT (...)`` wrapper whose
            ``EXCEPT`` name-set equals exactly the reserved columns.
        """
        # Arrange
        sql = (
            "SELECT * FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b"
        )

        # Act
        output = _generate_for_target(
            sql, tables_with_peaks_and_genes, DataFusionTarget()
        )

        # Assert
        assert _wrapper_except_names(output) == set(_RESERVED_FALLBACK_COLUMNS)

    def test_expand_nearest_should_not_wrap_when_projection_is_explicit_on_datafusion(
        self, tables_with_peaks_and_genes
    ):
        """Test an explicitly-projected fallback query is not wrapped.

        Given:
            A correlated NEAREST projecting named columns (``SELECT a.start,
            b.start``) on the DataFusion target — no reserved column surfaces.
        When:
            Transpiling.
        Then:
            No ``* EXCEPT`` wrapper should be added (wrapping absent columns would
            fail at engine runtime).
        """
        # Arrange
        sql = (
            "SELECT a.start, b.start FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b"
        )

        # Act
        output = _generate_for_target(
            sql, tables_with_peaks_and_genes, DataFusionTarget()
        )

        # Assert
        assert "EXCEPT" not in output

    def test_expand_nearest_should_not_wrap_when_star_nests_over_explicit_inner_on_datafusion(  # noqa: E501
        self, tables_with_peaks_and_genes
    ):
        """Test an outer ``SELECT *`` over an explicit inner fallback is not wrapped.

        Given:
            A correlated NEAREST whose reserved columns stay inside an inner
            subquery that projects only named columns, with an outer ``SELECT *``
            over it on the DataFusion target.
        When:
            Transpiling.
        Then:
            No wrapper should be added anywhere — the inner select does not surface
            the reserved columns, so the outer ``*`` is clean (the false-positive
            regression guard).
        """
        # Arrange
        sql = (
            "SELECT * FROM (SELECT a.start AS s FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b) sub"
        )

        # Act
        output = _generate_for_target(
            sql, tables_with_peaks_and_genes, DataFusionTarget()
        )

        # Assert
        assert "EXCEPT" not in output

    def test_expand_nearest_should_wrap_inner_select_when_star_nests_over_b_star_on_datafusion(  # noqa: E501
        self, tables_with_peaks_and_genes
    ):
        """Test the wrapper lands on the inner select when a star nests over ``b.*``.

        Given:
            A correlated NEAREST projected as ``SELECT b.*`` inside a subquery, with
            an outer ``SELECT *`` over it on the DataFusion target.
        When:
            Transpiling.
        Then:
            The ``* EXCEPT`` wrapper should land on the *inner* select (the one that
            surfaces the reserved columns), leaving the outer ``SELECT *`` plain.
        """
        # Arrange
        sql = (
            "SELECT * FROM (SELECT b.* FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b) sub"
        )

        # Act
        output = _generate_for_target(
            sql, tables_with_peaks_and_genes, DataFusionTarget()
        )

        # Assert
        assert not output.startswith("SELECT * EXCEPT (")
        inner = output[output.index("SELECT * EXCEPT (") :]
        assert _wrapper_except_names(inner) == set(_RESERVED_FALLBACK_COLUMNS)

    def test_expand_nearest_should_not_wrap_on_lateral_capable_target(
        self, tables_with_peaks_and_genes
    ):
        """Test the LATERAL form is never wrapped (no reserved columns leak).

        Given:
            The same correlated ``SELECT b.*`` NEAREST on the generic
            (lateral-capable) target, which emits the LATERAL form with no reserved
            columns.
        When:
            Transpiling.
        Then:
            No ``* EXCEPT`` wrapper should be added — the finalizer is registered
            only by the decorrelated fallback path.
        """
        # Arrange
        sql = (
            "SELECT b.* FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b"
        )

        # Act
        output = _generate(sql, tables_with_peaks_and_genes)

        # Assert
        assert "EXCEPT" not in output
        assert "LATERAL" in output.upper()

    def test_expand_nearest_should_except_strand_key_when_stranded_b_star_on_datafusion(
        self, tables_with_peaks_and_genes
    ):
        """Test a stranded correlated ``SELECT b.*`` fallback excepts the strand key.

        Given:
            A **stranded** correlated NEAREST projected as ``SELECT b.*`` on the
            DataFusion target — the fallback's reference key includes strand, so it
            exposes an extra reserved ``__giql_x_rk_strand`` column on the join.
        When:
            Transpiling.
        Then:
            The ``* EXCEPT`` wrapper's name-set should equal the base reserved
            columns plus ``__giql_x_rk_strand`` — the strand key is hidden too.
        """
        # Arrange
        sql = (
            "SELECT b.* FROM peaks a "
            "CROSS JOIN LATERAL "
            "NEAREST(genes, reference := a.interval, k := 1, stranded := true) b"
        )

        # Act
        output = _generate_for_target(
            sql, tables_with_peaks_and_genes, DataFusionTarget()
        )

        # Assert
        assert _wrapper_except_names(output) == set(_RESERVED_FALLBACK_COLUMNS) | {
            _STRAND_RESERVED_COLUMN
        }

    def test_expand_nearest_should_not_wrap_when_projection_is_outer_relation_star_on_datafusion(  # noqa: E501
        self, tables_with_peaks_and_genes
    ):
        """Test a star over the outer relation only is not wrapped.

        Given:
            A correlated NEAREST projecting ``SELECT a.*`` — a star over the outer
            relation, which carries none of the fallback's reserved columns — on the
            DataFusion target.
        When:
            Transpiling.
        Then:
            No ``* EXCEPT`` wrapper should be added; ``a.*`` surfaces no reserved
            columns, and wrapping absent columns would fail at engine runtime.
        """
        # Arrange
        sql = (
            "SELECT a.* FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b"
        )

        # Act
        output = _generate_for_target(
            sql, tables_with_peaks_and_genes, DataFusionTarget()
        )

        # Assert
        assert "EXCEPT" not in output

    def test_expand_nearest_should_not_double_wrap_when_projection_already_excepts_reserved_on_datafusion(  # noqa: E501
        self, tables_with_peaks_and_genes
    ):
        """Test a star already excepting the reserved columns is not re-wrapped.

        Given:
            A correlated NEAREST whose user projection is already
            ``SELECT * EXCEPT (<all reserved columns>)`` on the DataFusion target.
        When:
            Transpiling.
        Then:
            No second wrapper should be added — the star already excepts every
            reserved name, so the finalizer's idempotency guard treats it as not
            surfacing (exactly one ``EXCEPT`` clause in the output).
        """
        # Arrange
        excepted = ", ".join(f'"{name}"' for name in _RESERVED_FALLBACK_COLUMNS)
        sql = (
            f"SELECT * EXCEPT ({excepted}) FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b"
        )

        # Act
        output = _generate_for_target(
            sql, tables_with_peaks_and_genes, DataFusionTarget()
        )

        # Assert
        assert output.count("EXCEPT") == 1
        assert _wrapper_except_names(output) == set(_RESERVED_FALLBACK_COLUMNS)

    def test_expand_nearest_should_wrap_when_b_star_qualifier_is_mixed_case_on_datafusion(  # noqa: E501
        self, tables_with_peaks_and_genes
    ):
        """Test an unquoted mixed-case ``B.*`` over alias ``b`` is still wrapped.

        Given:
            A correlated NEAREST projected as ``SELECT B.*`` — an unquoted
            qualifier whose case differs from the lateral alias ``b`` — on the
            DataFusion target. Engines fold unquoted identifiers, so ``B`` binds to
            the same relation as ``b``.
        When:
            Transpiling.
        Then:
            The output should be a top-level ``SELECT * EXCEPT (...)`` wrapper whose
            ``EXCEPT`` name-set equals the reserved columns — a case-variant
            qualifier must not slip the reserved columns through.
        """
        # Arrange
        sql = (
            "SELECT B.* FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b"
        )

        # Act
        output = _generate_for_target(
            sql, tables_with_peaks_and_genes, DataFusionTarget()
        )

        # Assert
        assert _wrapper_except_names(output) == set(_RESERVED_FALLBACK_COLUMNS)

    def test_expand_nearest_should_not_mint_wrapper_alias_when_no_wrap_on_datafusion(
        self, tables_with_peaks_and_genes
    ):
        """Test the no-wrap path does not advance the run's alias sequence.

        Given:
            The same correlated NEAREST projected as a wrapping ``SELECT b.*`` and
            as a non-wrapping explicit projection on the DataFusion target — the
            wrapper alias is minted lazily, only when a wrapper is emitted.
        When:
            Transpiling both.
        Then:
            The explicit (no-wrap) projection should mint exactly one fewer
            ``__giql_x_<n>`` alias than the wrapped one — the wrapper's derived-table
            alias is never minted, so the shared alias sequence is not advanced.
        """
        # Arrange
        join = (
            "FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b"
        )

        # Act
        wrapped = _generate_for_target(
            f"SELECT b.* {join}", tables_with_peaks_and_genes, DataFusionTarget()
        )
        explicit = _generate_for_target(
            f"SELECT a.start, b.start {join}",
            tables_with_peaks_and_genes,
            DataFusionTarget(),
        )

        # Assert
        wrapped_aliases = set(re.findall(r"__giql_x_\d+", wrapped))
        explicit_aliases = set(re.findall(r"__giql_x_\d+", explicit))
        assert explicit_aliases < wrapped_aliases
        assert len(wrapped_aliases) - len(explicit_aliases) == 1

    def test_expand_cluster_over_nearest_should_not_crash_and_leaks_as_documented_residual_on_datafusion(  # noqa: E501
        self, tables_with_peaks_and_genes
    ):
        """Test a ``SELECT *`` CLUSTER wrapping a correlated NEAREST is a documented leak.

        Given:
            A correlated NEAREST fallback nested under a ``SELECT *``-projecting
            CLUSTER on the DataFusion target. The CLUSTER's copy+transplant detaches
            the join the finalizer captured (``parent_select`` becomes ``None``), so
            the finalizer no-ops and the outer ``SELECT *`` re-surfaces the reserved
            columns.
        When:
            Transpiling.
        Then:
            It should transpile without error and leak the reserved columns — the
            documented, unwrapped residual tracked by #172 — exercising the
            finalizer's detached-join (``select is None``) no-op branch. This pins
            the residual so a future change that alters it is caught.
        """
        # Arrange
        sql = (
            "SELECT *, CLUSTER(interval) AS cid FROM ("
            "SELECT b.* FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b"
            ") sub"
        )

        # Act
        output = _generate_for_target(
            sql, tables_with_peaks_and_genes, DataFusionTarget()
        )

        # Assert
        assert "EXCEPT" not in output
        assert "__giql_x_rk_chrom" in output  # residual leak (#172), not wrapped
