"""Transpilation tests for NEAREST operator SQL generation.

Tests verify that NEAREST() is correctly transpiled to SQL
(LATERAL joins for correlated queries, ORDER BY + LIMIT for standalone).
"""

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
