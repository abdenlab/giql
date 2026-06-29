"""Direct unit tests for the spatial / set predicate expanders (#141).

These call ``expand_intersects`` / ``expand_contains`` / ``expand_within`` /
``expand_spatial_set`` directly with a hand-built :class:`ExpansionContext`,
characterizing each dispatch branch (column-to-column vs literal range; CONTAINS
point vs range; ANY/OR vs ALL/AND) and pinning the chosen error messages on
invalid input. They sit outside ``tests/test_expander.py`` so they do not touch
that file's shared, operator-agnostic fixture/infra region.
"""

import pytest
from sqlglot import exp
from sqlglot import parse_one

from giql.dialect import GIQLDialect
from giql.expander import REGISTRY
from giql.expander import ExpansionContext
from giql.expanders.intersects import expand_contains
from giql.expanders.intersects import expand_intersects
from giql.expanders.intersects import expand_spatial_set
from giql.expanders.intersects import expand_within
from giql.expressions import Contains
from giql.expressions import Intersects
from giql.expressions import SpatialSetPredicate
from giql.expressions import Within
from giql.resolver import OperatorResolution
from giql.resolver import ResolvedColumn
from giql.table import Tables
from giql.targets import DataFusionTarget
from giql.targets import GenericTarget
from giql.transpile import transpile

_LEFT = ResolvedColumn(
    chrom='a."chrom"', start='a."start"', end='a."end"', strand=None, table=None
)
_RIGHT = ResolvedColumn(
    chrom='b."chrom"', start='b."start"', end='b."end"', strand=None, table=None
)
_OPERATOR_TYPES = (Intersects, Contains, Within, SpatialSetPredicate)


def _context(query: str, columns: dict[str, ResolvedColumn]) -> tuple:
    """Find the spatial operator in *query* and build a context with *columns*."""
    root = parse_one(query, dialect=GIQLDialect)
    node = next(n for n in root.walk() if isinstance(n, _OPERATOR_TYPES))
    resolution = OperatorResolution(
        operator=type(node).__name__, slots={}, columns=columns
    )
    ctx = ExpansionContext(node, resolution, GenericTarget(), Tables())
    return node, ctx


def _sql(expression: exp.Expression) -> str:
    """Serialize a built expression through the GIQL dialect."""
    return expression.sql(dialect=GIQLDialect)


class TestSpatialExpanders:
    """Direct expansion of the spatial / set predicate expanders (#141)."""

    def test_expand_intersects_should_build_overlap_predicate_when_literal_range(self):
        """Test that a literal-range INTERSECTS expands to the overlap predicate.

        Given:
            An INTERSECTS node whose right operand is a literal range and a
            context resolving only the left column.
        When:
            Expanding it.
        Then:
            It should build the overlap boolean (chrom = lit AND start < end2 AND
            end > start2) with the right operand as numeric literals.
        """
        # Arrange
        node, ctx = _context(
            "SELECT * FROM a WHERE interval INTERSECTS 'chr1:1000-2000'",
            {"this": _LEFT},
        )

        # Act
        result = expand_intersects(node, ctx)

        # Assert
        assert _sql(result) == (
            '(a."chrom" = \'chr1\' AND a."start" < 2000 AND a."end" > 1000)'
        )

    def test_expand_intersects_should_build_join_predicate_when_column_to_column(self):
        """Test that a column-to-column INTERSECTS expands to a join predicate.

        Given:
            An INTERSECTS node with a resolved right *column* (the dispatch keys on
            ctx.resolution.column("expression")).
        When:
            Expanding it.
        Then:
            It should compare the two columns' endpoints rather than literals.
        """
        # Arrange
        node, ctx = _context(
            "SELECT * FROM a JOIN b ON a.interval INTERSECTS b.interval",
            {"this": _LEFT, "expression": _RIGHT},
        )

        # Act
        result = expand_intersects(node, ctx)

        # Assert
        assert _sql(result) == (
            '(a."chrom" = b."chrom" AND a."start" < b."end" AND a."end" > b."start")'
        )

    def test_expand_contains_should_build_point_predicate_when_single_base(self):
        """Test that a single-base CONTAINS expands to the point-containment form.

        Given:
            A CONTAINS node whose literal range is a single base (end == start+1).
        When:
            Expanding it.
        Then:
            It should use the point form (start <= point AND end > point), not the
            range form.
        """
        # Arrange
        node, ctx = _context(
            "SELECT * FROM a WHERE interval CONTAINS 'chr1:1000'", {"this": _LEFT}
        )

        # Act
        result = expand_contains(node, ctx)

        # Assert
        assert _sql(result) == (
            '(a."chrom" = \'chr1\' AND a."start" <= 1000 AND a."end" > 1000)'
        )

    def test_expand_contains_should_build_range_predicate_when_multi_base(self):
        """Test that a multi-base CONTAINS expands to the range-containment form.

        Given:
            A CONTAINS node whose literal range spans more than one base.
        When:
            Expanding it.
        Then:
            It should use the range form (start <= start2 AND end >= end2).
        """
        # Arrange
        node, ctx = _context(
            "SELECT * FROM a WHERE interval CONTAINS 'chr1:1000-2000'",
            {"this": _LEFT},
        )

        # Act
        result = expand_contains(node, ctx)

        # Assert
        assert _sql(result) == (
            '(a."chrom" = \'chr1\' AND a."start" <= 1000 AND a."end" >= 2000)'
        )

    def test_expand_contains_should_build_join_predicate_when_column_to_column(self):
        """Test that a column-to-column CONTAINS expands to the containment join.

        Given:
            A CONTAINS node with a resolved right *column* (the dispatch keys on
            ctx.resolution.column("expression")).
        When:
            Expanding it.
        Then:
            It should build the left-contains-right predicate (start1 <= start2
            AND end1 >= end2) comparing the two columns' endpoints.
        """
        # Arrange
        node, ctx = _context(
            "SELECT * FROM a JOIN b ON a.interval CONTAINS b.interval",
            {"this": _LEFT, "expression": _RIGHT},
        )

        # Act
        result = expand_contains(node, ctx)

        # Assert
        assert _sql(result) == (
            '(a."chrom" = b."chrom" AND a."start" <= b."start" AND a."end" >= b."end")'
        )

    def test_expand_within_should_build_containment_predicate_when_literal_range(self):
        """Test that WITHIN expands to the left-within-right containment form.

        Given:
            A WITHIN node with a literal range.
        When:
            Expanding it.
        Then:
            It should build start >= start2 AND end <= end2.
        """
        # Arrange
        node, ctx = _context(
            "SELECT * FROM a WHERE interval WITHIN 'chr1:1000-2000'", {"this": _LEFT}
        )

        # Act
        result = expand_within(node, ctx)

        # Assert
        assert _sql(result) == (
            '(a."chrom" = \'chr1\' AND a."start" >= 1000 AND a."end" <= 2000)'
        )

    def test_expand_within_should_build_join_predicate_when_column_to_column(self):
        """Test that a column-to-column WITHIN expands to the within join predicate.

        Given:
            A WITHIN node with a resolved right *column* (the dispatch keys on
            ctx.resolution.column("expression")).
        When:
            Expanding it.
        Then:
            It should build the left-within-right predicate (start1 >= start2 AND
            end1 <= end2) comparing the two columns' endpoints.
        """
        # Arrange
        node, ctx = _context(
            "SELECT * FROM a JOIN b ON a.interval WITHIN b.interval",
            {"this": _LEFT, "expression": _RIGHT},
        )

        # Act
        result = expand_within(node, ctx)

        # Assert
        assert _sql(result) == (
            '(a."chrom" = b."chrom" AND a."start" >= b."start" AND a."end" <= b."end")'
        )

    def test_expand_spatial_set_should_or_combine_conditions_when_any(self):
        """Test that an ANY set predicate OR-combines its per-range conditions.

        Given:
            An INTERSECTS ANY node over two literal ranges.
        When:
            Expanding it.
        Then:
            The two per-range overlap predicates should be OR-combined inside one
            outer paren.
        """
        # Arrange
        node, ctx = _context(
            "SELECT * FROM a WHERE interval "
            "INTERSECTS ANY ('chr1:1-100', 'chr1:200-300')",
            {"this": _LEFT},
        )

        # Act
        result = expand_spatial_set(node, ctx)

        # Assert
        assert _sql(result) == (
            '((a."chrom" = \'chr1\' AND a."start" < 100 AND a."end" > 1) OR '
            '(a."chrom" = \'chr1\' AND a."start" < 300 AND a."end" > 200))'
        )

    def test_expand_spatial_set_should_and_combine_conditions_when_all(self):
        """Test that an ALL set predicate AND-combines its per-range conditions.

        Given:
            An INTERSECTS ALL node over two literal ranges.
        When:
            Expanding it.
        Then:
            The two per-range overlap predicates should be AND-combined inside one
            outer paren.
        """
        # Arrange
        node, ctx = _context(
            "SELECT * FROM a WHERE interval "
            "INTERSECTS ALL ('chr1:1-100', 'chr1:200-300')",
            {"this": _LEFT},
        )

        # Act
        result = expand_spatial_set(node, ctx)

        # Assert
        assert _sql(result) == (
            '((a."chrom" = \'chr1\' AND a."start" < 100 AND a."end" > 1) AND '
            '(a."chrom" = \'chr1\' AND a."start" < 300 AND a."end" > 200))'
        )


@pytest.fixture
def isolated_registry():
    """Snapshot/restore the process REGISTRY so a test can register an override."""
    saved = REGISTRY.snapshot()
    try:
        yield REGISTRY
    finally:
        REGISTRY.restore(saved)


class TestBinnedTargetOverrideDeferral:
    """A target-specific Intersects override defers the binned join rewrite (#141)."""

    def test_transpile_should_skip_binned_rewrite_when_target_override(
        self, isolated_registry
    ):
        """Test that a (target, Intersects) override bypasses the binned transformer.

        Given:
            A column-to-column INTERSECTS join on the generic binned path
            (dialect='datafusion') with a (DataFusionTarget, Intersects) override
            registered.
        When:
            Transpiling.
        Then:
            The override's sentinel reaches the SQL and no binned equi-join
            artifact is emitted — the override takes over the join rewrite that the
            built-in binned transformer would otherwise perform.
        """
        # Arrange
        isolated_registry.register(
            DataFusionTarget(),
            Intersects,
            lambda n, c: exp.column("BINNED_OVERRIDE_SENTINEL"),
        )
        query = (
            "SELECT a.start FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval"
        )

        # Act
        sql = transpile(query, tables=["peaks", "genes"], dialect="datafusion")

        # Assert
        assert "BINNED_OVERRIDE_SENTINEL" in sql
        assert "_bins" not in sql

    def test_transpile_should_reject_bin_size_when_target_override(
        self, isolated_registry
    ):
        """Test that bin size is rejected under a binned-target Intersects override.

        Given:
            A (DataFusionTarget, Intersects) override registered.
        When:
            Transpiling with intersects_bin_size set (which only configures the
            built-in binned transformer the override supersedes).
        Then:
            transpile() raises ValueError rather than silently dropping the bin
            size, parallel to the iejoin rejection.
        """
        # Arrange
        isolated_registry.register(
            DataFusionTarget(),
            Intersects,
            lambda n, c: exp.column("BINNED_OVERRIDE_SENTINEL"),
        )
        query = (
            "SELECT a.start FROM peaks a "
            "JOIN genes b ON a.interval INTERSECTS b.interval"
        )

        # Act & assert
        with pytest.raises(ValueError, match=r"intersects_bin_size has no effect"):
            transpile(
                query,
                tables=["peaks", "genes"],
                dialect="datafusion",
                intersects_bin_size=5000,
            )


class TestSpatialExpanderErrors:
    """Characterization tests pinning the chosen error messages on invalid input."""

    def test_expand_intersects_should_wrap_parse_error_when_invalid_range(self):
        """Test that an unparseable literal range raises the wrapped diagnostic.

        Given:
            An INTERSECTS node whose literal range string cannot be parsed.
        When:
            Expanding it.
        Then:
            It should raise ValueError with the historical "Could not parse
            genomic range" wrapper, chained from the underlying parse error.
        """
        # Arrange
        node, ctx = _context(
            "SELECT * FROM a WHERE interval INTERSECTS 'invalid'", {"this": _LEFT}
        )

        # Act & assert
        with pytest.raises(ValueError, match=r"Could not parse genomic range") as exc:
            expand_intersects(node, ctx)
        assert exc.value.__cause__ is not None

    def test_expand_intersects_should_raise_invariant_when_left_unresolved(self):
        """Test that a missing left-operand resolution raises the invariant error.

        Given:
            An INTERSECTS node whose context resolved no "this" column (pass 1 did
            not run).
        When:
            Expanding it.
        Then:
            It should raise ValueError naming the unresolved operand and pointing
            at the ResolveOperatorRefs pass.
        """
        # Arrange
        node, ctx = _context(
            "SELECT * FROM a WHERE interval INTERSECTS 'chr1:1-100'", {}
        )

        # Act & assert
        with pytest.raises(
            ValueError, match=r"Spatial predicate operand 'this' was not resolved"
        ):
            expand_intersects(node, ctx)

    def test_expand_spatial_set_should_not_wrap_parse_error_when_invalid_range(self):
        """Test that a set-predicate bad range surfaces the raw parser error.

        Given:
            An INTERSECTS ANY node with one unparseable range.
        When:
            Expanding it.
        Then:
            The raw RangeParser ValueError propagates *unwrapped* — the set-
            predicate path does NOT apply the "Could not parse genomic range"
            wrapper the single-operand path does. This pins the current
            (pre-existing) asymmetry so any future unification is a conscious
            change, not an accident.
        """
        # Arrange
        node, ctx = _context(
            "SELECT * FROM a WHERE interval INTERSECTS ANY ('bad', 'chr1:1-2')",
            {"this": _LEFT},
        )

        # Act & assert
        with pytest.raises(ValueError, match=r"Invalid genomic range format") as exc:
            expand_spatial_set(node, ctx)
        assert "Could not parse genomic range" not in str(exc.value)
