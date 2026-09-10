"""Unit tests for the datafusion-bio spatial-predicate overrides (#77).

These exercise :mod:`giql.expanders.intersects_datafusion_bio` directly with a
hand-built :class:`ExpansionContext` (the closed overlap form and the
literal-range deferral) and through :func:`transpile` (the join-kind gate, the
count-overlaps zero-fill rewrite, the bare ``COUNT(*)`` carried projection, and
the CONTAINS / WITHIN gating). They live in their own module so
``tests/expanders/test_intersects.py``'s shared fixture region stays untouched.
"""

import pytest
from sqlglot import exp
from sqlglot import parse_one

from giql.dialect import GIQLDialect
from giql.expander import REGISTRY
from giql.expander import ExpansionContext
from giql.expanders.intersects import SPATIAL_PREDICATE_META
from giql.expanders.intersects import expand_intersects
from giql.expanders.intersects_datafusion_bio import expand_intersects_datafusion_bio
from giql.expressions import Intersects
from giql.resolver import OperatorResolution
from giql.resolver import ResolvedColumn
from giql.table import Table
from giql.table import Tables
from giql.targets import DataFusionBioTarget
from giql.targets import GenericTarget
from giql.transpile import transpile

_LEFT = ResolvedColumn(
    chrom='a."chrom"', start='a."start"', end='a."end"', strand=None, table=None
)
_RIGHT = ResolvedColumn(
    chrom='b."chrom"', start='b."start"', end='b."end"', strand=None, table=None
)

_CLOSED_FORM = (
    '(a."chrom" = b."chrom" AND a."start" <= b."end" - 1 AND a."end" - 1 >= b."start")'
)
_DEFEATED_FORM = (
    '(a."chrom" = b."chrom" AND a."start" <= b."end" - 1 AND a."end" - 1 >= b."start"'
    ' AND a."start" < b."end")'
)
_STRICT_FORM = '(a."chrom" = b."chrom" AND a."start" < b."end" AND a."end" > b."start")'
_TABLES = ["peaks", "genes"]


def _context(query: str, columns: dict[str, ResolvedColumn]) -> tuple:
    """Find the INTERSECTS in *query* and build a datafusion-bio context for it."""
    root = parse_one(query, dialect=GIQLDialect)
    node = next(n for n in root.walk() if isinstance(n, Intersects))
    resolution = OperatorResolution(
        operator=type(node).__name__, slots={}, columns=columns
    )
    ctx = ExpansionContext(node, resolution, DataFusionBioTarget(), Tables())
    return node, ctx


def _sql(expression: exp.Expression) -> str:
    """Serialize a built expression through the GIQL dialect."""
    return expression.sql(dialect=GIQLDialect)


def _bio(query: str, tables=None) -> str:
    """Transpile *query* for the datafusion-bio target."""
    return transpile(query, tables=tables or _TABLES, dialect="datafusion-bio")


@pytest.fixture
def isolated_registry():
    """Snapshot/restore the process REGISTRY so a test can register an override."""
    saved = REGISTRY.snapshot()
    try:
        yield REGISTRY
    finally:
        REGISTRY.restore(saved)


class TestExpandIntersectsDatafusionBio:
    """Direct expansion through the datafusion-bio override."""

    def test_expand_intersects_datafusion_bio_should_emit_closed_form_when_inner_join(
        self,
    ):
        """Test that an inner-join column-to-column INTERSECTS emits the closed form.

        Given:
            An INTERSECTS node in a plain JOIN with both operands resolved to columns.
        When:
            Expanding it through the datafusion-bio override.
        Then:
            It should emit the non-strict ``start <= end - 1 AND end - 1 >= start``
            form polars-bio's interval-join rule accepts on Int64 coordinates, and
            tag it as a spatial predicate.
        """
        # Arrange
        node, ctx = _context(
            "SELECT a.start FROM a JOIN b ON a.interval INTERSECTS b.interval",
            {"this": _LEFT, "expression": _RIGHT},
        )

        # Act
        result = expand_intersects_datafusion_bio(node, ctx)

        # Assert
        assert _sql(result) == _CLOSED_FORM
        assert result.meta.get(SPATIAL_PREDICATE_META) is True

    def test_expand_intersects_datafusion_bio_should_defer_to_generic_when_literal_range(  # noqa: E501
        self,
    ):
        """Test that a literal-range INTERSECTS is left to the generic expander.

        Given:
            An INTERSECTS node whose right operand is a literal range.
        When:
            Expanding it through the datafusion-bio override.
        Then:
            It should produce exactly what the generic expander produces, since a
            literal-range predicate is a filter the interval-join rule never sees.
        """
        # Arrange
        query = "SELECT * FROM a WHERE interval INTERSECTS 'chr1:1000-2000'"
        node, ctx = _context(query, {"this": _LEFT})
        generic_node, generic_ctx = _context(query, {"this": _LEFT})

        # Act
        result = expand_intersects_datafusion_bio(node, ctx)
        expected = expand_intersects(generic_node, generic_ctx)

        # Assert
        assert _sql(result) == _sql(expected)

    def test_expand_intersects_datafusion_bio_should_honor_generic_override_when_literal_range(  # noqa: E501
        self, isolated_registry
    ):
        """Test that the literal-range deferral goes through the registry.

        Given:
            A user override registered on the (GenericTarget, Intersects) slot.
        When:
            A literal-range INTERSECTS is expanded through the datafusion-bio override.
        Then:
            It should dispatch to the user override rather than the built-in, so
            the generic slot governs this target's declined shapes too.
        """
        # Arrange
        isolated_registry.register(
            GenericTarget(),
            Intersects,
            lambda n, c: exp.column("GENERIC_OVERRIDE_SENTINEL"),
        )
        node, ctx = _context(
            "SELECT * FROM a WHERE interval INTERSECTS 'chr1:1000-2000'",
            {"this": _LEFT},
        )

        # Act
        result = expand_intersects_datafusion_bio(node, ctx)

        # Assert
        assert _sql(result) == "GENERIC_OVERRIDE_SENTINEL"


class TestJoinKindGate:
    """Which join kinds keep the rule-eligible shape and which defeat the rule."""

    @pytest.mark.parametrize(
        "query",
        [
            "SELECT a.start FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval",
            "SELECT a.start FROM peaks a INNER JOIN genes b "
            "ON a.interval INTERSECTS b.interval",
            "SELECT a.start FROM peaks a, genes b "
            "WHERE a.interval INTERSECTS b.interval",
            "WITH x AS (SELECT a.start FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval) SELECT * FROM x",
            "SELECT * FROM (SELECT a.start FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval) t",
        ],
        ids=["join", "inner-join", "comma-where", "cte", "from-subquery"],
    )
    def test_transpile_should_emit_closed_form_when_inner_semantics(self, query):
        """Test that inner-semantics shapes keep the rule-eligible closed form.

        Given:
            A column-to-column INTERSECTS in a plain JOIN, a comma-join WHERE, a
            CTE body, or a FROM-clause subquery.
        When:
            Transpiling for datafusion-bio.
        Then:
            It should emit the closed overlap form with no defeating conjunct, so
            polars-bio plans the join through IntervalJoinExec.
        """
        # Act
        sql = _bio(query)

        # Assert
        assert _CLOSED_FORM in sql
        assert _DEFEATED_FORM not in sql

    @pytest.mark.parametrize(
        "query",
        [
            "SELECT a.start FROM peaks a LEFT JOIN genes b "
            "ON a.interval INTERSECTS b.interval",
            "SELECT a.start FROM peaks a RIGHT JOIN genes b "
            "ON a.interval INTERSECTS b.interval",
            "SELECT a.start FROM peaks a FULL OUTER JOIN genes b "
            "ON a.interval INTERSECTS b.interval",
            "SELECT a.start FROM peaks a WHERE EXISTS "
            "(SELECT 1 FROM genes b WHERE a.interval INTERSECTS b.interval)",
            "SELECT a.start FROM peaks a WHERE NOT EXISTS "
            "(SELECT 1 FROM genes b WHERE a.interval INTERSECTS b.interval)",
            "SELECT a.start, (SELECT COUNT(*) FROM genes b "
            "WHERE a.interval INTERSECTS b.interval) AS n FROM peaks a",
        ],
        ids=["left", "right", "full", "exists", "not-exists", "scalar-subquery"],
    )
    def test_transpile_should_emit_defeated_form_when_non_inner_semantics(self, query):
        """Test that outer, semi, anti and correlated shapes defeat the rule.

        Given:
            A column-to-column INTERSECTS under an outer join, an EXISTS / NOT
            EXISTS operand, or a correlated scalar subquery.
        When:
            Transpiling for datafusion-bio.
        Then:
            It should append the implied ``a.start < b.end`` comparison so the
            interval-join rule declines, since IntervalJoinExec returns the inner
            pairs for every join type upstream.
        """
        # Act
        sql = _bio(query)

        # Assert
        assert _DEFEATED_FORM in sql
        assert _STRICT_FORM not in sql

    @pytest.mark.parametrize(
        ("operator", "generic_form", "defeated_form"),
        [
            (
                "CONTAINS",
                '(a."chrom" = b."chrom" AND a."start" <= b."start" '
                'AND a."end" >= b."end")',
                '(a."chrom" = b."chrom" AND a."start" <= b."start" '
                'AND a."end" >= b."end" AND a."start" - 1 < b."start")',
            ),
            (
                "WITHIN",
                '(a."chrom" = b."chrom" AND a."start" >= b."start" '
                'AND a."end" <= b."end")',
                '(a."chrom" = b."chrom" AND a."start" >= b."start" '
                'AND a."end" <= b."end" AND a."start" + 1 > b."start")',
            ),
        ],
    )
    def test_transpile_should_gate_containment_by_join_kind(
        self, operator, generic_form, defeated_form
    ):
        """Test CONTAINS / WITHIN keep the generic form on inner joins only.

        Given:
            A column-to-column CONTAINS or WITHIN in an inner join and in a LEFT
            join.
        When:
            Transpiling both for datafusion and for datafusion-bio.
        Then:
            The inner join should match the vanilla DataFusion output exactly (the
            generic non-strict form is already rule-eligible), while the LEFT join
            should carry the implied third comparison that defeats the rule.
        """
        # Arrange
        inner = f"SELECT a.start FROM peaks a JOIN genes b ON a.interval {operator} b.interval"  # noqa: E501
        left = f"SELECT a.start FROM peaks a LEFT JOIN genes b ON a.interval {operator} b.interval"  # noqa: E501

        # Act
        inner_bio = _bio(inner)
        inner_vanilla = transpile(inner, tables=_TABLES, dialect="datafusion")
        left_bio = _bio(left)

        # Assert
        assert inner_bio == inner_vanilla
        assert generic_form in inner_bio
        assert defeated_form in left_bio


class TestCountOverlapsZeroFill:
    """The LEFT JOIN + COUNT(right col) + GROUP BY left-keys rewrite."""

    def test_transpile_should_zero_fill_inner_count_when_script_shape(self):
        """Test the cCRE insertion-matrix shape is rewritten to an inner count.

        Given:
            The downstream script's LEFT JOIN + COUNT(r.start) + GROUP BY query
            over a custom-column insertions table.
        When:
            Transpiling for datafusion-bio.
        Then:
            It should emit an inner interval-join count under a LEFT-joined
            zero-fill against the distinct left keys, with the closed overlap form
            on the inner join and the user's output names preserved.
        """
        # Arrange
        query = (
            "SELECT c.ccre_id, COUNT(r.start) AS n FROM ccres c "
            "LEFT JOIN insertions r ON c.interval INTERSECTS r.interval "
            "GROUP BY c.ccre_id"
        )
        tables = [
            "ccres",
            Table(
                "insertions",
                chrom_col="rname",
                start_col="start",
                end_col="end",
                strand_col=None,
            ),
        ]

        # Act
        sql = _bio(query, tables)

        # Assert
        assert sql == (
            "SELECT __giql_x_0.ccre_id AS ccre_id, COALESCE(__giql_x_1.n, 0) AS n "
            "FROM (SELECT DISTINCT c.ccre_id AS ccre_id FROM ccres AS c) AS __giql_x_0 "
            "LEFT JOIN (SELECT c.ccre_id AS ccre_id, COUNT(r.start) AS n "
            "FROM ccres AS c JOIN insertions AS r ON "
            '(c."chrom" = r."rname" AND c."start" <= r."end" - 1 '
            'AND c."end" - 1 >= r."start") GROUP BY c.ccre_id) AS __giql_x_1 '
            "ON __giql_x_0.ccre_id = __giql_x_1.ccre_id"
        )

    def test_transpile_should_preserve_projection_order_and_quoting_when_zero_filling(
        self,
    ):
        """Test the zero-fill keeps the SELECT list's order, aliases and quoting.

        Given:
            The count shape with the aggregate first, two keys including a quoted
            reserved word aliased by the user.
        When:
            Transpiling for datafusion-bio.
        Then:
            It should project the count first, then the keys under their aliases,
            and join on every key.
        """
        # Arrange
        query = (
            'SELECT COUNT(r.start) AS n, c.chrom, c."end" AS e FROM ccres c '
            "LEFT JOIN insertions r ON c.interval INTERSECTS r.interval "
            'GROUP BY c.chrom, c."end"'
        )

        # Act
        sql = _bio(query, ["ccres", "insertions"])

        # Assert
        assert sql.startswith(
            "SELECT COALESCE(__giql_x_1.n, 0) AS n, __giql_x_0.chrom AS chrom, "
            "__giql_x_0.e AS e FROM "
        )
        assert 'SELECT DISTINCT c.chrom AS chrom, c."end" AS e FROM ccres AS c' in sql
        assert sql.endswith(
            "ON __giql_x_0.chrom = __giql_x_1.chrom AND __giql_x_0.e = __giql_x_1.e"
        )

    @pytest.mark.parametrize(
        "query",
        [
            # COUNT(*) is not a right-column count.
            "SELECT c.chrom, COUNT(*) AS n FROM ccres c LEFT JOIN insertions r "
            "ON c.interval INTERSECTS r.interval GROUP BY c.chrom",
            # A WHERE residual.
            "SELECT c.chrom, COUNT(r.start) AS n FROM ccres c LEFT JOIN insertions r "
            "ON c.interval INTERSECTS r.interval WHERE c.chrom = 'chr1' "
            "GROUP BY c.chrom",
            # An inner join needs no zero-fill.
            "SELECT c.chrom, COUNT(r.start) AS n FROM ccres c JOIN insertions r "
            "ON c.interval INTERSECTS r.interval GROUP BY c.chrom",
        ],
        ids=["count-star", "where-residual", "inner"],
    )
    def test_transpile_should_not_zero_fill_when_shape_differs(self, query):
        """Test that only the exact count-overlaps shape is rewritten.

        Given:
            A near-miss of the count shape: a COUNT(*), a WHERE residual, or an
            inner join.
        When:
            Transpiling for datafusion-bio.
        Then:
            It should not introduce the zero-fill wrapper.
        """
        # Act
        sql = _bio(query, ["ccres", "insertions"])

        # Assert
        assert "COALESCE" not in sql
        assert "SELECT DISTINCT" not in sql


class TestBareCountStar:
    """The carried projection for ``SELECT COUNT(*)`` over an inner interval join."""

    @pytest.mark.parametrize(
        "query",
        [
            "SELECT COUNT(*) FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval",
            "SELECT COUNT(*) AS n FROM peaks a INNER JOIN genes b "
            "ON a.interval INTERSECTS b.interval",
            "SELECT COUNT(*) FROM peaks a JOIN genes b "
            "ON b.interval INTERSECTS a.interval",
        ],
        ids=["plain", "aliased", "operands-swapped"],
    )
    def test_transpile_should_carry_from_side_start_when_bare_count_star(self, query):
        """Test that a bare COUNT(*) over an inner join counts the FROM side's start.

        Given:
            ``SELECT COUNT(*)`` over an inner column-to-column INTERSECTS join with
            no other projected column, whichever way the operands are written.
        When:
            Transpiling for datafusion-bio.
        Then:
            It should count the FROM table's start column so the interval join
            projects a column, since an empty projection fails upstream.
        """
        # Act
        sql = _bio(query)

        # Assert
        assert 'COUNT(a."start")' in sql
        assert "COUNT(*)" not in sql

    @pytest.mark.parametrize(
        "query",
        [
            "SELECT a.chrom, COUNT(*) AS n FROM peaks a JOIN genes b "
            "ON a.interval INTERSECTS b.interval GROUP BY a.chrom",
            "SELECT COUNT(*) FROM peaks a LEFT JOIN genes b "
            "ON a.interval INTERSECTS b.interval",
            "SELECT COUNT(*) FROM peaks a JOIN genes b ON a.interval INTERSECTS "
            "b.interval JOIN other c ON a.chrom = c.chrom",
        ],
        ids=["grouped", "left-join", "two-joins"],
    )
    def test_transpile_should_keep_count_star_when_shape_is_not_bare_inner(self, query):
        """Test that COUNT(*) is left alone outside the bare inner shape.

        Given:
            A COUNT(*) that shares its SELECT with a grouped column, sits over a
            LEFT join (steered to the hash join, which handles an empty projection),
            or over more than one join.
        When:
            Transpiling for datafusion-bio.
        Then:
            It should keep ``COUNT(*)``.
        """
        # Act
        sql = _bio(query, ["peaks", "genes", "other"])

        # Assert
        assert "COUNT(*)" in sql
