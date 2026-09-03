"""The DuckDB IEJoin column-to-column INTERSECTS join override (epic #137, #169).

The registered ``(DuckDBTarget, Intersects)`` override expander for DuckDB.
:class:`IntersectsDuckDBIEJoinTransformer` rewrites a column-to-column INTERSECTS
*join* into a per-chromosome dynamic-SQL pattern DuckDB plans through its
range-join family (``IE_JOIN`` / ``PIECEWISE_MERGE_JOIN``);
:func:`expand_intersects_duckdb` wires it into the operator-expander registry so
pass 3 (:class:`giql.expander.ExpandOperators`) dispatches DuckDB's ``Intersects``
nodes here. For the shapes it declines, and for every literal-range or residual
``Intersects`` predicate, the override defers to whichever expander serves
``(GenericTarget, Intersects)`` — by default the shared naive overlap predicate
that ``dialect=None`` / ``"datafusion"`` emit on every target, a plain ``ON``
condition the engine plans as a range join (#167). Resolving that through the
registry rather than calling the built-in directly is what lets a user's
``(GenericTarget, Intersects)`` override reach the shapes this target declines,
instead of being honoured under ``dialect=None`` and silently ignored here.

Because the IEJoin rewrite emits a *whole-query* multi-statement string
(``SET VARIABLE ...; SELECT ... FROM query(getvariable(...))``) rather than a
node-local expression, :func:`expand_intersects_duckdb` produces it through the
query-level :meth:`~giql.expander.ExpansionContext.add_statement_finalizer` seam,
wrapping the built string in an :class:`sqlglot.expressions.Command` so it
serializes verbatim through the stock per-target serializer. This relocates the
former capability-gated pre-pass IEJoin transformer into the registry, mirroring
the CLUSTER / MERGE relocation (#144); INTERSECTS join strategy is now fully
registry-driven and target-overridable. The generic binned equi-join transformer
was dropped in favor of the naive predicate earlier (#167).
"""

from dataclasses import dataclass
from typing import ClassVar
from typing import Literal

from sqlglot import exp

from giql.canonical import canonical_end
from giql.canonical import canonical_start
from giql.constants import DEFAULT_CHROM_COL
from giql.constants import DEFAULT_END_COL
from giql.constants import DEFAULT_START_COL
from giql.expander import REGISTRY
from giql.expander import ExpansionContext
from giql.expander import register
from giql.expanders._per_chrom import CHROM_LITERAL
from giql.expanders._per_chrom import dynamic_relation
from giql.expanders._per_chrom import set_variable_statement
from giql.expanders.intersects import SPATIAL_PREDICATE_META
from giql.expanders.intersects import _expand_spatial_op
from giql.expressions import Contains
from giql.expressions import Intersects
from giql.expressions import SpatialSetPredicate
from giql.expressions import Within
from giql.table import Tables
from giql.targets import DuckDBTarget
from giql.targets import GenericTarget

#: Which side of an outer join keeps its unmatched rows. Narrowed by explicit
#: per-branch assignment rather than a cast: nothing type-checks this package, so
#: a cast asserts without verifying, and the value decides whether the emission
#: swaps the FROM and joined tables.
JoinSide = Literal["left", "right"]


class _UnqualifiedProjectionError(Exception):
    """Signal an unrecoverable projection-resolution failure in the DuckDB IEJoin path.

    Raised by :class:`IntersectsDuckDBIEJoinTransformer` when a query engages
    the dialect path but contains a SELECT-list shape the rebuild cannot
    attribute to a join side: an unqualified column, an unknown table
    qualifier, an aggregate over an unqualified column, or a right-side
    reference under a SEMI / ANTI join, plus an extras predicate that
    references an unqualified or unknown-aliased column. The unknown-table,
    aggregate-of-unqualified, and right-side cases are also rejected by the
    naive plan; a *bare* unqualified column may in fact be naive-valid when it
    is unambiguous, but the dialect has no live schema to attribute it to a
    side, so it surfaces a clean transpile-time error rather than guessing.
    Caught and re-raised as :class:`ValueError` at the dialect's public
    boundary.

    Shapes the naive plan *accepts* but the IEJoin cannot rebuild do NOT raise
    — star projections (#202), and expressions / window aggregates / ``FILTER``
    clauses / scalar subqueries / aggregate-nested stars (#204, #205) — they
    signal :class:`_DeclineIEJoin` and fall back to the naive plan instead.
    """


class _DeclineIEJoin(Exception):
    """Signal that the IEJoin should decline a valid query to the naive plan.

    Distinct from :class:`_UnqualifiedProjectionError`: this is not a user
    error but a shape the naive-predicate plan compiles for free while the
    IEJoin projection rebuild cannot express it (an expression, a window
    aggregate, a ``FILTER`` clause, a scalar subquery, an aggregate nested in
    an expression, or a star nested in an aggregate argument — #204, #205).
    Caught in :meth:`IntersectsDuckDBIEJoinTransformer.transform_to_parts`,
    which returns ``None`` so the generic naive-predicate plan handles the query,
    keeping ``dialect="duckdb"`` consistent with every other backend.
    """


# These predicates support :class:`IntersectsDuckDBIEJoinTransformer`; they live
# at module scope rather than on the class so a future target-specific INTERSECTS
# override can reuse them without importing the transformer. Style §2's strict
# "functions after classes" ordering is relaxed here under the topic-interleaved
# exception.


def _is_column_intersects(node: Intersects) -> bool:
    """Return True if *node* is column-to-column :class:`Intersects`.

    Both operands must be table-qualified column references.
    """
    return (
        isinstance(node.this, exp.Column)
        and bool(node.this.table)
        and isinstance(node.expression, exp.Column)
        and bool(node.expression.table)
    )


def _find_column_intersects_in(expr: exp.Expression) -> Intersects | None:
    """Return the first column-to-column :class:`Intersects` in *expr*, or ``None``."""
    return next(
        (n for n in expr.find_all(Intersects) if _is_column_intersects(n)),
        None,
    )


def _fold_identifier(name: str) -> str:
    """Return *name* case-folded, matching how DuckDB compares identifiers.

    DuckDB resolves identifiers case-insensitively, and unlike the SQL standard
    it does so even when they are quoted: given ``SELECT 1 AS "x", 2 AS "X"``,
    both ``"x"`` and ``"X"`` select the first column. Comparing identifiers with
    raw Python equality would therefore reject mixed-case queries the engine
    itself accepts.

    Folding quoted identifiers too is what makes the duplicate-name gate in
    :func:`_match_outer_join_decomposition` correct: that gate exists because
    ``AS x`` alongside ``AS X`` would otherwise bind both positions to whichever
    column came first.

    Named for identifiers rather than aliases because it is applied to all of
    them -- table aliases, column qualifiers, and projection output names.
    """
    return name.casefold()


def _has_outer_join_intersects(query: exp.Select) -> bool:
    """Return True if any outer JOIN has a column-to-column :class:`Intersects`.

    Handles two surfaces the user can write the INTERSECTS under an outer
    join: directly in the join's ``ON`` clause, or in the top-level
    ``WHERE`` while the join keeps its ``LEFT``/``RIGHT``/``FULL`` side
    modifier (``LEFT JOIN ... ON TRUE WHERE a.interval INTERSECTS
    b.interval``). The dialect's inner-join IEJoin rewrite would silently
    drop the outer-join's unmatched-row guarantee in the second shape,
    so we fall back whenever an outer-join is present and an
    INTERSECTS sits anywhere it could influence the join result.

    Since #95 this is no longer the first gate an outer join meets:
    :func:`_match_outer_join_decomposition` runs ahead of it and claims the
    LEFT / RIGHT shapes it can decompose. What reaches here is therefore
    ``FULL OUTER``, the ``WHERE``-INTERSECTS shape above (where the filter
    would discard the very rows the join preserves), and every LEFT / RIGHT
    shape the decomposition declined.
    """
    outer_joins = [
        join for join in query.args.get("joins") or [] if join.args.get("side")
    ]
    if not outer_joins:
        return False
    for join in outer_joins:
        on = join.args.get("on")
        if on is not None and _find_column_intersects_in(on):
            return True
    where = query.args.get("where")
    if where and _find_column_intersects_in(where.this):
        return True
    return False


def _has_left_only_join_where_intersects(query: exp.Select) -> bool:
    """Return True if a SEMI/ANTI join carries its INTERSECTS in the WHERE.

    A left-only (``SEMI`` / ``ANTI``) join's right table is out of scope in
    the top-level ``WHERE`` — the reference plans (``dialect=None`` /
    ``"datafusion"``) reject a right-referencing ``WHERE`` predicate with a
    binder error. The dialect's ``WHERE``-branch (which exists to support the
    legitimate implicit comma/cross-join idiom, an ``INNER`` join whose right
    table *is* in scope in the ``WHERE``) would instead relocate that
    ``INTERSECTS`` into the join ``ON`` and invent anti/semi-overlap results,
    diverging from every other backend (#201). Decline whenever an explicit
    ``SEMI`` / ``ANTI`` join is present and a column-to-column ``INTERSECTS``
    sits in the top-level ``WHERE`` (necessarily out of scope), so the naive
    plan surfaces the same binder error the reference plans raise. The
    ``INTERSECTS`` search recurses (via ``find_all``), so a column-to-column
    ``INTERSECTS`` nested inside a ``WHERE`` subquery also declines here;
    that is a harmless over-decline (it only ever routes to the correct
    naive plan) and such shapes already decline via other gates.
    """
    has_left_only_join = any(
        (join.args.get("kind") or "").upper() in ("SEMI", "ANTI")
        for join in query.args.get("joins") or []
    )
    if not has_left_only_join:
        return False
    where = query.args.get("where")
    return bool(where and _find_column_intersects_in(where.this))


def _unwrap_projection_target(sel: exp.Expression) -> exp.Expression:
    """Return the expression a SELECT item projects, past its alias and parens.

    Peeling both means a paren-wrapped projection such as ``(SUM(a.score)) AS s``
    is classified identically to its bare form, rather than falling through to
    whatever generic branch handles an unrecognized node.
    """
    target = sel.this if isinstance(sel, exp.Alias) else sel
    while isinstance(target, exp.Paren):
        target = target.this
    return target


def _has_star_projection(query: exp.Select) -> bool:
    """Return True if the SELECT list contains any star projection.

    The IEJoin rewrite restructures the query into per-chromosome subqueries
    that project **explicitly-named** ``__giql_p<n>`` columns up through a
    wrapper relation, so building that dynamic SQL requires enumerating the
    star's columns at transpile time. But GIQL only knows the genomic columns
    declared in the :class:`~giql.table.Table` config (chrom / start / end /
    strand) — it has no view of the live table schema, so it cannot see
    arbitrary base-table columns (``name`` / ``score`` / ...). Expanding
    ``a.*`` here would silently drop them, narrowing the result relative to
    every other backend (#202). Decline any star projection (bare ``*`` or
    qualified ``t.*``, with or without a user alias) so the naive-predicate
    plan expands it against DuckDB's live schema at bind time, keeping the
    projection identical across backends.
    """
    for sel in query.expressions:
        target = _unwrap_projection_target(sel)
        if isinstance(target, exp.Star):
            return True
        if isinstance(target, exp.Column) and isinstance(target.this, exp.Star):
            return True
    return False


def _projection_declines_to_naive(
    target: exp.Expression,
    l_alias: str,
    r_alias: str,
    is_left_only: bool,
) -> bool:
    """Return True if *target* is a projection the IEJoin should decline to naive.

    The IEJoin projection rebuild only handles a qualified column (``a.col``)
    or a plain aggregate over qualified columns (``SUM(a.score)`` / ``COUNT(*)``);
    it emits explicitly-named ``__giql_p<n>`` inner columns and reconstructs the
    outer SELECT from them. Any other shape the naive-predicate plan compiles
    for free — an expression (``a.start + 1``), a window aggregate, a ``FILTER``
    clause, a scalar subquery, an aggregate nested in an expression
    (``COUNT(*) * 2``), or a star nested in an aggregate argument
    (``COUNT(a.*)`` / ``MIN(COLUMNS(*))``) — should decline so ``dialect="duckdb"``
    stays consistent with every other backend instead of hard-erroring or, for
    the aggregate-nested star, miscompiling (#204, #205).

    A projection that references a column **in its own scope** the rebuild
    cannot attribute to a join side (unqualified, an unknown table, or the
    right side under a SEMI / ANTI join) does NOT decline — the dispatch raises
    a clean :class:`_UnqualifiedProjectionError` instead. The unknown-table and
    right-side cases are rejected by the naive plan too; a bare unqualified
    column may be naive-valid when unambiguous, but the dialect has no live
    schema to pick a side, so a transpile-time error is more actionable than
    guessing. Columns inside a nested subquery resolve against that subquery,
    not the join, so they are ignored by this test.
    """
    # A qualified, non-star column is rebuilt directly — never declines here.
    if (
        isinstance(target, exp.Column)
        and target.table
        and not isinstance(target.this, exp.Star)
    ):
        return False
    # An aggregate is decided here in full (it never falls through to the
    # column scan below — an unqualified aggregate argument like ``SUM(score)``
    # is diagnosed in ``render_aggregate``, not treated as a decline). A star
    # nested in an aggregate argument (``COUNT(a.*)`` / ``MIN(COLUMNS(*))`` /
    # ``COUNT(a.* EXCLUDE (name))``) cannot be enumerated, so it declines
    # regardless of any EXCLUDE / REPLACE modifier columns (#204); a plain
    # aggregate over qualified columns (``COUNT(*)`` / ``SUM(a.score)`` /
    # ``COUNT(DISTINCT a.col)``) is rebuilt directly.
    if isinstance(target, exp.AggFunc):
        has_star_arg = any(
            isinstance(col.this, exp.Star) for col in target.find_all(exp.Column)
        )
        return has_star_arg or target.find(exp.Columns) is not None
    # Unsupported by the rebuild. Decline to naive unless an own-scope column
    # is out of scope (a user error the naive plan rejects too — let the
    # dispatch raise instead).
    for col in target.find_all(exp.Column):
        if isinstance(col.this, exp.Star):
            continue  # qualified star — declinable, not a plain column ref
        # Skip a column inside a nested subquery within *target* (it resolves
        # against that subquery, not the join). *target* is an ancestor of
        # *col*, so walking up from ``col.parent`` reaches it; a Select /
        # Subquery strictly between them is a nested scope. When ``col`` IS
        # ``target`` (a bare column projection), there is no interior to walk.
        in_subquery = False
        node = col.parent if col is not target else target
        while node is not None and node is not target:
            if isinstance(node, (exp.Select, exp.Subquery)):
                in_subquery = True
                break
            node = node.parent
        if in_subquery:
            continue
        tbl = _fold_identifier(col.table) if col.table else ""
        if tbl not in (l_alias, r_alias) or (is_left_only and tbl == r_alias):
            return False
    return True


def _strip_intersects(
    expr: exp.Expression | None, intersects: Intersects
) -> exp.Expression | None:
    """Return the parts of an AND tree that aren't the given :class:`Intersects`."""
    if expr is None or expr is intersects:
        return None
    if isinstance(expr, exp.And):
        if expr.this is intersects:
            return expr.expression
        if expr.expression is intersects:
            return expr.this
        left = _strip_intersects(expr.this, intersects)
        right = _strip_intersects(expr.expression, intersects)
        if left is None:
            return right
        if right is None:
            return left
        return exp.And(this=left, expression=right)
    return expr


def _count_column_intersects(query: exp.Select) -> int:
    """Count column-to-column :class:`Intersects` predicates anywhere in *query*.

    Walks the full AST so the dialect's "exactly one INTERSECTS" gate cannot
    be circumvented by an INTERSECTS hidden inside a subquery (e.g. inside
    a SELECT-list scalar subquery or inside ``EXISTS (SELECT ...)``). The
    SELECT-list / modifier-clause subquery gates in ``transform_to_parts``
    already reject those shapes, but this widened count is defense-in-depth
    against future gate relaxations.
    """
    return sum(1 for n in query.find_all(Intersects) if _is_column_intersects(n))


@dataclass(frozen=True)
class _IEJoinSides:
    """Resolved two-sided join shape for the DuckDB IEJoin dialect path.

    Bundles the left/right user aliases (normalized via :func:`_fold_identifier`)
    with the un-aliased :class:`~sqlglot.expressions.Table` AST nodes and a
    ``kind`` discriminator identifying the join variant. The hardcoded inner
    aliases ``a`` and ``b`` used inside each per-chromosome subquery are
    exposed as :data:`left_inner_alias` and :data:`right_inner_alias`.

    The :attr:`is_left_only_join` property is True for SEMI and ANTI shapes,
    which restrict the outer SELECT to left-side projections only.
    """

    left_user_alias: str
    right_user_alias: str
    left_table: exp.Table
    right_table: exp.Table
    kind: Literal["INNER", "SEMI", "ANTI"]

    left_inner_alias: ClassVar[str] = "a"
    right_inner_alias: ClassVar[str] = "b"

    @property
    def is_left_only_join(self) -> bool:
        """True when the join's output only references the left side (SEMI / ANTI)."""
        return self.kind in ("SEMI", "ANTI")

    @classmethod
    def from_intersects(
        cls,
        from_table: exp.Table,
        join: exp.Join,
        intersects: Intersects,
    ) -> "_IEJoinSides | None":
        """Resolve sides, run the FROM-side orientation swap, normalize aliases.

        Returns ``None`` when the join's table operand isn't a base
        :class:`~sqlglot.expressions.Table`, when neither side of the
        :class:`Intersects` predicate references the FROM table's alias, or
        when the join's ``kind`` is outside the supported set
        (``INNER`` / ``CROSS`` / ``SEMI`` / ``ANTI``; ``CROSS`` folds to
        ``INNER`` because per-chromosome partitioning makes them equivalent
        inside the inner subquery).
        """
        join_table = join.this
        if not isinstance(join_table, exp.Table):
            return None
        from_alias = _fold_identifier(from_table.alias_or_name)
        left_alias = _fold_identifier(intersects.this.table)
        right_alias = _fold_identifier(intersects.expression.table)
        if left_alias == from_alias:
            l_user_alias = from_alias
            l_table = from_table
            r_user_alias = right_alias
            r_table = join_table
        elif right_alias == from_alias:
            # Inner subquery template hardcodes "a" to the FROM side; swap so
            # that contract holds regardless of which operand the user put
            # first in INTERSECTS.
            l_user_alias = from_alias
            l_table = from_table
            r_user_alias = left_alias
            r_table = join_table
        else:
            return None
        raw_kind = (join.args.get("kind") or "").upper()
        kind: Literal["INNER", "SEMI", "ANTI"]
        if raw_kind in ("", "CROSS"):
            kind = "INNER"
        elif raw_kind == "SEMI":
            kind = "SEMI"
        elif raw_kind == "ANTI":
            kind = "ANTI"
        else:
            return None
        return cls(
            left_user_alias=l_user_alias,
            right_user_alias=r_user_alias,
            left_table=l_table,
            right_table=r_table,
            kind=kind,
        )


def _resolve_two_base_tables(
    query: exp.Select,
) -> tuple[exp.Join, exp.Table, exp.Table] | None:
    """Return *query*'s single join and its two base-table operands, or ``None``.

    Accepts only a one-join query whose FROM and joined operands are both named
    base tables. A GIQL table-function operand such as ``DISJOIN(genes)`` parses
    as an :class:`~sqlglot.expressions.Table` with an empty ``name``, which would
    interpolate an empty identifier into broken SQL, so an unnamed operand
    declines here rather than downstream. (A CTE-named operand is caught by the
    ``with_`` guard; an unregistered base table is fine -- the dialect uses
    default columns exactly as the naive-predicate path does.)

    Every rewrite in this module has to name both relations before it can emit
    anything, so these gates precede all of them.

    What they deliberately do not cover is alias distinctness. The callers
    disagree about which aliases to compare -- syntactic FROM/JOIN position here
    versus INTERSECTS-operand orientation in :meth:`transform_to_parts` -- and
    those differ for a query like ``FROM peaks a JOIN genes b ON a.interval
    INTERSECTS c.interval``, so folding the comparison in would change which
    shapes decline.
    """
    joins = query.args.get("joins") or []
    if len(joins) != 1:
        return None
    the_join = joins[0]

    from_node = query.args.get("from_")
    from_table = from_node.this if from_node else None
    join_table = the_join.this
    if not isinstance(from_table, exp.Table) or not from_table.name:
        return None
    if not isinstance(join_table, exp.Table) or not join_table.name:
        return None
    return the_join, from_table, join_table


@dataclass(frozen=True)
class _IntersectsJoin:
    """A two-base-table join whose ``ON`` is exactly one column-to-column INTERSECTS.

    The resolved form the two shape matchers share. ``from_table`` and
    ``join_table`` keep their aliases, since the emissions render them into
    ``FROM`` clauses whose columns resolve through those aliases;
    ``from_alias`` and ``join_alias`` are the case-folded forms used to
    attribute each projected column to a side.
    """

    side: JoinSide
    from_table: exp.Table
    join_table: exp.Table
    from_alias: str
    join_alias: str

    @classmethod
    def from_query(
        cls, query: exp.Select, allowed_sides: tuple[JoinSide, ...]
    ) -> "_IntersectsJoin | None":
        """Return *query*'s single resolved INTERSECTS join, or ``None``.

        Builds on :func:`_resolve_two_base_tables` and adds the gates specific to
        a bare outer-join ``ON``: the accepted side, a plain join kind, exactly
        one column-to-column INTERSECTS with no residual, and two distinct
        aliases. Shared by the two shape matchers; the INNER / SEMI / ANTI path
        in :meth:`IntersectsDuckDBIEJoinTransformer.transform_to_parts` accepts
        shapes these gates reject and so builds on the base-table helper
        directly.

        *allowed_sides* entries must already be case-folded, since the join's own
        side is folded before the membership test.
        """
        resolved = _resolve_two_base_tables(query)
        if resolved is None:
            return None
        the_join, from_table, join_table = resolved

        # Assigned per branch rather than cast, so an unexpected side declines
        # here instead of reaching the emission as an unhandled value.
        raw_side = _fold_identifier(the_join.args.get("side") or "")
        side: JoinSide
        if raw_side == "left":
            side = "left"
        elif raw_side == "right":
            side = "right"
        else:
            return None
        if side not in allowed_sides:
            return None
        if (the_join.args.get("kind") or "").upper() not in ("", "OUTER"):
            return None

        on = the_join.args.get("on")
        if on is None:
            return None
        intersects = _find_column_intersects_in(on)
        if intersects is None or _count_column_intersects(query) != 1:
            return None
        # v1 supports only a bare INTERSECTS ON. A residual would have to be
        # applied identically to every relation a composed rewrite emits -- the
        # INNER half's join and the unmatched half's NOT EXISTS correlation --
        # which none of them coordinate yet.
        if _strip_intersects(on, intersects) is not None:
            return None

        from_alias = _fold_identifier(from_table.alias_or_name)
        join_alias = _fold_identifier(join_table.alias_or_name)
        if from_alias == join_alias:
            return None

        return cls(
            side=side,
            from_table=from_table,
            join_table=join_table,
            from_alias=from_alias,
            join_alias=join_alias,
        )


@dataclass(frozen=True)
class _CountOverlapsShape:
    """The matched ``count_overlaps`` LEFT-join-with-COUNT shape (#209).

    Bundles the left (FROM) base table, the left-side group-key SELECT items,
    and the single ``COUNT`` aggregate SELECT item. The emission reuses the
    INNER path to build the per-chromosome IEJoin count and wraps it in a
    zero-fill LEFT join against the distinct left keys.

    ``left_table`` keeps its alias: the distinct-key relation renders it
    directly and the group-key projections resolve through that alias.
    """

    left_table: exp.Table
    key_items: tuple[exp.Expression, ...]
    agg_item: exp.Alias


def _match_count_overlaps(query: exp.Select) -> "_CountOverlapsShape | None":
    """Return the count_overlaps shape when *query* matches it, else ``None``.

    Matches exactly ``SELECT <left cols>, COUNT(<right col>) AS <alias> FROM
    <base> a LEFT JOIN <base> b ON <one column-to-column INTERSECTS> GROUP BY
    <left cols>`` -- the only right-side reference is the COUNT argument, and the
    projected left columns are the group keys. Every other shape (``COUNT(*)``,
    non-COUNT or multiple aggregates, a WHERE / HAVING / ORDER BY / LIMIT, a
    non-LEFT outer join, extra ON residuals) returns ``None`` so the caller keeps
    the existing behaviour (the naive-predicate plan).
    """
    if query.args.get("with_") or query.args.get("distinct"):
        return None
    if any(
        query.args.get(k) is not None for k in ("having", "order", "limit", "offset")
    ):
        return None
    if query.args.get("group") is None or query.args.get("where") is not None:
        return None

    join = _IntersectsJoin.from_query(query, ("left",))
    if join is None:
        return None
    from_table = join.from_table
    left_alias = join.from_alias
    right_alias = join.join_alias

    key_items: list[exp.Expression] = []
    agg_item: exp.Alias | None = None
    for sel in query.expressions:
        target = _unwrap_projection_target(sel)
        if isinstance(target, exp.AggFunc):
            # Exactly one COUNT(<right col>) / COUNT(DISTINCT <right col>) aggregate,
            # aliased so the wrapper can reference and re-project it.
            if agg_item is not None or not isinstance(target, exp.Count):
                return None
            if not isinstance(sel, exp.Alias) or not sel.output_name:
                return None
            arg = target.this
            if isinstance(arg, exp.Distinct):
                if len(arg.expressions) != 1:
                    return None
                arg = arg.expressions[0]
            if not (
                isinstance(arg, exp.Column)
                and not isinstance(arg.this, exp.Star)
                and _fold_identifier(arg.table) == right_alias
            ):
                return None
            agg_item = sel
        elif (
            isinstance(target, exp.Column)
            and target.table
            and not isinstance(target.this, exp.Star)
        ):
            if _fold_identifier(target.table) != left_alias:
                return None
            key_items.append(sel)
        else:
            return None

    if agg_item is None or not key_items:
        return None
    # The GROUP BY keys must be exactly the projected left columns. An extra
    # group key produces more rows than the distinct projected keys (which the
    # zero-fill base relation reconstructs and so cannot reproduce), and a
    # projected key absent from the group would be an invalid aggregate query.
    projected_key_exprs = {
        (item.this if isinstance(item, exp.Alias) else item).sql(dialect="duckdb")
        for item in key_items
    }
    group_exprs = {e.sql(dialect="duckdb") for e in query.args["group"].expressions}
    if group_exprs != projected_key_exprs:
        return None
    # The zero-fill wrapper joins on the user-facing output names and re-projects
    # the count under its alias, so those names must be distinct — otherwise a
    # key column and the count alias collide and USING / COALESCE bind to the
    # wrong column.
    output_names = [item.output_name for item in key_items] + [agg_item.output_name]
    if len(set(output_names)) != len(output_names):
        return None

    return _CountOverlapsShape(
        left_table=from_table,
        key_items=tuple(key_items),
        agg_item=agg_item,
    )


@dataclass(frozen=True)
class _OuterJoinShape:
    """A LEFT / RIGHT outer-join ``INTERSECTS`` decomposable into two halves (#95).

    Carries the resolved join plus the SELECT-list positions belonging to each
    side. ``side`` is ``"left"`` or ``"right"``; for ``"right"`` the emission
    first swaps the FROM and joined tables so one LEFT-shaped code path serves
    both orientations.

    The row-preserving side is derived from the join rather than copied out of
    it, so the table node and the alias cannot drift apart. ``preserved_table``
    keeps its alias intact -- the NULL-chromosome branch selects from it
    directly and references its columns through that alias, so stripping the
    alias would leave that branch unable to bind.

    ``preserved_indexes`` and ``null_indexes`` partition the projection by
    position: the preserved side's columns survive into both halves, while the
    other side's become NULL fills in the unmatched half.
    """

    join: _IntersectsJoin
    preserved_indexes: tuple[int, ...]
    null_indexes: tuple[int, ...]

    @property
    def side(self) -> JoinSide:
        """Which side keeps its unmatched rows."""
        return self.join.side

    @property
    def preserved_table(self) -> exp.Table:
        """The row-preserving side's table node, alias intact."""
        return self.join.from_table if self.side == "left" else self.join.join_table

    @property
    def preserved_alias(self) -> str:
        """The row-preserving side's case-folded alias."""
        return self.join.from_alias if self.side == "left" else self.join.join_alias


def _match_outer_join_decomposition(query: exp.Select) -> "_OuterJoinShape | None":
    """Return the decomposable outer-join shape for *query*, else ``None``.

    Matches ``SELECT <side-attributable columns> FROM <base> a LEFT|RIGHT JOIN
    <base> b ON <one column-to-column INTERSECTS>`` -- the shape a ``UNION ALL``
    of an INNER half and an unmatched-rows half reproduces exactly (#95).

    Everything else returns ``None`` so the caller keeps the naive-predicate
    plan. In particular:

    - ``FULL OUTER`` preserves rows on both sides and has no two-half form here.
    - An ``INTERSECTS`` in the top-level ``WHERE`` under an outer join is a
      different query -- the filter drops the very unmatched rows the outer join
      preserves -- so it must keep declining, as must any other ``WHERE``
      residual this rewrite does not yet thread through both halves.
    - Top-level modifiers (``GROUP BY`` / ``HAVING`` / ``ORDER BY`` / ``LIMIT``
      / ``OFFSET`` / ``DISTINCT``) belong on the union, but the recursive
      transpile of each half would apply them per branch instead.
    - A projection that is not a side-attributable qualified column (a star, an
      expression, an aggregate) cannot be split into "survives" and "becomes
      NULL", which is exactly what the unmatched half needs. This is also why
      the ``bedtools -wao`` ``CASE`` projection still declines until #109 lands.
    """
    # Whitelist, not a blocklist: the rewrite re-emits the query as a UNION ALL
    # of two independently transpiled halves, so ANY top-level clause it does
    # not itself read would be applied per-half instead of over the union --
    # silently changing the answer. Enumerating clauses to *reject* leaves the
    # rewrite exposed to every clause sqlglot grows later, which is exactly how
    # ``QUALIFY`` slipped onto the fast path and had its filter dropped. Reject
    # anything outside the set ``_build_outer_join_parts`` actually consumes.
    if any(
        value is not None and value != []
        for key, value in query.args.items()
        if key not in ("expressions", "from_", "joins")
    ):
        return None

    join = _IntersectsJoin.from_query(query, ("left", "right"))
    if join is None:
        return None
    # A TABLESAMPLE hangs off the table node rather than the top-level SELECT,
    # so the clause whitelist above cannot see it. Each half would sample the
    # relation independently, making the union neither a superset nor a subset
    # of a real outer join over one sample -- and the sample clause lands in a
    # position the per-chromosome template cannot even parse.
    if join.from_table.args.get("sample") is not None:
        return None
    if join.join_table.args.get("sample") is not None:
        return None

    if join.side == "left":
        preserved_alias, other_alias = join.from_alias, join.join_alias
    else:
        preserved_alias, other_alias = join.join_alias, join.from_alias

    if not query.expressions:
        return None
    preserved_indexes: list[int] = []
    null_indexes: list[int] = []
    for idx, sel in enumerate(query.expressions):
        target = _unwrap_projection_target(sel)
        if (
            not isinstance(target, exp.Column)
            or isinstance(target.this, exp.Star)
            or not target.table
        ):
            return None
        tbl = _fold_identifier(target.table)
        if tbl == preserved_alias:
            preserved_indexes.append(idx)
        elif tbl == other_alias:
            null_indexes.append(idx)
        else:
            return None

    # The unmatched half re-projects the preserved columns out of a subquery by
    # their output names, so duplicates there would bind ambiguously. The
    # comparison must be case-INSENSITIVE: DuckDB resolves identifiers that way
    # even when they are quoted, so ``AS x`` alongside ``AS X`` would otherwise
    # pass this gate and then bind both positions to whichever column came
    # first -- returning the wrong value and, because the branch's type then
    # wins the UNION ALL unification, corrupting the matched rows too.
    # Deliberately preserved-side only: duplicate names among the NULL-filled
    # positions are never re-projected by name, so they stay supported.
    preserved_names = [
        _fold_identifier(query.expressions[i].output_name) for i in preserved_indexes
    ]
    if len(set(preserved_names)) != len(preserved_names):
        return None

    return _OuterJoinShape(
        join=join,
        preserved_indexes=tuple(preserved_indexes),
        null_indexes=tuple(null_indexes),
    )


class IntersectsDuckDBIEJoinTransformer:
    """Transform column-to-column INTERSECTS joins into a DuckDB IEJoin pattern.

    Emits dynamic SQL in two steps per declared session variable — one
    variable for the INNER / SEMI / ANTI shapes, two for a decomposed outer
    join:

    1. ``SET VARIABLE __giql_iejoin_<token> = COALESCE(<dynamic-sql>, <empty-schema>);``
       where ``<dynamic-sql>`` is a ``string_agg`` over the chromosome
       partition that emits one per-chromosome overlap query joined by
       ``UNION ALL``. For INNER, that is ``SELECT ... FROM (... WHERE chrom =
       '<c>') a JOIN (... WHERE chrom = '<c>') b ON <range-overlap>``; for
       SEMI / ANTI it is ``SELECT ... FROM (... WHERE chrom = '<c>') a WHERE
       [NOT] EXISTS (SELECT 1 FROM (... WHERE chrom = '<c>') b WHERE
       <range-overlap>)``.
    2. ``SELECT <projections> FROM query(getvariable('__giql_iejoin_<token>'))``

    Each per-chromosome subquery uses pure inequality predicates which DuckDB
    plans through its range-join family (``IE_JOIN`` or ``PIECEWISE_MERGE_JOIN``).
    The dialect supports INNER, SEMI, and ANTI variants, plus LEFT / RIGHT
    outer joins. Crucially, DuckDB only selects ``IE_JOIN`` for INNER
    pure-inequality joins — a bare ``SEMI JOIN`` / ``ANTI JOIN`` inequality join
    is planned as a ``BLOCKWISE_NL_JOIN`` (a nested loop, quadratic). So the
    left-only shapes are emitted as a correlated
    ``WHERE EXISTS`` (SEMI) / ``WHERE NOT EXISTS`` (ANTI) subquery instead of a
    semi/anti join keyword, which reaches ``IE_JOIN`` while preserving exact
    semantics (#208). SEMI / ANTI restrict the outer SELECT to left-side
    projections only.

    See :doc:`/transpilation/performance` for the full join-shape support
    matrix, fallback enumeration, and hard-error projection shapes.

    Notes:

    - Star projections (bare ``*`` and ``a.*`` / ``b.*``) decline to the
      generic naive-predicate plan, which delegates star expansion to DuckDB
      and projects every base-table column against the live schema. The
      dialect cannot enumerate a star at transpile time — it only knows the
      :class:`~giql.table.Table` config's genomic columns (chrom / start /
      end / strand), not arbitrary user columns — so declining keeps the
      projection identical across backends rather than silently narrowing
      it (#202).
    - A projection the rebuild cannot express but the naive plan handles —
      an expression (``a.start + 1``), a window aggregate, a ``FILTER``
      clause, a scalar subquery, an aggregate nested in an expression
      (``COUNT(*) * 2``), or a star nested in an aggregate argument
      (``COUNT(a.*)`` / ``MIN(COLUMNS(*))``) — also declines to the naive
      plan, keeping ``dialect="duckdb"`` consistent with every backend
      rather than hard-erroring or miscompiling (#204, #205). Only a
      genuine out-of-scope column reference (unqualified / unknown-table /
      right-side) still raises.
    - ``SEMI JOIN`` and ``ANTI JOIN`` outputs are left-side only. Any
      ``b.col`` reference in the outer SELECT, aggregate arguments, GROUP BY,
      HAVING, or ORDER BY raises :class:`_UnqualifiedProjectionError`
      (wrapped to :class:`ValueError` at the public boundary). A right-side
      ``b.*`` under SEMI / ANTI declines to the naive plan with every other
      star (which then rejects the out-of-scope right table at bind time).
    - ``ANTI JOIN`` partitions on the chromosome INTERSECT like INNER and
      SEMI, then unions one non-partitioned branch carrying the left rows
      whose chromosome the right table lacks. Those rows cannot match, so a
      branch apiece would scan both tables to prove an emptiness the
      partition already knows (#95).
    - ``LEFT`` / ``RIGHT`` outer joins are not emitted as outer joins at all:
      DuckDB's ``IE_JOIN`` is INNER-only, so the query is decomposed into an
      INNER half and a ``NOT EXISTS`` unmatched half combined with
      ``UNION ALL``, plus a third branch carrying the preserved side's
      NULL-chromosome rows, which no per-chromosome partition can reach
      (#95). ``RIGHT`` is served by swapping the FROM and joined tables so a
      single LEFT-shaped path covers both orientations.
    """

    def __init__(self, tables: Tables):
        """Initialize transformer.

        :param tables:
            Table configurations for column mapping.
        """
        self.tables = tables

    def transform_to_sql(self, query: exp.Expression) -> str | None:
        """Return DuckDB IEJoin SQL for *query* if applicable, else ``None``.

        Renders what :meth:`transform_to_parts` returns as one script. Both
        reach the same primitive: a rewrite composing on top of this one calls
        that method instead of splitting this string, which cannot be taken
        apart safely -- the separator appears inside the script itself whenever
        an identifier contains it.

        The result is a multi-statement script whose ``SET VARIABLE`` setup and
        final ``SELECT`` must run on one connection, in order.

        *query* is borrowed, not owned: it is never mutated. The registry entry
        point passes the live statement AST that pass 3 goes on transforming,
        and the naive-predicate fallback rewrites that same node in place when
        this path declines, so that has to hold.

        Raises
        ------
        ValueError
            If a projection names a column the rewrite cannot attribute to a
            side.
        """
        parts = self.transform_to_parts(query)
        if parts is None:
            return None
        setup, select = parts
        return f"{setup};\n{select}"

    def transform_to_parts(self, query: exp.Expression) -> tuple[str, str] | None:
        """Return ``(setup statements, final SELECT)`` for *query*, else ``None``.

        The dialect rewrites the *whole* query as
        ``SET VARIABLE __giql_iejoin_<token> = ...; SELECT ... FROM
        query(getvariable(...))``, so any top-level clause :meth:`_build_sql`
        does not read would be silently dropped. The check below is a
        whitelist: engage the dialect path only for a query that is
        exactly ``SELECT <qualified columns / plain aggregates> FROM
        <registered base table> [JOIN <registered base table>] {ON|WHERE}
        <one column-to-column INTERSECTS>`` with no other top-level clause.
        Star projections (schema-less enumeration would narrow them, #202),
        projections the rebuild cannot express but the naive plan handles —
        expressions / window aggregates / FILTER clauses / scalar subqueries /
        aggregate-nested stars (#204, #205) — and a ``SEMI`` / ``ANTI`` join
        whose ``INTERSECTS`` sits out of scope in the ``WHERE`` (#201) also
        decline. Everything else returns ``None`` so the generic
        naive-predicate plan handles it.
        """
        if not isinstance(query, exp.Select):
            return None

        # Top-level WITH (CTEs) still falls back to the naive-predicate plan.
        # ORDER BY / LIMIT / OFFSET / DISTINCT / GROUP BY / HAVING ride
        # on the outer SELECT wrapper. ``DISTINCT ON (...)`` is the one
        # DISTINCT variant we don't try to translate — it would interact
        # with the outer alias rewrite in ways the naive-predicate plan handles
        # for free.
        if query.args.get("with_"):
            return None
        distinct_node = query.args.get("distinct")
        if distinct_node is not None and distinct_node.args.get("on"):
            return None

        # count_overlaps fast path (#209): a LEFT JOIN whose only right-side use
        # is a COUNT(b.col) with a GROUP BY on the left keys reaches IE_JOIN via
        # the INNER path plus a zero-fill LEFT join, rather than declining to the
        # naive HASH_JOIN + inequality filter that every other outer join takes.
        # Checked before the outer-join decline below since it *is* an outer join.
        count_shape = _match_count_overlaps(query)
        if count_shape is not None:
            return self._build_count_overlaps_parts(query, count_shape)

        # Outer-join decomposition (#95): a LEFT / RIGHT JOIN whose projections
        # are side-attributable columns is emitted as the INNER pairs UNION ALL
        # the unmatched preserved-side rows with NULL fills, so both halves reach
        # IE_JOIN. Checked after count_overlaps -- also a LEFT JOIN, with its own
        # faster zero-fill path -- and before the decline below, which still
        # catches FULL OUTER and the WHERE-INTERSECTS shape.
        outer_shape = _match_outer_join_decomposition(query)
        if outer_shape is not None:
            decomposed = self._build_outer_join_parts(query, outer_shape)
            if decomposed is not None:
                return decomposed

        if _has_outer_join_intersects(query):
            return None

        # A SEMI / ANTI join whose INTERSECTS lives in the top-level WHERE
        # (out of scope for the right table) diverges from the reference
        # plans, which reject it with a binder error. Decline so the naive
        # plan surfaces that error instead of inventing results (#201).
        if _has_left_only_join_where_intersects(query):
            return None

        # Star projections (bare ``*`` or ``a.*`` / ``b.*``) need column
        # enumeration the schema-less transpile cannot do without silently
        # narrowing to the configured genomic columns. Decline so the naive
        # plan expands the real star against DuckDB's live schema (#202).
        if _has_star_projection(query):
            return None

        if _count_column_intersects(query) != 1:
            return None

        # The wrapper-relation rewriter is not scope-aware: any subquery in
        # GROUP BY / HAVING / ORDER BY (or inside a SELECT-list aggregate
        # argument) would have its inner column refs substituted with
        # wrapper-relation aliases that don't exist in the subquery's own
        # scope. Reject those shapes here. ``EXISTS (SELECT ...)`` parses
        # as ``Exists(Select(...))`` without an ``exp.Subquery`` wrapper,
        # so we also check ``exp.Select`` directly. A bare scalar subquery
        # in the SELECT list (``SELECT (SELECT SUM(x) FROM ...)``) is not
        # gated here — it declines to the naive plan later via the
        # :meth:`_resolve_projections` projection pre-scan (#205), since its
        # columns resolve against its own scope.
        for clause_key in ("group", "having", "order"):
            modifier = query.args.get(clause_key)
            if modifier is None:
                continue
            if (
                modifier.find(exp.Subquery) is not None
                or modifier.find(exp.Select) is not None
            ):
                return None
        for sel in query.expressions:
            target = _unwrap_projection_target(sel)
            if isinstance(target, exp.AggFunc) and (
                target.find(exp.Subquery) is not None
                or target.find(exp.Select) is not None
            ):
                return None

        # The supported shape connects exactly two named base tables via one
        # INTERSECTS (either an explicit JOIN ... ON or an implicit comma-join
        # in WHERE). A third comma/CROSS-joined table would otherwise be
        # silently dropped by _build_sql.
        resolved = _resolve_two_base_tables(query)
        if resolved is None:
            return None
        the_join, from_table, join_table = resolved

        # NATURAL needs schema introspection we don't have at transpile time
        # (Table configs only expose chrom/start/end/strand; user columns
        # like name/score are invisible). The naive-predicate plan defers to
        # DuckDB, which does see the schema, so falling back IS the optimal plan.
        if the_join.args.get("method"):
            return None

        # USING admits only the single-column form whose column matches both
        # tables' ``chrom_col`` (the per-chromosome partition is exactly the
        # equi-join that ``USING(chrom)`` requests). The chrom-equivalence
        # check itself happens in ``_build_sql`` after ``l_table`` /
        # ``r_table`` are resolved. Multi-column USING and USING(non-chrom)
        # fall back; inline support for those is a documented follow-up.
        using_cols = the_join.args.get("using") or []
        if using_cols and len(using_cols) != 1:
            return None

        # INNER (the default), CROSS (functionally equivalent to INNER once
        # the chromosome partition is in place), SEMI, and ANTI all engage.
        # Everything else falls back so the naive-predicate plan can preserve
        # its semantics -- FULL OUTER, and the LEFT / RIGHT shapes the outer-join
        # decomposition declined before handing the query here (#95), which the
        # naive overlap predicate expresses as a plain outer ``ON`` condition.
        kind_str = the_join.args.get("kind")
        if kind_str and kind_str.upper() not in ("INNER", "CROSS", "SEMI", "ANTI"):
            return None

        # Locates which of the two relations carries the INTERSECTS, and where
        # (an explicit ON, or the top-level WHERE of a comma join). The single
        # join gate above means it can only return the join already resolved.
        intersects, _ = self._find_target_join(query)
        if intersects is None:
            return None

        sides = _IEJoinSides.from_intersects(from_table, the_join, intersects)
        if sides is None:
            return None

        if sides.left_user_alias == sides.right_user_alias:
            # Same alias on both sides of the INTERSECTS — the inner subquery
            # would alias the same underlying relation as both "a" and "b",
            # which the naive-predicate plan renders correctly in place (no
            # inner subquery, so no ambiguous double-alias).
            return None

        # Compare fully-qualified identifiers (catalog/db/name) so two
        # distinct same-named tables in different schemas are not treated
        # as a self-join. Same qualified identifier still falls back so the
        # naive-predicate plan handles the legitimate self-join shape.
        if self._qualified_table_ident(sides.left_table) == self._qualified_table_ident(
            sides.right_table
        ):
            return None

        # Extra predicates are inlined into each per-chromosome subquery's
        # join ON, except WHERE residuals under a left-only (SEMI / ANTI) join,
        # which :meth:`_build_sql` applies as an outer wrapper filter (#200).
        # Soft-fallback residuals (subquery,
        # aggregate, window, or an INTERSECTS still wrapped in OR/NOT/Paren)
        # cause ``_build_sql`` to return None. User-mistake residuals
        # (unqualified or unknown-aliased columns, or columns with a
        # catalog/schema qualifier) raise :class:`_UnqualifiedProjectionError`,
        # which we wrap to :class:`ValueError` here so the public surface is
        # uniform.
        try:
            return self._build_sql(query, intersects, sides, the_join)
        except _DeclineIEJoin:
            # A projection the naive plan handles but the IEJoin cannot rebuild
            # (expression / window / FILTER / scalar subquery / aggregate-nested
            # star) — fall back to the naive-predicate plan (#204, #205).
            return None
        except _UnqualifiedProjectionError as exc:
            raise ValueError(str(exc)) from exc

    @staticmethod
    def _sql_escape(s: str) -> str:
        """Return *s* with single quotes doubled for safe SQL string literal use."""
        return s.replace("'", "''")

    @staticmethod
    def _sql_quote_ident(name: str) -> str:
        """Render *name* as a DuckDB-quoted identifier.

        Delegates to sqlglot's identifier emitter so embedded double quotes
        are doubled correctly and DuckDB-specific quoting rules (case,
        reserved words) are honored. Used at every identifier-interpolation
        site in the dialect's dynamic SQL.
        """
        return exp.to_identifier(name, quoted=True).sql(dialect="duckdb")

    @staticmethod
    def _qualified_table_ident(table: exp.Table) -> str:
        """Render *table*'s catalog/db/name as a qualified DuckDB identifier.

        Drops the user-supplied alias (``peaks AS a``) so the rendered string
        is suitable both as a bare ``FROM`` operand (where the caller adds
        its own alias) and as a subquery source. Catalog and schema qualifiers
        survive — sqlglot's :meth:`~sqlglot.expressions.Expression.sql` emitter
        handles identifier quoting per DuckDB's rules (reserved words and
        special characters get quoted; safe unquoted identifiers stay
        unquoted).
        """
        no_alias = table.copy()
        no_alias.set("alias", None)
        return no_alias.sql(dialect="duckdb")

    @staticmethod
    def _classify_extras(
        extras: list[exp.Expression],
    ) -> Literal["inline", "fallback"]:
        """Return ``"fallback"`` if any extra needs the naive-predicate plan.

        Soft-fallback shapes: an INTERSECTS still embedded in the residual
        (because :func:`_strip_intersects` couldn't peel an OR/NOT/Paren
        wrapper), a subquery, an aggregate function, or a window function.
        Returns ``"inline"`` when every extra is safe to inline into the
        per-chromosome subquery's join ON clause.
        """
        for predicate in extras:
            if predicate.find(Intersects) is not None:
                return "fallback"  # INTERSECTS wrapped in OR/NOT/Paren
            if predicate.find(exp.Subquery) is not None:
                return "fallback"
            if predicate.find(exp.Select) is not None:
                return "fallback"
            if predicate.find(exp.AggFunc) is not None:
                return "fallback"
            if predicate.find(exp.Window) is not None:
                return "fallback"
        return "inline"

    @staticmethod
    def _validate_extra_qualifiers(
        extras: list[exp.Expression],
        sides: "_IEJoinSides",
    ) -> None:
        """Raise for unqualified, unknown-aliased, or schema-qualified columns.

        Surfaces user mistakes at transpile time rather than deferring to a
        downstream binder error. Three failure modes raise:

        - **Unqualified column:** no ``.table`` qualifier.
        - **Unknown alias:** ``.table`` doesn't match either side of the join.
        - **Catalog/schema-qualified column:** ``mycat.myschema.a.score``
          would silently mis-bind once the alias-rewriter remaps ``a`` to
          the inner subquery's hardcoded ``a`` alias (the ``mycat.myschema``
          qualifier survives the rewrite and addresses a different relation
          than intended).
        """
        valid_aliases = {sides.left_user_alias, sides.right_user_alias}
        for predicate in extras:
            for col in predicate.find_all(exp.Column):
                col_table = _fold_identifier(col.table) if col.table else ""
                if not col_table:
                    raise _UnqualifiedProjectionError(
                        f"dialect='duckdb' cannot inline extra predicate "
                        f"{predicate.sql()!r}: column {col.sql()!r} must be "
                        f"qualified with {sides.left_user_alias!r} or "
                        f"{sides.right_user_alias!r}."
                    )
                if col_table not in valid_aliases:
                    raise _UnqualifiedProjectionError(
                        f"dialect='duckdb' cannot inline extra predicate "
                        f"{predicate.sql()!r}: column references unknown "
                        f"table qualifier {col.table!r}; expected "
                        f"{sides.left_user_alias!r} or "
                        f"{sides.right_user_alias!r}."
                    )
                if col.args.get("db") or col.args.get("catalog"):
                    raise _UnqualifiedProjectionError(
                        f"dialect='duckdb' cannot inline extra predicate "
                        f"{predicate.sql()!r}: column {col.sql()!r} carries "
                        f"a catalog/schema qualifier; use the alias-only "
                        f"form ({sides.left_user_alias!r}.<col> or "
                        f"{sides.right_user_alias!r}.<col>)."
                    )

    @staticmethod
    def _rewrite_refs_for_per_chrom_subquery(
        node: exp.Expression,
        sides: "_IEJoinSides",
    ) -> str:
        """Render a residual predicate against the inner subquery's ``a``/``b`` scope.

        The user's AST references their original aliases (``peaks a`` /
        ``genes b`` or ``peaks p`` / ``genes g``), but the per-chromosome
        inner subquery is emitted with hardcoded aliases ``a`` and ``b``
        regardless of what the user wrote. Rewrite each column reference
        accordingly before emitting the residual SQL into the inner join's
        ON clause. Alias matching is case-folded so mixed-case identifiers
        match the DuckDB-equivalent contract.

        Assumes residuals have already been validated to reference only
        ``sides.left_user_alias`` or ``sides.right_user_alias``; columns
        without a table qualifier or pointing at an unknown alias must
        have been rejected upstream.
        """

        def replace(n: exp.Expression) -> exp.Expression:
            if not isinstance(n, exp.Column):
                return n
            col_table = _fold_identifier(n.table) if n.table else ""
            if col_table == sides.left_user_alias:
                new = n.copy()
                new.set("table", exp.to_identifier(sides.left_inner_alias))
                return new
            if col_table == sides.right_user_alias:
                new = n.copy()
                new.set("table", exp.to_identifier(sides.right_inner_alias))
                return new
            return n

        rewritten = node.transform(replace, copy=True)
        return rewritten.sql(dialect="duckdb")

    @staticmethod
    def _rewrite_refs_for_wrapper_relation(
        node: exp.Expression,
        projection_map: dict[tuple[str, str], str],
        sides: "_IEJoinSides",
        clause_label: str,
    ) -> str:
        """Rewrite table-qualified column refs to the wrapper relation's inner alias.

        The dialect's outer ``SELECT`` reads from ``query(getvariable(...)) AS
        <alias>``. The wrapper relation's columns are whatever the inner
        per-chromosome subqueries projected — i.e., the ``__giql_p<n>`` names
        on the left-hand side of each ``AS`` in the inner SELECT. Any user
        ``ORDER BY`` / ``GROUP BY`` / ``HAVING`` (and any aggregate-function
        argument) written against the original table aliases (``a.start``,
        ``b.score``) must be rewritten to those inner alias names so DuckDB
        can resolve them against the wrapper relation. Alias matching is
        case-folded for DuckDB-equivalent semantics.

        An upstream pre-allocation step ensures every referenced column has
        an entry in ``projection_map``; missing keys raise
        :class:`_UnqualifiedProjectionError` as a guard against future scope
        relaxations that would otherwise produce silent miscompiles.

        Scope caveat: this walker is not scope-aware — it rewrites every
        :class:`~sqlglot.expressions.Column` whose ``.table`` matches either
        side's alias anywhere inside ``node``, including correlated or
        non-correlated subqueries that may rebind those aliases to a
        different table. The dispatcher gates already reject any subquery
        in GROUP BY / HAVING / ORDER BY / SELECT list so the unsafe shape
        can't reach here today. A future relaxation of those gates must
        add a scope-aware walker here.
        """
        valid_aliases = (sides.left_user_alias, sides.right_user_alias)

        def replace(n: exp.Expression) -> exp.Expression:
            if not isinstance(n, exp.Column):
                return n
            col_table = _fold_identifier(n.table) if n.table else ""
            if col_table not in valid_aliases:
                return n
            key = (col_table, n.name)
            if key not in projection_map:
                raise _UnqualifiedProjectionError(
                    f"dialect='duckdb' cannot resolve {n.sql()!r} in "
                    f"{clause_label}: ensure the column is referenced "
                    f"with a table qualifier matching the join's left "
                    f"or right alias."
                )
            return exp.column(projection_map[key])

        rewritten = node.transform(replace, copy=True)
        return rewritten.sql(dialect="duckdb")

    def _iejoin_variant_parts(
        self, query: exp.Select, *, kind: str | None
    ) -> tuple[str, str] | None:
        """Transpile a copy of *query* with its single join retargeted to *kind*.

        The copy is taken here rather than by the caller, so no caller ever hands
        its own node to the transform. Whether the emission mutates what it is
        given stops being something each composed rewrite has to know and get
        right; before this, one of them defended the invariant with a comment
        and a second defensive copy.
        """
        variant = query.copy()
        join = variant.args["joins"][0]
        join.set("side", None)
        join.set("kind", kind)
        return self.transform_to_parts(variant)

    def _build_count_overlaps_parts(
        self, query: exp.Select, shape: "_CountOverlapsShape"
    ) -> tuple[str, str] | None:
        """Render the count_overlaps zero-fill SQL for a matched LEFT-join shape (#209).

        Reuses the INNER IEJoin path: an INNER copy of the query (same
        projections and GROUP BY, LEFT downgraded to INNER) is transpiled by
        :meth:`transform_to_parts` into a setup statement and a count
        ``SELECT``, which DuckDB plans through ``IE_JOIN`` + hash aggregate. The
        INNER count drops left rows with no overlap, so its ``SELECT`` is
        wrapped as a CTE and LEFT-joined back onto the distinct left keys,
        zero-filling the missing counts. Returns ``None`` (fall back to the
        naive plan) if the INNER copy itself declines.
        """
        inner_parts = self._iejoin_variant_parts(query, kind=None)
        if inner_parts is None:
            return None
        set_var_stmt, inner_select = inner_parts

        q = self._sql_quote_ident
        key_names = [item.output_name for item in shape.key_items]
        key_projections = [item.sql(dialect="duckdb") for item in shape.key_items]
        group_exprs = [
            (item.this if isinstance(item, exp.Alias) else item).sql(dialect="duckdb")
            for item in shape.key_items
        ]
        agg_name = shape.agg_item.output_name

        # Distinct left-key relation, rendering the FROM table with its own alias
        # so the key projections (``a.chrom`` ...) resolve against it.
        base_cte = (
            f"SELECT {', '.join(key_projections)} "
            f"FROM {shape.left_table.sql(dialect='duckdb')} "
            f"GROUP BY {', '.join(group_exprs)}"
        )
        using_cols = ", ".join(q(name) for name in key_names)
        base_out = ", ".join(f"base.{q(name)}" for name in key_names)
        wrapper = (
            f"WITH __giql_counts AS ({inner_select}), "
            f"__giql_base AS ({base_cte}) "
            f"SELECT {base_out}, COALESCE(c.{q(agg_name)}, 0) AS {q(agg_name)} "
            f"FROM __giql_base AS base "
            f"LEFT JOIN __giql_counts AS c USING ({using_cols})"
        )
        return set_var_stmt, wrapper

    def _build_outer_join_parts(
        self, query: exp.Select, shape: "_OuterJoinShape"
    ) -> tuple[str, str] | None:
        """Render an outer join as INNER pairs UNION ALL unmatched rows (#95).

        Reuses the two IEJoin paths that already exist rather than emitting an
        outer join: an INNER copy supplies the matched pairs, and an ANTI copy
        (the #208 correlated ``NOT EXISTS`` rewrite) supplies the preserved-side
        rows with no overlap. Both reach ``IE_JOIN``. A per-chromosome ``LEFT
        JOIN`` does not, even though ``EXPLAIN`` reports ``IE_JOIN`` for it.
        Measured, that plan runs in about 0.04s at 100,000 rows per side and
        does not finish within 300s at 4,194,304, where this decomposition
        returns in about two seconds. Plan inspection alone does not separate
        them; only execution at scale does.

        Both halves partition on the chromosome ``INTERSECT``, since a
        chromosome present on one side only yields no pairs. What completes the
        union is what the ANTI half adds on top: one non-partitioned branch for
        the preserved side's chromosomes the other side lacks, so preserved rows
        there still surface with NULLs.

        Returns ``None`` (fall back to the naive-predicate plan) when either
        half declines, so the query falls back as one unit.
        """
        # RIGHT JOIN is the mirror of LEFT: swapping the FROM and joined tables
        # rewrites ``a RIGHT JOIN b`` as the equivalent ``b LEFT JOIN a``, so a
        # single code path serves both. The projections are alias-qualified and
        # resolve unchanged across the swap.
        work = query.copy()
        if shape.side == "right":
            work_join = work.args["joins"][0]
            from_node = work.args["from_"]
            swapped_from = work_join.this
            work_join.set("this", from_node.this)
            from_node.set("this", swapped_from)

        inner_parts = self._iejoin_variant_parts(work, kind=None)
        if inner_parts is None:
            return None
        inner_setup, inner_select = inner_parts

        # The unmatched half projects the preserved side only: the other side is
        # out of scope under ANTI, and a NULL literal in the SELECT list would
        # itself decline the projection rebuild. The NULLs are re-inserted by the
        # wrapper below instead.
        preserved_config = self.tables.get(shape.preserved_table.name)
        chrom_col = preserved_config.chrom_col if preserved_config else DEFAULT_CHROM_COL

        preserved_sels = [query.expressions[i].copy() for i in shape.preserved_indexes]
        if not preserved_sels:
            # Nothing from the preserved side is projected (``SELECT b.chrom FROM
            # a LEFT JOIN b ...``). The half still needs a column to select, so
            # carry the preserved side's chromosome column and drop it in the
            # wrapper, which projects NULL for every original position anyway.
            preserved_sels = [exp.column(chrom_col, table=shape.preserved_alias)]
        # Safe to mutate in place: the INNER half above transpiled a copy, so
        # nothing downstream holds a reference to ``work``.
        work.set("expressions", preserved_sels)
        anti_parts = self._iejoin_variant_parts(work, kind="ANTI")
        if anti_parts is None:
            return None
        anti_setup, anti_select = anti_parts

        # Both remaining branches re-project the original SELECT order, filling
        # the other side's positions with NULL. That NULL-fill rule is what makes
        # the union reproduce outer-join semantics, so it is expressed once here
        # rather than once per branch. DuckDB resolves each NULL's type from the
        # first UNION ALL branch, so the dialect needs no schema knowledge.
        #
        # They differ only in how a preserved position is spelled: the unmatched
        # half reads it back out of a subquery by output name, while the
        # NULL-chromosome branch selects from the base table and so renders the
        # original qualified column.
        q = self._sql_quote_ident
        preserved_names = [
            query.expressions[i].output_name for i in shape.preserved_indexes
        ]
        wrapper_items: list[str] = []
        null_chrom_items: list[str] = []
        preserved_seen = 0
        for idx, sel in enumerate(query.expressions):
            output = q(sel.output_name)
            if idx in shape.null_indexes:
                null_fill = f"NULL AS {output}"
                wrapper_items.append(null_fill)
                null_chrom_items.append(null_fill)
                continue
            wrapper_items.append(f"{q(preserved_names[preserved_seen])} AS {output}")
            preserved_seen += 1
            column = _unwrap_projection_target(sel)
            null_chrom_items.append(f"{column.sql(dialect='duckdb')} AS {output}")

        unmatched_select = (
            f"SELECT {', '.join(wrapper_items)} FROM ({anti_select}) AS __giql_unmatched"
        )

        # A preserved-side row whose chromosome is NULL can never match: the
        # per-chromosome partition is an equality on chrom and ``NULL = x`` is
        # never true. Neither half can surface those rows on its own -- both
        # partitions come from ``SELECT DISTINCT <chrom>``, where a NULL renders
        # as a NULL literal that ``string_agg`` skips, so no branch is ever
        # emitted for it. Union them in directly, or the decomposition would
        # silently drop the very rows an outer join exists to preserve.
        # Qualify the WHERE with the same identifier node the FROM clause renders
        # rather than with the case-folded alias. The two bind to one relation
        # today only because DuckDB folds identifiers -- quoted ones included --
        # so ``FROM peaks AS "A"`` happened to accept ``a.chrom``. Deriving both
        # from one node removes that dependency instead of resting on it.
        #
        # The sibling above, which carries the preserved chromosome into the
        # ANTI half, deliberately keeps the folded alias: that column is fed
        # back through transform_to_parts, which matches aliases
        # case-insensitively and re-emits it against the per-chromosome
        # subquery's own ``a``, never against the user's alias.
        alias_node = shape.preserved_table.args.get("alias")
        preserved_qualifier = (
            alias_node.this if alias_node is not None else shape.preserved_table.this
        )
        preserved_chrom = exp.column(chrom_col, table=preserved_qualifier)
        null_chrom_select = (
            f"SELECT {', '.join(null_chrom_items)} "
            f"FROM {shape.preserved_table.sql(dialect='duckdb')} "
            f"WHERE {preserved_chrom.sql(dialect='duckdb')} IS NULL"
        )

        # The matched half is used as a UNION ALL branch directly. Wrapping it in
        # ``SELECT * FROM (...)`` would make DuckDB de-duplicate any repeated
        # output name (``chrom`` becoming ``chrom_1``), and that renaming would
        # become the union's schema -- so a flag documented as a pure performance
        # opt-in would silently relabel the result columns, and precisely for the
        # ``a.chrom, ..., b.chrom`` projection bedtools users write.
        union = (
            f"{inner_select} UNION ALL {unmatched_select} UNION ALL {null_chrom_select}"
        )
        return ";\n".join([inner_setup, anti_setup]), union

    def _extract_extra_predicates(
        self,
        query: exp.Select,
        intersects: Intersects,
    ) -> tuple[list[exp.Expression], list[exp.Expression]]:
        """Return the non-INTERSECTS residuals split by provenance.

        Walks the AND tree under each JOIN ON and the WHERE, strips the
        target INTERSECTS node, and returns ``(on_residuals, where_residuals)``
        in document order. ``on_residuals`` are the predicates that sat in a
        ``JOIN ... ON`` alongside the INTERSECTS; ``where_residuals`` are the
        predicates that sat in the top-level ``WHERE``. Both lists empty means
        the query has no extras beyond the target INTERSECTS.

        The split matters because the two sources are *not* interchangeable
        for every join kind. An ON residual is a genuine join condition and is
        correct to inline into the per-chromosome ON for INNER / SEMI / ANTI
        alike. A WHERE residual is a post-join filter: inlining it into the ON
        is only sound for INNER (where ``ON overlap WHERE r`` ≡ ``ON overlap
        AND r``) and SEMI (a left-only ``r`` is invariant over the right side);
        for ANTI it inverts the anti-join for rows failing ``r`` (#200). See
        :meth:`_build_sql`, which routes WHERE residuals to an outer filter for
        the left-only (SEMI / ANTI) shapes.

        Note: :func:`_strip_intersects` only descends ``exp.And``. Predicates
        whose AND tree wraps the INTERSECTS in ``exp.Or`` / ``exp.Not`` /
        ``exp.Paren`` surface here with the INTERSECTS still embedded;
        :meth:`_classify_extras` then routes them to the naive-predicate plan, while
        :meth:`_validate_extra_qualifiers` enforces qualifier rules for the
        remaining inlinable residuals.
        """
        on_residuals: list[exp.Expression] = []
        for join in query.args.get("joins") or []:
            on = join.args.get("on")
            if on is None:
                continue
            residual = _strip_intersects(on, intersects)
            if residual is not None:
                on_residuals.append(residual)
        where_residuals: list[exp.Expression] = []
        where = query.args.get("where")
        if where:
            residual = _strip_intersects(where.this, intersects)
            if residual is not None:
                where_residuals.append(residual)
        return on_residuals, where_residuals

    def _find_target_join(
        self, query: exp.Select
    ) -> tuple[Intersects | None, exp.Join | None]:
        """Locate the single INTERSECTS join target.

        Returns ``(intersects_node, join)`` where ``join`` is the
        :class:`~sqlglot.expressions.Join` AST node carrying the target
        INTERSECTS — either directly via its ``ON`` clause, or via an
        implicit cross-join whose ``INTERSECTS`` predicate lives in the
        top-level ``WHERE``. Alias matching on the WHERE-branch is
        case-folded for DuckDB-equivalent semantics.
        """
        for join in query.args.get("joins") or []:
            on = join.args.get("on")
            if on:
                intersects = _find_column_intersects_in(on)
                if intersects:
                    if not isinstance(join.this, exp.Table):
                        return None, None
                    return intersects, join

        where = query.args.get("where")
        if where:
            intersects = _find_column_intersects_in(where.this)
            if intersects:
                from_table = query.args["from_"].this
                if not isinstance(from_table, exp.Table):
                    return None, None
                from_alias = _fold_identifier(from_table.alias_or_name)
                left_alias = _fold_identifier(intersects.this.table)
                right_alias = _fold_identifier(intersects.expression.table)
                target_alias = right_alias if left_alias == from_alias else left_alias
                for join in query.args.get("joins") or []:
                    if isinstance(join.this, exp.Table):
                        if _fold_identifier(join.this.alias_or_name) == target_alias:
                            return intersects, join

        return None, None

    def _build_sql(
        self,
        query: exp.Select,
        intersects: Intersects,
        sides: "_IEJoinSides",
        the_join: exp.Join,
    ) -> tuple[str, str] | None:
        """Render the DuckDB IEJoin ``(setup statements, final SELECT)`` for *query*.

        The two halves are returned separately rather than pre-joined so a
        rewrite composed on top never has to recover them by splitting rendered
        SQL, which would break on an identifier containing ``";\\n"``.

        Returns ``None`` when a soft-fallback condition fires (the residual
        of the join carries a shape the naive-predicate plan can handle but the
        dialect cannot inline, or a single-column USING references a column
        that is not both tables' ``chrom_col``). May propagate
        :class:`_DeclineIEJoin` from the :meth:`_resolve_projections` projection
        pre-scan (a naive-valid projection the rebuild cannot express — #204,
        #205); :meth:`transform_to_parts` catches it and declines to the naive
        plan. Raises :class:`_UnqualifiedProjectionError` when a user-mistake
        condition fires; callers translating to the public surface wrap the
        raise to :class:`ValueError`.
        """
        on_residuals, where_residuals = self._extract_extra_predicates(query, intersects)
        all_residuals = on_residuals + where_residuals
        if all_residuals and self._classify_extras(all_residuals) == "fallback":
            return None
        if all_residuals:
            self._validate_extra_qualifiers(all_residuals, sides)

        # ON residuals are genuine join conditions and inline into the
        # per-chromosome ON for every kind. WHERE residuals are post-join
        # filters: inlining them into the ON is only sound for INNER (and the
        # equivalent SEMI case), but inverts the anti-join for rows failing the
        # residual under ANTI (#200). For the left-only shapes (SEMI / ANTI) we
        # therefore layer WHERE residuals as an outer filter on the wrapper
        # relation instead of folding them into the join. INNER keeps inlining
        # them (byte-identical output, provably equivalent).
        if sides.is_left_only_join:
            inline_residuals = on_residuals
            outer_where_residuals = where_residuals
        else:
            inline_residuals = on_residuals + where_residuals
            outer_where_residuals = []

        # Tables registry is keyed by the bare table name; an unregistered
        # table uses default column names like the naive-predicate plan does.
        l_table = self.tables.get(sides.left_table.name)
        r_table = self.tables.get(sides.right_table.name)
        l_chrom = l_table.chrom_col if l_table else DEFAULT_CHROM_COL
        l_start = l_table.start_col if l_table else DEFAULT_START_COL
        l_end = l_table.end_col if l_table else DEFAULT_END_COL
        r_chrom = r_table.chrom_col if r_table else DEFAULT_CHROM_COL
        r_start = r_table.start_col if r_table else DEFAULT_START_COL
        r_end = r_table.end_col if r_table else DEFAULT_END_COL

        # USING(<col>) admission check (the gate already restricted to
        # single-column USING). The dialect's per-chromosome partition IS
        # the equi-join that USING(<chrom_col>) requests; if the USING
        # column doesn't match both tables' chrom_col (or the chrom_cols
        # diverge), fall back so the naive-predicate plan preserves the
        # USING equi-join for the engine to resolve.
        using_cols = the_join.args.get("using") or []
        if using_cols:
            using_name = _fold_identifier(using_cols[0].name)
            if (
                _fold_identifier(l_chrom) != _fold_identifier(r_chrom)
                or _fold_identifier(l_chrom) != using_name
            ):
                return None

        (
            inner_projections,
            outer_projections,
            projection_map,
        ) = self._resolve_projections(
            query,
            sides,
            outer_where_residuals,
        )

        # Canonicalize endpoints to 0-based half-open and use strict
        # operators. Mixing inclusive operators with raw closed-closed
        # coordinates would silently report non-overlapping touching
        # intervals as overlapping.
        q = self._sql_quote_ident
        l_start_expr = canonical_start(f"a.{q(l_start)}", l_table)
        l_end_expr = canonical_end(f"a.{q(l_end)}", l_table)
        r_start_expr = canonical_start(f"b.{q(r_start)}", r_table)
        r_end_expr = canonical_end(f"b.{q(r_end)}", r_table)

        # Render each table operand qualified (catalog/db/name) without the
        # user's alias so the empty-schema fallback can supply its own
        # ``a``/``b`` alias without producing ``FROM peaks AS a a`` shapes.
        l_table_ident = self._qualified_table_ident(sides.left_table)
        r_table_ident = self._qualified_table_ident(sides.right_table)

        # The overlap predicate (canonical, strict) plus any inline residuals,
        # shared by both the INNER join-ON form and the SEMI / ANTI EXISTS form.
        overlap_predicates = [
            f"{l_start_expr} < {r_end_expr}",
            f"{l_end_expr} > {r_start_expr}",
        ]
        for residual in inline_residuals:
            overlap_predicates.append(
                self._rewrite_refs_for_per_chrom_subquery(residual, sides)
            )
        predicate_sql = " AND ".join(overlap_predicates)

        # The dynamic SQL builder, expressed as a SQL string-concat that
        # ``string_agg`` aggregates per-chromosome. Single quotes are
        # doubled because the result is itself an SQL string literal that
        # contains SQL. The chromosome literal is interpolated twice — once
        # into the left partition filter, once into the right.
        chrom_literal = CHROM_LITERAL
        esc = self._sql_escape
        inner_select_list = ", ".join(inner_projections)

        if sides.is_left_only_join:
            # SEMI / ANTI: a bare per-chromosome ``SEMI JOIN`` / ``ANTI JOIN``
            # with an inequality overlap predicate is planned by DuckDB as a
            # ``BLOCKWISE_NL_JOIN`` (a nested loop), not its fast ``IE_JOIN``
            # range join — DuckDB selects IE_JOIN only for INNER pure-inequality
            # joins (#208). A correlated ``WHERE EXISTS`` / ``WHERE NOT EXISTS``
            # over the same per-chromosome right partition reaches IE_JOIN and
            # preserves exact SEMI / ANTI semantics (one row per qualifying left
            # row, no dedup-on-projection hazard). Inline (ON) residuals sit in
            # the subquery WHERE beside the overlap predicate; WHERE residuals
            # are still layered as an outer wrapper filter (``outer_where_residuals``
            # below, #200).
            exists_keyword = "EXISTS" if sides.kind == "SEMI" else "NOT EXISTS"
            select_from_left = (
                f"SELECT {inner_select_list} "
                f"FROM (SELECT * FROM {l_table_ident} WHERE {q(l_chrom)} = "
            )
            exists_open = (
                f") a WHERE {exists_keyword} (SELECT 1 FROM "
                f"(SELECT * FROM {r_table_ident} WHERE {q(r_chrom)} = "
            )
            exists_close = f") b WHERE {predicate_sql})"
            per_chrom_sql_expr = (
                f"'{esc(select_from_left)}' "
                f"|| {chrom_literal} "
                f"|| '{esc(exists_open)}' "
                f"|| {chrom_literal} "
                f"|| '{esc(exists_close)}'"
            )
        else:
            # INNER (CROSS folded to INNER inside ``_IEJoinSides``): a
            # per-chromosome INNER join whose pure-inequality ON predicate DuckDB
            # plans through ``IE_JOIN``.
            per_chrom_select_prefix = f"SELECT {inner_select_list} "
            from_left = f"FROM (SELECT * FROM {l_table_ident} WHERE {q(l_chrom)} = "
            from_right_open = (
                f") a JOIN (SELECT * FROM {r_table_ident} WHERE {q(r_chrom)} = "
            )
            on_clause = f") b ON {predicate_sql}"
            per_chrom_sql_expr = (
                f"'{esc(per_chrom_select_prefix)}{esc(from_left)}' "
                f"|| {chrom_literal} "
                f"|| '{esc(from_right_open)}' "
                f"|| {chrom_literal} "
                f"|| '{esc(on_clause)}'"
            )

        # Empty-schema fallback when no chromosomes pass the partition
        # filter. Project from the source tables under WHERE FALSE so
        # DuckDB resolves every output column to its real declared type.
        # SEMI / ANTI outputs reference only the left side, so we drop
        # the right operand from the empty-schema FROM clause.
        empty_select_list = ", ".join(inner_projections)
        if sides.is_left_only_join:
            empty_schema = (
                f"SELECT {empty_select_list} FROM {l_table_ident} a WHERE FALSE"
            )
        else:
            empty_schema = (
                f"SELECT {empty_select_list} "
                f"FROM {l_table_ident} a, {r_table_ident} b WHERE FALSE"
            )

        # Partition source: the chromosomes both sides share. A chromosome
        # present on only one side can never satisfy the overlap predicate, so a
        # per-chromosome branch for it is wasted work -- it scans both tables to
        # prove an emptiness the partition already knows. ANTI must still emit
        # its left-only chromosomes' rows, but as one non-partitioned tail
        # branch rather than one branch apiece (#95).
        chrom_partition_subquery = (
            f"SELECT DISTINCT {q(l_chrom)} AS chrom FROM {l_table_ident} "
            f"INTERSECT SELECT DISTINCT {q(r_chrom)} AS chrom "
            f"FROM {r_table_ident}"
        )

        # ANTI's left-only chromosomes, in one branch. Both NULL guards are
        # load-bearing: ``x NOT IN (<set containing NULL>)`` is never true, and a
        # left row whose chromosome is NULL must stay out of this branch to match
        # the per-chromosome path, which drops it because ``string_agg`` skips
        # the NULL partition row (#224). The outer-join decomposition depends on
        # that exclusion -- it emits NULL-chromosome rows from its own third
        # branch and would otherwise double-count them.
        tail_branch: str | None = None
        if sides.kind == "ANTI":
            tail_branch = (
                f"SELECT {inner_select_list} FROM {l_table_ident} a "
                f"WHERE a.{q(l_chrom)} IS NOT NULL "
                f"AND a.{q(l_chrom)} NOT IN ("
                f"SELECT DISTINCT {q(r_chrom)} FROM {r_table_ident} "
                f"WHERE {q(r_chrom)} IS NOT NULL)"
            )

        # The same result without the per-chromosome partition, selected at
        # execution time when the chromosome count exceeds the ceiling. The
        # chromosome equality moves from the partition filter into the join or
        # correlation predicate; ANTI keeps the IS NOT NULL guard so both forms
        # agree about NULL-chromosome rows.
        chrom_equality = f"a.{q(l_chrom)} = b.{q(r_chrom)}"
        if sides.is_left_only_join:
            null_guard = ""
            if sides.kind == "ANTI":
                null_guard = f"a.{q(l_chrom)} IS NOT NULL AND "
            unpartitioned_query = (
                f"SELECT {inner_select_list} FROM {l_table_ident} a "
                f"WHERE {null_guard}{exists_keyword} ("
                f"SELECT 1 FROM {r_table_ident} b "
                f"WHERE {chrom_equality} AND {predicate_sql})"
            )
        else:
            unpartitioned_query = (
                f"SELECT {inner_select_list} FROM {l_table_ident} a "
                f"JOIN {r_table_ident} b "
                f"ON {chrom_equality} AND {predicate_sql}"
            )

        var_name, set_var_stmt = set_variable_statement(
            per_chrom_sql_expr,
            chrom_partition_subquery,
            empty_schema,
            unpartitioned_query,
            tail_branch=tail_branch,
        )

        # Constant wrapper alias — uniqueness is provided by the
        # token-bearing session variable; the outer SELECT's relation
        # alias is never user-visible and doesn't need to vary per call.
        wrapper_alias = "__giql_iejoin_wrapper"
        outer_projection = ", ".join(outer_projections)
        distinct_kw = "DISTINCT " if query.args.get("distinct") else ""
        outer_select_parts = [
            f"SELECT {distinct_kw}{outer_projection} "
            f"FROM {dynamic_relation(var_name, wrapper_alias)}"
        ]

        # WHERE residuals for the left-only (SEMI / ANTI) shapes are applied
        # here as a post-join filter on the wrapper relation, rewritten to the
        # wrapper's inner-alias column names — never folded into the per-chrom
        # ON, which would invert the anti-join for rows failing the residual
        # (#200). Every referenced column was pre-allocated into the wrapper's
        # projection list by :meth:`_resolve_projections`.
        if outer_where_residuals:
            where_sql = " AND ".join(
                self._rewrite_refs_for_wrapper_relation(
                    residual, projection_map, sides, "WHERE"
                )
                for residual in outer_where_residuals
            )
            outer_select_parts.append(f"WHERE {where_sql}")

        # The rewriter translates any ``a.col`` / ``b.col`` references in
        # the GROUP BY / HAVING / ORDER BY clauses to the wrapper
        # relation's inner-alias column names (``__giql_p<n>``), since
        # the user's original aliases are not in scope in the outer SELECT.
        group_node = query.args.get("group")
        if group_node is not None:
            outer_select_parts.append(
                self._rewrite_refs_for_wrapper_relation(
                    group_node, projection_map, sides, "GROUP BY"
                )
            )

        having_node = query.args.get("having")
        if having_node is not None:
            outer_select_parts.append(
                self._rewrite_refs_for_wrapper_relation(
                    having_node, projection_map, sides, "HAVING"
                )
            )

        order_node = query.args.get("order")
        if order_node is not None:
            outer_select_parts.append(
                self._rewrite_refs_for_wrapper_relation(
                    order_node, projection_map, sides, "ORDER BY"
                )
            )

        limit_node = query.args.get("limit")
        if limit_node is not None:
            outer_select_parts.append(limit_node.sql(dialect="duckdb"))

        offset_node = query.args.get("offset")
        if offset_node is not None:
            outer_select_parts.append(offset_node.sql(dialect="duckdb"))

        outer_select = " ".join(outer_select_parts)

        return set_var_stmt, outer_select

    def _resolve_projections(
        self,
        query: exp.Select,
        sides: "_IEJoinSides",
        outer_where_residuals: list[exp.Expression] | None = None,
    ) -> tuple[list[str], list[str], dict[tuple[str, str], str]]:
        """Resolve the SELECT list (and any modifier-referenced columns).

        ``outer_where_residuals`` are the WHERE residuals that :meth:`_build_sql`
        applies as an outer filter on the wrapper relation (the left-only
        SEMI / ANTI shapes — #200). Their columns are pre-allocated into
        ``inner_projections`` here so the wrapper exposes them even when they
        appear only in the filter and not in the SELECT list. For INNER the
        list is empty because those residuals inline into the per-chromosome ON.

        Returns ``(inner_projections, outer_projections, projection_map)``:

        * ``inner_projections`` is the per-chromosome subquery's SELECT
          list — every column the query touches in any clause, projected
          once each under a unique ``__giql_p<n>`` alias.
        * ``outer_projections`` is the outer SELECT's projection list,
          rendered using the inner aliases (qualified columns and plain
          aggregate calls). Star projections never reach here — they decline
          upstream (#202) — nor do the naive-valid shapes the projection
          pre-scan declines (expressions, window aggregates, FILTER clauses,
          scalar subqueries, aggregate-nested stars — #204, #205).
        * ``projection_map`` binds each ``(normalized_alias, column_name)``
          referenced anywhere in the query to its inner alias name so the
          wrapper-relation rewriter and aggregate-argument rewriting can
          translate user-written references to the wrapper relation.

        Raises :class:`_DeclineIEJoin` (caught upstream, declining to the naive
        plan) for a SELECT-list shape the naive plan handles but the rebuild
        cannot express. Raises :class:`_UnqualifiedProjectionError` only for a
        genuine user error the naive plan also rejects — an unqualified column,
        an unknown table qualifier, an aggregate over an unqualified column, or
        (when ``sides.is_left_only_join``) a reference to the right side.
        """
        q = self._sql_quote_ident
        l_alias = sides.left_user_alias
        r_alias = sides.right_user_alias

        # Decline (fall back to the naive-predicate plan) for any projection the
        # IEJoin cannot rebuild but the naive plan handles — an expression, a
        # window aggregate, a FILTER clause, a scalar subquery, an aggregate
        # nested in an expression, or a star nested in an aggregate argument
        # (#204, #205). Run before any allocation so an aggregate-nested star
        # (COUNT(a.*)) declines rather than allocating a bogus ``a."*"`` column.
        # A projection with an out-of-scope column reference is a user error the
        # naive plan also rejects, so it is left to the dispatch below to raise.
        for sel in query.expressions:
            target = _unwrap_projection_target(sel)
            if _projection_declines_to_naive(
                target, l_alias, r_alias, sides.is_left_only_join
            ):
                raise _DeclineIEJoin(
                    f"dialect='duckdb' cannot rebuild projection {target.sql()!r}; "
                    "deferring to the naive-predicate plan."
                )

        inner_projections: list[str] = []
        projection_map: dict[tuple[str, str], str] = {}

        def reject_right_side_if_left_only(tbl_normalized: str, source_sql: str) -> None:
            """Raise when a SEMI / ANTI join projects from the right side."""
            if sides.is_left_only_join and tbl_normalized == r_alias:
                raise _UnqualifiedProjectionError(
                    f"dialect='duckdb' with kind={sides.kind!r} only "
                    f"projects left-side columns; got {source_sql!r} "
                    f"which references the right side "
                    f"({sides.right_user_alias!r})."
                )

        def allocate(tbl: str, col: str) -> str:
            """Return the inner alias for ``(tbl, col)``, allocating if new.

            Allocation order is deterministic per query but depends on
            (a) the fixed pre-walk order in :meth:`_resolve_projections` —
            modifier columns first (``("group", "having", "order")``),
            then aggregate-argument columns, then the SELECT list — and
            (b) sqlglot's ``find_all`` traversal order within each
            sub-clause. Tests therefore should NOT couple to a specific
            ``__giql_p<n>`` ↔ source-column mapping (e.g.
            ``assert projection_map[(a, start)] == "__giql_p3"``).
            """
            tbl = _fold_identifier(tbl)
            key = (tbl, col)
            existing = projection_map.get(key)
            if existing is not None:
                return existing
            if tbl == l_alias:
                side = sides.left_inner_alias
            elif tbl == r_alias:
                side = sides.right_inner_alias
            else:
                raise _UnqualifiedProjectionError(
                    f"Unknown table qualifier {tbl!r} in projection; "
                    f"expected {l_alias!r} or {r_alias!r}."
                )
            alias = f"__giql_p{len(inner_projections)}"
            inner_projections.append(f"{side}.{q(col)} AS {alias}")
            projection_map[key] = alias
            return alias

        def render_aggregate(agg_node: exp.Expression, user_alias: str | None) -> None:
            """Validate aggregate-argument qualifiers and emit the outer projection."""
            for col in agg_node.find_all(exp.Column):
                col_table = _fold_identifier(col.table) if col.table else ""
                if col_table not in (l_alias, r_alias):
                    raise _UnqualifiedProjectionError(
                        "dialect='duckdb' requires qualified aggregate "
                        f"arguments; got {agg_node.sql()!r}. Use a "
                        "qualified column reference like "
                        f"'{type(agg_node).__name__.upper()}(a.col)' or "
                        f"'{type(agg_node).__name__.upper()}(b.col)'."
                    )
                reject_right_side_if_left_only(col_table, agg_node.sql())
            rewritten_agg = self._rewrite_refs_for_wrapper_relation(
                agg_node, projection_map, sides, "aggregate argument"
            )
            outer_name = (
                user_alias
                if user_alias is not None
                else f"__giql_agg_{len(outer_projections)}"
            )
            outer_projections.append(f"{rewritten_agg} AS {q(outer_name)}")

        # Pre-allocate inner aliases for every column referenced in
        # GROUP BY / HAVING / ORDER BY so the wrapper relation exposes
        # them even if the user didn't put them in the SELECT list.
        for clause_key in ("group", "having", "order"):
            modifier_node = query.args.get(clause_key)
            if modifier_node is None:
                continue
            for col in modifier_node.find_all(exp.Column):
                col_table = _fold_identifier(col.table) if col.table else ""
                if col_table in (l_alias, r_alias):
                    reject_right_side_if_left_only(col_table, col.sql())
                    allocate(col_table, col.name)

        # Pre-allocate inner aliases for columns referenced in the WHERE
        # residuals routed to the outer wrapper filter (SEMI / ANTI, #200), so
        # the wrapper exposes them even when they appear only in the filter.
        # A right-side reference is rejected here (the left-only outer SELECT
        # cannot resolve it), matching the SELECT-list / modifier behavior.
        for residual in outer_where_residuals or []:
            for col in residual.find_all(exp.Column):
                col_table = _fold_identifier(col.table) if col.table else ""
                if col_table in (l_alias, r_alias):
                    reject_right_side_if_left_only(col_table, col.sql())
                    allocate(col_table, col.name)

        # Pre-allocate inner aliases for columns referenced inside any
        # aggregate function in the SELECT list (we'll rewrite the
        # aggregate's argument to use the inner alias when rendering the
        # outer projection below).
        for sel in query.expressions:
            target = _unwrap_projection_target(sel)
            if isinstance(target, exp.AggFunc):
                for col in target.find_all(exp.Column):
                    col_table = _fold_identifier(col.table) if col.table else ""
                    if col_table in (l_alias, r_alias):
                        reject_right_side_if_left_only(col_table, col.sql())
                        allocate(col_table, col.name)

        outer_projections: list[str] = []

        for sel in query.expressions:
            target = _unwrap_projection_target(sel)
            user_alias = sel.alias if isinstance(sel, exp.Alias) else None

            # Star projections (bare ``*`` and ``a.*`` / ``b.*``) never reach
            # here — :func:`_has_star_projection` declines them upstream (#202) —
            # and neither do the naive-valid shapes the projection pre-scan
            # above declined (expressions / window aggregates / FILTER clauses /
            # scalar subqueries / aggregate-nested stars — #204, #205).

            if isinstance(target, exp.AggFunc):
                render_aggregate(target, user_alias)
                continue

            if isinstance(target, exp.Column) and target.table:
                tbl = _fold_identifier(target.table)
                reject_right_side_if_left_only(tbl, target.sql())
                col = target.name
                alias = allocate(tbl, col)
                outer_name = user_alias if user_alias is not None else col
                outer_projections.append(f"{alias} AS {q(outer_name)}")
                continue

            # Everything reaching here is an unsupported wrapper (an
            # expression, a window aggregate, a FILTER clause, or an aggregate
            # nested in an expression) that the pre-scan did NOT decline because
            # it references an out-of-scope column — an unqualified column, an
            # unknown table, or the right side under a SEMI / ANTI join (e.g.
            # ``SUM(score) OVER (...)``, ``score + 1``). A plain aggregate with
            # an unqualified argument (``SUM(score)``) is diagnosed earlier in
            # ``render_aggregate``.
            #
            # A right-side column wrapped in an unsupported shape under a
            # left-only join gets the dedicated left-side-only message rather
            # than the generic "qualify the column" one, which would steer the
            # user toward another invalid form (``b.col`` is never valid here).
            for col in target.find_all(exp.Column):
                col_table = _fold_identifier(col.table) if col.table else ""
                reject_right_side_if_left_only(col_table, target.sql())
            # Otherwise the fault is an unqualified / unknown-table column; the
            # accurate guidance is to qualify it (once qualified, the wrapper
            # itself declines to the naive plan — #204, #205 — rather than
            # erroring).
            raise _UnqualifiedProjectionError(
                "dialect='duckdb' requires qualified projections; got "
                f"{target.sql()!r}. Use a qualified column reference like "
                "'a.col' or 'b.col [AS x]'."
            )

        if not outer_projections:
            raise _UnqualifiedProjectionError(
                "dialect='duckdb' requires at least one qualified projection."
            )

        # An aggregate without any source-column argument (e.g.
        # ``COUNT(*)``) leaves ``inner_projections`` empty in queries that
        # have no other column references. The per-chromosome subquery
        # still needs at least one column to project, so emit a constant
        # placeholder that DuckDB optimizes away. (Without this the
        # generated SQL ``SELECT  FROM (...) a JOIN (...) b ON ...`` would
        # fail to parse.)
        if not inner_projections:
            inner_projections.append("1 AS __giql_placeholder")

        return inner_projections, outer_projections, projection_map


def _has_sibling_spatial_predicate(node: Intersects, root: exp.Expression) -> bool:
    """Return True if *root* holds a spatial predicate other than the join *node*.

    The former pre-pass ran on the raw parsed AST, so a residual GIQL spatial
    predicate (``CONTAINS`` / ``WITHIN`` / ``ANY`` / ``ALL`` / a second
    ``INTERSECTS``) sharing the query was always still an un-expanded GIQL node —
    the transformer detected it and declined its IEJoin rewrite. In pass 3 the
    ``ExpandOperators`` walk expands nodes deepest-first, so a sibling spatial
    predicate may already have been rewritten to plain comparison SQL by the time
    this override runs, which the transformer's ``_classify_extras`` can no longer
    recognize (it would wrongly inline the residual into a per-chromosome ``ON``,
    corrupting SEMI/ANTI results — #169). Detect the sibling either way: as an
    un-expanded GIQL spatial node, or via the :data:`SPATIAL_PREDICATE_META` tag the
    generic spatial expanders stamp on their output. Either signal means the join
    must defer to the naive predicate, which composes correctly with any residual.
    """
    spatial_types = (Intersects, Contains, Within, SpatialSetPredicate)
    for candidate in root.walk():
        if candidate is node:
            continue
        if isinstance(candidate, spatial_types):
            return True
        if candidate.meta.get(SPATIAL_PREDICATE_META):
            return True
    return False


def _decline_to_generic(node: exp.Expression, ctx: ExpansionContext) -> exp.Expression:
    """Expand *node* with whatever expander serves ``(GenericTarget, Intersects)``.

    Going through the registry keeps this target's fallback and the generic
    target's expansion the same thing by construction. Calling the built-in
    directly made them the same only by coincidence: a user override registered
    on ``(GenericTarget, Intersects)`` took effect under ``dialect=None`` and was
    silently ignored under ``dialect="duckdb"`` for every shape this override
    declines.

    Only the generic slot is consulted. Resolving this target's own slot would
    return this function's caller and recurse without bound, since
    :meth:`~giql.expander.ExpanderRegistry.resolve` tries the exact
    ``(target, operator)`` entry first.

    Falls back to the built-in when the generic slot is empty, which a cleared
    registry produces.
    """
    generic = REGISTRY.resolve(GenericTarget(), Intersects)
    if generic is not None:
        return generic(node, ctx)
    return _expand_spatial_op(node, ctx, "intersects")


@register(DuckDBTarget, Intersects)
def expand_intersects_duckdb(
    node: exp.Expression, ctx: ExpansionContext
) -> exp.Expression:
    """Expand a DuckDB ``Intersects`` node — IEJoin join rewrite or naive predicate.

    The registered ``(DuckDBTarget, Intersects)`` override. For a column-to-column
    INTERSECTS *join* whose whole-query shape the IEJoin path supports,
    :meth:`IntersectsDuckDBIEJoinTransformer.transform_to_sql` builds the
    per-chromosome ``SET VARIABLE ...; SELECT ... FROM query(getvariable(...))``
    string; since that is a whole-query multi-statement rewrite, not a node-local
    expression, this registers a query-level finalizer
    (:meth:`~giql.expander.ExpansionContext.add_statement_finalizer`) that replaces
    the statement root with a :class:`~sqlglot.expressions.Command` wrapping the
    built SQL (which serializes verbatim), and returns the ``Intersects`` node
    unchanged (mirroring the CLUSTER / MERGE whole-query expanders).

    Every other shape — a join the IEJoin path declines (FULL OUTER, the LEFT /
    RIGHT shapes the decomposition declines, self-join, multiple INTERSECTS,
    extra predicates, non-base operands, 3+ tables), a literal-range predicate,
    or a residual column-to-column predicate — defers to whichever expander
    serves ``(GenericTarget, Intersects)`` (see :func:`_decline_to_generic`),
    which by default is the same naive overlap predicate every other target
    emits (#167). A query carrying any *sibling* spatial predicate also defers (see
    :func:`_has_sibling_spatial_predicate`), so the IEJoin never inlines a
    residual it cannot see. ``transform_to_sql`` then runs on the pass-3 AST, which
    for the shapes the IEJoin path accepts is structurally identical to the raw
    parse — those shapes provably carry no other GIQL operator, and pass 1 only
    attaches resolution metadata while pass 2 synthesizes no CTE for a
    column-to-column INTERSECTS (its operands are column slots, not reference
    slots) — so the rewrite reads the same query the former pre-pass did.
    """
    if isinstance(node, Intersects) and _is_column_intersects(node):
        root = node.root()
        if isinstance(root, exp.Select) and not _has_sibling_spatial_predicate(
            node, root
        ):
            transformer = IntersectsDuckDBIEJoinTransformer(ctx.tables)
            iejoin_sql = transformer.transform_to_sql(root)
            if iejoin_sql is not None:
                ctx.add_statement_finalizer(lambda _root: exp.Command(this=iejoin_sql))
                return node
    return _decline_to_generic(node, ctx)
