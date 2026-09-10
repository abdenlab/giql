"""The polars-bio (``dialect="datafusion-bio"``) spatial-predicate overrides (#77).

The registered ``(DataFusionBioTarget, Intersects / Contains / Within)`` expanders.
polars-bio ships DataFusion with an optimizer rule that replaces any hash /
nested-loop join whose filter is exactly two range comparisons (plus any equality
keys) with its native ``IntervalJoinExec``. The generic overlap predicate
``a.chrom = b.chrom AND a.start < b.end AND a.end > b.start`` already has that
shape, so vanilla-DataFusion output *plans* as an interval join -- and then meets
two upstream defects this module steers around. Everything it emits is valid
vanilla-DataFusion SQL; only the predicate shape changes.

**Defect 1 -- strict comparisons fail on Int64.** The rule's strict-comparison
adjustment builds an ``Int32`` literal without coercion, so the strict form fails
at execution on ``Int64`` coordinates (the pyarrow / Parquet default) with
``Invalid arithmetic operation: Int64 - Int32``. A column-to-column INTERSECTS
therefore emits the **closed, non-strict** form
``(l.chrom = r.chrom AND l.start <= r.end - 1 AND l.end - 1 >= r.start)``, which the
rule accepts unadjusted and which executes correctly on ``Int32``, ``Int64`` and
mixed-width sides, book-ended pairs excluded. CONTAINS / WITHIN already emit
non-strict comparisons, so their predicate text is unchanged.

**Defect 2 -- ``IntervalJoinExec`` ignores the join type.** Whatever join the rule
replaces, the operator produces the *inner* pairs: a LEFT / RIGHT / FULL join drops
its unmatched rows, a semi join (``WHERE EXISTS``) emits one row per match, and an
anti join (``WHERE NOT EXISTS``) returns the matched rows. Measured on every
``bio.interval_join_algorithm``. Row identity across engines is GIQL's contract,
so a predicate that would reach the operator through anything but an inner join is
emitted in a **rule-defeating** spelling: the same predicate plus one more
comparison the others already imply (``l.start < r.end`` for the closed overlap,
``l.start - 1 < r.start`` for CONTAINS, ``l.start + 1 > r.start`` for WITHIN). A
third range comparison makes the rule decline, and DataFusion plans the join as a
hash join on ``chrom`` with a residual filter -- correct, merely not accelerated.

One outer shape keeps the native join anyway: the count-overlaps shape
``SELECT <left keys>, COUNT(<right col>) AS n FROM l LEFT JOIN r ON <INTERSECTS>
GROUP BY <left keys>`` (the downstream cCRE insertion-matrix query, and the shape
DuckDB's IEJoin path accelerates as #209). A statement finalizer rewrites it as an
*inner* interval-join count zero-filled through a hash LEFT join against the
distinct left keys, so the interval work reaches ``IntervalJoinExec`` and the
unmatched keys still count 0.

**Defect 3 -- an empty projection fails.** ``SELECT COUNT(*)`` over an inner
interval join with no other projected column plans an ``IntervalJoinExec`` with
``projection=[]`` and fails (``must either specify a row count or at least one
column``). A finalizer rewrites that ``COUNT(*)`` to ``COUNT(<FROM-side start>)``;
GIQL assumes non-null coordinates throughout, so the count is unchanged.

All three workarounds retire when the upstream fixes land (tracked under #77). A
literal-range predicate (``interval INTERSECTS 'chr1:1000-2000'``) is a filter, not
a join, so it defers to whichever expander serves the ``GenericTarget`` slot
through the registry -- the same deferral the DuckDB override uses, so a user's
generic override reaches this target's declined shapes too.
"""

from __future__ import annotations

from typing import Callable

from sqlglot import exp

from giql.expander import REGISTRY
from giql.expander import ExpansionContext
from giql.expander import register
from giql.expanders.intersects import _column_join
from giql.expanders.intersects import _expand_spatial_op
from giql.expanders.intersects import _fragment
from giql.expanders.intersects import _predicate_column
from giql.expanders.intersects import _tag_spatial
from giql.expanders.intersects_duckdb import _CountOverlapsShape
from giql.expanders.intersects_duckdb import _match_count_overlaps
from giql.expressions import Contains
from giql.expressions import Intersects
from giql.expressions import Within
from giql.resolver import ResolvedColumn
from giql.targets import DataFusionBioTarget
from giql.targets import GenericTarget

#: Join kinds that keep inner semantics, the only ones ``IntervalJoinExec``
#: executes correctly. ``CROSS`` folds to inner once the predicate is applied.
_INNER_KINDS = frozenset({"", "INNER", "CROSS"})


def _one() -> exp.Expression:
    return exp.Literal.number(1)


def _closed_overlap(left: ResolvedColumn, right: ResolvedColumn) -> exp.Expression:
    """Build the non-strict ``l.start <= r.end - 1 AND l.end - 1 >= r.start`` overlap.

    Both operands' fragments are pre-canonicalized to 0-based half-open (pass 2),
    so subtracting one from the exclusive end yields the inclusive last base and
    the ``<=`` / ``>=`` comparisons exclude book-ended pairs exactly as the strict
    form does. Returns the unparenthesized conjunction.
    """
    return exp.and_(
        exp.EQ(this=_fragment(left.chrom), expression=_fragment(right.chrom)),
        exp.LTE(
            this=_fragment(left.start),
            expression=exp.Sub(this=_fragment(right.end), expression=_one()),
        ),
        exp.GTE(
            this=exp.Sub(this=_fragment(left.end), expression=_one()),
            expression=_fragment(right.start),
        ),
    )


def _implied_comparison(
    left: ResolvedColumn, right: ResolvedColumn, op_type: str
) -> exp.Expression:
    """Build the redundant third comparison that makes the interval-join rule decline.

    Each is implied by the operator's own predicate on integer coordinates, so
    adding it never changes the rows -- only the predicate's shape.
    """
    l_start, r_start = _fragment(left.start), _fragment(right.start)
    if op_type == "intersects":
        # closed: l.start <= r.end - 1  =>  l.start < r.end
        return exp.LT(this=l_start, expression=_fragment(right.end))
    if op_type == "contains":
        # l.start <= r.start  =>  l.start - 1 < r.start
        return exp.LT(this=exp.Sub(this=l_start, expression=_one()), expression=r_start)
    if op_type == "within":
        # l.start >= r.start  =>  l.start + 1 > r.start
        return exp.GT(this=exp.Add(this=l_start, expression=_one()), expression=r_start)
    raise ValueError(f"Unknown spatial op_type: {op_type!r}")


def _column_predicate(
    left: ResolvedColumn, right: ResolvedColumn, op_type: str, accelerate: bool
) -> exp.Expression:
    """Build the column-to-column predicate for *op_type* in the requested shape.

    *accelerate* selects the rule-eligible spelling (closed overlap for
    INTERSECTS, the generic form for CONTAINS / WITHIN); otherwise the implied
    third comparison is appended so the rule declines. Returns a parenthesized,
    spatial-tagged boolean.
    """
    if op_type == "intersects":
        cond = _closed_overlap(left, right)
    else:
        # The generic form is parenthesized; unwrap so the conjunct appends flat.
        cond = _column_join(left, right, op_type).this
    if not accelerate:
        conjuncts = list(cond.flatten())
        cond = exp.and_(*conjuncts, _implied_comparison(left, right, op_type))
    return _tag_spatial(exp.paren(cond))


def _decline_to_generic(
    node: exp.Expression, ctx: ExpansionContext, operator: type, op_type: str
) -> exp.Expression:
    """Expand *node* with whatever expander serves ``(GenericTarget, operator)``.

    Only the generic slot is consulted (resolving this target's own slot would
    recurse into this module); the built-in is the last resort for a cleared
    registry.
    """
    generic = REGISTRY.resolve(GenericTarget(), operator)
    if generic is not None:
        return generic(node, ctx)
    return _expand_spatial_op(node, ctx, op_type)


def _is_inner_join(join: exp.Join) -> bool:
    """Return True if *join* has inner semantics (no side, an inner kind)."""
    if join.args.get("side"):
        return False
    return (join.args.get("kind") or "").upper() in _INNER_KINDS


def _reaches_operator_as_inner(node: exp.Expression) -> bool:
    """Return True if *node*'s predicate would reach the join as an inner join.

    A predicate in a join's ``ON`` follows that join's kind. A predicate in a
    ``WHERE`` is an inner filter unless the enclosing SELECT is itself a
    correlated subquery -- an ``EXISTS`` / ``IN`` operand or a scalar subquery --
    which DataFusion decorrelates into a semi / anti / outer join before the rule
    sees it. A FROM-clause subquery or a CTE body keeps inner semantics.
    """
    join = node.find_ancestor(exp.Join)
    if join is not None:
        return _is_inner_join(join)
    for ancestor in _ancestors(node):
        if isinstance(ancestor, (exp.Exists, exp.In)):
            return False
        if isinstance(ancestor, exp.Subquery) and not isinstance(
            ancestor.parent, (exp.From, exp.Join)
        ):
            return False
    return True


def _ancestors(node: exp.Expression):
    parent = node.parent
    while parent is not None:
        yield parent
        parent = parent.parent


def _is_bare_count_star_select(select: exp.Select) -> bool:
    """Return True if *select* projects only ``COUNT(*)``-style items and no column.

    The shape whose join projection DataFusion prunes to nothing: no ``GROUP BY``
    and no column reference anywhere in the SELECT list, with at least one
    ``COUNT(*)`` to rewrite.
    """
    if select.args.get("group") is not None:
        return False
    has_count_star = False
    for item in select.expressions:
        for sub in item.walk():
            if isinstance(sub, exp.Column):
                return False
            if isinstance(sub, exp.Count) and isinstance(sub.this, exp.Star):
                has_count_star = True
    return has_count_star


def _from_side_column(
    select: exp.Select, join: exp.Join, left: ResolvedColumn, right: ResolvedColumn
) -> ResolvedColumn | None:
    """Return the operand living on the FROM side of the sole inner join, or None."""
    joins = select.args.get("joins") or []
    if len(joins) != 1 or joins[0] is not join or not _is_inner_join(join):
        return None
    from_ = select.args.get("from_")
    if from_ is None or not isinstance(from_.this, exp.Table):
        return None
    from_alias = from_.this.alias_or_name.casefold()
    for column in (left, right):
        start = _fragment(column.start)
        if isinstance(start, exp.Column) and start.table.casefold() == from_alias:
            return column
    return None


def _count_star_finalizer(start_fragment: str) -> Callable:
    """Build the finalizer that rewrites every root-level ``COUNT(*)`` to a column."""

    def _finalize(root: exp.Expression) -> exp.Expression:
        if isinstance(root, exp.Select):
            for item in root.expressions:
                for count in item.find_all(exp.Count):
                    if isinstance(count.this, exp.Star):
                        count.set("this", _fragment(start_fragment))
        return root

    return _finalize


def _output_identifier(item: exp.Expression) -> exp.Identifier:
    """Return the identifier a SELECT item is output under (alias or column name)."""
    if isinstance(item, exp.Alias):
        return item.args["alias"].copy()
    if isinstance(item, exp.Column):
        return item.this.copy()
    return exp.to_identifier(item.output_name)


def _zero_fill_finalizer(
    shape: _CountOverlapsShape, keys_alias: str, counts_alias: str
) -> Callable:
    """Build the finalizer for the count-overlaps LEFT-join shape.

    Receives the root after node replacement (so the join ``ON`` already holds the
    closed overlap form) and rebuilds it as::

        SELECT k.<key>..., COALESCE(x.<n>, 0) AS <n>
        FROM (SELECT DISTINCT <key exprs> FROM <left>) AS k
        LEFT JOIN (<root with the LEFT join turned INNER>) AS x
          ON k.<key> = x.<key> AND ...

    The inner count reaches ``IntervalJoinExec``; the zero-fill is an equi hash
    join on the user's own key names, so the output columns and their order match
    the original SELECT list. Keys are compared with ``=`` (``IS NOT DISTINCT
    FROM`` does not type-coerce on DataFusion 53), so a NULL key zero-fills rather
    than matching -- the same behavior as DuckDB's ``USING`` zero-fill (#209).
    """

    def _finalize(root: exp.Expression) -> exp.Expression:
        if not isinstance(root, exp.Select):
            return root
        inner = root.copy()
        join = inner.args["joins"][0]
        join.set("side", None)
        join.set("kind", None)
        # Pin every inner output name so the wrapper can address it verbatim.
        inner.set(
            "expressions",
            [
                exp.alias_(
                    (item.this if isinstance(item, exp.Alias) else item).copy(),
                    _output_identifier(item),
                )
                for item in inner.expressions
            ],
        )

        key_items = list(shape.key_items)
        key_idents = [_output_identifier(item) for item in key_items]
        keys = (
            exp.select(
                *[
                    exp.alias_(
                        (item.this if isinstance(item, exp.Alias) else item).copy(),
                        ident.copy(),
                    )
                    for item, ident in zip(key_items, key_idents)
                ]
            )
            .distinct()
            .from_(shape.left_table.copy())
        )
        on = exp.and_(
            *[
                exp.EQ(
                    this=exp.column(ident.copy(), table=keys_alias),
                    expression=exp.column(ident.copy(), table=counts_alias),
                )
                for ident in key_idents
            ]
        )

        outer_items = []
        for item in root.expressions:
            ident = _output_identifier(item)
            if item is shape.agg_item:
                value = exp.Coalesce(
                    this=exp.column(ident.copy(), table=counts_alias),
                    expressions=[exp.Literal.number(0)],
                )
            else:
                value = exp.column(ident.copy(), table=keys_alias)
            outer_items.append(exp.alias_(value, ident.copy()))

        return (
            exp.select(*outer_items)
            .from_(keys.subquery(keys_alias))
            .join(inner.subquery(counts_alias), on=on, join_type="LEFT")
        )

    return _finalize


def _expand_column_predicate(
    node: exp.Expression, ctx: ExpansionContext, operator: type, op_type: str
) -> exp.Expression:
    """Shared body of the three overrides; see the module docstring for the rules."""
    resolution = ctx.resolution
    right = resolution.column("expression") if resolution is not None else None
    if right is None:
        return _decline_to_generic(node, ctx, operator, op_type)
    left = _predicate_column(ctx, "this")

    root = node.root()
    join = node.find_ancestor(exp.Join)

    if op_type == "intersects" and isinstance(root, exp.Select):
        shape = _match_count_overlaps(root)
        if shape is not None and join is not None and join.args.get("on") is not None:
            ctx.add_statement_finalizer(
                _zero_fill_finalizer(shape, ctx.alias(), ctx.alias())
            )
            return _column_predicate(left, right, op_type, accelerate=True)

    accelerate = _reaches_operator_as_inner(node)
    if (
        accelerate
        and join is not None
        and isinstance(root, exp.Select)
        and _is_bare_count_star_select(root)
    ):
        from_side = _from_side_column(root, join, left, right)
        if from_side is not None:
            ctx.add_statement_finalizer(_count_star_finalizer(from_side.start))
    return _column_predicate(left, right, op_type, accelerate=accelerate)


@register(DataFusionBioTarget, Intersects)
def expand_intersects_datafusion_bio(
    node: exp.Expression, ctx: ExpansionContext
) -> exp.Expression:
    """Expand a datafusion-bio ``Intersects`` node (closed form, gated by join kind)."""
    return _expand_column_predicate(node, ctx, Intersects, "intersects")


@register(DataFusionBioTarget, Contains)
def expand_contains_datafusion_bio(
    node: exp.Expression, ctx: ExpansionContext
) -> exp.Expression:
    """Expand a datafusion-bio ``Contains`` node (generic form, gated by join kind)."""
    return _expand_column_predicate(node, ctx, Contains, "contains")


@register(DataFusionBioTarget, Within)
def expand_within_datafusion_bio(
    node: exp.Expression, ctx: ExpansionContext
) -> exp.Expression:
    """Expand a datafusion-bio ``Within`` node (generic form, gated by join kind)."""
    return _expand_column_predicate(node, ctx, Within, "within")
