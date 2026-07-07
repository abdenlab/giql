"""Generic expanders for the spatial predicates and set predicates (epic #137).

Migrates INTERSECTS / CONTAINS / WITHIN and the ``ANY`` / ``ALL`` set predicates
off the legacy ``*_sql`` generator emitters and onto the operator-expander
registry. Each expander turns one predicate node
into standard sqlglot AST built from the pass-1 :class:`~giql.resolver.ResolvedColumn`
metadata (already canonicalized to 0-based half-open by pass 2), so the emitted SQL
is byte-identical to the strings the legacy emitter produced.

These are *node-local* predicate rewrites: an INTERSECTS / CONTAINS / WITHIN node
expands to a boolean ``(chrom = ... AND start < ... AND end > ...)`` expression that
replaces it in place. For a column-to-column INTERSECTS *join* this in-place
expansion IS the join plan on every target but DuckDB: the boolean overlap
predicate is a plain ``ON`` condition the engine plans as a range join (a hash join
keyed on ``chrom`` with the position inequalities as a residual filter), correct for
both inner and outer joins. The one whole-query pre-pass rewrite that survives is
the DuckDB IEJoin transformer in :mod:`giql.transformer`, gated on
``capabilities.range_join_strategy == "iejoin"``; it consumes a column-to-column
INTERSECTS *join* before this pass runs for the shapes it supports, and every shape
it declines falls through to this expander's naive predicate. (The generic binned
equi-join was dropped in favor of the naive predicate — #167.) A literal-range or
residual column-to-column INTERSECTS *predicate* (e.g. inside an ``OR``) reaches the
expander and is rendered the same way.

Only :class:`~giql.targets.GenericTarget` expanders are registered: spatial-predicate
*emission* is portable SQL-92 and does not vary by engine, so one generic expander
covers every target via the registry's ``(generic, op)`` fallback.
"""

from __future__ import annotations

from sqlglot import exp
from sqlglot import maybe_parse

from giql.dialect import GIQLDialect
from giql.expander import ExpansionContext
from giql.expander import register
from giql.expressions import Contains
from giql.expressions import Intersects
from giql.expressions import SpatialSetPredicate
from giql.expressions import Within
from giql.range_parser import ParsedRange
from giql.range_parser import RangeParser
from giql.resolver import ResolvedColumn
from giql.targets import GenericTarget


def _fragment(fragment: str) -> exp.Expression:
    """Parse a resolved SQL fragment (e.g. ``a."end"`` / ``'chr1'``) into AST.

    The pass-1 :class:`~giql.resolver.ResolvedColumn` carries column references as
    pre-canonicalized SQL string fragments; parse them through the GIQL dialect so
    the rebuilt predicate reserializes identically to the legacy emitter's string.
    """
    parsed = maybe_parse(fragment, dialect=GIQLDialect)
    if parsed is None:
        # maybe_parse returns None only for an empty/None input; a ResolvedColumn
        # fragment is never empty, so this is an internal invariant violation.
        raise ValueError(f"Could not parse resolved column fragment: {fragment!r}")
    return parsed


def _predicate_column(ctx: ExpansionContext, arg: str) -> ResolvedColumn:
    """Return the :class:`ResolvedColumn` for predicate operand *arg*.

    Mirrors the legacy ``_predicate_operand`` emitter helper: the
    expander consumes only the pass-1 resolution; a missing column means pass 1 did
    not run (an internal invariant violation), so raise the historical message.
    """
    resolution = ctx.resolution
    if resolution is not None:
        resolved = resolution.column(arg)
        if resolved is not None:
            return resolved
    raise ValueError(
        f"Spatial predicate operand {arg!r} was not resolved; run the "
        "ResolveOperatorRefs pass (transpile pipeline) before generation."
    )


def _range_predicate(
    column: ResolvedColumn, parsed: ParsedRange, op_type: str
) -> exp.Expression:
    """Build the boolean AST for ``column <op_type> <literal range>``.

    Reproduces the legacy ``_generate_range_predicate`` emitter as AST. The
    column fragments are already canonical 0-based half-open (pass 2); the parsed
    range is canonicalized by the caller. Returns a parenthesized boolean.
    """
    chrom = _fragment(column.chrom)
    start = _fragment(column.start)
    end = _fragment(column.end)
    chrom_lit = exp.Literal.string(parsed.chromosome)
    r_start = exp.Literal.number(parsed.start)
    r_end = exp.Literal.number(parsed.end)

    if op_type == "intersects":
        # Ranges overlap if: start1 < end2 AND end1 > start2
        cond = exp.and_(
            exp.EQ(this=chrom, expression=chrom_lit),
            exp.LT(this=start, expression=r_end),
            exp.GT(this=end, expression=r_start),
        )
    elif op_type == "contains":
        if parsed.end == parsed.start + 1:
            # Point query: start1 <= point < end1
            cond = exp.and_(
                exp.EQ(this=chrom, expression=chrom_lit),
                exp.LTE(this=start, expression=r_start),
                exp.GT(this=end, expression=r_start),
            )
        else:
            # Range query: start1 <= start2 AND end1 >= end2
            cond = exp.and_(
                exp.EQ(this=chrom, expression=chrom_lit),
                exp.LTE(this=start, expression=r_start),
                exp.GTE(this=end, expression=r_end),
            )
    elif op_type == "within":
        # left within right: start1 >= start2 AND end1 <= end2
        cond = exp.and_(
            exp.EQ(this=chrom, expression=chrom_lit),
            exp.GTE(this=start, expression=r_start),
            exp.LTE(this=end, expression=r_end),
        )
    else:
        raise ValueError(f"Unknown spatial op_type: {op_type!r}")

    return exp.paren(cond)


def _column_join(
    left: ResolvedColumn, right: ResolvedColumn, op_type: str
) -> exp.Expression:
    """Build the boolean AST for a column-to-column spatial predicate.

    Reproduces the legacy ``_generate_column_join`` emitter as AST. Both
    operands' fragments are pre-canonicalized (pass 2). Returns a parenthesized
    boolean.
    """
    l_chrom, r_chrom = _fragment(left.chrom), _fragment(right.chrom)
    l_start, r_start = _fragment(left.start), _fragment(right.start)
    l_end, r_end = _fragment(left.end), _fragment(right.end)

    if op_type == "intersects":
        cond = exp.and_(
            exp.EQ(this=l_chrom, expression=r_chrom),
            exp.LT(this=l_start, expression=r_end),
            exp.GT(this=l_end, expression=r_start),
        )
    elif op_type == "contains":
        cond = exp.and_(
            exp.EQ(this=l_chrom, expression=r_chrom),
            exp.LTE(this=l_start, expression=r_start),
            exp.GTE(this=l_end, expression=r_end),
        )
    elif op_type == "within":
        cond = exp.and_(
            exp.EQ(this=l_chrom, expression=r_chrom),
            exp.GTE(this=l_start, expression=r_start),
            exp.LTE(this=l_end, expression=r_end),
        )
    else:
        raise ValueError(f"Unknown spatial op_type: {op_type!r}")

    return exp.paren(cond)


def _expand_spatial_op(
    node: exp.Expression, ctx: ExpansionContext, op_type: str
) -> exp.Expression:
    """Expand one INTERSECTS / CONTAINS / WITHIN node to a boolean predicate.

    Dispatches on the right operand exactly as the legacy emitter did: the
    presence of a resolved right *column* — keyed off
    ``ctx.resolution.column("expression")``, the slot pass 1 attaches a
    :class:`ResolvedColumn` to when the right operand is a column reference —
    selects the column-to-column path; its absence means the right operand is a
    literal range, parsed in place.
    """
    resolution = ctx.resolution
    right_column = resolution.column("expression") if resolution is not None else None
    left = _predicate_column(ctx, "this")

    if right_column is not None:
        return _column_join(left, right_column, op_type)

    # Literal range string (e.g. interval INTERSECTS 'chr1:1000-2000'). Reproduce
    # the legacy emitter's parse-and-wrap-error behavior verbatim: any parse
    # failure (including the RangeParser's own ValueError) is wrapped in the
    # historical "Could not parse genomic range" message.
    right_expr = node.args.get("expression")
    raw = right_expr.sql(dialect=GIQLDialect) if right_expr is not None else ""
    try:
        range_str = raw.strip("'\"")
        parsed = RangeParser.parse(range_str).to_zero_based_half_open()
        return _range_predicate(left, parsed, op_type)
    except Exception as e:
        raise ValueError(f"Could not parse genomic range: {raw}. Error: {e}") from e


@register(GenericTarget, Intersects)
def expand_intersects(node: exp.Expression, ctx: ExpansionContext) -> exp.Expression:
    """Expand an INTERSECTS predicate to standard boolean SQL AST."""
    return _expand_spatial_op(node, ctx, "intersects")


@register(GenericTarget, Contains)
def expand_contains(node: exp.Expression, ctx: ExpansionContext) -> exp.Expression:
    """Expand a CONTAINS predicate to standard boolean SQL AST."""
    return _expand_spatial_op(node, ctx, "contains")


@register(GenericTarget, Within)
def expand_within(node: exp.Expression, ctx: ExpansionContext) -> exp.Expression:
    """Expand a WITHIN predicate to standard boolean SQL AST."""
    return _expand_spatial_op(node, ctx, "within")


@register(GenericTarget, SpatialSetPredicate)
def expand_spatial_set(node: exp.Expression, ctx: ExpansionContext) -> exp.Expression:
    """Expand a quantified set predicate (``ANY`` / ``ALL``) to boolean SQL AST.

    Reproduces the legacy ``_generate_spatial_set`` emitter: the single left
    column is compared against every literal range, and the per-range conditions
    are OR-combined for ``ANY`` / AND-combined for ``ALL``, all wrapped in one
    outer paren.
    """
    operator = node.args["operator"]
    quantifier = node.args["quantifier"]
    ranges = node.args["ranges"]

    column = _predicate_column(ctx, "this")
    op_type = operator.lower()

    conditions: list[exp.Expression] = []
    for range_expr in ranges:
        range_str = range_expr.sql(dialect=GIQLDialect).strip("'\"")
        parsed = RangeParser.parse(range_str).to_zero_based_half_open()
        conditions.append(_range_predicate(column, parsed, op_type))

    if quantifier.upper() == "ANY":
        combined = exp.or_(*conditions)
    else:
        combined = exp.and_(*conditions)

    return exp.paren(combined)
