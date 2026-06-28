"""The generic DISTANCE operator expander (epic #137, step / issue #140).

DISTANCE is the proof-of-concept that validates the expander protocol, the
registry dispatch, and the cross-target result-oracle workflow before the
harder operators migrate. It is the simplest operator — a single ``CASE``
expression with no joins, CTEs, or per-target divergence — so a single
*generic* expander registered for :class:`~giql.targets.GenericTarget` serves
every target. DISTANCE emits identical SQL on DuckDB, DataFusion, and the
generic baseline, so no per-target override is needed.

The CASE this expander builds matches
:meth:`giql.generators.base.BaseGIQLGenerator._generate_distance_case` exactly
(bedtools ``closest -d`` semantics): overlapping intervals report ``0``,
book-ended (adjacent) intervals report ``1``, and a raw half-open gap of ``N``
bases reports ``N + 1``. The ``+ 1`` is applied to the absolute gap magnitude
before any directional sign, so a downstream book-ended pair reports ``+1`` and
an upstream one ``-1`` in signed mode. There are four shapes — the cartesian
product of unsigned/signed and non-stranded/stranded — preserved verbatim from
the legacy emitter.

Because the returned CASE is reserialized by the active target's serializer
(rather than spliced in as a raw string, as the legacy ``giqldistance_sql``
emitter did), the emitted text changes cosmetically — most visibly ``!=``
renders as the SQL-standard ``<>``. The two are semantically identical.
"""

from __future__ import annotations

from sqlglot import exp
from sqlglot import parse_one

from giql.dialect import GIQLDialect
from giql.expander import ExpansionContext
from giql.expander import register
from giql.expressions import GIQLDistance
from giql.resolver import ResolvedColumn
from giql.targets import GenericTarget

__all__ = ["expand_distance"]


def _frag(fragment: str) -> exp.Expression:
    """Parse one canonicalized SQL fragment into AST.

    The pass-1 :class:`~giql.resolver.ResolvedColumn` endpoints are SQL string
    fragments (e.g. ``a."end"``, or ``'chr1'`` for a literal range) already
    canonicalized in place by pass 2, so they are parsed — not rebuilt — back
    into AST under the GIQL dialect.

    Parameters
    ----------
    fragment : str
        A canonicalized SQL fragment (a column reference or literal).

    Returns
    -------
    exp.Expression
        The parsed fragment as a sqlglot AST node.
    """
    return parse_one(fragment, dialect=GIQLDialect)


def _gap(minuend: str, subtrahend: str) -> exp.Expression:
    """Build ``(minuend - subtrahend + 1)`` — the bedtools-parity gap magnitude.

    Mirrors the legacy ``({start} - {end} + 1)`` fragment: the ``+ 1`` lifts a
    book-ended (adjacent) pair from a raw half-open gap of ``0`` to a reported
    distance of ``1``.

    Parameters
    ----------
    minuend : str
        The SQL fragment subtracted *from* (the left operand of the ``-``).
    subtrahend : str
        The SQL fragment subtracted (the right operand of the ``-``).

    Returns
    -------
    exp.Expression
        The parenthesized ``(minuend - subtrahend + 1)`` AST.
    """
    diff = exp.Sub(this=_frag(minuend), expression=_frag(subtrahend))
    return exp.paren(exp.Add(this=diff, expression=exp.Literal.number(1)))


def _bool_param(param: exp.Expression | None) -> bool:
    """Coerce an optional DISTANCE boolean argument to a Python ``bool``.

    Mirrors ``BaseGIQLGenerator._extract_bool_param`` so the expander reads the
    ``stranded`` / ``signed`` keyword arguments identically to the legacy
    emitter.

    Parameters
    ----------
    param : exp.Expression | None
        The ``stranded`` or ``signed`` argument node, or ``None`` if absent.

    Returns
    -------
    bool
        The coerced Python boolean (``False`` when the argument is absent).
    """
    if not param:
        return False
    if isinstance(param, exp.Boolean):
        return bool(param.this)
    return str(param).upper() in ("TRUE", "1", "YES")


def _operand(ctx: ExpansionContext, arg: str, position: str) -> ResolvedColumn:
    """Return the resolved column for one DISTANCE interval operand.

    Reads the pass-1 metadata attached to the node. A deferred operand (a
    literal range, or an unqualified column the resolver could not resolve) has
    no column attached; this raises the historical literal-range diagnostic so
    the public error contract is preserved.

    Parameters
    ----------
    ctx : ExpansionContext
        The expansion context carrying the node's pass-1 resolution.
    arg : str
        The operand slot key (``"this"`` or ``"expression"``).
    position : str
        Human-readable operand position (``"first"`` / ``"second"``) for the
        diagnostic message.

    Returns
    -------
    ResolvedColumn
        The resolved column metadata for the operand.

    Raises
    ------
    ValueError
        If the operand was deferred (a literal range or unresolved column).
    """
    # TODO(#146): this read-required-column-or-raise pattern is duplicated across
    # expanders; hoist it to a shared ``ExpansionContext.require_column`` helper
    # once a second expander needs it.
    resolution = ctx.resolution
    if resolution is not None:
        resolved = resolution.column(arg)
        if resolved is not None:
            return resolved
    raise ValueError(f"Literal range as {position} argument not yet supported")


def _unsigned_distance(col_a: ResolvedColumn, col_b: ResolvedColumn) -> exp.Expression:
    """Branch 1: unsigned (absolute) non-stranded distance, returning ``|gap| + 1``."""
    return _wrap_overlap_case(
        col_a,
        col_b,
        downstream=_gap(col_b.start, col_a.end),
        upstream=_gap(col_a.start, col_b.end),
    )


def _signed_distance(col_a: ResolvedColumn, col_b: ResolvedColumn) -> exp.Expression:
    """Branch 2: signed non-stranded distance.

    ``+`` downstream (B after A), ``-`` upstream (B before A).
    """
    return _wrap_overlap_case(
        col_a,
        col_b,
        downstream=_gap(col_b.start, col_a.end),
        upstream=exp.Neg(this=_gap(col_a.start, col_b.end)),
    )


def _stranded_distance(
    col_a: ResolvedColumn, col_b: ResolvedColumn, signed: bool
) -> exp.Expression:
    """Branches 3 & 4: stranded distance, flipping sign on A's ``-`` strand.

    The downstream and upstream gaps each become a nested ``CASE`` keyed on
    ``strand_a``. The ``signed`` flag additionally layers the directional sign
    on top of the strand flip, exactly as the legacy emitter's two stranded
    branches do.
    """
    strand_a = col_a.strand
    strand_b = col_b.strand
    assert strand_a is not None and strand_b is not None  # gated by caller

    down_gap = _gap(col_b.start, col_a.end)
    up_gap = _gap(col_a.start, col_b.end)

    # Downstream (B after A) is identical across the signed and unsigned arms:
    # positive by default, flipped negative on A's '-' strand. Only *upstream*
    # differs between the arms, so hoist downstream and compute upstream per arm.
    downstream = _strand_flip_case(
        strand_a, neg=exp.Neg(this=down_gap), pos=down_gap.copy()
    )
    if signed:
        # Stranded + signed: upstream (B before A) is negative by default but
        # flips positive on A's '-' strand (the directional sign layered on top
        # of the strand flip).
        upstream = _strand_flip_case(
            strand_a, neg=up_gap, pos=exp.Neg(this=up_gap.copy())
        )
    else:
        # Stranded but not signed: upstream carries the strand flip only.
        upstream = _strand_flip_case(
            strand_a, neg=exp.Neg(this=up_gap), pos=up_gap.copy()
        )

    case = _wrap_overlap_case(col_a, col_b, downstream=downstream, upstream=upstream)
    # Prepend the strand-validity guards ahead of the overlap guards. The WHEN
    # order matters: chrom mismatch, then strand NULL/'.'/'?', then overlap.
    return _prepend_strand_guards(case, strand_a, strand_b)


def _strand_flip_case(
    strand_a: str, neg: exp.Expression, pos: exp.Expression
) -> exp.Expression:
    """Build ``CASE WHEN strand_a = '-' THEN neg ELSE pos END``."""
    return (
        exp.Case()
        .when(exp.EQ(this=_frag(strand_a), expression=exp.Literal.string("-")), neg)
        .else_(pos)
    )


def _wrap_overlap_case(
    col_a: ResolvedColumn,
    col_b: ResolvedColumn,
    downstream: exp.Expression,
    upstream: exp.Expression,
) -> exp.Expression:
    """Build the shared distance CASE skeleton common to all four branches.

    ``CASE WHEN chrom_a != chrom_b THEN NULL WHEN <overlap> THEN 0 WHEN
    end_a <= start_b THEN <downstream> ELSE <upstream> END``.

    Parameters
    ----------
    col_a, col_b : ResolvedColumn
        The resolved A and B interval operands.
    downstream : exp.Expression
        The B-after-A gap expression (already carrying any sign / strand flip).
    upstream : exp.Expression
        The B-before-A gap expression (already carrying any sign / strand flip).

    Returns
    -------
    exp.Expression
        The assembled distance ``CASE`` expression.
    """
    overlap = exp.and_(
        exp.LT(this=_frag(col_a.start), expression=_frag(col_b.end)),
        exp.GT(this=_frag(col_a.end), expression=_frag(col_b.start)),
    )
    chrom_mismatch = exp.NEQ(this=_frag(col_a.chrom), expression=_frag(col_b.chrom))
    end_a_le_start_b = exp.LTE(this=_frag(col_a.end), expression=_frag(col_b.start))
    return (
        exp.Case()
        .when(chrom_mismatch, exp.Null())
        .when(overlap, exp.Literal.number(0))
        .when(end_a_le_start_b, downstream)
        .else_(upstream)
    )


def _prepend_strand_guards(
    case: exp.Case, strand_a: str, strand_b: str
) -> exp.Case:
    """Insert the strand-validity WHEN guards after the chrom guard.

    Distance is undefined for an unstranded (``'.'``/``'?'``) or missing strand,
    matching the legacy emitter.

    Parameters
    ----------
    case : exp.Case
        The overlap/gap ``CASE`` to prepend the strand guards onto (mutated in
        place).
    strand_a, strand_b : str
        The A and B strand-column SQL fragments.

    Returns
    -------
    exp.Case
        The same *case*, whose WHEN order is now: chrom mismatch -> NULL, either
        strand NULL -> NULL, ``strand_a`` is ``'.'``/``'?'`` -> NULL,
        ``strand_b`` is ``'.'``/``'?'`` -> NULL, then the original overlap/gap
        branches.
    """
    sa = _frag(strand_a)
    sb = _frag(strand_b)
    null_guard = exp.condition(exp.Is(this=sa, expression=exp.Null())).or_(
        exp.Is(this=sb.copy(), expression=exp.Null())
    )
    a_unstranded = exp.condition(
        exp.EQ(this=sa.copy(), expression=exp.Literal.string("."))
    ).or_(exp.EQ(this=sa.copy(), expression=exp.Literal.string("?")))
    b_unstranded = exp.condition(
        exp.EQ(this=sb.copy(), expression=exp.Literal.string("."))
    ).or_(exp.EQ(this=sb.copy(), expression=exp.Literal.string("?")))

    guards = [
        exp.If(this=null_guard, true=exp.Null()),
        exp.If(this=a_unstranded, true=exp.Null()),
        exp.If(this=b_unstranded, true=exp.Null()),
    ]
    # The existing WHENs keep their order; the strand guards slot in right after
    # the leading chrom-mismatch guard.
    existing = case.args["ifs"]
    case.set("ifs", existing[:1] + guards + existing[1:])
    return case


# KEEP IN SYNC: this expander and
# ``BaseGIQLGenerator._generate_distance_case`` (base.py) build the *same*
# distance CASE by two routes. The legacy method is retained only because
# NEAREST still calls it for its ORDER BY / filter math; once NEAREST migrates
# to the expander path that method can be deleted and this duplication retired.
# Until then, any change to the distance math here must be mirrored there (and
# vice versa). The parity test in tests/test_distance_udf.py guards the drift.
@register(GenericTarget, GIQLDistance)
def expand_distance(node: exp.Expression, ctx: ExpansionContext) -> exp.Expression:
    """Expand a ``GIQLDistance`` node into a standard SQL ``CASE`` expression.

    The portable expander registered for every target. Reads ``stranded`` /
    ``signed`` from the node and the pass-1 resolved interval operands, then
    builds the matching one of the four CASE shapes
    (unsigned/signed x non-stranded/stranded).

    Parameters
    ----------
    node : exp.Expression
        The ``GIQLDistance`` operator node being expanded.
    ctx : ExpansionContext
        The expansion context carrying the node's pass-1 resolution.

    Returns
    -------
    exp.Expression
        The distance ``CASE`` that replaces the operator node and is rendered by
        the active target's serializer.
    """
    stranded = _bool_param(node.args.get("stranded"))
    signed = _bool_param(node.args.get("signed"))

    col_a = _operand(ctx, "this", "first")
    col_b = _operand(ctx, "expression", "second")

    # Strand columns are consumed only in stranded mode, and only when both
    # operands actually carry a strand fragment — mirroring the legacy emitter's
    # `strand_a is None or strand_b is None` fall-through to the unstranded path.
    if stranded and col_a.strand is not None and col_b.strand is not None:
        return _stranded_distance(col_a, col_b, signed=signed)
    if signed:
        return _signed_distance(col_a, col_b)
    return _unsigned_distance(col_a, col_b)
