"""Pass 1 of the GIQL normalization pipeline: ``ResolveOperatorRefs``.

This module implements the first normalization pass described in epic #114. It
runs between the transformer chain and the generator and attaches resolution
metadata to every GIQL operator node so that later epic steps can turn the
generator into a pure formatter.

The pass is built on ``sqlglot``'s scope machinery
(:func:`sqlglot.optimizer.scope.traverse_scope`), which yields scopes
leaves-first so CTE and derived-table sources are registered before their
consumers, and whose ``cte_sources`` natively distinguishes an in-query CTE from
a registered table — exactly the registered-table / CTE / subquery trichotomy
each operator reference slot must resolve against.

For every GIQL operator node the pass resolves each declared *reference slot*
(see :class:`giql.expressions.SlotSpec`) against the enclosing scope plus the
registered :class:`giql.table.Tables` container and attaches a frozen
:class:`ResolvedRef` to the node through a single namespaced
``node.meta["giql"]`` key — ``sqlglot``'s established metadata channel, which is
deep-copied by ``Expression.copy()`` so resolution survives tree copies. The
pass closes with a validation boundary (:func:`validate_operator_refs`) that
asserts every operator slot carries well-formed resolution metadata, mirroring
``sqlglot``'s ``validate_qualify_columns`` and Spark's ``CheckAnalysis``.

Scope note (epic #114 / #137)
-----------------------------
The pass is behavior-preserving. Every operator in the ``_OPERATORS`` roster now
takes the ``ExpandOperators`` path (epic #137): the DISJOIN, NEAREST, DISTANCE,
and spatial-predicate (INTERSECTS / CONTAINS / WITHIN / ``SpatialSetPredicate``)
expanders in :mod:`giql.expanders` consume the metadata attached here, resolving
their reference slots and column operands through it. CLUSTER and MERGE (#144)
declare no reference slots, so they resolve to an empty-but-valid
:class:`OperatorResolution` and their expanders derive columns from the enclosing
FROM table rather than from resolution. The resolution semantics computed here
mirror the generator's historical ``_resolve_target_table`` /
``_resolve_disjoin_reference`` / ``_enclosing_cte_names`` (DISJOIN) and
``_resolve_nearest_reference`` / ``_find_outer_table_in_lateral_join`` (NEAREST)
behavior exactly; all of those helpers now live only here.

Two consequences of the zero-behavior-change constraint shape the
implementation:

* The pass never raises a *user-facing* diagnostic. When a slot cannot be
  resolved (unregistered target, unknown reference name, reserved ``__giql_dj_``
  prefix, unsupported reference shape) the pass simply leaves that slot
  unresolved and the generator raises its existing error exactly as before.
  Promoting these into GIQL-level diagnostics is a later epic step.
* Table-shaped *reference slots* (DISJOIN ``this``/``reference`` and NEAREST
  ``this``) resolve to a :class:`ResolvedRef`. NEAREST's ``reference`` slot —
  whose accepted shapes are the non-table literal-range / column /
  implicit-outer forms — resolves to a :class:`ResolvedInterval` (epic #114,
  step 3). DISTANCE's two *column* operands (``this`` / ``expression``) are
  resolved to a :class:`ResolvedColumn` and attached through the separate
  :attr:`OperatorResolution.columns` channel (epic #114, step 4). The spatial
  predicates' column operands (:class:`~giql.expressions.Intersects` /
  ``Contains`` / ``Within`` ``this`` and column-shaped ``expression``, and
  ``SpatialSetPredicate`` ``this``) resolve onto the same
  :class:`ResolvedColumn` type through that columns channel (epic #114,
  step #120); a literal-range ``expression`` slot stays deferred and the emitter
  parses it on its existing path.
"""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field
from dataclasses import replace
from typing import Literal

from sqlglot import exp
from sqlglot.optimizer.scope import Scope
from sqlglot.optimizer.scope import traverse_scope

from giql.constants import DEFAULT_CHROM_COL
from giql.constants import DEFAULT_END_COL
from giql.constants import DEFAULT_START_COL
from giql.constants import DEFAULT_STRAND_COL
from giql.constants import DJ_PREFIX
from giql.expressions import Contains
from giql.expressions import GIQLCluster
from giql.expressions import GIQLDisjoin
from giql.expressions import GIQLDistance
from giql.expressions import GIQLMerge
from giql.expressions import GIQLNearest
from giql.expressions import Intersects
from giql.expressions import SlotSpec
from giql.expressions import SpatialSetPredicate
from giql.expressions import Within
from giql.range_parser import RangeParser
from giql.table import Table
from giql.table import Tables

__all__ = [
    "META_KEY",
    "RefKind",
    "IntervalKind",
    "ResolvedRef",
    "ResolvedInterval",
    "ResolvedColumn",
    "SlotDeferral",
    "OperatorResolution",
    "ResolutionError",
    "resolve_operator_refs",
    "validate_operator_refs",
    "validate_projection_contracts",
]

#: The single namespaced key under which resolution metadata is stored on a
#: node's ``.meta`` mapping. All GIQL metadata lives under this one key so the
#: channel stays a single, greppable namespace.
META_KEY = "giql"

#: The kinds a resolved reference slot can take — the registered-table / CTE /
#: subquery trichotomy that ``sqlglot``'s ``scope.sources`` distinguishes.
RefKind = Literal["registered_table", "cte", "subquery"]

#: The kinds a resolved *interval* slot can take. Unlike a :data:`RefKind`, an
#: interval slot does not resolve to a whole relation but to a column-qualified
#: genomic interval (``column`` / ``implicit_outer``) or a parsed literal range
#: (``literal_range``). These are the non-table shapes NEAREST's ``reference``
#: slot accepts; DISTANCE and the spatial predicates (#119/#120) resolve the
#: same shapes onto the same :class:`ResolvedInterval` type.
IntervalKind = Literal["literal_range", "column", "implicit_outer"]

#: Canonical default genomic column names. A CTE or subquery reference is
#: assumed (in the existing generator and here) to expose canonical
#: ``chrom`` / ``start`` / ``end`` columns; validating that contract is a later
#: epic step.
_DEFAULT_COLS: tuple[str, str, str] = (
    DEFAULT_CHROM_COL,
    DEFAULT_START_COL,
    DEFAULT_END_COL,
)

#: The GIQL operator expression classes the pass inspects. CLUSTER and MERGE
#: (#144) declare no reference slots, so they resolve to an empty-but-valid
#: ``OperatorResolution``; the ``ExpandOperators`` pass requires that metadata to
#: be present, and their expanders derive columns from the enclosing FROM table
#: rather than from resolution.
_OPERATORS: tuple[type[exp.Expression], ...] = (
    GIQLDisjoin,
    GIQLNearest,
    GIQLDistance,
    GIQLCluster,
    GIQLMerge,
    Intersects,
    Contains,
    Within,
    SpatialSetPredicate,
)


class ResolutionError(Exception):
    """Raised when attached resolution metadata is malformed.

    This is an internal invariant check fired by :func:`validate_operator_refs`,
    not a user-facing GIQL diagnostic. It indicates a bug in the resolver pass —
    a slot carrying metadata that does not satisfy its declared
    :class:`~giql.expressions.SlotSpec` — rather than an invalid user query.
    """


@dataclass(frozen=True, slots=True)
class ResolvedRef:
    """Resolved metadata for one table-shaped operator reference slot.

    Attributes
    ----------
    kind : RefKind
        Whether the slot resolved to a registered table, an in-query CTE, or a
        subquery.
    name : str | None
        The relation name for a registered-table or CTE reference; ``None`` for
        an (anonymous) subquery reference.
    cols : tuple[str, str, str]
        The ``(chrom, start, end)`` physical column names. For a registered
        table these come from its :class:`~giql.table.Table` config; for a CTE
        or subquery they are the canonical defaults the relation is assumed to
        expose.
    table : Table | None
        The :class:`~giql.table.Table` config backing a registered-table
        reference (carrying its coordinate system); ``None`` for CTE and
        subquery references, which are assumed canonical.
    coverage_skippable : bool
        Whether this reference is provably the same relation as the operator's
        target set — i.e. a self-reference — which lets DISJOIN omit its
        coverage ``EXISTS`` filter. Always ``False`` for non-reference (target)
        slots.
    """

    kind: RefKind
    name: str | None
    cols: tuple[str, str, str]
    table: Table | None
    coverage_skippable: bool


@dataclass(frozen=True, slots=True)
class ResolvedInterval:
    """Resolved metadata for one non-table interval operator slot.

    Where :class:`ResolvedRef` resolves a slot to a whole relation, a
    ``ResolvedInterval`` resolves a slot to a single genomic interval expressed
    as SQL fragments — either column references qualified by a relation alias
    (``column`` / ``implicit_outer``) or the literal endpoints of a parsed
    genomic range (``literal_range``). NEAREST's ``reference`` slot is the first
    consumer; DISTANCE's two operands and the spatial predicates' column / range
    slots (#119, #120) resolve onto this same type.

    Coordinate canonicalization deliberately stays in the emitter for now
    (epic #114 step 8 / #123). The ``start`` / ``end`` fields therefore carry
    the *raw, un-canonicalized* column SQL for the ``column`` /
    ``implicit_outer`` kinds — the emitter wraps them with
    :func:`giql.canonical.canonical_start` / ``canonical_end`` using
    :attr:`table`. For the ``literal_range`` kind the endpoints are already
    canonical 0-based half-open integer literals (parsed via
    :meth:`~giql.range_parser.ParsedRange.to_zero_based_half_open`) and
    :attr:`table` is ``None``; the emitter uses them verbatim.

    Attributes
    ----------
    kind : IntervalKind
        Whether the interval came from a literal range, an explicit column
        reference, or an implicit LATERAL outer relation.
    chrom : str
        SQL fragment for the chromosome — a quoted literal (``'chr1'``) for a
        ``literal_range`` or a qualified column (``alias."chrom"``) otherwise.
    start : str
        SQL fragment for the start endpoint. Canonical literal for a
        ``literal_range``; raw column SQL (to be canonicalized by the emitter)
        otherwise.
    end : str
        SQL fragment for the end endpoint, mirroring :attr:`start`.
    strand : str | None
        SQL fragment for the strand — a quoted literal for a ``literal_range``
        (``None`` when the range carries no strand), a qualified column for an
        explicit ``column`` reference (always present), or a qualified column
        for an ``implicit_outer`` reference only when the outer table declares a
        strand column (``None`` otherwise). The emitter reads it only in
        stranded mode.
    table : Table | None
        The :class:`~giql.table.Table` config backing a ``column`` /
        ``implicit_outer`` interval, carried so the emitter can canonicalize the
        endpoints; ``None`` for a ``literal_range`` (already canonical) or when
        the qualifier resolves to no registered table.
    """

    kind: IntervalKind
    chrom: str
    start: str
    end: str
    strand: str | None
    table: Table | None


@dataclass(frozen=True, slots=True)
class ResolvedColumn:
    """Resolved metadata for one column-shaped interval operand.

    Models the resolution of a DISTANCE operand (``a.interval``) — an
    ``exp.Column`` qualified by a table alias — into the physical genomic
    columns it references, qualified by that alias. Unlike a
    :class:`ResolvedRef` (which names a whole relation), a column operand
    resolves to concrete SQL column expressions ready to drop into the
    emitter's distance arithmetic.

    Coordinate canonicalization stays in the emitter (epic #114 defers its
    removal to step 8 / issue #123); this metadata therefore carries the
    backing :class:`~giql.table.Table` so the emitter keeps wrapping the
    endpoints with :func:`giql.canonical.canonical_start` /
    :func:`giql.canonical.canonical_end` exactly as before.

    Attributes
    ----------
    chrom : str
        The chromosome column qualified by the operand's alias, e.g.
        ``a."chrom"``.
    start : str
        The start column qualified by the operand's alias, e.g. ``a."start"``.
    end : str
        The end column qualified by the operand's alias, e.g. ``a."end"``.
    strand : str | None
        The strand column qualified by the operand's alias, e.g.
        ``a."strand"``. Always resolved when the operand is a column (mirroring
        the generator's ``_get_column_refs(..., include_strand=True)``); the
        emitter consumes it only in stranded mode.
    table : Table | None
        The :class:`~giql.table.Table` config backing the operand's relation
        (carrying its coordinate system), or ``None`` when the alias does not
        resolve to a registered table (an unregistered relation is assumed
        canonical, exactly as the generator's ``_resolve_table`` returns
        ``None``).
    """

    chrom: str
    start: str
    end: str
    strand: str | None
    table: Table | None


@dataclass(frozen=True, slots=True)
class SlotDeferral:
    """Why a slot was left unresolved, so the emitter raises the right error.

    The resolver is behavior-preserving during the epic #114 migration: when a
    slot cannot be resolved it is *deferred* and the generator re-raises its
    historical diagnostic. Some of those diagnostics need information the
    emitter can no longer recompute on its own (the resolver moved the relevant
    scope/ancestor walk out of the generator). A ``SlotDeferral`` carries that
    information forward.

    NEAREST's implicit-outer reference is the motivating case: distinguishing
    "no outer relation found in the LATERAL context" from "outer relation found
    but not registered" requires the ancestor walk that now lives in this pass,
    so the resolver records which failure occurred and, for the latter, the
    offending relation label. DISTANCE and the spatial predicates (#119/#120)
    reuse this channel for their own deferred-shape diagnostics.

    Attributes
    ----------
    reason : str
        A stable machine token naming the deferral cause (e.g.
        ``"implicit_outer_missing"``, ``"implicit_outer_unregistered"``).
    detail : str | None
        Optional supporting datum for the emitter's message (e.g. the outer
        relation label for ``"implicit_outer_unregistered"``).
    """

    reason: str
    detail: str | None = None


@dataclass(frozen=True, slots=True)
class OperatorResolution:
    """Resolution metadata attached to a single GIQL operator node.

    Attributes
    ----------
    operator : str
        The operator expression class name (e.g. ``"GIQLDisjoin"``).
    slots : dict[str, ResolvedRef | ResolvedInterval]
        Mapping from a slot's :attr:`~giql.expressions.SlotSpec.arg` key to its
        resolved metadata — a :class:`ResolvedRef` for a table-shaped reference
        slot, or a :class:`ResolvedInterval` for an interval slot. Only
        successfully resolved slots appear; a slot left out was either not a
        resolvable slot or could not be resolved (and the generator raises its
        existing error, possibly aided by :attr:`deferrals`).
    deferrals : dict[str, SlotDeferral]
        Mapping from a slot key to a :class:`SlotDeferral` recording why the
        slot was deferred, when the emitter needs that context to raise the
        historical diagnostic verbatim. Empty for slots that resolved or whose
        deferral the emitter can reclassify unaided.
    columns : dict[str, ResolvedColumn]
        Mapping from a *column* operand's arg key to its resolved
        :class:`ResolvedColumn`. Carries DISTANCE's two interval operands and
        the spatial predicates' column operands; an operand the pass could not
        resolve (a literal range, or an unqualified column outside a
        current-table context) is left out and the generator raises its existing
        error.
    output_tables : dict[str, Table]
        Mapping from a slot's arg key to the *original* :class:`~giql.table.Table`
        whose declared encoding that slot carried before
        :func:`giql.canonicalizer.canonicalize_coordinates` (pass 2) wrapped it
        in a canonical CTE and blanked its ``ResolvedRef.table``. The pass
        populates this for every slot it wraps so the operator's emitter can
        round-trip *synthesized* output columns (DISJOIN's ``disjoin_*`` and its
        passed-through interval) back into that encoding — columns that do not
        exist as AST in pass 2 and that a ``SELECT *`` consumer hides from any
        outer-projection rewrite. Empty until pass 2 wraps a slot.
    """

    operator: str
    slots: dict[str, ResolvedRef | ResolvedInterval]
    deferrals: dict[str, SlotDeferral] = field(default_factory=dict)
    columns: dict[str, ResolvedColumn] = field(default_factory=dict)
    output_tables: dict[str, Table] = field(default_factory=dict)

    def slot(self, arg: str) -> ResolvedRef | ResolvedInterval | None:
        """Return the resolved metadata for slot *arg*, or ``None``."""
        return self.slots.get(arg)

    def deferral(self, arg: str) -> SlotDeferral | None:
        """Return the deferral recorded for slot *arg*, or ``None``."""
        return self.deferrals.get(arg)

    def column(self, arg: str) -> ResolvedColumn | None:
        """Return the resolved column operand for *arg*, or ``None``."""
        return self.columns.get(arg)


def resolve_operator_refs(expression: exp.Expression, tables: Tables) -> exp.Expression:
    """Attach resolution metadata to every GIQL operator in *expression*.

    Walks the tree's scopes leaves-first, resolves each operator's declared
    reference slots against the enclosing scope's visible CTEs plus *tables*,
    attaches an :class:`OperatorResolution` under ``node.meta["giql"]``, and
    closes with :func:`validate_operator_refs`.

    The pass mutates and returns *expression* in place. It is behavior-
    preserving: the attached metadata is additive and the generator ignores it.

    Parameters
    ----------
    expression : exp.Expression
        The (post-transformer) AST to annotate.
    tables : Tables
        The registered table configurations.

    Returns
    -------
    exp.Expression
        The same *expression*, with metadata attached.
    """
    seen: set[int] = set()

    for scope in _safe_traverse_scope(expression):
        cte_names = frozenset(scope.cte_sources)
        for node in scope.walk():
            if isinstance(node, _OPERATORS) and id(node) not in seen:
                seen.add(id(node))
                _resolve_operator(node, tables, cte_names)

    # Fallback for any operator a scope walk did not reach (e.g. if scope
    # construction failed). Resolving with no visible CTE names keeps the pass
    # behavior-preserving: a missed CTE reference or column operand simply stays
    # unresolved and the generator handles it on its existing path.
    for node in expression.walk():
        if isinstance(node, _OPERATORS) and id(node) not in seen:
            seen.add(id(node))
            _resolve_operator(node, tables, frozenset())

    validate_operator_refs(expression)
    validate_projection_contracts(expression)
    return expression


def _safe_traverse_scope(expression: exp.Expression) -> list[Scope]:
    """Return the leaves-first scopes of *expression*, or ``[]`` on failure.

    ``traverse_scope`` raises on a handful of exotic constructs. Swallowing the
    failure keeps the pass behavior-preserving: the fallback walk in
    :func:`resolve_operator_refs` still annotates every operator, just without
    scoped CTE visibility.
    """
    try:
        return traverse_scope(expression)
    except Exception:
        return []


def _resolve_operator(
    node: exp.Expression, tables: Tables, cte_names: frozenset[str]
) -> None:
    """Resolve *node*'s slots and attach an :class:`OperatorResolution`."""
    slots: dict[str, ResolvedRef | ResolvedInterval] = {}
    deferrals: dict[str, SlotDeferral] = {}
    columns: dict[str, ResolvedColumn] = {}

    if isinstance(node, GIQLDisjoin):
        target_ref = _resolve_target(node.this, tables)
        if target_ref is not None:
            slots["this"] = target_ref
            reference = node.args.get("reference")
            ref = _resolve_disjoin_reference(reference, target_ref, cte_names, tables)
            if ref is not None:
                slots["reference"] = ref
    elif isinstance(node, GIQLNearest):
        target_ref = _resolve_target(node.this, tables)
        if target_ref is not None:
            slots["this"] = target_ref
        # The reference slot resolves independently of the target — a literal
        # range parses without any registered table, and an implicit-outer
        # reference resolves against the enclosing scope. The interval and the
        # deferral are mutually exclusive (one is always ``None``).
        interval, deferral = _resolve_nearest_reference(node, tables)
        if interval is not None:
            slots["reference"] = interval
        if deferral is not None:
            deferrals["reference"] = deferral
    elif isinstance(node, GIQLDistance):
        columns = _resolve_distance_columns(node, tables)
    elif isinstance(node, (Intersects, Contains, Within, SpatialSetPredicate)):
        columns = _resolve_predicate_columns(node, tables)

    node.meta[META_KEY] = OperatorResolution(
        type(node).__name__, slots, deferrals, columns
    )


def _target_name(target: exp.Expression) -> str:
    """Extract the target table name from an operator's ``this`` slot.

    Mirrors the legacy generator's ``_resolve_target_table`` exactly.
    """
    if isinstance(target, exp.Table):
        return target.name
    if isinstance(target, exp.Column):
        return target.table if target.table else str(target.this)
    return str(target)


def _resolve_target(target: exp.Expression, tables: Tables) -> ResolvedRef | None:
    """Resolve a target (``this``) slot to a registered-table :class:`ResolvedRef`.

    Returns ``None`` when the target is not a registered table, in which case the
    generator raises its existing "not found in tables" error.
    """
    name = _target_name(target)
    table = tables.get(name)
    if table is None:
        return None
    return ResolvedRef(
        kind="registered_table",
        name=name,
        cols=(table.chrom_col, table.start_col, table.end_col),
        table=table,
        coverage_skippable=False,
    )


def _resolve_distance_columns(
    node: GIQLDistance, tables: Tables
) -> dict[str, ResolvedColumn]:
    """Resolve DISTANCE's two interval operands to :class:`ResolvedColumn`\\s.

    DISTANCE's ``this`` and ``expression`` operands are both column refs (per
    its :attr:`~giql.expressions.GIQLDistance.GIQL_SLOTS`). Each is resolved
    against the operand's table alias, mirroring the generator's historical
    ``_get_column_refs`` / ``_resolve_table`` path: the alias's physical
    genomic columns (from the registered :class:`~giql.table.Table` config, or
    the canonical defaults) qualified by that alias, plus the backing table
    config for the emitter's canonicalization wrapping.

    An operand the pass cannot resolve — a literal range, or an unqualified
    column — is omitted; the generator then raises its historical
    "Literal range as ... argument not yet supported" error.
    """
    alias_map, current_table = _enclosing_alias_map(node)
    columns: dict[str, ResolvedColumn] = {}
    for arg in ("this", "expression"):
        operand = node.args.get(arg)
        resolved = _resolve_column_operand(operand, tables, alias_map, current_table)
        if resolved is not None:
            columns[arg] = resolved
    return columns


def _resolve_predicate_columns(
    node: exp.Expression, tables: Tables
) -> dict[str, ResolvedColumn]:
    """Resolve a spatial predicate's column operands to :class:`ResolvedColumn`.

    Mirrors the generator's historical ``_generate_spatial_op`` operand handling
    exactly so the emitted SQL is byte-identical:

    * A literal-range predicate (``column INTERSECTS 'chr1:...'``) and every
      ``SpatialSetPredicate`` resolve only their left ``this`` operand, against
      the FROM-clause current table — matching ``_generate_range_predicate``,
      which passes ``self._current_table`` as the resolution context.
    * A column-to-column predicate (``a.interval CONTAINS b.interval``) resolves
      both ``this`` and ``expression`` with *no* current-table context, so each
      operand resolves through the alias map — matching ``_generate_column_join``,
      which passes ``None`` as the context for both sides.

    Most column-to-column ``INTERSECTS`` joins never reach here: the
    :class:`~giql.transformer.IntersectsBinnedJoinTransformer` rewrites them into
    binned equi-joins before this pass runs, deleting the ``Intersects`` node.
    Column-to-column ``CONTAINS`` / ``WITHIN`` (which that transformer leaves
    untouched) and non-join ``INTERSECTS`` shapes still reach the emitter, so the
    column-join branch is exercised.

    Reuses the shared ``_enclosing_alias_map`` (FROM/JOIN alias derivation) and
    ``_physical_cols`` helpers; the predicate-specific bit lives in
    :func:`_resolve_predicate_column`, whose current-table-vs-alias-only
    precedence differs from DISTANCE's ``_resolve_column_operand``.
    """
    alias_map, current_table = _enclosing_alias_map(node)
    this_node = node.this

    if isinstance(node, SpatialSetPredicate):
        if not isinstance(this_node, exp.Column):
            return {}
        return {
            "this": _resolve_predicate_column(
                this_node, current_table, alias_map, current_table, tables
            )
        }

    right = node.args.get("expression")
    if isinstance(right, exp.Column) and right.table:
        # Column-to-column: resolve both operands with no current-table context,
        # so each resolves through the alias map.
        columns: dict[str, ResolvedColumn] = {}
        if isinstance(this_node, exp.Column):
            columns["this"] = _resolve_predicate_column(
                this_node, None, alias_map, current_table, tables
            )
        columns["expression"] = _resolve_predicate_column(
            right, None, alias_map, current_table, tables
        )
        return columns

    # Literal-range predicate: resolve only the left operand, anchored to the
    # FROM-clause current table.
    if not isinstance(this_node, exp.Column):
        return {}
    return {
        "this": _resolve_predicate_column(
            this_node, current_table, alias_map, current_table, tables
        )
    }


def _resolve_predicate_column(
    column: exp.Column,
    ctx_table: str | None,
    alias_map: dict[str, str],
    current_table: str | None,
    tables: Tables,
) -> ResolvedColumn:
    """Resolve one spatial-predicate column operand to a :class:`ResolvedColumn`.

    Replicates the generator's ``_get_column_refs`` / ``_resolve_table`` (and the
    ``_resolve_table_name`` precedence underneath them) exactly:

    * the operand's table qualifier is kept verbatim for output formatting;
    * the *config* table name is resolved by precedence — an explicit
      ``ctx_table`` (the caller-supplied current table) wins; otherwise a
      qualified operand resolves through the alias map (falling back to
      ``current_table``), and an unqualified operand resolves to no table;
    * physical column names come from the resolved :class:`~giql.table.Table`
      config via the shared ``_physical_cols`` helper, or the canonical defaults
      when no table backs the operand.

    This is the minimal predicate-specific layer DISTANCE's
    ``_resolve_column_operand`` cannot serve: it formats unqualified operands
    (bare ``interval``) and lets the literal-range path anchor to ``ctx_table``
    directly, where ``_resolve_column_operand`` requires a qualifier and always
    routes through the alias map.
    """
    alias = column.table or ""
    qualified = bool(alias)

    if ctx_table:
        table_name: str | None = ctx_table
    elif qualified:
        table_name = alias_map.get(alias, current_table)
    else:
        table_name = None

    table = tables.get(table_name) if table_name else None
    chrom_col, start_col, end_col, strand_col = _physical_cols(table)

    if qualified:
        return ResolvedColumn(
            chrom=f'{alias}."{chrom_col}"',
            start=f'{alias}."{start_col}"',
            end=f'{alias}."{end_col}"',
            strand=f'{alias}."{strand_col}"',
            table=table,
        )
    return ResolvedColumn(
        chrom=f'"{chrom_col}"',
        start=f'"{start_col}"',
        end=f'"{end_col}"',
        strand=f'"{strand_col}"',
        table=table,
    )


def _resolve_column_operand(
    operand: exp.Expression | None,
    tables: Tables,
    alias_map: dict[str, str],
    current_table: str | None,
) -> ResolvedColumn | None:
    """Resolve a single column operand, or ``None`` if it is not a column ref.

    Mirrors the generator's ``_get_column_refs`` / ``_resolve_table`` exactly:
    the operand's alias is resolved to an underlying table name via the alias
    map (with the current FROM table as fallback), the physical column names
    come from that table's config (or the canonical defaults), and every column
    is qualified by the operand's alias.

    Returns ``None`` for an operand that is not a qualified column (a literal
    range or an unaliased column), deferring the diagnostic to the generator.
    """
    if not isinstance(operand, exp.Column):
        return None
    alias = operand.table
    if not alias:
        # An unqualified column has no alias to qualify by; the generator
        # treats it as a literal range and raises its existing error.
        return None

    table_name = alias_map.get(alias, current_table)
    table = tables.get(table_name) if table_name else None
    chrom_col, start_col, end_col, strand_col = _physical_cols(table)
    return ResolvedColumn(
        chrom=f'{alias}."{chrom_col}"',
        start=f'{alias}."{start_col}"',
        end=f'{alias}."{end_col}"',
        strand=f'{alias}."{strand_col}"',
        table=table,
    )


def _physical_cols(table: Table | None) -> tuple[str, str, str, str]:
    """Return the ``(chrom, start, end, strand)`` physical column names.

    Mirrors ``_get_column_refs``: the canonical defaults unless the registered
    :class:`~giql.table.Table` overrides them, with the strand column falling
    back to the default when the table declares none.
    """
    if table is None:
        return (
            DEFAULT_CHROM_COL,
            DEFAULT_START_COL,
            DEFAULT_END_COL,
            DEFAULT_STRAND_COL,
        )
    return (
        table.chrom_col,
        table.start_col,
        table.end_col,
        table.strand_col or DEFAULT_STRAND_COL,
    )


def _resolve_disjoin_reference(
    reference: exp.Expression | None,
    target_ref: ResolvedRef,
    cte_names: frozenset[str],
    tables: Tables,
) -> ResolvedRef | None:
    """Resolve a DISJOIN ``reference`` slot.

    Mirrors the generator's historical ``_resolve_disjoin_reference`` exactly
    (ported here in epic #114, step 2), with one deliberate difference: CTE
    visibility comes from the operator's scope (``scope.cte_sources``) rather
    than the hand-rolled ``_enclosing_cte_names`` ancestor walk, as epic #114
    specifies. The two agree on every supported shape.

    Returns ``None`` for the unresolvable shapes (reserved prefix, unknown
    name, unsupported node type); the generator raises its existing error.
    """
    # No reference: default to the target set — a self-reference whose coverage
    # EXISTS is provably always true.
    if reference is None:
        return replace(target_ref, coverage_skippable=True)

    # Subquery reference: a distinct relation assumed to expose canonical
    # default columns. Proving equivalence with the target is out of scope, so
    # it is never a self-reference.
    if isinstance(reference, (exp.Subquery, exp.Select, exp.Union)):
        return ResolvedRef(
            kind="subquery",
            name=None,
            cols=_DEFAULT_COLS,
            table=None,
            coverage_skippable=False,
        )

    # Anything that is not a bare table/CTE name is unsupported.
    if not isinstance(reference, (exp.Table, exp.Column, exp.Identifier)):
        return None

    ref_name = reference.name

    # The __giql_dj_ prefix names the operator's internal CTEs.
    if ref_name.startswith(DJ_PREFIX):
        return None

    # A CTE from an enclosing WITH shadows a registered table of the same name
    # and is assumed to expose canonical default columns. A CTE may hold rows
    # distinct from a same-named registered table, so it is never a
    # self-reference.
    if ref_name in cte_names:
        return ResolvedRef(
            kind="cte",
            name=ref_name,
            cols=_DEFAULT_COLS,
            table=None,
            coverage_skippable=False,
        )

    ref_table = tables.get(ref_name)
    if ref_table is not None:
        is_self = ref_name == target_ref.name and ref_table is target_ref.table
        return ResolvedRef(
            kind="registered_table",
            name=ref_name,
            cols=(ref_table.chrom_col, ref_table.start_col, ref_table.end_col),
            table=ref_table,
            coverage_skippable=is_self,
        )

    return None


def _resolve_nearest_reference(
    node: GIQLNearest, tables: Tables
) -> tuple[ResolvedInterval | None, SlotDeferral | None]:
    """Resolve a NEAREST ``reference`` slot to a :class:`ResolvedInterval`.

    Ported here (epic #114, step 3) from the generator's historical
    ``_resolve_nearest_reference`` / ``_find_outer_table_in_lateral_join``. The
    accepted shapes are declarative on :class:`giql.expressions.GIQLNearest`:

    * **literal range** — a quoted genomic-range string, parsed to canonical
      0-based half-open literal endpoints.
    * **column reference** — an explicit ``alias.col`` interval, resolved
      against the enclosing SELECT's ``FROM`` / ``JOIN`` aliases exactly as the
      generator's ``select_sql`` builds ``_alias_to_table`` / ``_current_table``.
    * **implicit LATERAL outer** — no reference given, resolved to the outer
      relation of the enclosing ``CROSS JOIN LATERAL`` via the same ancestor
      walk the generator used.

    Returns ``(interval, None)`` on success or ``(None, deferral)`` when the
    slot cannot be resolved; the generator then re-raises its historical error
    (a ``None`` deferral means the emitter reclassifies the failure unaided, as
    for a literal range that fails to parse).
    """
    reference = node.args.get("reference")

    if reference is None:
        return _resolve_implicit_outer(node, tables)

    # Literal genomic-range string (e.g. 'chr1:1000-2000'). sqlglot parses it as
    # a string Literal; the generator detected it by the leading quote on the
    # rendered SQL.
    if isinstance(reference, exp.Literal) and reference.is_string:
        return _resolve_literal_range(reference.name)

    # Explicit column reference (e.g. peaks.interval).
    if isinstance(reference, exp.Column):
        return _resolve_column_interval(node, reference.table or None, tables)

    # Fallback for any other shape: mirror the generator's string-level
    # classification on the rendered reference (quote prefix => literal,
    # otherwise a possibly-qualified column).
    rendered = reference.sql()
    if rendered.startswith("'") or rendered.startswith('"'):
        return _resolve_literal_range(rendered.strip("'\""))
    qualifier = rendered.rsplit(".", 1)[0] if "." in rendered else None
    return _resolve_column_interval(node, qualifier, tables)


def _resolve_literal_range(
    range_str: str,
) -> tuple[ResolvedInterval | None, SlotDeferral | None]:
    """Resolve a literal genomic-range reference to canonical literal endpoints.

    Defers (``(None, None)``) when the range fails to parse; the generator
    re-parses and raises the historical "Could not parse reference genomic
    range" diagnostic with the original exception text.
    """
    try:
        parsed = RangeParser.parse(range_str).to_zero_based_half_open()
    except Exception:
        return None, None
    strand = f"'{parsed.strand}'" if parsed.strand else None
    return (
        ResolvedInterval(
            kind="literal_range",
            chrom=f"'{parsed.chromosome}'",
            start=str(parsed.start),
            end=str(parsed.end),
            strand=strand,
            table=None,
        ),
        None,
    )


def _resolve_column_interval(
    node: GIQLNearest, qualifier: str | None, tables: Tables
) -> tuple[ResolvedInterval, None]:
    """Resolve an explicit column reference to a column :class:`ResolvedInterval`.

    Mirrors the generator's ``_get_column_refs`` / ``_resolve_table`` for an
    explicit reference: the *qualifier* (the column's relation alias) is kept
    verbatim for output, while the backing :class:`~giql.table.Table` — and thus
    the physical column names and coordinate system — is looked up by resolving
    that alias through the enclosing SELECT's alias map. An unresolved qualifier
    falls back to default column names and no table (so the emitter applies no
    canonicalization), exactly as the generator did. The strand column is always
    emitted for an explicit reference, matching ``include_strand=True``.
    """
    table = _lookup_aliased_table(node, qualifier, tables)
    chrom_col = table.chrom_col if table else DEFAULT_CHROM_COL
    start_col = table.start_col if table else DEFAULT_START_COL
    end_col = table.end_col if table else DEFAULT_END_COL
    strand_col = table.strand_col if (table and table.strand_col) else DEFAULT_STRAND_COL
    prefix = f"{qualifier}." if qualifier else ""
    return (
        ResolvedInterval(
            kind="column",
            chrom=f'{prefix}"{chrom_col}"',
            start=f'{prefix}"{start_col}"',
            end=f'{prefix}"{end_col}"',
            strand=f'{prefix}"{strand_col}"',
            table=table,
        ),
        None,
    )


def _resolve_implicit_outer(
    node: GIQLNearest, tables: Tables
) -> tuple[ResolvedInterval | None, SlotDeferral | None]:
    """Resolve an implicit-outer reference from the enclosing LATERAL join.

    Mirrors the generator's ``_find_outer_table_in_lateral_join`` followed by
    its outer-table column resolution: the outer relation's label (alias or
    name) qualifies the emitted columns, while the backing table is resolved by
    mapping that label through the enclosing SELECT's alias map. Unlike an
    explicit reference, the strand column is emitted only when the outer table
    declares one — preserving the generator's divergent strand handling between
    the two paths.

    Defers with a :class:`SlotDeferral` when no outer relation is found
    (``implicit_outer_missing``) or the outer relation is unregistered
    (``implicit_outer_unregistered``, carrying the offending label); the
    generator raises the matching historical diagnostic.
    """
    outer_label = _find_outer_table(node)
    if outer_label is None:
        return None, SlotDeferral("implicit_outer_missing")

    alias_map, _current = _enclosing_alias_map(node)
    actual_name = alias_map.get(outer_label, outer_label)
    table = tables.get(actual_name)
    if table is None:
        return None, SlotDeferral("implicit_outer_unregistered", outer_label)

    strand = f'{outer_label}."{table.strand_col}"' if table.strand_col else None
    return (
        ResolvedInterval(
            kind="implicit_outer",
            chrom=f'{outer_label}."{table.chrom_col}"',
            start=f'{outer_label}."{table.start_col}"',
            end=f'{outer_label}."{table.end_col}"',
            strand=strand,
            table=table,
        ),
        None,
    )


def _lookup_aliased_table(
    node: GIQLNearest, qualifier: str | None, tables: Tables
) -> Table | None:
    """Resolve the :class:`~giql.table.Table` backing a qualified column.

    Replicates the generator's ``_resolve_table`` / ``_resolve_table_name``: a
    dotted reference's *qualifier* is mapped to a registered-table name through
    the enclosing SELECT's alias map (with ``_current_table`` as the FROM-clause
    fallback); an unqualified reference resolves to no table.
    """
    if not qualifier:
        return None
    alias_map, current_table = _enclosing_alias_map(node)
    table_name = alias_map.get(qualifier, current_table)
    return tables.get(table_name) if table_name else None


def _enclosing_alias_map(node: exp.Expression) -> tuple[dict[str, str], str | None]:
    """Build the alias->table map of *node*'s enclosing SELECT.

    Mirrors the legacy generator's ``select_sql`` exactly: the FROM-clause table and
    every JOIN-clause table contribute an ``(alias or name) -> name`` entry, and
    the FROM-clause table name is the ``_current_table`` fallback. Only direct
    ``exp.Table`` operands participate (a LATERAL/derived join contributes
    nothing), matching the generator. CTE and derived-table aliasing is out of
    scope here because NEAREST column references resolve against physical
    relations, exactly as the generator's hand-rolled map did.
    """
    select = node.parent_select
    alias_to_table: dict[str, str] = {}
    current_table: str | None = None
    if isinstance(select, exp.Select):
        from_ = select.args.get("from_")
        if from_ is not None and isinstance(from_.this, exp.Table):
            current_table = from_.this.name
            if from_.this.alias:
                alias_to_table[from_.this.alias] = current_table
            else:
                alias_to_table[current_table] = current_table
        for join in select.args.get("joins") or []:
            if isinstance(join.this, exp.Table):
                name = join.this.name
                if join.this.alias:
                    alias_to_table[join.this.alias] = name
                else:
                    alias_to_table[name] = name
    return alias_to_table, current_table


def _find_outer_table(node: exp.Expression) -> str | None:
    """Find the outer relation label of *node*'s enclosing LATERAL join.

    Walks up from the NEAREST node through the enclosing ``exp.Lateral`` to the
    ``exp.Join``, then reads the join-owning SELECT's FROM-clause relation,
    returning its alias or name. A byte-for-byte port of the generator's
    ``_find_outer_table_in_lateral_join``; returns ``None`` when no such outer
    relation exists (e.g. a standalone NEAREST).
    """
    current: exp.Expression | None = node
    while current is not None:
        parent = current.parent
        if parent is None:
            break
        if isinstance(parent, exp.Lateral):
            current = parent
            continue
        if isinstance(parent, exp.Join):
            select = parent.parent
            if isinstance(select, exp.Select):
                from_ = select.args.get("from_")
                if from_ is not None:
                    table_expr = from_.this
                    if isinstance(table_expr, exp.Table):
                        return table_expr.alias or table_expr.name
                    elif isinstance(table_expr, exp.Alias):
                        return table_expr.alias
            break
        current = parent
    return None


def validate_operator_refs(expression: exp.Expression) -> None:
    """Assert every GIQL operator carries well-formed resolution metadata.

    The closing validation boundary of pass 1, modeled on ``sqlglot``'s
    ``validate_qualify_columns`` and Spark's ``CheckAnalysis``: a pure check
    over the annotated tree. For each GIQL operator it asserts that the node
    carries an :class:`OperatorResolution`, and that every reference slot which
    *was* resolved carries a :class:`ResolvedRef` whose ``kind`` is permitted by
    the slot's declared :class:`~giql.expressions.SlotSpec`.

    Consistent with the step-1 zero-behavior-change constraint, an unresolved
    reference slot is **not** an error here — that diagnostic is deferred to the
    generator (and, in later epic steps, promoted to a GIQL-level message). This
    function raises only :class:`ResolutionError`, signalling an internal bug in
    the resolver rather than an invalid user query.

    Parameters
    ----------
    expression : exp.Expression
        The annotated AST to validate.

    Raises
    ------
    ResolutionError
        If an operator is missing its resolution metadata or a resolved slot is
        malformed.
    """
    for node in expression.walk():
        if not isinstance(node, _OPERATORS):
            continue

        resolution = node.meta.get(META_KEY)
        if not isinstance(resolution, OperatorResolution):
            raise ResolutionError(
                f"{type(node).__name__} is missing resolution metadata; "
                "the ResolveOperatorRefs pass did not annotate it."
            )

        specs: tuple[SlotSpec, ...] = getattr(node, "GIQL_SLOTS", ())
        for spec in specs:
            resolved = resolution.slots.get(spec.arg)
            if resolved is not None:
                if spec.is_ref_slot:
                    _validate_ref(resolved, spec, type(node).__name__)
                else:
                    _validate_interval(resolved, spec, type(node).__name__)
                continue
            # A column operand (DISTANCE) resolves through the separate columns
            # channel rather than the slots map; validate it there.
            column = resolution.columns.get(spec.arg)
            if column is not None:
                _validate_column(column, spec, type(node).__name__)
            # Otherwise deferred: an unresolved slot is handled by the generator
            # on its existing path (and may carry a SlotDeferral).


def _validate_interval(interval: object, spec: SlotSpec, operator: str) -> None:
    """Assert a single resolved interval is well-formed against its slot spec."""
    if not isinstance(interval, ResolvedInterval):
        raise ResolutionError(
            f"{operator} slot {spec.arg!r} carries {type(interval).__name__}, "
            "expected ResolvedInterval."
        )
    if interval.kind not in spec.accepts:
        raise ResolutionError(
            f"{operator} slot {spec.arg!r} resolved to kind {interval.kind!r}, "
            f"which is not accepted by the slot (accepts {sorted(spec.accepts)})."
        )
    if not all(
        isinstance(part, str) for part in (interval.chrom, interval.start, interval.end)
    ):
        raise ResolutionError(
            f"{operator} slot {spec.arg!r} has malformed interval endpoints; "
            "expected SQL fragment strings for chrom/start/end."
        )
    if interval.kind == "literal_range" and interval.table is not None:
        raise ResolutionError(
            f"{operator} slot {spec.arg!r} resolved to a literal_range but "
            "carries a Table config; literal ranges are already canonical."
        )


def _validate_column(column: object, spec: SlotSpec, operator: str) -> None:
    """Assert a single resolved column operand is well-formed against its slot.

    DISTANCE's two interval operands resolve to a :class:`ResolvedColumn` in the
    :attr:`OperatorResolution.columns` channel. A column operand is only ever
    attached for a slot whose declared shapes include ``"column"``; the endpoint
    fragments must be SQL strings ready to drop into the emitter's arithmetic.
    """
    if not isinstance(column, ResolvedColumn):
        raise ResolutionError(
            f"{operator} slot {spec.arg!r} carries {type(column).__name__}, "
            "expected ResolvedColumn."
        )
    if "column" not in spec.accepts:
        raise ResolutionError(
            f"{operator} slot {spec.arg!r} resolved to a column operand, which is "
            f"not accepted by the slot (accepts {sorted(spec.accepts)})."
        )
    if not all(
        isinstance(part, str) for part in (column.chrom, column.start, column.end)
    ):
        raise ResolutionError(
            f"{operator} slot {spec.arg!r} has malformed column endpoints; "
            "expected SQL fragment strings for chrom/start/end."
        )


def _validate_ref(ref: object, spec: SlotSpec, operator: str) -> None:
    """Assert a single resolved reference is well-formed against its slot spec."""
    if not isinstance(ref, ResolvedRef):
        raise ResolutionError(
            f"{operator} slot {spec.arg!r} carries {type(ref).__name__}, "
            "expected ResolvedRef."
        )
    if ref.kind not in spec.accepts:
        raise ResolutionError(
            f"{operator} slot {spec.arg!r} resolved to kind {ref.kind!r}, "
            f"which is not accepted by the slot (accepts {sorted(spec.accepts)})."
        )
    if not (
        isinstance(ref.cols, tuple)
        and len(ref.cols) == 3
        and all(isinstance(col, str) for col in ref.cols)
    ):
        raise ResolutionError(
            f"{operator} slot {spec.arg!r} has malformed cols {ref.cols!r}; "
            "expected a 3-tuple of column names."
        )
    if ref.kind == "registered_table" and ref.table is None:
        raise ResolutionError(
            f"{operator} slot {spec.arg!r} resolved to a registered table but "
            "carries no Table config."
        )
    if ref.kind in ("cte", "subquery") and ref.table is not None:
        raise ResolutionError(
            f"{operator} slot {spec.arg!r} resolved to a {ref.kind} but carries "
            "a Table config; CTE and subquery references are assumed canonical."
        )


def validate_projection_contracts(expression: exp.Expression) -> None:
    """Assert every CTE / subquery operator operand projects canonical columns.

    A user-facing GIQL diagnostic boundary (epic #114, step 9). When an operator
    reference slot resolves to an in-query CTE or a ``(SELECT ...)`` subquery
    (:class:`ResolvedRef` kind ``"cte"`` / ``"subquery"``), the resolver and the
    canonicalizer both *assume* it exposes the canonical ``chrom`` / ``start`` /
    ``end`` columns the operator will reference — an assumption that, when wrong,
    historically surfaced only as an opaque engine ``column not found`` error at
    execution time. This pass inspects the operand's projection at transpile time
    and raises a clear :class:`ValueError` naming the operator, the offending
    relation, and the missing column(s) instead.

    Unlike :func:`validate_operator_refs` (an internal invariant check that raises
    :class:`ResolutionError`), this is a *user-facing* check that raises
    :class:`ValueError` so it surfaces verbatim through ``transpile()``'s
    "Resolution error" boundary.

    False-positive guard
    --------------------
    A projection's output column names are only statically known when **every**
    projection expression is an intentional named column — a bare ``exp.Column``
    (whose output name is the column) or an aliased expression ``exp.Alias``
    (whose output name is the explicit alias). If any projection is a ``*`` star
    (whose columns expand from an unknown source schema) or an unnamed expression
    (a bare literal like ``SELECT 1`` or a computed ``chrom + 1`` with no alias),
    the relation's schema is not statically determinable and the operand is
    **passed** without a diagnostic — a documented limitation that avoids
    false-positives on star and placeholder projections. The check fires only
    when the schema is fully known *and* a canonical column is provably absent.

    Parameters
    ----------
    expression : exp.Expression
        The pass-1-annotated AST to validate.

    Raises
    ------
    ValueError
        If a CTE or subquery operand's statically-known projection omits one or
        more of the canonical ``chrom`` / ``start`` / ``end`` columns the
        operator references.
    """
    cte_bodies = _cte_body_index(expression)
    for node in expression.walk():
        if not isinstance(node, _OPERATORS):
            continue
        resolution = node.meta.get(META_KEY)
        if not isinstance(resolution, OperatorResolution):
            continue
        operator = resolution.operator
        for arg, ref in resolution.slots.items():
            if not isinstance(ref, ResolvedRef):
                continue
            if ref.kind not in ("cte", "subquery"):
                continue
            select = _projection_select(node, arg, ref, cte_bodies)
            if select is None:
                continue
            _check_projection(node, arg, ref, select, operator)


def _cte_body_index(expression: exp.Expression) -> dict[str, exp.Expression]:
    """Map each in-query CTE alias to the query body it defines.

    Built once per validation so a CTE-kind reference can locate its ``WITH``
    definition and read its output projection. The body is the CTE's ``this`` —
    an ``exp.Select`` or set operation (``exp.Union`` / ``exp.Intersect`` /
    ``exp.Except``). When two CTEs share an alias (inner-WITH redeclaration), the
    last one wins; the resolver's scope-based reference resolution already
    selected which one the operator sees, and the projection contract is
    identical across redeclarations for the columns we check.
    """
    index: dict[str, exp.Expression] = {}
    for cte in expression.find_all(exp.CTE):
        body = cte.this
        if body is not None:
            index[cte.alias] = body
    return index


def _projection_select(
    node: exp.Expression,
    arg: str,
    ref: ResolvedRef,
    cte_bodies: dict[str, exp.Expression],
) -> exp.Expression | None:
    """Return the query body whose projection backs reference slot *arg*.

    For a CTE reference the body is looked up by name in *cte_bodies*; for a
    subquery reference it is the operand AST node itself (an ``exp.Subquery``
    unwrapped to its inner query, or a bare ``exp.Select`` / set operation).
    Returns ``None`` when no projection-bearing body can be located — e.g. a CTE
    whose definition is not in the index (out of scope) — so the caller skips it.
    """
    if ref.kind == "cte":
        return cte_bodies.get(ref.name) if ref.name else None
    operand = node.args.get(arg)
    if isinstance(operand, exp.Subquery):
        operand = operand.this
    if isinstance(operand, (exp.Select, exp.SetOperation)):
        return operand
    return None


def _static_output_columns(select: exp.Expression) -> frozenset[str] | None:
    """Return a body's statically-known output column names, or ``None``.

    ``None`` signals an *unresolvable* projection — one whose output schema
    cannot be determined statically (a ``*`` star or any unnamed expression) —
    which the contract check treats as a pass (the false-positive guard). A
    non-``None`` result is the complete set of intentional output column names,
    safe to test for the canonical columns.

    A set operation (``UNION`` / ``INTERSECT`` / ``EXCEPT``) defers to its left
    branch, whose projection determines the result's column names in SQL.
    """
    if isinstance(select, exp.SetOperation):
        return _static_output_columns(select.this)
    if not isinstance(select, exp.Select):
        return None

    names: set[str] = set()
    for projection in select.selects:
        if not isinstance(projection, (exp.Column, exp.Alias)):
            # A star or an unnamed expression (bare literal, computed column):
            # the output schema is not statically determinable, so the whole
            # projection is unresolvable and the operand passes.
            return None
        if isinstance(projection, exp.Column) and isinstance(projection.this, exp.Star):
            # A qualified star (``a.*``) also expands from an unknown schema.
            return None
        names.add(projection.output_name)
    return frozenset(names)


def _check_projection(
    node: exp.Expression,
    arg: str,
    ref: ResolvedRef,
    select: exp.Expression,
    operator: str,
) -> None:
    """Raise a GIQL diagnostic if *select* omits a canonical column.

    Compares the operand's statically-known output columns against the canonical
    column names the operator references (``ref.cols``). A star or otherwise
    unresolvable projection (``_static_output_columns`` returns ``None``) is a
    pass. When the schema is known and any required column is absent, raises a
    :class:`ValueError` naming the operator, the offending relation, the missing
    column(s), and a source position when one is available.
    """
    columns = _static_output_columns(select)
    if columns is None:
        return
    missing = [col for col in ref.cols if col not in columns]
    if not missing:
        return

    relation = f"CTE {ref.name!r}" if ref.kind == "cte" else "subquery"
    missing_str = ", ".join(repr(col) for col in missing)
    have_str = ", ".join(repr(col) for col in sorted(columns)) or "no columns"
    position = _source_position(node.args.get(arg))
    raise ValueError(
        f"{operator} {arg} {relation} does not project the canonical "
        f"genomic column(s) {missing_str} that the operator references; "
        f"it projects {have_str}. A CTE or subquery operand must project "
        f"canonical 0-based half-open 'chrom'/'start'/'end' columns."
        f"{position}"
    )


def _source_position(operand: exp.Expression | None) -> str:
    """Return a `` (line L, col C)`` suffix for *operand*, or ``""``.

    Mirrors ``sqlglot``'s ``validate_qualify_columns`` precedent of attaching the
    parser-populated source position to a validation diagnostic. The position is
    best-effort: when the parser did not record one (``meta`` carries no
    ``line``), an empty string is returned and the textual diagnostic stands on
    its own.
    """
    if operand is None:
        return ""
    line = operand.meta.get("line")
    col = operand.meta.get("col")
    if line is None:
        return ""
    if col is None:
        return f" (line {line})"
    return f" (line {line}, col {col})"
