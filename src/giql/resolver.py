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

Scope note (epic #114, steps 1-2)
---------------------------------
The pass is behavior-preserving. DISJOIN's emitter
(``BaseGIQLGenerator.giqldisjoin_sql``) consumes the attached metadata (step 2);
the remaining operators still use the generator's legacy resolver paths and
ignore everything attached here until their port issues land. The resolution
semantics computed for the table-shaped reference slots mirror the generator's
historical ``_resolve_target_table`` / ``_resolve_disjoin_reference`` /
``_enclosing_cte_names`` behavior exactly (the latter two now live only here).

Two consequences of the zero-behavior-change constraint shape the
implementation:

* The pass never raises a *user-facing* diagnostic. When a slot cannot be
  resolved (unregistered target, unknown reference name, reserved ``__giql_dj_``
  prefix, unsupported reference shape) the pass simply leaves that slot
  unresolved and the generator raises its existing error exactly as before.
  Promoting these into GIQL-level diagnostics is a later epic step.
* Only *reference slots* — slots whose accepted shapes are the table
  trichotomy (DISJOIN ``this``/``reference`` and NEAREST ``this``) — are
  resolved to a :class:`ResolvedRef` here. The column / literal /
  implicit-outer slots of NEAREST, DISTANCE, and the spatial predicates are
  declared on the expression classes but their resolution metadata type is
  designed by the per-operator port issues (#118, #119, #120).
"""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import replace
from typing import Literal

from sqlglot import exp
from sqlglot.optimizer.scope import Scope
from sqlglot.optimizer.scope import traverse_scope

from giql.constants import DEFAULT_CHROM_COL
from giql.constants import DEFAULT_END_COL
from giql.constants import DEFAULT_START_COL
from giql.expressions import Contains
from giql.expressions import GIQLDisjoin
from giql.expressions import GIQLDistance
from giql.expressions import GIQLNearest
from giql.expressions import Intersects
from giql.expressions import SlotSpec
from giql.expressions import SpatialSetPredicate
from giql.expressions import Within
from giql.table import Table
from giql.table import Tables

__all__ = [
    "META_KEY",
    "RefKind",
    "ResolvedRef",
    "OperatorResolution",
    "ResolutionError",
    "resolve_operator_refs",
    "validate_operator_refs",
]

#: The single namespaced key under which resolution metadata is stored on a
#: node's ``.meta`` mapping. All GIQL metadata lives under this one key so the
#: channel stays a single, greppable namespace.
META_KEY = "giql"

#: The kinds a resolved reference slot can take — the registered-table / CTE /
#: subquery trichotomy that ``sqlglot``'s ``scope.sources`` distinguishes.
RefKind = Literal["registered_table", "cte", "subquery"]

#: Canonical default genomic column names. A CTE or subquery reference is
#: assumed (in the existing generator and here) to expose canonical
#: ``chrom`` / ``start`` / ``end`` columns; validating that contract is a later
#: epic step.
_DEFAULT_COLS: tuple[str, str, str] = (
    DEFAULT_CHROM_COL,
    DEFAULT_START_COL,
    DEFAULT_END_COL,
)

#: The GIQL operator expression classes the pass inspects.
_OPERATORS: tuple[type[exp.Expression], ...] = (
    GIQLDisjoin,
    GIQLNearest,
    GIQLDistance,
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
class OperatorResolution:
    """Resolution metadata attached to a single GIQL operator node.

    Attributes
    ----------
    operator : str
        The operator expression class name (e.g. ``"GIQLDisjoin"``).
    slots : dict[str, ResolvedRef]
        Mapping from a slot's :attr:`~giql.expressions.SlotSpec.arg` key to its
        resolved :class:`ResolvedRef`. Only successfully resolved reference
        slots appear; a slot left out was either not a reference slot or could
        not be resolved (and the generator will raise its existing error).
    """

    operator: str
    slots: dict[str, ResolvedRef]

    def slot(self, arg: str) -> ResolvedRef | None:
        """Return the resolved reference for slot *arg*, or ``None``."""
        return self.slots.get(arg)


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
    # behavior-preserving: a missed CTE reference simply stays unresolved and
    # the generator handles it on its existing path.
    for node in expression.walk():
        if isinstance(node, _OPERATORS) and id(node) not in seen:
            seen.add(id(node))
            _resolve_operator(node, tables, frozenset())

    validate_operator_refs(expression)
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
    """Resolve *node*'s reference slots and attach an :class:`OperatorResolution`."""
    slots: dict[str, ResolvedRef] = {}

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
    # DISTANCE and the spatial predicates declare only column / literal slots,
    # whose resolution metadata is designed by their port issues; the pass
    # attaches an (empty-slot) resolution so every operator carries metadata.

    node.meta[META_KEY] = OperatorResolution(type(node).__name__, slots)


def _target_name(target: exp.Expression) -> str:
    """Extract the target table name from an operator's ``this`` slot.

    Mirrors ``BaseGIQLGenerator._resolve_target_table`` exactly.
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
    if ref_name.startswith("__giql_dj_"):
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
            if not spec.is_ref_slot:
                continue
            ref = resolution.slots.get(spec.arg)
            if ref is None:
                # Deferred: unresolved reference slots are handled by the
                # generator on its existing path in step 1.
                continue
            _validate_ref(ref, spec, type(node).__name__)


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
