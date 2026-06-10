"""Pass 2 of the GIQL normalization pipeline: ``CanonicalizeCoordinates``.

This module implements the second normalization pass described in epic #114. It
runs after :func:`giql.resolver.resolve_operator_refs` (pass 1) and consumes the
:class:`~giql.resolver.ResolvedRef` metadata that pass attached to every GIQL
operator slot. For every interval-bearing relation reachable from an operator
slot whose ``ResolvedRef`` carries a :class:`~giql.table.Table` config declaring
a *non-canonical* encoding (``coordinate_system`` / ``interval_type`` other than
``0based`` / ``half_open``), the pass synthesizes a hidden ``__giql_canon_*``
WITH-CTE that projects the relation to canonical 0-based half-open coordinates
(replicating the :mod:`giql.canonical` arithmetic as SQL in the CTE projection)
and rewrites the operator slot — both the AST node and its ``ResolvedRef``
metadata — to point at the canonical CTE.

The synthesis follows ``sqlglot``'s ``eliminate_subqueries`` mechanics
(:mod:`sqlglot.optimizer.eliminate_subqueries`):

* **Collision-safe aliases.** Candidate names come from
  :func:`sqlglot.helper.name_sequence` with the ``__giql_canon_`` prefix, checked
  against a taken-set collected across every scope's table sources and CTE
  aliases; :func:`sqlglot.helper.find_new_name` resolves any residual collision
  against a user identifier.
* **Wrapper-body deduplication.** Two operator slots referencing the same source
  under the same encoding share a single canonical CTE — one source is wrapped at
  most once per query, keyed on the rendered wrapper body.
* **DAG-ordered insertion.** Canonical CTEs depend only on base (registered)
  tables, never on other CTEs, so they are prepended to the outermost ``WITH``;
  any existing CTE or main-body operator that now references a canonical CTE
  therefore sees it defined first.

Identity-encoded relations pass through **unwrapped**: a relation already in
0-based half-open form, or one with no ``Table`` config (a CTE or subquery
reference, assumed canonical), is left untouched — the readability-vs-volume
tradeoff the epic calls out (only synthesize a wrapper when canonicalization
actually changes columns).

Gating (epic #114, step 6)
--------------------------
The pass is gated per operator by a ``GIQL_CANONICALIZE`` class attribute on the
operator's expression class. An operator opts in by setting
``GIQL_CANONICALIZE = True``; absent or ``False`` (the default for every operator
as of this issue) the pass ignores it entirely. The operator port issues — #122
(DISJOIN) and #123 (NEAREST / DISTANCE / predicates) — flip these flags as each
operator's emitter is moved off in-emitter canonicalization
(:mod:`giql.canonical`) and onto this pass's output. **With every flag off the
pass is a strict no-op and the emitted SQL is byte-identical**, so the existing
suite is the migration oracle.

De-canonicalization hook
-------------------------
The outermost ``SELECT`` projection receives a de-canonicalization rewrite for
any output column that a migrated operator emitted in canonical form but that
must land in the user's preferred encoding. With no operator migrated in this
issue that rewrite has nothing to act on; :func:`_decanonicalize_outputs` is the
designed-but-inert hook the port issues will fill in.
"""

from __future__ import annotations

from dataclasses import replace

from sqlglot import exp
from sqlglot.helper import find_new_name
from sqlglot.helper import name_sequence
from sqlglot.optimizer.scope import traverse_scope

from giql.expressions import Contains
from giql.expressions import GIQLDisjoin
from giql.expressions import GIQLDistance
from giql.expressions import GIQLNearest
from giql.expressions import Intersects
from giql.expressions import SlotSpec
from giql.expressions import SpatialSetPredicate
from giql.expressions import Within
from giql.resolver import META_KEY
from giql.resolver import OperatorResolution
from giql.resolver import ResolvedRef
from giql.table import Table

__all__ = [
    "CANON_PREFIX",
    "canonicalize_coordinates",
]

#: The alias prefix for every synthesized canonical wrapper CTE. Mirrors the
#: ``__giql_dj_`` reserved prefix the DISJOIN emitter already uses; the leading
#: double underscore keeps the namespace clear of user identifiers.
CANON_PREFIX = "__giql_canon_"

#: The GIQL operator expression classes the pass inspects — the same set pass 1
#: annotates. Membership alone does not opt an operator in; the per-class
#: ``GIQL_CANONICALIZE`` flag does (see module docstring).
_OPERATORS: tuple[type[exp.Expression], ...] = (
    GIQLDisjoin,
    GIQLNearest,
    GIQLDistance,
    Intersects,
    Contains,
    Within,
    SpatialSetPredicate,
)


def canonicalize_coordinates(expression: exp.Expression) -> exp.Expression:
    """Synthesize canonical wrapper CTEs for non-canonical operator operands.

    Walks *expression* for opted-in GIQL operators (those whose expression class
    sets ``GIQL_CANONICALIZE = True``), and for each reference slot that resolved
    to a registered table with a non-canonical encoding, synthesizes a
    ``__giql_canon_*`` WITH-CTE projecting that table to canonical 0-based
    half-open coordinates and rewrites the slot (AST node + ``ResolvedRef``
    metadata) to point at the canonical CTE.

    The pass mutates and returns *expression* in place. **When no operator opts
    in — the state as of issue #121 — it is a strict no-op: no node is touched
    and the emitted SQL is byte-identical.**

    Parameters
    ----------
    expression : exp.Expression
        The pass-1-annotated AST.

    Returns
    -------
    exp.Expression
        The same *expression*, with canonical wrapper CTEs inserted and migrated
        operator slots rewritten (none, while every flag is off).
    """
    targets = _collect_targets(expression)
    if not targets:
        # No opted-in operator carries a non-canonical operand: strict no-op.
        return expression

    taken = _collect_taken_names(expression)
    next_name = name_sequence(CANON_PREFIX)
    body_to_name: dict[str, str] = {}
    new_ctes: list[exp.CTE] = []

    for node, arg, ref in targets:
        body = _canonical_projection(ref)
        body_sql = body.sql()
        name = body_to_name.get(body_sql)
        if name is None:
            # First time this exact wrapper body is needed: mint a fresh,
            # collision-safe alias and emit one CTE for it.
            name = _fresh_name(next_name, taken)
            taken.add(name)
            body_to_name[body_sql] = name
            new_ctes.append(_make_cte(name, body))
        # Subsequent slots with an identical body reuse the same CTE (dedup).
        _rewrite_slot(node, arg, ref, name)

    _insert_ctes(expression, new_ctes)
    _decanonicalize_outputs(expression, targets)
    return expression


def _collect_targets(
    expression: exp.Expression,
) -> list[tuple[exp.Expression, str, ResolvedRef]]:
    """Return ``(node, slot_arg, ref)`` triples needing canonicalization.

    A slot qualifies when its operator opts in (``GIQL_CANONICALIZE``), the slot
    is a reference slot that resolved to a registered table, and that table's
    declared encoding is non-canonical. CTE / subquery references (assumed
    canonical) and identity-encoded tables are skipped — they pass through
    unwrapped.
    """
    targets: list[tuple[exp.Expression, str, ResolvedRef]] = []
    for node in expression.walk():
        if not isinstance(node, _OPERATORS):
            continue
        if not getattr(node, "GIQL_CANONICALIZE", False):
            continue
        resolution = node.meta.get(META_KEY)
        if not isinstance(resolution, OperatorResolution):
            continue
        specs: tuple[SlotSpec, ...] = getattr(node, "GIQL_SLOTS", ())
        for spec in specs:
            if not spec.is_ref_slot:
                continue
            ref = resolution.slots.get(spec.arg)
            if ref is None:
                continue
            # Only registered tables carry a declared encoding; CTE / subquery
            # references are assumed canonical and pass through unwrapped.
            if ref.kind != "registered_table":
                continue
            if _is_canonical(ref.table):
                continue
            targets.append((node, spec.arg, ref))
    return targets


def _is_canonical(table: Table | None) -> bool:
    """Whether *table*'s declared encoding is already canonical 0-based half-open.

    ``None`` (an unregistered relation) is treated as canonical, matching the
    :mod:`giql.canonical` convention.
    """
    if table is None:
        return True
    return table.coordinate_system == "0based" and table.interval_type == "half_open"


def _collect_taken_names(expression: exp.Expression) -> set[str]:
    """Collect every relation name already in use across all scopes.

    Mirrors ``eliminate_subqueries``'s taken-set: every registered-table source
    name plus every CTE alias, so a freshly minted ``__giql_canon_*`` alias can
    be proven collision-free. A scope-traversal failure falls back to a plain
    tree walk so the set is never under-populated.
    """
    taken: set[str] = set()
    try:
        scopes = traverse_scope(expression)
    except Exception:
        scopes = []
    for scope in scopes:
        for source in scope.sources.values():
            if isinstance(source, exp.Table):
                taken.add(source.name)
        taken.update(scope.cte_sources)
    # Belt-and-suspenders: a raw walk catches any name the scope machinery did
    # not surface (e.g. when traverse_scope bailed on an exotic construct).
    for table in expression.find_all(exp.Table):
        taken.add(table.name)
    for cte in expression.find_all(exp.CTE):
        taken.add(cte.alias)
    return taken


def _fresh_name(next_name, taken: set[str]) -> str:
    """Return a ``__giql_canon_*`` alias not present in *taken*.

    Draws sequential candidates from the ``name_sequence`` generator and, on the
    rare collision with a user identifier, defers to ``find_new_name`` for a
    suffixed variant — the same two-step ``eliminate_subqueries`` uses.
    """
    candidate = next_name()
    if candidate in taken:
        candidate = find_new_name(taken, candidate)
    return candidate


def _canonical_projection(ref: ResolvedRef) -> exp.Select:
    """Build the ``SELECT`` body that projects *ref*'s table to canonical form.

    The projection exposes the canonical ``chrom`` / ``start`` / ``end`` columns
    under their original physical names, with ``start`` / ``end`` rewritten by
    the :mod:`giql.canonical` arithmetic for the table's declared encoding. This
    is the interval contract every CTE / subquery reference is assumed to satisfy
    (canonical 0-based half-open ``chrom`` / ``start`` / ``end``); operator port
    issues #122 / #123 may extend it with pass-through columns as their emitters
    require.
    """
    chrom, start, end = ref.cols
    table = ref.table
    relation = ref.name
    return exp.select(
        exp.alias_(exp.column(chrom), chrom),
        exp.alias_(_canonical_start_expr(start, table), start),
        exp.alias_(_canonical_end_expr(end, table), end),
    ).from_(exp.to_table(relation))


def _canonical_start_expr(start: str, table: Table | None) -> exp.Expression:
    """Canonical 0-based half-open start expression for a raw start column.

    SQL analog of :func:`giql.canonical.canonical_start`:

    - ``0based``: ``start`` (identity)
    - ``1based``: ``start - 1``
    """
    col = exp.column(start)
    if table is None or table.coordinate_system == "0based":
        return col
    return exp.paren(exp.Sub(this=col, expression=exp.Literal.number(1)))


def _canonical_end_expr(end: str, table: Table | None) -> exp.Expression:
    """Canonical 0-based half-open end expression for a raw end column.

    SQL analog of :func:`giql.canonical.canonical_end`:

    - ``0based`` / ``half_open``: ``end`` (identity)
    - ``0based`` / ``closed``:    ``end + 1``
    - ``1based`` / ``half_open``: ``end - 1``
    - ``1based`` / ``closed``:    ``end`` (identity)
    """
    col = exp.column(end)
    if table is None:
        return col
    key = (table.coordinate_system, table.interval_type)
    if key == ("0based", "closed"):
        return exp.paren(exp.Add(this=col, expression=exp.Literal.number(1)))
    if key == ("1based", "half_open"):
        return exp.paren(exp.Sub(this=col, expression=exp.Literal.number(1)))
    return col


def _make_cte(name: str, body: exp.Select) -> exp.CTE:
    """Wrap a projection *body* as a named CTE."""
    return exp.CTE(
        this=body,
        alias=exp.TableAlias(this=exp.to_identifier(name)),
    )


def _rewrite_slot(node: exp.Expression, arg: str, ref: ResolvedRef, name: str) -> None:
    """Point operator slot *arg* at the canonical CTE *name*.

    Rewrites both the AST slot (to a bare table reference naming the CTE) and the
    attached :class:`~giql.resolver.OperatorResolution` (a fresh
    :class:`~giql.resolver.ResolvedRef` of kind ``"cte"`` carrying no ``Table``
    config — the CTE is canonical by construction). The physical column names are
    preserved because the CTE re-exposes them under their original names.
    """
    new_ref = replace(ref, kind="cte", name=name, table=None)
    resolution = node.meta[META_KEY]
    resolution.slots[arg] = new_ref
    node.set(arg, exp.to_table(name))


def _insert_ctes(expression: exp.Expression, new_ctes: list[exp.CTE]) -> None:
    """Prepend *new_ctes* to the outermost ``WITH`` in DAG order.

    Canonical CTEs depend only on base tables, so prepending them guarantees any
    existing CTE or main-body reference resolves them. An existing ``WITH`` is
    extended in place; otherwise a fresh one is attached to the root query.
    """
    if not new_ctes:
        return
    query = expression.expression if isinstance(expression, exp.DDL) else expression
    existing = query.args.get("with_")
    if existing is not None:
        existing.set("expressions", new_ctes + list(existing.expressions))
    else:
        query.set("with_", exp.With(expressions=new_ctes))


def _decanonicalize_outputs(
    expression: exp.Expression,
    targets: list[tuple[exp.Expression, str, ResolvedRef]],
) -> None:
    """De-canonicalize migrated operator outputs in the outermost projection.

    Inert hook (epic #114, step 6). The outermost ``SELECT`` projection list
    should rewrite any output column a migrated operator emitted in canonical
    form back into the user's preferred encoding. No operator is migrated in
    issue #121, so there is nothing to rewrite; the operator port issues (#122,
    #123) fill this in alongside flipping their ``GIQL_CANONICALIZE`` flags.
    """
    # Intentionally empty until an operator opts in (see module docstring).
    return None
