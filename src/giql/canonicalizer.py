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

Engine portability (capability-driven, issue #145)
---------------------------------------------------
The wrapper projection canonicalizes the interval columns while passing every
other source column through untouched (the registry declares only the genomic
columns, so an explicit full-column projection is not available). The emit
strategy is chosen from the active target's :class:`~giql.targets.Capabilities`,
following the DISJOIN passthrough's precedent (#143) — the same two emit forms,
though this wrapper additionally accepts ``capabilities=None`` (a direct caller
default, see :func:`canonicalize_coordinates`):

* ``SELECT * REPLACE (...)`` when ``capabilities.supports_star_replace`` holds
  (DuckDB / BigQuery / Snowflake / ClickHouse) — substitutes start/end in place,
  preserving source column order;
* the portable ``SELECT * EXCEPT (start, end), <start>, <end>`` form otherwise
  (the generic baseline / DataFusion family), which every ``* EXCEPT``-capable
  engine plans. This form is row-equivalent but **not column-order-equivalent**:
  ``* EXCEPT`` drops the interval columns and re-appends the recomputed ones at
  the end of the projection. It is also not SQL-92 and not DuckDB-runnable.

Identity-encoded (default 0-based half-open) relations are unaffected either way:
they skip wrapping entirely and emit portable SQL. The capability is threaded in
from :func:`giql.transpile.transpile` via the active target; a direct caller that
passes no capabilities defaults to the ``* REPLACE`` form (the historical
behavior). This finalizes the dialect-aware emit strategy formerly tracked by
https://github.com/abdenlab/giql/issues/132.

Gating (epic #114, step 6)
--------------------------
The pass is gated per operator by a ``GIQL_CANONICALIZE`` class attribute on the
operator's expression class. An operator opts in by setting
``GIQL_CANONICALIZE = True``; absent or ``False`` the pass ignores it entirely.
The operator port issues — #122 (DISJOIN) and #123 (NEAREST / DISTANCE /
predicates) — flipped these flags as each operator's emitter moved off in-emitter
canonicalization (:mod:`giql.canonical`) and onto this pass's output. As of those
ports every migrated operator opts in by default, so the pass actively
synthesizes wrappers; an operator can still toggle its flag off (a test or a
not-yet-migrated operator), in which case the pass leaves it untouched and the
emitted SQL is byte-identical for it.

De-canonicalization hook
-------------------------
A migrated operator's *output* columns must land back in the target relation's
declared encoding. Epic #114 step 6 envisioned a rewrite of the outermost
``SELECT`` projection, but that placement is wrong for a table function: DISJOIN
and NEAREST synthesize their output columns and pass whole source rows through at
*generation* time, so those columns do not exist as AST in this pass, and a
``SELECT *`` consumer hides them from any outer-projection rewrite. So
:func:`_decanonicalize_outputs` instead records each wrapped slot's *original*
:class:`~giql.table.Table` on the operator's
:class:`~giql.resolver.OperatorResolution`, and the operator's emitter reads it
to de-canonicalize those synthesized columns where it generates them (DISJOIN,
issue #122; NEAREST's target row passthrough, issue #123).

Column / interval operands (epic #114, step 8 / issue #123)
-----------------------------------------------------------
A reference slot (DISJOIN/NEAREST target, DISJOIN reference) *owns* the relation
it names, so the pass can wrap that whole relation in a ``__giql_canon_*`` CTE and
redirect the slot's AST node to it. A *column* operand cannot be wrapped that way:
``DISTANCE(a.interval, b.interval)`` and ``a.interval INTERSECTS b.interval``
reference an alias bound in the *enclosing* query's ``FROM`` / ``JOIN``, shared
with the user's own projection (``SELECT a.start, DISTANCE(...)``). Rewriting that
``FROM`` source to a canonical CTE would silently canonicalize the user's own
``a.start`` too — a behavior change. NEAREST's column / implicit-outer
``reference`` slot is the same shape (an alias from the outer LATERAL relation).

For those operands the pass therefore canonicalizes the resolution metadata *in
place* rather than synthesizing a CTE: :func:`_canonicalize_column_operands`
rewrites each :class:`~giql.resolver.ResolvedColumn` (DISTANCE's two operands and
the spatial predicates' column operands) and each non-table
:class:`~giql.resolver.ResolvedInterval` (NEAREST's ``column`` /
``implicit_outer`` reference) so its ``start`` / ``end`` fragments carry the
canonical 0-based half-open arithmetic and its ``table`` is blanked. The emitter
then consumes the fragments verbatim — no in-emitter
:func:`giql.canonical.canonical_start` / ``canonical_end``. The arithmetic is the
same the emitter used to emit inline, so the SQL stays byte-identical for these
operands; only the *owner* of the arithmetic moves from the generator to the pass.
A ``literal_range`` interval is already canonical and is left untouched, and an
operand whose ``table`` is already canonical (or ``None``) is a no-op.
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
from giql.resolver import ResolvedColumn
from giql.resolver import ResolvedInterval
from giql.resolver import ResolvedRef
from giql.table import Table
from giql.targets import Capabilities

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


def canonicalize_coordinates(
    expression: exp.Expression, capabilities: Capabilities | None = None
) -> exp.Expression:
    """Synthesize canonical wrapper CTEs for non-canonical operator operands.

    Walks *expression* for opted-in GIQL operators (those whose expression class
    sets ``GIQL_CANONICALIZE = True``), and for each reference slot that resolved
    to a registered table with a non-canonical encoding, synthesizes a
    ``__giql_canon_*`` WITH-CTE projecting that table to canonical 0-based
    half-open coordinates and rewrites the slot (AST node + ``ResolvedRef``
    metadata) to point at the canonical CTE.

    The pass mutates and returns *expression* in place. For an operator whose
    ``GIQL_CANONICALIZE`` flag is off, or whose operands are already in the
    canonical 0-based half-open encoding, it touches nothing and leaves the
    emitted SQL byte-identical.

    Parameters
    ----------
    expression : exp.Expression
        The pass-1-annotated AST.
    capabilities : Capabilities | None
        The active target's capabilities, used to choose the wrapper projection's
        emit strategy (``* REPLACE`` vs the portable ``* EXCEPT`` form — see the
        module docstring). :func:`giql.transpile.transpile` passes the active
        target's capabilities; ``None`` (a direct caller) defaults to the
        ``* REPLACE`` form, preserving the historical behavior.

    Returns
    -------
    exp.Expression
        The same *expression*, with canonical wrapper CTEs inserted and the
        opted-in operator slots that reference non-canonical tables rewritten to
        point at them.
    """
    # Column / interval operands (DISTANCE, predicates, NEAREST's non-table
    # reference) canonicalize their metadata in place; this is independent of the
    # ref-slot CTE synthesis below and runs for every opted-in operator.
    _canonicalize_column_operands(expression)

    targets = _collect_targets(expression)
    if not targets:
        # No opted-in operator carries a non-canonical *reference-slot* operand:
        # no wrapper CTE is synthesized (the in-place column canonicalization
        # above already ran).
        return expression

    taken = _collect_taken_names(expression)
    next_name = name_sequence(CANON_PREFIX)
    body_to_name: dict[str, str] = {}
    new_ctes: list[exp.CTE] = []

    for node, arg, ref in targets:
        body = _canonical_projection(ref, capabilities)
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


def _canonicalize_column_operands(expression: exp.Expression) -> None:
    """Canonicalize column / interval operand metadata in place for opted-in ops.

    For every opted-in operator (``GIQL_CANONICALIZE``) carrying column operands
    (DISTANCE's two operands and the spatial predicates' column operands, in the
    :attr:`OperatorResolution.columns` channel) or a non-table interval reference
    (NEAREST's ``column`` / ``implicit_outer`` slot, in
    :attr:`OperatorResolution.slots`), rewrite the operand's ``start`` / ``end``
    SQL fragments to carry the canonical 0-based half-open arithmetic for its
    declared encoding and blank its ``table``.

    This replaces the in-emitter :func:`giql.canonical.canonical_start` /
    ``canonical_end`` wrapping for those operands (epic #114, step 8). Unlike the
    ref-slot CTE synthesis, no relation is wrapped: the operand references an
    alias bound in the enclosing query's ``FROM`` shared with the user's own
    projection, so canonicalizing the whole relation would change unrelated
    columns. Operands already canonical (``table`` is ``None`` or 0-based
    half-open) are left untouched, keeping their SQL byte-identical; a
    ``literal_range`` interval is already canonical and is skipped.
    """
    for node in expression.walk():
        if not isinstance(node, _OPERATORS):
            continue
        if not getattr(node, "GIQL_CANONICALIZE", False):
            continue
        resolution = node.meta.get(META_KEY)
        if not isinstance(resolution, OperatorResolution):
            continue
        for arg, column in list(resolution.columns.items()):
            resolution.columns[arg] = _canonicalize_column(column)
        for arg, slot in list(resolution.slots.items()):
            if isinstance(slot, ResolvedInterval):
                resolution.slots[arg] = _canonicalize_interval(slot)


def _canonicalize_column(column: ResolvedColumn) -> ResolvedColumn:
    """Return *column* with canonical start/end fragments and a blanked table.

    A no-op (returns *column* unchanged) when its backing table is already
    canonical or ``None``.
    """
    if _is_canonical(column.table):
        return column
    return replace(
        column,
        start=_canonical_start_sql(column.start, column.table),
        end=_canonical_end_sql(column.end, column.table),
        table=None,
    )


def _canonicalize_interval(interval: ResolvedInterval) -> ResolvedInterval:
    """Return *interval* with canonical start/end fragments and a blanked table.

    A no-op for a ``literal_range`` (already canonical, no table) and for any
    interval whose backing table is already canonical or ``None``.
    """
    if interval.kind == "literal_range" or _is_canonical(interval.table):
        return interval
    return replace(
        interval,
        start=_canonical_start_sql(interval.start, interval.table),
        end=_canonical_end_sql(interval.end, interval.table),
        table=None,
    )


def _canonical_start_sql(start: str, table: Table | None) -> str:
    """SQL-fragment analog of :func:`giql.canonical.canonical_start`.

    - ``0based``: ``start`` (identity)
    - ``1based``: ``(start - 1)``
    """
    if table is None or table.coordinate_system == "0based":
        return start
    return f"({start} - 1)"


def _canonical_end_sql(end: str, table: Table | None) -> str:
    """SQL-fragment analog of :func:`giql.canonical.canonical_end`.

    - ``0based`` / ``half_open``: ``end`` (identity)
    - ``0based`` / ``closed``:    ``(end + 1)``
    - ``1based`` / ``half_open``: ``(end - 1)``
    - ``1based`` / ``closed``:    ``end`` (identity)
    """
    if table is None:
        return end
    key = (table.coordinate_system, table.interval_type)
    if key == ("0based", "closed"):
        return f"({end} + 1)"
    if key == ("1based", "half_open"):
        return f"({end} - 1)"
    return end


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


def _canonical_projection(
    ref: ResolvedRef, capabilities: Capabilities | None
) -> exp.Select:
    """Build the ``SELECT`` body that projects *ref*'s table to canonical form.

    The projection is a **full-row passthrough**: ``SELECT *`` keeps every
    physical column of the source relation, and only the two interval columns —
    ``start`` / ``end``, under their original physical names — are rewritten with
    the :mod:`giql.canonical` arithmetic for the table's declared encoding.
    ``chrom`` and every non-interval column flow through the star untouched.

    The emit strategy is chosen from *capabilities* (issue #145), following the
    precedent of :func:`giql.expanders.disjoin._disjoin_passthrough` — the same
    two emit forms, with an added ``capabilities is None`` arm for direct callers
    (the passthrough always receives a concrete ``ctx.capabilities``):

    * ``SELECT * REPLACE (...)`` when ``supports_star_replace`` holds (or no
      capabilities are supplied) — substitutes the interval columns in place,
      preserving source column order;
    * the portable ``SELECT * EXCEPT (start, end), <start>, <end>`` form otherwise
      — drops the interval columns from the star and re-appends them recomputed.
      Row-equivalent but not column-order-equivalent, and not DuckDB-runnable.

    The full row (rather than a bare ``chrom`` / ``start`` / ``end`` triple) is
    required by table-function operators whose final projection passes the whole
    source row through — DISJOIN's ``SELECT t.*`` (#122) — and by their join-back
    semantics, which key on the source's physical columns. A CTE / subquery
    reference that only needs the canonical interval triple still reads those
    three columns from the same wrapper.
    """
    _chrom, start, end = ref.cols
    table = ref.table
    relation = ref.name
    # Quote the interval identifiers: the canonical column names are physical and
    # routinely reserved words (the default genomic layout's ``start`` / ``end``),
    # so the executed wrapper must quote them.
    start_id = exp.to_identifier(start, quoted=True)
    end_id = exp.to_identifier(end, quoted=True)
    start_proj = exp.alias_(_canonical_start_expr(start, table), start_id)
    end_proj = exp.alias_(_canonical_end_expr(end, table), end_id)
    if capabilities is None or capabilities.supports_star_replace:
        star = exp.Star(replace=[start_proj, end_proj])
        return exp.Select(expressions=[star]).from_(exp.to_table(relation))
    # Portable form: drop the interval columns from the star and re-project them
    # recomputed under their own names. EXCEPT removes them from the row; the
    # trailing projections add them back in canonical form.
    star = exp.Star(except_=[exp.column(start_id), exp.column(end_id)])
    return exp.Select(expressions=[star, start_proj, end_proj]).from_(
        exp.to_table(relation)
    )


def _canonical_start_expr(start: str, table: Table | None) -> exp.Expression:
    """Canonical 0-based half-open start expression for a raw start column.

    SQL analog of :func:`giql.canonical.canonical_start`:

    - ``0based``: ``start`` (identity)
    - ``1based``: ``start - 1``
    """
    col = exp.column(exp.to_identifier(start, quoted=True))
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
    col = exp.column(exp.to_identifier(end, quoted=True))
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
    """Preserve each wrapped slot's original encoding for the emitter's output.

    A wrapped slot's :class:`~giql.resolver.ResolvedRef` is rewritten to a
    ``Table``-free canonical-CTE ref, which would otherwise lose the
    (non-canonical) encoding the operator's *output* must round-trip back into.

    The de-canonicalization itself cannot be applied on the AST in this pass for
    a table-function operator: DISJOIN synthesizes its ``disjoin_*`` columns and
    its passed-through interval at *generation* time, so those columns do not
    exist as AST here, and a ``SELECT *`` consumer hides them from any
    outer-projection rewrite. The originally-envisioned outermost-projection
    rewrite (epic #114, step 6) is therefore wrong for projected
    table-function columns; instead this hook records the per-slot original
    :class:`~giql.table.Table` on the :class:`~giql.resolver.OperatorResolution`,
    and the operator's emitter reads it to de-canonicalize those synthesized
    columns where it generates them (see :issue:`122`).

    *targets* carries the original (pre-rewrite) refs, so ``ref.table`` is the
    source relation's declared encoding.
    """
    for node, arg, ref in targets:
        resolution = node.meta.get(META_KEY)
        if isinstance(resolution, OperatorResolution):
            resolution.output_tables[arg] = ref.table
