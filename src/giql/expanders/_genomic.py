"""Operator-neutral CLUSTER/MERGE expansion toolkit (epic #137, issue #163).

The CLUSTER (:mod:`giql.expanders.cluster`) and MERGE
(:mod:`giql.expanders.merge`) expanders are both **whole-query rewrites** built
on the same generic plumbing: resolve the genomic columns of an enclosing FROM
source/relation, locate the operator in a top-level projection, guard the shapes
that have no coherent rewrite, and transplant a freshly-built SELECT onto the
original in place. None of that plumbing is specific to either operator.

This module is that shared surface — "operator-neutral" in that no helper is
specific to one operator. The lone exception proves the rule: the CLUSTER/MERGE
co-occurrence guard (:func:`reject_cluster_merge_mix`) references both concrete
operator types, but it is *symmetric* between them rather than agnostic of them —
a guard owned by neither expander belongs here, not in either one. It was
extracted from ``cluster.py`` (#163) so the shared toolkit has an explicit home
rather than living in one expander and being imported by its sibling — an
arrangement that made the dependency direction between the two expanders look
accidental and invited further expander-to-expander imports. The extraction is a
pure code-move: the emitted SQL is unchanged.

The name is ``_``-prefixed so the :mod:`giql.expanders` auto-discovery import
skips it — it registers no expander; it is a helper module the expanders import.
The one genuinely CLUSTER-specific reuse point, ``expand_cluster_query`` (MERGE
composes the CLUSTER rewrite to build its clustered subquery), deliberately stays
in ``cluster.py`` and is imported by ``merge.py`` as a documented composition.
"""

from __future__ import annotations

from typing import NamedTuple
from typing import TypeVar

from sqlglot import exp

from giql.constants import DEFAULT_CHROM_COL
from giql.constants import DEFAULT_END_COL
from giql.constants import DEFAULT_START_COL
from giql.constants import DEFAULT_STRAND_COL
from giql.expressions import GIQLCluster
from giql.expressions import GIQLMerge
from giql.table import Table
from giql.table import Tables

_T = TypeVar("_T", bound=exp.Expression)


class GenomicColumns(NamedTuple):
    """The resolved physical column names CLUSTER / MERGE operate over.

    Derived from the enclosing FROM table by :func:`genomic_columns`. A
    :class:`~typing.NamedTuple` so it still unpacks and indexes positionally
    (``chrom, start, end, strand = columns``) while giving the four fields names.
    """

    chrom: str
    start: str
    end: str
    strand: str


_CANONICAL_COLUMNS = GenomicColumns(
    DEFAULT_CHROM_COL, DEFAULT_START_COL, DEFAULT_END_COL, DEFAULT_STRAND_COL
)


def genomic_columns(select: exp.Select, tables: Tables) -> GenomicColumns:
    """Return the ``(chrom, start, end, strand)`` columns for *select*'s FROM source.

    Part of the shared, operator-neutral CLUSTER/MERGE expansion toolkit. Resolves
    the physical genomic columns of the enclosing FROM relation: a bare registered
    table yields its configured mapping (and the default strand column when the table
    declares none); a derived table or CTE is resolved *through* to the underlying
    table when the subquery passes its columns through unchanged (``SELECT *`` /
    ``SELECT alias.*``). An unregistered table, an empty FROM, or a projection that
    already exposes the canonical ``chrom`` / ``start`` / ``end`` columns falls back
    to those defaults.

    Unlike ``validate_projection_contracts`` / ``_resolve_disjoin_reference`` in
    ``resolver.py``, which require canonical columns from an operator's reference
    *operands*, this deliberately traces a star-passthrough FROM *source* to its
    physical columns — CLUSTER/MERGE operate over the enclosing relation's own
    coordinates, a different operator semantics. That FROM-source resolution lives in
    the CLUSTER/MERGE toolkit rather than the pass-1 resolver because it serves the
    whole-query rewrite's need to name the physical window/grouping columns, not the
    resolver's reference-slot contract validation.

    :raises ValueError: When the FROM source is an untraceable relation — a derived
        table, CTE, ``VALUES``, table-valued function, or ``LATERAL`` — whose genomic
        columns cannot be determined (it neither exposes the canonical columns nor
        resolves to a single registered mapping), or when the FROM combines multiple
        sources through a join. Silently defaulting to the canonical columns in these
        cases would emit SQL referencing columns that do not exist, so the ambiguity
        is surfaced instead (#164).
    """
    if select.args.get("joins"):
        # Only the query being rewritten carries a join that ``transplant`` would
        # drop; a join nested inside a derived-table/CTE body is preserved intact,
        # so the rejection guards the top-level FROM, not the recursive descent.
        raise ValueError(
            "CLUSTER/MERGE over a join is not supported; run it over a single "
            "relation, or over a subquery."
        )
    mapping, opaque = _resolve_from_columns(select, tables, frozenset())
    if mapping is not None:
        return mapping
    if opaque:
        raise ValueError(
            "Cannot resolve the genomic columns for CLUSTER/MERGE over this FROM "
            "clause: it is a derived table or CTE that neither exposes the canonical "
            f"{DEFAULT_CHROM_COL!r} / {DEFAULT_START_COL!r} / {DEFAULT_END_COL!r} "
            "columns nor resolves to a single registered table. Alias the genomic "
            "columns to their canonical names inside the subquery, or run "
            "CLUSTER/MERGE directly over a registered table."
        )
    return _CANONICAL_COLUMNS


def _resolve_from_columns(
    select: exp.Select, tables: Tables, seen: frozenset[str]
) -> tuple[GenomicColumns | None, bool]:
    """Resolve *select*'s driving FROM source to a ``(mapping, opaque)`` result.

    *mapping* is the resolved :class:`GenomicColumns` when a registered mapping backs
    the source, else ``None``. *opaque* is ``True`` when the source is a derived
    table/CTE whose columns could not be traced — a signal the caller turns into a
    hard error rather than a silent canonical default.

    Only the driving ``from_`` relation is resolved; any joins are ignored. A join on
    the top-level rewrite target is rejected upstream in :func:`genomic_columns`,
    while a join inside a traced-through derived-table/CTE body is preserved by the
    rewrite, so its column list is taken from the driving relation as before (#164).
    """
    from_clause = select.args.get("from_")
    if from_clause is None:
        return None, False
    return _resolve_source(from_clause.this, select, tables, seen)


def _resolve_source(
    source: exp.Expression,
    enclosing: exp.Select,
    tables: Tables,
    seen: frozenset[str],
) -> tuple[GenomicColumns | None, bool]:
    """Resolve one FROM/JOIN *source* node to a ``(mapping, opaque)`` result."""
    if isinstance(source, exp.Table):
        name = source.name
        cte_query = _find_cte(enclosing, name)
        if cte_query is not None:
            if name in seen:
                # A CTE that (transitively) selects from itself: stop recursing and
                # report it untraceable rather than looping forever.
                return None, True
            return _resolve_query(cte_query, tables, seen | {name})
        table = tables.get(name)
        if table is not None:
            return _mapping_for(table), False
        # An unregistered bare table is assumed to expose the canonical columns,
        # matching how the resolver treats unknown relations elsewhere.
        return None, False
    if isinstance(source, exp.Subquery):
        return _resolve_query(source.this, tables, seen)
    # A LATERAL, VALUES, or table-valued function has no traceable table mapping.
    return None, True


def _resolve_query(
    query: exp.Expression, tables: Tables, seen: frozenset[str]
) -> tuple[GenomicColumns | None, bool]:
    """Resolve the columns a derived-table/CTE *query* exposes to ``(mapping, opaque)``.

    Any projection *containing* a star — a bare ``*``, a qualified ``alias.*``, or a
    star alongside other items (``SELECT r.*, 1 AS extra``) — passes its FROM
    source's columns through, so the mapping is resolved from that source. An explicit
    projection instead names the output columns directly: it is transparent only when
    it re-exposes the canonical genomic columns, otherwise its genomic columns cannot
    be named and the source is reported opaque.
    """
    if isinstance(query, exp.SetOperation):
        # A set operation's arms share one column list; the left arm determines it.
        return _resolve_query(query.left, tables, seen)
    if not isinstance(query, exp.Select):
        # Defensive: a parenthesized subquery / set-operation arm the grammar can
        # reach is always a SELECT or a SetOperation (caught above, including
        # INTERSECT / EXCEPT), so this total-function guard for any other body (e.g. a
        # bare VALUES, which instead parses as an un-wrapped FROM source handled in
        # _resolve_source) is not reachable via a GIQL query.
        return None, True

    projections = query.expressions
    passes_through = any(
        isinstance(projection, exp.Star)
        or (isinstance(projection, exp.Column) and isinstance(projection.this, exp.Star))
        for projection in projections
    )
    if passes_through:
        return _resolve_from_columns(query, tables, seen)

    output_names = {
        projection.alias if isinstance(projection, exp.Alias) else projection.name
        for projection in projections
        if isinstance(projection, (exp.Alias, exp.Column))
    }
    canonical = {DEFAULT_CHROM_COL, DEFAULT_START_COL, DEFAULT_END_COL}
    if canonical <= output_names:
        return None, False
    return None, True


def _find_cte(select: exp.Select, name: str) -> exp.Expression | None:
    """Return the query body of an enclosing CTE named *name*, or ``None``."""
    node: exp.Expression | None = select
    while node is not None:
        with_clause = node.args.get("with_")
        if with_clause is not None:
            for cte in with_clause.expressions:
                if cte.alias.lower() == name.lower():
                    return cte.this
        node = node.parent
    return None


def _mapping_for(table: Table) -> GenomicColumns:
    """Return the genomic column mapping configured on a registered *table*."""
    return GenomicColumns(
        table.chrom_col,
        table.start_col,
        table.end_col,
        table.strand_col or DEFAULT_STRAND_COL,
    )


def extract_stranded(stranded_expr: exp.Expression | None) -> bool:
    """Coerce a CLUSTER/MERGE ``stranded`` operand to a bool. Shared toolkit.

    Mirrors the legacy per-transformer coercion exactly: a missing operand is
    ``False``; an ``exp.Boolean`` yields its raw ``.this``; an ``exp.Literal``
    compares case-folded to ``TRUE``. The final two arms (``exp.Literal`` and the
    string-truthiness fallback) are defensive — the GIQL grammar only ever produces
    ``exp.Boolean`` for ``stranded := <bool>`` — retained for parity with the
    legacy port.
    """
    if stranded_expr is None:
        return False
    if isinstance(stranded_expr, exp.Boolean):
        return stranded_expr.this
    if isinstance(stranded_expr, exp.Literal):
        return str(stranded_expr.this).upper() == "TRUE"
    return str(stranded_expr).upper() in ("TRUE", "1", "YES")


def find_projected(select: exp.Select, op_type: type[_T]) -> list[_T]:
    """Return *select*'s projected operators of *op_type* (bare or aliased). Toolkit.

    Shared CLUSTER/MERGE primitive: both expanders and the co-occurrence guard
    locate their operator the same way — a top-level SELECT projection item that
    either *is* the operator or is an ``exp.Alias`` wrapping it.
    """
    found: list[_T] = []
    for expression in select.expressions:
        if isinstance(expression, op_type):
            found.append(expression)
        elif isinstance(expression, exp.Alias) and isinstance(expression.this, op_type):
            found.append(expression.this)
    return found


def require_top_level_projection(
    select: exp.Select, node: exp.Expression, op_type: type
) -> None:
    """Raise if *node* is buried inside a projection expression. Shared toolkit.

    A CLUSTER / MERGE is only expandable as a *top-level* projection item — bare
    or directly aliased — because the whole-query rewrite restructures the SELECT
    around it. One nested inside a larger projection expression such as
    ``ABS(CLUSTER(interval))`` has no coherent rewrite and would otherwise leak an
    unexpanded operator to the generator, so fail loudly here (#144 A16). An
    operator that parses *outside* the projection entirely (e.g. in WHERE /
    ORDER BY) is not under any projection item and is left for the expander's
    existing no-op path.
    """
    operator = op_type.__name__.removeprefix("GIQL").upper()
    for projection in select.expressions:
        inner = projection.this if isinstance(projection, exp.Alias) else projection
        if inner is node:
            return
        if any(descendant is node for descendant in inner.walk()):
            raise ValueError(
                f"{operator} must be a top-level projection item; it cannot be "
                "nested inside another expression (e.g. a function call or "
                "arithmetic)."
            )


def reject_cluster_merge_mix(select: exp.Select) -> None:
    """Raise if *select* projects both a CLUSTER and a MERGE. Shared toolkit.

    The two are mutually incompatible in one SELECT: MERGE aggregates the rows
    away (it rewrites the query into a ``GROUP BY`` over a clustered subquery)
    while CLUSTER is a per-row window over those same rows, so no coherent single
    query expresses both. The legacy pre-pass chained the two transformers and
    emitted *non-executable* SQL for this shape (a window over ``GROUP BY``-
    aggregated rows — a DuckDB ``BinderException``, never a leaked operator). The
    new in-place ``transplant`` cannot express both at all: whichever expander
    runs first rewrites the shared SELECT and strands the sibling as an unexpanded
    node. Fail loudly here — mirroring the ``Multiple MERGE expressions not yet
    supported`` guard — so the combination raises a clear diagnostic rather than
    emitting broken SQL.
    """
    if find_projected(select, GIQLCluster) and find_projected(select, GIQLMerge):
        raise ValueError(
            "CLUSTER and MERGE cannot be combined in a single SELECT; MERGE "
            "aggregates rows while CLUSTER is a per-row window. Use them in "
            "separate queries (e.g. CLUSTER over a subquery, or MERGE over one)."
        )


def transplant(select: exp.Select, new: exp.Select) -> None:
    """Replace *select*'s contents with *new*'s, preserving *select*'s identity.

    Part of the shared CLUSTER/MERGE expansion toolkit. Clears every argument of
    *select* and re-installs *new*'s, so *select* keeps its position in the
    surrounding tree (and its identity as the object the pass returns) while taking
    on the rewritten structure. This is how a whole-query rewrite is applied to a
    *root* SELECT, which has no parent to ``replace`` through.

    Precondition: *new* MUST be a detached throwaway — a freshly-built
    ``exp.Select``, as the ``_transform_for_*`` helpers return. Its children are
    re-parented onto *select*, so passing a node still attached elsewhere would
    corrupt that other tree.
    """
    assert new.parent is None, "transplant() requires a detached `new` subtree"
    select.args.clear()
    for key, value in list(new.args.items()):
        select.set(key, value)
