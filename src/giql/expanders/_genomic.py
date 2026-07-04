"""Operator-neutral CLUSTER/MERGE expansion toolkit (epic #137, issue #163).

The CLUSTER (:mod:`giql.expanders.cluster`) and MERGE
(:mod:`giql.expanders.merge`) expanders are both **whole-query rewrites** built
on the same generic plumbing: resolve the genomic columns of an enclosing FROM
table, locate the operator in a top-level projection, guard the shapes that have
no coherent rewrite, and transplant a freshly-built SELECT onto the original in
place. None of that plumbing is specific to either operator.

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


def genomic_columns(select: exp.Select, tables: Tables) -> GenomicColumns:
    """Return the ``(chrom, start, end, strand)`` columns for *select*'s FROM table.

    Part of the shared, operator-neutral CLUSTER/MERGE expansion toolkit. Mirrors
    the legacy ``ClusterTransformer._get_genomic_columns`` / ``_get_table_name``:
    read the FROM-clause table name, look it up in *tables*, and use its configured
    column names, falling back to the canonical defaults (and to the default strand
    column when the table declares none).
    """
    table_name: str | None = None
    from_clause = select.args.get("from_")
    if from_clause is not None and isinstance(from_clause.this, exp.Table):
        table_name = from_clause.this.name

    chrom_col = DEFAULT_CHROM_COL
    start_col = DEFAULT_START_COL
    end_col = DEFAULT_END_COL
    strand_col = DEFAULT_STRAND_COL

    if table_name:
        table = tables.get(table_name)
        if table:
            chrom_col = table.chrom_col
            start_col = table.start_col
            end_col = table.end_col
            if table.strand_col:
                strand_col = table.strand_col

    return GenomicColumns(chrom_col, start_col, end_col, strand_col)


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
