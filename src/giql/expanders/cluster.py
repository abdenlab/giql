"""The CLUSTER operator expander (epic #137, issue #144).

CLUSTER assigns a cluster id to every interval, grouping each run of mutually
adjacent intervals (optionally within a maximum gap, strand-aware, and gated by a
pairwise predicate) under one id. It cannot be a single window function because it
needs a window over a window (``SUM`` over a ``LAG``-derived flag), so the
expansion restructures the enclosing query into a two-level form::

    SELECT *, CLUSTER(interval) AS cluster_id FROM features

becomes::

    SELECT *, SUM(is_new_cluster) OVER (PARTITION BY chrom ORDER BY start) AS cluster_id
    FROM (
        SELECT *,
               CASE WHEN LAG(end) OVER (...) >= start THEN 0 ELSE 1 END AS is_new_cluster
        FROM features
    ) AS lag_calc

(The ``becomes::`` form is simplified for readability; the emitted SQL quotes
identifiers and appends ``NULLS LAST`` to the window ORDER BY.)

This module is the AST-expansion replacement for the legacy
:class:`giql.transformer.ClusterTransformer`, which ran as a pre-pass transformer
on the raw parsed AST. It produces the same SQL; the existing CLUSTER
transpilation and bedtools oracle tests are the migration oracle.

Unlike the node-local expanders for DISTANCE / NEAREST / DISJOIN — whose result
replaces a single node — CLUSTER is a **whole-query rewrite**. It shares only the
*shape* of :func:`giql.expanders.nearest._fallback_form` — return the operator
node unchanged so the pass's ``node.replace`` is a no-op — but the *mechanism*
differs. ``nearest`` does an in-place child ``.replace`` because its LATERAL has a
parent to replace through; CLUSTER instead copies the enclosing
:class:`~sqlglot.expressions.Select`, restructures the copy, and *transplants* its
contents back onto the original SELECT (see :func:`transplant`). That root-
preserving transplant is required because the canonical
``SELECT *, CLUSTER(...) FROM t`` puts CLUSTER at the *root* ``SELECT``, which has
no parent to replace it through; transplanting preserves the root's identity (the
pass returns the same expression object it was handed).

The pass walks the tree and collects every operator node, expanding deepest-first,
so this expander never recurses into CTEs / subqueries: a nested CLUSTER is
collected and expanded on its own. The shared cluster-restructure helper
(:func:`expand_cluster_query`) is reused by :mod:`giql.expanders.merge`, since
MERGE is built on CLUSTER.
"""

from __future__ import annotations

from typing import NamedTuple
from typing import TypeVar

from sqlglot import exp

from giql.constants import DEFAULT_CHROM_COL
from giql.constants import DEFAULT_END_COL
from giql.constants import DEFAULT_START_COL
from giql.constants import DEFAULT_STRAND_COL
from giql.expander import ExpansionContext
from giql.expander import expand_operators
from giql.expander import register
from giql.expressions import GIQLCluster
from giql.expressions import GIQLMerge
from giql.table import Tables
from giql.targets import GenericTarget

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


@register(GenericTarget, GIQLCluster)
def expand_cluster(node: GIQLCluster, ctx: ExpansionContext) -> exp.Expression:
    """Expand a CLUSTER node by restructuring its enclosing SELECT in place.

    Registered for :class:`~giql.targets.GenericTarget`, so every target resolves
    to it through the registry's generic chain (CLUSTER emits identical SQL across
    targets). Locates the enclosing :class:`~sqlglot.expressions.Select`, derives
    the genomic columns from the FROM table, and rewrites that SELECT in place into
    the two-level ``lag_calc`` form. Returns the original CLUSTER node (now
    unreachable from the rewritten root SELECT) so the pass's ``node.replace`` is a
    no-op.

    Parameters
    ----------
    node : GIQLCluster
        The CLUSTER node being expanded.
    ctx : ExpansionContext
        The expansion context; only its :attr:`~ExpansionContext.tables` is read
        (CLUSTER derives its columns from the enclosing FROM table, not from
        resolution metadata).

    Returns
    -------
    exp.Expression
        The CLUSTER node, unchanged — the surrounding SELECT was mutated in place.
    """
    select = node.find_ancestor(exp.Select)
    if select is None:
        # Defensive: CLUSTER only ever parses inside a SELECT projection. With no
        # enclosing SELECT there is nothing to restructure; leave the node for the
        # generator to error on, exactly as the legacy transformer's
        # non-Select guard did.
        return node
    require_top_level_projection(select, node, GIQLCluster)
    reject_cluster_merge_mix(select)
    columns = genomic_columns(select, ctx.tables)
    # Build the transformed query from a detached copy so the intermediate never
    # aliases live nodes, then transplant its args onto the original SELECT to
    # preserve its identity (and so a root SELECT is rewritten without a parent to
    # replace it through).
    source = select.copy()
    transformed = expand_cluster_query(source, columns)
    if transformed is None:
        # No-op for a CLUSTER that parses outside the SELECT projection (e.g. in
        # WHERE / ORDER BY): find_projected finds none in the copy, so there is
        # nothing to restructure and the operator leaks to the generator exactly
        # as on `main`. require_top_level_projection above already rejects the
        # in-projection-expression case; this guards the out-of-projection case.
        return node
    transplant(select, transformed)
    # copy()+transplant duplicated the enclosing WHERE / HAVING into the new
    # lag_calc subquery; the originals the pass collected are now unreachable, so
    # re-run the pass over the restructured SELECT to expand any sibling pass-3
    # operators (spatial predicates, DISTANCE) carried into it. Safe from
    # recursion: the CLUSTER node is already replaced by its SUM window. (#144 B1)
    #
    # expand_operators may return a new root (a registered statement finalizer can
    # wrap it), and its contract requires callers to use the return value rather
    # than assume in-place mutation, so reinstall a new root in place of `select`.
    # Today the branch is never taken here — the sole finalizer-registering operator
    # (a correlated NEAREST fallback) is already a plain join by this deepest-first
    # re-walk, and its wrapper targets an inner SELECT rather than this re-walk root,
    # so `result is select` in practice — but honoring the contract keeps the seam
    # future-proof. A NEAREST fallback whose reserved columns are re-surfaced by this
    # enclosing CLUSTER `SELECT *` remains a documented residual (#172), not a lost
    # root.
    result = expand_operators(select, ctx.target, ctx.tables, ctx.registry)
    if result is not select:
        select.replace(result)
    return node


def genomic_columns(select: exp.Select, tables: Tables) -> GenomicColumns:
    """Return the ``(chrom, start, end, strand)`` columns for *select*'s FROM table.

    Part of the shared CLUSTER/MERGE expansion toolkit (reused by
    :mod:`giql.expanders.merge`). Mirrors the legacy
    ``ClusterTransformer._get_genomic_columns`` / ``_get_table_name``: read the
    FROM-clause table name, look it up in *tables*, and use its configured column
    names, falling back to the canonical defaults (and to the default strand
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


def expand_cluster_query(
    query: exp.Select, columns: GenomicColumns
) -> exp.Select | None:
    """Restructure *query* for every CLUSTER in its projection; ``None`` if none.

    Finds CLUSTER expressions in *query*'s SELECT list and rewrites the query into
    the two-level ``lag_calc`` form once per CLUSTER (chaining, as the legacy
    transformer did). Operates on *query* directly and returns the rewritten
    query. Returns ``None`` when *query* projects no CLUSTER, so callers can treat
    that as a no-op.

    Reused by :mod:`giql.expanders.merge`, whose intermediate clustered query is
    restructured through this same helper.
    """
    cluster_exprs = find_projected(query, GIQLCluster)
    if not cluster_exprs:
        return None
    if len(cluster_exprs) > 1:
        # Chaining the rewrite per CLUSTER yields a duplicate ``lag_calc`` alias
        # and an ``is_new_cluster`` binder error — non-executable SQL. Fail loudly,
        # mirroring the multiple-MERGE guard, rather than emitting it (#144 A15).
        raise ValueError("Multiple CLUSTER expressions not yet supported")
    for cluster_expr in cluster_exprs:
        query = _transform_for_cluster(query, cluster_expr, columns)
    return query


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


def _transform_for_cluster(
    query: exp.Select, cluster_expr: GIQLCluster, columns: GenomicColumns
) -> exp.Select:
    """Rewrite *query* into the two-level ``lag_calc`` form for one CLUSTER.

    A byte-for-byte port of the legacy
    ``ClusterTransformer._transform_for_cluster``, with the genomic columns passed
    in (the legacy method re-derived them from the FROM table) rather than read off
    ``self``. Builds an inner ``lag_calc`` subquery that materializes an
    ``is_new_cluster`` flag from a ``LAG`` window, then an outer query whose
    ``SUM(is_new_cluster) OVER (...)`` window replaces the CLUSTER projection.
    """
    chrom_col, start_col, end_col, strand_col = columns

    # Extract CLUSTER parameters
    distance_expr = cluster_expr.args.get("distance")

    # Handle distance parameter - could be int literal or None
    if distance_expr:
        if isinstance(distance_expr, exp.Literal):
            distance = int(distance_expr.this)
        else:
            # Defensive: the grammar only yields exp.Literal for a distance
            # argument, so this non-Literal fallback is unreachable via the parser;
            # retained for parity with the legacy port.
            try:
                distance = int(str(distance_expr.this))
            except (ValueError, AttributeError):
                distance = 0
    else:
        distance = 0

    stranded = extract_stranded(cluster_expr.args.get("stranded"))

    # Build partition clause
    partition_cols = [exp.column(chrom_col, quoted=True)]
    if stranded:
        partition_cols.append(exp.column(strand_col, quoted=True))

    # Build ORDER BY for window
    order_by = [exp.Ordered(this=exp.column(start_col, quoted=True))]

    # Create LAG window spec
    lag_window = exp.Window(
        this=exp.Anonymous(this="LAG", expressions=[exp.column(end_col, quoted=True)]),
        partition_by=partition_cols,
        order=exp.Order(expressions=order_by),
    )

    # Add distance offset if specified
    if distance > 0:
        lag_with_distance = exp.Add(
            this=lag_window, expression=exp.Literal.number(distance)
        )
    else:
        lag_with_distance = lag_window

    # Build the adjacency condition (predecessor end >= current start).
    adjacency = exp.GTE(
        this=lag_with_distance,
        expression=exp.column(start_col, quoted=True),
    )

    # An optional predicate further restricts which adjacent intervals
    # are kept together: a row stays in the current cluster only when it
    # is adjacent to its predecessor AND the predicate holds between them.
    # ``PREV(col)`` references in the predicate resolve to the predecessor
    # row via LAG over the same partition/order as the adjacency window.
    predicate_expr = cluster_expr.args.get("predicate")
    if predicate_expr is not None:
        rewritten_predicate = _rewrite_predecessor_refs(
            predicate_expr, partition_cols, order_by
        )
        keep_together = exp.And(
            this=adjacency,
            expression=exp.Paren(this=rewritten_predicate),
        )
    else:
        keep_together = adjacency

    # Create CASE expression for is_new_cluster
    case_expr = exp.Case(
        ifs=[
            exp.If(
                this=keep_together,
                true=exp.Literal.number(0),
            )
        ],
        default=exp.Literal.number(1),
    )

    # Build CTE SELECT expressions (all original except CLUSTER, plus is_new_cluster)
    cte_expressions = []
    for expression in query.expressions:
        # Skip CLUSTER expressions
        if isinstance(expression, GIQLCluster):
            continue
        elif isinstance(expression, exp.Alias) and isinstance(
            expression.this, GIQLCluster
        ):
            continue
        else:
            cte_expressions.append(expression)

    # Ensure required columns for window functions are included
    required_cols = {chrom_col, start_col, end_col}
    if stranded:
        required_cols.add(strand_col)

    # The predicate is evaluated inside the lag_calc CTE, so every column
    # it references (current-row columns and PREV() arguments alike) must
    # be projected into that CTE. Folding them into required_cols makes the
    # scope dependency explicit and keeps the columns available even when a
    # later operator wraps this query in a further subquery.
    if predicate_expr is not None:
        required_cols |= {col.name for col in predicate_expr.find_all(exp.Column)}

    # Check if required columns are already in the select list
    selected_cols = set()
    for expr in cte_expressions:
        if isinstance(expr, exp.Column):
            selected_cols.add(expr.name)
        elif isinstance(expr, exp.Alias):
            # Don't count aliases as the source column
            pass
        elif isinstance(expr, exp.Star):
            # SELECT * includes all columns
            selected_cols = required_cols  # Assume all are covered
            break

    # Add missing required columns
    # Sort the residual so the injected-column order is deterministic across runs
    # (set-difference iteration order is PYTHONHASHSEED-dependent; the legacy port
    # left it nondeterministic — see review #144 A2).
    for col in sorted(required_cols - selected_cols):
        cte_expressions.append(exp.column(col, quoted=True))

    # Add is_new_cluster calculation
    # NOTE: synthesized name; the missing __giql_ reserved prefix is left to #161.
    cte_expressions.append(exp.alias_(case_expr, "is_new_cluster", quoted=False))

    # Build CTE query
    cte_select = exp.Select()
    cte_select.select(*cte_expressions, copy=False)

    # Copy FROM, WHERE, GROUP BY, HAVING from original (but not ORDER BY)
    # Use copy() to avoid sharing references between queries
    if query.args.get("from_"):
        from_clause = query.args["from_"].copy()
        cte_select.set("from_", from_clause)
    if query.args.get("where"):
        cte_select.set("where", query.args["where"].copy())
    if query.args.get("group"):
        cte_select.set("group", query.args["group"].copy())
    if query.args.get("having"):
        cte_select.set("having", query.args["having"].copy())

    # Create outer query with SUM over is_new_cluster
    sum_window = exp.Window(
        this=exp.Sum(this=exp.column("is_new_cluster")),
        partition_by=partition_cols,
        order=exp.Order(expressions=order_by),
    )

    # Build outer SELECT expressions (replace CLUSTER with SUM)
    new_expressions = []
    for expression in query.expressions:
        if isinstance(expression, GIQLCluster):
            new_expressions.append(sum_window)
        elif isinstance(expression, exp.Alias) and isinstance(
            expression.this, GIQLCluster
        ):
            # Keep the alias but replace the expression
            new_expressions.append(
                exp.alias_(sum_window, expression.alias, quoted=False)
            )
        else:
            new_expressions.append(expression)

    # Build new query
    new_query = exp.Select()
    new_query.select(*new_expressions, copy=False)

    # Wrap CTE in subquery and set as FROM clause
    # NOTE: synthesized name; the missing __giql_ reserved prefix is left to #161.
    subquery = exp.Subquery(
        this=cte_select,
        alias=exp.TableAlias(this=exp.Identifier(this="lag_calc")),
    )
    new_query.from_(subquery, copy=False)

    # Copy ORDER BY from original to outer query
    if query.args.get("order"):
        new_query.order_by(*query.args["order"].expressions, copy=False)

    return new_query


def _rewrite_predecessor_refs(
    predicate: exp.Expression,
    partition_cols: list[exp.Expression],
    order_by: list[exp.Ordered],
) -> exp.Expression:
    """Rewrite ``PREV(col)`` calls in a predicate to LAG windows.

    A byte-for-byte port of the legacy
    ``ClusterTransformer._rewrite_predecessor_refs``. Bare column references in the
    predicate denote the current interval. Each ``PREV(col)`` call denotes the
    sorted predecessor's value of that column and is rewritten to ``LAG(col) OVER
    (...)`` using the same partition/order as the cluster's adjacency window, so
    the predicate is evaluated pairwise against the immediately preceding row.
    Every column identifier (current-row columns and LAG arguments alike) is quoted
    so that reserved-word genomic columns such as ``start`` / ``end`` are emitted as
    valid SQL.

    :param predicate:
        Boolean predicate expression to rewrite (not mutated).
    :param partition_cols:
        Window partition columns (chromosome, optionally strand).
    :param order_by:
        Window ORDER BY terms (start position).
    :return:
        A copy of the predicate with every ``PREV(...)`` call replaced by an
        equivalent LAG window and all column identifiers quoted.
    :raises ValueError:
        If a ``PREV()`` call does not take exactly one argument, or if a
        ``PREV()`` call is nested inside another (predicates compare only the
        immediate predecessor).
    """

    def _is_prev(node: exp.Expression) -> bool:
        return isinstance(node, exp.Anonymous) and node.name.upper() == "PREV"

    def _replace(node: exp.Expression) -> exp.Expression:
        if _is_prev(node):
            args = node.expressions
            if len(args) != 1:
                raise ValueError(
                    f"PREV() takes exactly one column argument; got {len(args)}."
                )
            if any(_is_prev(inner) for inner in args[0].find_all(exp.Anonymous)):
                raise ValueError(
                    "PREV() cannot be nested; a CLUSTER/MERGE predicate "
                    "compares only the immediate predecessor."
                )
            return exp.Window(
                this=exp.Anonymous(this="LAG", expressions=[args[0].copy()]),
                partition_by=[col.copy() for col in partition_cols],
                order=exp.Order(expressions=[term.copy() for term in order_by]),
            )
        return node

    rewritten = predicate.copy().transform(_replace)
    for column in rewritten.find_all(exp.Column):
        column.this.set("quoted", True)
    return rewritten


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
