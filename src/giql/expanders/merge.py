"""The MERGE operator expander (epic #137, issue #144).

MERGE combines overlapping (and, with parameters, adjacent / strand-matched /
predicate-gated) intervals into single intervals. It is built on CLUSTER: assign
a cluster id, then aggregate ``MIN(start)`` / ``MAX(end)`` per cluster::

    SELECT MERGE(interval) FROM features

becomes::

    SELECT chrom, MIN(start) AS start, MAX(end) AS end
    FROM (SELECT *, CLUSTER(interval) AS __giql_cluster_id FROM features) AS clustered
    GROUP BY chrom, __giql_cluster_id
    ORDER BY chrom, start

(The ``becomes::`` form is simplified for readability; the emitted SQL quotes
identifiers, appends ``NULLS LAST``, and the inner ``CLUSTER(...)`` is itself
expanded into the two-level ``lag_calc`` form.)

This module is the AST-expansion replacement for the legacy
:class:`giql.transformer.MergeTransformer`; it produces the same SQL (the existing
MERGE transpilation and bedtools oracle tests are the migration oracle).

Like CLUSTER (:mod:`giql.expanders.cluster`), MERGE is a **whole-query rewrite**:
it navigates to the enclosing :class:`~sqlglot.expressions.Select`, mutates it in
place, and returns the operator node unchanged so the pass's ``node.replace`` is a
no-op. The intermediate clustered subquery is restructured through CLUSTER's shared
:func:`giql.expanders.cluster.expand_cluster_query`, mirroring how the legacy
``MergeTransformer`` composed ``ClusterTransformer``.
"""

from __future__ import annotations

from sqlglot import exp

from giql.expander import ExpansionContext
from giql.expander import expand_operators
from giql.expander import register

# The shared, operator-neutral expansion toolkit lives in the neutral ``_genomic``
# module (#163). ``expand_cluster_query`` is the one deliberate CLUSTER-reuse
# point — MERGE composes the CLUSTER rewrite to build its clustered subquery — so
# it is imported from the CLUSTER expander as a documented composition rather than
# incidental coupling.
from giql.expanders._genomic import GenomicColumns
from giql.expanders._genomic import extract_stranded
from giql.expanders._genomic import find_projected
from giql.expanders._genomic import genomic_columns
from giql.expanders._genomic import reject_cluster_merge_mix
from giql.expanders._genomic import require_top_level_projection
from giql.expanders._genomic import transplant
from giql.expanders.cluster import expand_cluster_query
from giql.expressions import GIQLCluster
from giql.expressions import GIQLMerge
from giql.targets import GenericTarget


@register(GenericTarget, GIQLMerge)
def expand_merge(node: GIQLMerge, ctx: ExpansionContext) -> exp.Expression:
    """Expand a MERGE node by restructuring its enclosing SELECT in place.

    Registered for :class:`~giql.targets.GenericTarget` (MERGE emits identical SQL
    across targets). Locates the enclosing
    :class:`~sqlglot.expressions.Select`, derives the genomic columns from the FROM
    table, and rewrites that SELECT in place into the clustered-aggregation form.
    Returns the original MERGE node (now unreachable from the rewritten root SELECT)
    so the pass's ``node.replace`` is a no-op.

    Parameters
    ----------
    node : GIQLMerge
        The MERGE node being expanded.
    ctx : ExpansionContext
        The expansion context; only its :attr:`~ExpansionContext.tables` is read.

    Returns
    -------
    exp.Expression
        The MERGE node, unchanged — the surrounding SELECT was mutated in place.
    """
    select = node.find_ancestor(exp.Select)
    if select is None:
        return node
    require_top_level_projection(select, node, GIQLMerge)
    reject_cluster_merge_mix(select)
    columns = genomic_columns(select, ctx.tables)
    # Build from a detached copy, then transplant onto the original SELECT to
    # preserve its identity (and rewrite a root SELECT without a parent to replace
    # through), exactly as the CLUSTER expander does.
    source = select.copy()
    transformed = _transform_select_merge(source, columns)
    if transformed is None:
        # No-op for a MERGE that parses outside the SELECT projection (e.g. in
        # WHERE / ORDER BY): find_projected finds none in the copy, so there is
        # nothing to restructure and the operator leaks to the generator exactly
        # as on `main`. require_top_level_projection above already rejects the
        # in-projection-expression case; this guards the out-of-projection case.
        return node
    transplant(select, transformed)
    # copy()+transplant duplicated the enclosing WHERE into the new clustered
    # subquery; the originals the pass collected are now unreachable, so re-run the
    # pass over the restructured SELECT to expand any sibling pass-3 operators
    # carried into it. Safe from recursion: the MERGE is already gone. (#144 B1)
    #
    # expand_operators may return a new root (a registered statement finalizer can
    # wrap it), and its contract requires callers to use the return value rather
    # than assume in-place mutation, so reinstall a new root in place of `select`.
    # Today the branch is never taken here — the sole finalizer-registering operator
    # (a correlated NEAREST fallback) is already a plain join by this deepest-first
    # re-walk, and MERGE's final projection is explicit, so a nested NEAREST fallback
    # never surfaces its reserved columns here anyway — but honoring the contract
    # keeps the seam future-proof.
    result = expand_operators(select, ctx.target, ctx.tables, ctx.registry)
    if result is not select:
        select.replace(result)
    return node


def _transform_select_merge(
    query: exp.Select, columns: GenomicColumns
) -> exp.Select | None:
    """Rewrite *query* for the MERGE in its projection; ``None`` if none.

    Mirrors the legacy ``MergeTransformer.transform`` dispatch: find the MERGE
    expressions, reject more than one, and delegate the single supported case to
    :func:`_transform_for_merge`.
    """
    merge_exprs = find_projected(query, GIQLMerge)
    if not merge_exprs:
        return None
    # For now, support only one MERGE expression
    if len(merge_exprs) > 1:
        raise ValueError("Multiple MERGE expressions not yet supported")
    return _transform_for_merge(query, merge_exprs[0], columns)


def _transform_for_merge(
    query: exp.Select, merge_expr: GIQLMerge, columns: GenomicColumns
) -> exp.Select:
    """Rewrite *query* into the clustered-aggregation form for one MERGE.

    A byte-for-byte port of the legacy ``MergeTransformer._transform_for_merge``,
    with the genomic columns passed in and the intermediate clustered query
    restructured through CLUSTER's shared
    :func:`giql.expanders.cluster.expand_cluster_query` (the legacy method called
    ``ClusterTransformer.transform``). Builds an inner ``clustered`` subquery that
    appends ``__giql_cluster_id``, then an outer query that aggregates
    ``MIN(start)`` / ``MAX(end)`` per cluster.
    """
    chrom_col, start_col, end_col, strand_col = columns

    # Extract MERGE parameters (same as CLUSTER)
    distance_expr = merge_expr.args.get("distance")
    stranded_expr = merge_expr.args.get("stranded")
    predicate_expr = merge_expr.args.get("predicate")

    # Build CLUSTER expression with same parameters
    cluster_kwargs = {"this": merge_expr.this}
    if distance_expr:
        cluster_kwargs["distance"] = distance_expr
    if stranded_expr:
        cluster_kwargs["stranded"] = stranded_expr
    if predicate_expr is not None:
        cluster_kwargs["predicate"] = predicate_expr

    cluster_expr = GIQLCluster(**cluster_kwargs)

    # Create intermediate query with CLUSTER
    # Start with original query's FROM/WHERE/etc
    cluster_query = exp.Select()
    cluster_query.select(exp.Star(), copy=False)
    # NOTE: __giql_cluster_id carries the reserved prefix; the un-prefixed
    # `clustered` / `lag_calc` / `is_new_cluster` siblings are left to #161.
    cluster_query.select(
        exp.alias_(cluster_expr, "__giql_cluster_id", quoted=False),
        append=True,
        copy=False,
    )

    # Copy FROM, WHERE from original
    # Use copy() to avoid sharing references between queries
    if query.args.get("from_"):
        cluster_query.set("from_", query.args["from_"].copy())
    if query.args.get("where"):
        cluster_query.set("where", query.args["where"].copy())

    # Apply CLUSTER transformation to get the CTE-based query
    clustered = expand_cluster_query(cluster_query, columns)
    assert clustered is not None, "intermediate MERGE cluster query has no CLUSTER"
    cluster_query = clustered

    # Build GROUP BY columns
    # Quote the chrom column (as every other chrom reference is) so a reserved-word
    # custom chrom column emits valid SQL (#144 A13).
    group_by_cols = [exp.column(chrom_col, quoted=True)]

    # Handle stranded parameter
    stranded = extract_stranded(stranded_expr)

    if stranded:
        group_by_cols.append(exp.column(strand_col, quoted=True))

    group_by_cols.append(exp.column("__giql_cluster_id"))

    # Build SELECT expressions for merged intervals
    select_exprs = []

    # Add group-by columns (non-aggregated)
    select_exprs.append(exp.column(chrom_col, quoted=True))
    if stranded:
        select_exprs.append(exp.column(strand_col, quoted=True))

    # Add merged interval bounds
    select_exprs.append(
        exp.alias_(
            exp.Min(this=exp.column(start_col, quoted=True)), start_col, quoted=False
        )
    )
    select_exprs.append(
        exp.alias_(exp.Max(this=exp.column(end_col, quoted=True)), end_col, quoted=False)
    )

    # Process other columns from original SELECT
    for expression in query.expressions:
        # Skip the MERGE expression itself
        if isinstance(expression, GIQLMerge):
            continue
        elif isinstance(expression, exp.Alias) and isinstance(
            expression.this, GIQLMerge
        ):
            continue
        # Include other columns (they should be aggregates or in GROUP BY)
        else:
            select_exprs.append(expression)

    # Build final query
    final_query = exp.Select()
    final_query.select(*select_exprs, copy=False)

    # FROM the clustered subquery
    subquery = exp.Subquery(
        this=cluster_query,
        alias=exp.TableAlias(this=exp.Identifier(this="clustered")),
    )
    final_query.from_(subquery, copy=False)

    # Add GROUP BY
    final_query.group_by(*group_by_cols, copy=False)

    # Add ORDER BY (chromosome, start)
    final_query.order_by(
        exp.Ordered(this=exp.column(chrom_col, quoted=True)), copy=False
    )
    final_query.order_by(
        exp.Ordered(this=exp.column(start_col, quoted=True)), append=True, copy=False
    )

    return final_query
