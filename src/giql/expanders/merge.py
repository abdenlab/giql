"""The MERGE operator expander (epic #137, issue #144).

MERGE combines overlapping (and, with parameters, adjacent / strand-matched /
predicate-gated) intervals into single intervals. It is built on CLUSTER: assign
a cluster id, then aggregate ``MIN(start)`` / ``MAX(end)`` per cluster::

    SELECT MERGE(interval) FROM features

becomes::

    SELECT chrom, MIN(start) AS start, MAX(end) AS end
    FROM (SELECT *, CLUSTER(interval) AS __giql_cluster_id FROM features)
        AS __giql_clustered
    GROUP BY chrom, __giql_cluster_id
    ORDER BY chrom, start

(The ``becomes::`` form is simplified for readability; the emitted SQL quotes
identifiers, appends ``NULLS LAST``, and the inner ``CLUSTER(...)`` is itself
expanded into the two-level ``__giql_lag_calc`` form.)

This module is the AST-expansion replacement for the legacy ``MergeTransformer``
(in the former ``giql.transformer`` module, removed once every join rewrite became
registry-driven — #144, #169); it produces the same SQL (the existing MERGE
transpilation and bedtools oracle tests are the migration oracle).

Like CLUSTER (:mod:`giql.expanders.cluster`), MERGE is a **whole-query rewrite**:
it navigates to the enclosing :class:`~sqlglot.expressions.Select`, mutates it in
place, and returns the operator node unchanged so the pass's ``node.replace`` is a
no-op. The intermediate clustered subquery is restructured through CLUSTER's shared
:func:`giql.expanders.cluster.expand_cluster_query`, mirroring how the legacy
``MergeTransformer`` composed ``ClusterTransformer``.
"""

from __future__ import annotations

from sqlglot import exp

from giql.constants import CLUSTER_ID_COL
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
from giql.expanders._genomic import is_star_projection
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
    _reject_star_projection(select)
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
    # copy()+transplant duplicated the enclosing WHERE into the new __giql_clustered
    # subquery; the originals the pass collected are now unreachable, so re-run the
    # pass over the restructured SELECT to expand any sibling pass-3 operators
    # carried into it. Safe from recursion: the MERGE is already gone. (#144 B1)
    #
    # expand_operators may return a new root (a registered statement finalizer can
    # wrap it), and its contract requires callers to use the return value rather
    # than assume in-place mutation, so reinstall a new root in place of `select`.
    # Today the branch is never taken here — the sole finalizer-registering operator
    # (a correlated NEAREST fallback) is already a plain join by this deepest-first
    # re-walk, so this nested run registers no finalizer — but honoring the contract
    # keeps the seam future-proof. A NEAREST fallback nested under this MERGE
    # registered its finalizer in the outer expand_operators run; that finalizer
    # re-locates the transplanted join by its reserved `meta` tag and would wrap a
    # surfacing star (#172). MERGE's final projection is explicit, so it never
    # surfaces the reserved columns — they are covered both by that (no-op) wrapper
    # check and by the explicit projection.
    result = expand_operators(select, ctx.target, ctx.tables, ctx.registry)
    if result is not select:
        select.replace(result)
    return node


def _reject_star_projection(select: exp.Select) -> None:
    """Raise if *select* projects a star alongside its MERGE. MERGE-specific guard.

    MERGE is an aggregate whole-query rewrite: its output is one row per merged
    region (``chrom``, ``MIN(start)``, ``MAX(end)``), so a star projected in the
    same SELECT has no coherent meaning — it names the *pre-aggregation* input
    columns of a relation the rewrite has already collapsed into ``__giql_clustered``
    and grouped away. Left unguarded, a bare ``*`` re-surfaces those non-grouped,
    non-aggregated columns under the synthesized ``GROUP BY`` and a qualified
    ``rel.*`` dangles a table alias the outer aggregation no longer exposes — both
    non-executable SQL (#189). Fail loudly here, mirroring the sibling MERGE guards
    (multiple-MERGE, CLUSTER/MERGE mix), so the shape raises a clear diagnostic
    instead.

    This is a genuinely MERGE-specific guard — unlike CLUSTER, a per-row window over
    which a star is meaningful and supported (#184, #185), MERGE collapses rows, so
    no star projection alongside it is well-formed — hence it lives here rather than
    in the operator-neutral toolkit. The check is gated on there being a *projected*
    MERGE (``find_projected``) so a MERGE parsing outside the projection (WHERE /
    ORDER BY) is left on the expander's existing out-of-projection no-op path
    rather than rejected. Legitimate sibling items — extra aggregates such as
    ``COUNT(*)`` (a ``Star`` wrapped in an aggregate, not a bare or qualified star)
    — are untouched; only a star is rejected.
    """
    if not find_projected(select, GIQLMerge):
        return
    if any(is_star_projection(projection) for projection in select.expressions):
        raise ValueError(
            "MERGE cannot be combined with a star projection (e.g. SELECT *, "
            "MERGE(...) or SELECT rel.*, MERGE(...)); MERGE aggregates rows into one "
            "row per merged region, so a star has no coherent per-row meaning. Drop "
            "the star, or project only grouping columns and aggregates (e.g. "
            "COUNT(*)) alongside MERGE — not raw input columns, which are neither "
            "grouped nor aggregated."
        )


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
    ``ClusterTransformer.transform``). Builds an inner ``__giql_clustered`` subquery
    that appends ``__giql_cluster_id``, then an outer query that aggregates
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
    # __giql_cluster_id and its ``__giql_clustered`` / ``__giql_lag_calc`` /
    # ``__giql_is_new_cluster`` siblings all carry the reserved prefix (#161).
    cluster_query.select(
        exp.alias_(cluster_expr, CLUSTER_ID_COL, quoted=False),
        append=True,
        copy=False,
    )

    # Copy FROM, WHERE from original
    # Use copy() to avoid sharing references between queries
    if query.args.get("from_"):
        cluster_query.set("from_", query.args["from_"].copy())
    if query.args.get("where"):
        cluster_query.set("where", query.args["where"].copy())

    # Apply CLUSTER transformation to get the CTE-based query. hide_reserved=False:
    # MERGE opts out of CLUSTER's flag-hiding (#184) because its explicit outer
    # projection (chrom, MIN/MAX) never surfaces __giql_is_new_cluster, so hiding it
    # would only add a needless ``* EXCEPT`` to this intermediate clustered subquery.
    clustered = expand_cluster_query(cluster_query, columns, hide_reserved=False)
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

    group_by_cols.append(exp.column(CLUSTER_ID_COL))

    # The columns the synthesized grouping exposes to the user's projection and GROUP
    # BY: chrom, plus strand when stranded (the per-merge __giql_cluster_id is internal
    # and never user-referenceable). Used below to reconcile the user's own projection
    # and GROUP BY against MERGE's aggregation (#192).
    grouping_keys = {chrom_col}
    if stranded:
        grouping_keys.add(strand_col)
    # Every identifier comparison below folds case. SQL binds unquoted identifiers
    # case-insensitively, so ``AS Start`` collides with the synthesized ``start`` and
    # ``GROUP BY CHROM`` names the ``chrom`` grouping key. Folding both sides matches a
    # case-variant rather than silently missing it — the safe direction, since a missed
    # collision is silently-wrong output while an over-match merely fails loud (#192).
    grouping_keys_folded = {key.lower() for key in grouping_keys}

    # Reconcile the user's GROUP BY with MERGE's synthesized grouping. MERGE groups by
    # chrom (and strand when stranded) plus its per-merge __giql_cluster_id, so a user
    # GROUP BY over exactly those grouping-key columns is subsumed by it (the documented
    # "merge by chromosome" shape). Anything else — another column, or an expression /
    # ordinal MERGE cannot map to its grouping — was previously dropped silently,
    # changing the result, so reject it loudly. Each item must be a *bare column*
    # reference to a grouping key: an expression (GROUP BY UPPER(chrom)) or an ordinal
    # (GROUP BY 1) is not an exp.Column, so it cannot be verified to name a grouping key
    # and is rejected rather than honored-then-discarded (#192).
    user_group = query.args.get("group")
    if user_group is not None:
        for item in user_group.expressions:
            if (
                not isinstance(item, exp.Column)
                or item.name.lower() not in grouping_keys_folded
            ):
                raise ValueError(
                    f"MERGE cannot honor GROUP BY {item.sql()}: MERGE groups rows by "
                    "chrom (and strand when stranded) and its synthesized per-merge "
                    "cluster id, so group only by a bare chrom (or strand) column "
                    "reference alongside MERGE — not another column, an expression, or "
                    "an ordinal, which it cannot faithfully honor."
                )

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

    # The output-column names MERGE synthesizes: the grouping keys plus the aggregated
    # ``start`` / ``end`` bounds. A user item may not silently collide with these. Folded
    # for case-insensitive comparison, as grouping_keys_folded above.
    synthesized_names_folded = {
        name.lower() for name in (grouping_keys | {start_col, end_col})
    }

    # Process the other projection items. MERGE already projects the grouping keys and
    # the MIN/MAX bounds and aggregates rows into one row per merged region, so each
    # remaining item is reconciled against that grouping rather than copied verbatim —
    # copying verbatim duplicated a grouping key the user also projected and, for a raw
    # non-grouped column, emitted SQL the engine rejects (#192).
    for expression in query.expressions:
        # Skip the MERGE expression itself
        if isinstance(expression, GIQLMerge):
            continue
        elif isinstance(expression, exp.Alias) and isinstance(
            expression.this, GIQLMerge
        ):
            continue

        # A projected item is well-formed under MERGE's per-cluster GROUP BY only if
        # every column it references *outside a grouping aggregate* is a grouping key.
        # A column counts as ungrouped when no aggregate wraps it (a raw ``score``, or
        # ``start`` in ``MAX(score) + start``) OR when the aggregate that wraps it is a
        # *window* aggregate (``SUM(score) OVER (...)``): a window is evaluated after the
        # GROUP BY and does not collapse the group, so its argument still has to be
        # grouped. Fail loudly rather than emitting SQL the engine rejects (or that a
        # lenient engine runs with an arbitrary-row value) (#192).
        ungrouped = {
            col.name.lower()
            for col in expression.find_all(exp.Column)
            if col.find_ancestor(exp.AggFunc) is None
            or col.find_ancestor(exp.Window) is not None
        }
        if not ungrouped <= grouping_keys_folded:
            raise ValueError(
                "MERGE cannot project the non-aggregated column "
                f"{', '.join(repr(c) for c in sorted(ungrouped - grouping_keys_folded))}"
                f" in '{expression.sql()}': MERGE aggregates rows into one row per "
                "merged region, so project only grouping columns (chrom, and "
                "strand when stranded) and plain aggregates (e.g. COUNT(*)) "
                "alongside MERGE."
            )

        out_name = expression.output_name
        inner = expression.this if isinstance(expression, exp.Alias) else expression

        # A grouping-key column MERGE already projects, under its own name (bare
        # ``chrom`` or a no-op ``chrom AS chrom`` of the "merge by chromosome" shape):
        # drop the duplicate rather than surfacing the column twice (#192).
        if (
            out_name.lower() in grouping_keys_folded
            and isinstance(inner, exp.Column)
            and inner.name.lower() == out_name.lower()
        ):
            continue

        # Reject an item whose *output name* collides with a synthesized column (the
        # aggregated ``start`` / ``end`` bounds, or a value aliased onto a grouping key
        # such as ``chrom AS start``): emitting two columns of that name silently
        # duplicates one, and for ``start`` / ``end`` hijacks the default
        # ``ORDER BY "chrom", "start"``, returning rows in the wrong order. Apply this
        # only to an explicit alias or a bare column, whose output name the engine emits
        # verbatim; for any other unaliased expression (``CAST(chrom AS VARCHAR)``,
        # ``chrom || ''``), ``output_name`` is a sqlglot heuristic that does NOT match
        # the engine's generated label, so checking it there would falsely reject a
        # well-formed projection over a grouping key. Fold case so an unquoted ``Start``
        # is caught (it binds to the synthesized ``start``) (#192).
        names_reliable = isinstance(expression, (exp.Alias, exp.Column))
        if names_reliable and out_name.lower() in synthesized_names_folded:
            raise ValueError(
                f"MERGE cannot project a column named {out_name!r} alongside MERGE: it "
                "collides with the chrom/start/end columns MERGE synthesizes for each "
                "merged region — alias it to a different name."
            )

        # An aggregate over the merged rows (COUNT(*), AVG(score), ...) or a
        # non-aggregate expression over only grouping-key columns (UPPER(chrom), a
        # rename ``chrom AS c``): keep it — well-formed under the GROUP BY, no collision.
        select_exprs.append(expression)

    # Build final query
    final_query = exp.Select()
    final_query.select(*select_exprs, copy=False)

    # FROM the clustered subquery. The alias carries the reserved __giql_ prefix for
    # namespace consistency with the other synthesized names (#161); a derived-table
    # alias never actually collides with a user relation (it lives in the enclosing
    # query's scope), so this is hygiene, not a fix.
    subquery = exp.Subquery(
        this=cluster_query,
        alias=exp.TableAlias(this=exp.Identifier(this="__giql_clustered")),
    )
    final_query.from_(subquery, copy=False)

    # Add GROUP BY
    final_query.group_by(*group_by_cols, copy=False)

    # Honor the user's HAVING over the merged output. MERGE aggregates rows into one row
    # per merged region, so a HAVING filters those regions — ``HAVING COUNT(*) > 1``
    # keeps only regions built from more than one input interval, resolving against the
    # per-region GROUP BY above. The rewrite builds ``final_query`` fresh and transplant
    # clears the original HAVING (a rebuilt root arg), so it must be re-attached here;
    # omitting it silently dropped the clause and returned unfiltered rows (#192). An
    # invalid HAVING (over a raw non-grouped column) surfaces a loud engine binder error,
    # exactly as an unsupported ORDER BY does — never silently-wrong output.
    user_having = query.args.get("having")
    if user_having is not None:
        final_query.set("having", user_having.copy())

    # ORDER BY: honor the user's ORDER BY over the merged output when present,
    # otherwise default to (chromosome, start) for deterministic output. MERGE
    # aggregates rows away, so the user's ORDER BY resolves against the merged
    # columns (e.g. ``ORDER BY "end" DESC`` over ``MAX("end") AS end``). Overriding
    # it with the fixed default silently returned the wrong rows once an outer LIMIT
    # was preserved (#181), so the user's ordering takes precedence.
    if query.args.get("order"):
        final_query.set("order", query.args["order"].copy())
    else:
        final_query.order_by(
            exp.Ordered(this=exp.column(chrom_col, quoted=True)), copy=False
        )
        final_query.order_by(
            exp.Ordered(this=exp.column(start_col, quoted=True)),
            append=True,
            copy=False,
        )

    return final_query
