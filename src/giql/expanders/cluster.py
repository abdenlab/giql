"""The CLUSTER operator expander (epic #137, issue #144).

CLUSTER assigns a cluster id to every interval, grouping each run of mutually
adjacent intervals (optionally within a maximum gap, strand-aware, and gated by a
pairwise predicate) under one id. It cannot be a single window function because it
needs a window over a window (``SUM`` over a running-edge-derived boundary flag),
so the expansion restructures the enclosing query into a two-level form::

    SELECT *, CLUSTER(interval) AS cluster_id FROM features

becomes::

    SELECT *,
           SUM(__giql_is_new_cluster) OVER (PARTITION BY chrom ORDER BY start)
               AS cluster_id
    FROM (
        SELECT *,
               CASE WHEN MAX(end) OVER (...) >= start THEN 0 ELSE 1 END
                    AS __giql_is_new_cluster
        FROM features
    ) AS __giql_lag_calc

The boundary flag keys off the running maximum end of the preceding rows (the
cluster's right edge so far), not the immediately preceding row's end, so a later
interval contained within an earlier wider one is not spuriously split (#214).

(The ``becomes::`` form is simplified for readability; the emitted SQL quotes
identifiers, appends ``NULLS LAST`` to the window ORDER BY and a ``ROWS BETWEEN
UNBOUNDED PRECEDING AND 1 PRECEDING`` frame, and — for a bare or qualified star
projection — emits the outer star as ``* EXCEPT (__giql_is_new_cluster)`` to hide
the synthesized flag from the output (#184, #185).)

This module is the AST-expansion replacement for the legacy
``ClusterTransformer``, which ran as a pre-pass transformer on the raw parsed AST
(in the former ``giql.transformer`` module, removed once every join rewrite became
registry-driven — #144, #169). It produces the same SQL; the existing CLUSTER
transpilation and bedtools oracle tests are the migration oracle.

Unlike the node-local expanders for DISTANCE / NEAREST / DISJOIN — whose result
replaces a single node — CLUSTER is a **whole-query rewrite**. It shares only the
*shape* of :func:`giql.expanders.nearest._fallback_form` — return the operator
node unchanged so the pass's ``node.replace`` is a no-op — but the *mechanism*
differs. ``nearest`` does an in-place child ``.replace`` because its LATERAL has a
parent to replace through; CLUSTER instead copies the enclosing
:class:`~sqlglot.expressions.Select`, restructures the copy, and *transplants* its
contents back onto the original SELECT (see
:func:`giql.expanders._genomic.transplant`). That root-
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

from sqlglot import exp

from giql.constants import CLUSTER_FLAG_COL
from giql.constants import CLUSTER_SIBLING_PREFIX
from giql.expander import ExpansionContext
from giql.expander import expand_operators
from giql.expander import register
from giql.expanders._genomic import GenomicColumns
from giql.expanders._genomic import extract_stranded
from giql.expanders._genomic import find_projected
from giql.expanders._genomic import genomic_columns
from giql.expanders._genomic import is_star_projection
from giql.expanders._genomic import reject_cluster_merge_mix
from giql.expanders._genomic import require_top_level_projection
from giql.expanders._genomic import transplant
from giql.expressions import GIQLCluster
from giql.targets import GenericTarget


@register(GenericTarget, GIQLCluster)
def expand_cluster(node: GIQLCluster, ctx: ExpansionContext) -> exp.Expression:
    """Expand a CLUSTER node by restructuring its enclosing SELECT in place.

    Registered for :class:`~giql.targets.GenericTarget`, so every target resolves
    to it through the registry's generic chain (CLUSTER emits identical SQL across
    targets). Locates the enclosing :class:`~sqlglot.expressions.Select`, derives
    the genomic columns from the FROM table, and rewrites that SELECT in place into
    the two-level ``__giql_lag_calc`` form. Returns the original CLUSTER node (now
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
    # Hiding the synthesized __giql_is_new_cluster flag from a surfacing star is the
    # safe default (#184; MERGE opts out — see expand_cluster_query). This runs for
    # every CLUSTER the pass expands, root or nested: harmless when the flag would
    # not surface, and it also stops a nested CLUSTER's flag leaking through an
    # enclosing star.
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
    # __giql_lag_calc subquery; the originals the pass collected are now unreachable, so
    # re-run the pass over the restructured SELECT to expand any sibling pass-3
    # operators (spatial predicates, DISTANCE) carried into it. Safe from
    # recursion: the CLUSTER node is already replaced by its SUM window. (#144 B1)
    #
    # expand_operators may return a new root (a registered statement finalizer can
    # wrap it), and its contract requires callers to use the return value rather
    # than assume in-place mutation, so reinstall a new root in place of `select`.
    # Today the branch is never taken here — the sole finalizer-registering operator
    # (a correlated NEAREST fallback) is already a plain join by this deepest-first
    # re-walk, so this nested run registers no finalizer and `result is select` in
    # practice — but honoring the contract keeps the seam future-proof. A NEAREST
    # fallback nested under this CLUSTER registered its finalizer in the *outer*
    # expand_operators run; that finalizer re-locates the transplanted join by its
    # reserved `meta` tag and wraps its enclosing SELECT after this rewrite, so a
    # `SELECT *`/`b.*` that would re-surface the fallback's reserved columns is
    # projected away rather than leaking (#172).
    result = expand_operators(select, ctx.target, ctx.tables, ctx.registry)
    if result is not select:
        select.replace(result)
    return node


def expand_cluster_query(
    query: exp.Select, columns: GenomicColumns, hide_reserved: bool = True
) -> exp.Select | None:
    """Restructure *query* for every CLUSTER in its projection; ``None`` if none.

    Finds CLUSTER expressions in *query*'s SELECT list and rewrites the query into
    the two-level ``__giql_lag_calc`` form once per CLUSTER (chaining, as the legacy
    transformer did). Operates on *query* directly and returns the rewritten
    query. Returns ``None`` when *query* projects no CLUSTER, so callers can treat
    that as a no-op.

    *hide_reserved* (default ``True``) rewrites a star in the outer projection to hide
    the synthesized ``__giql_is_new_cluster`` flag the star would otherwise surface
    (#184). Hiding is the default so a caller cannot silently reintroduce the leak by
    omitting it. MERGE — which reuses this helper to build its intermediate clustered
    query — is the one caller that opts out (``hide_reserved=False``): its explicit
    outer projection never surfaces the flag, so hiding it would only add a needless
    ``* EXCEPT`` to the intermediate ``__giql_clustered`` subquery.
    """
    cluster_exprs = find_projected(query, GIQLCluster)
    if not cluster_exprs:
        return None
    if len(cluster_exprs) > 1:
        # Chaining the rewrite per CLUSTER yields a duplicate ``__giql_lag_calc``
        # alias and an ``__giql_is_new_cluster`` binder error — non-executable SQL.
        # Fail loudly, mirroring the multiple-MERGE guard, rather than emitting the
        # broken SQL (#144 A15).
        raise ValueError("Multiple CLUSTER expressions not yet supported")
    if sum(is_star_projection(item) for item in query.expressions) > 1:
        # The two-level rewrite copies the whole projection into the inner
        # ``__giql_lag_calc`` subquery and re-expands each star at the outer level, so
        # N stars multiply the base columns (N inner expansions × N outer), silently
        # producing the wrong column set. A single star is the supported shape; fail
        # loudly on multiple, mirroring the multiple-CLUSTER guard above, rather than
        # emitting a silently multiplied projection (#194). (MERGE reuses this helper
        # with a synthesized single-star intermediate, so this never fires there.)
        raise ValueError(
            "CLUSTER does not support multiple star projections "
            "(e.g. SELECT *, *, CLUSTER(...)); project a single star"
        )
    for cluster_expr in cluster_exprs:
        query = _transform_for_cluster(query, cluster_expr, columns, hide_reserved)
    return query


def _transform_for_cluster(
    query: exp.Select,
    cluster_expr: GIQLCluster,
    columns: GenomicColumns,
    hide_reserved: bool = True,
) -> exp.Select:
    """Rewrite *query* into the two-level ``__giql_lag_calc`` form for one CLUSTER.

    A byte-for-byte port of the legacy
    ``ClusterTransformer._transform_for_cluster``, with the genomic columns passed
    in (the legacy method re-derived them from the FROM table) rather than read off
    ``self``. Builds an inner ``__giql_lag_calc`` subquery that materializes an
    ``__giql_is_new_cluster`` flag from a ``LAG`` window, then an outer query whose
    ``SUM(__giql_is_new_cluster) OVER (...)`` window replaces the CLUSTER projection.
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

    # Cluster boundaries key off the running maximum end of the preceding rows --
    # the cluster's right edge so far -- NOT the immediately preceding row's end.
    # LAG(end) would spuriously split a later interval that still overlaps a
    # cluster whose previous row is a contained interval ending early (#214).
    running_max_end = exp.Window(
        this=exp.Max(this=exp.column(end_col, quoted=True)),
        partition_by=partition_cols,
        order=exp.Order(expressions=order_by),
        spec=exp.WindowSpec(
            kind="ROWS",
            start="UNBOUNDED",
            start_side="PRECEDING",
            end="1",
            end_side="PRECEDING",
        ),
    )

    # Add distance offset if specified
    if distance > 0:
        edge_with_distance = exp.Add(
            this=running_max_end, expression=exp.Literal.number(distance)
        )
    else:
        edge_with_distance = running_max_end

    # Build the adjacency condition (running cluster edge >= current start).
    adjacency = exp.GTE(
        this=edge_with_distance,
        expression=exp.column(start_col, quoted=True),
    )

    # An optional predicate further restricts which adjacent intervals are kept
    # together: a row stays in the current cluster only when it is adjacent to the
    # cluster's running-max edge AND the predicate holds against its immediate
    # sorted predecessor. ``PREV(col)`` references in the predicate resolve to that
    # predecessor row via LAG over the same partition/order. Note the two notions
    # can reference different rows under containment -- adjacency comes from an
    # earlier, wider interval while the predicate compares to the immediate
    # predecessor -- matching the documented immediate-predecessor semantics (#214).
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

    # Create CASE expression for __giql_is_new_cluster
    case_expr = exp.Case(
        ifs=[
            exp.If(
                this=keep_together,
                true=exp.Literal.number(0),
            )
        ],
        default=exp.Literal.number(1),
    )

    # When the projection carries a star, every sibling item (each non-star, non-CLUSTER
    # projection alongside it) is re-projected at the outer level by *reference* to the
    # column the inner subquery materializes (see the outer loop below), pinned to its
    # written position. Referencing it by the user's own name is unsafe: the name can
    # collide with a base column the star also surfaces (``SELECT *, expr AS score`` —
    # the outer ``EXCEPT`` would strip the base ``score`` and the reference would bind
    # ambiguously) or with another sibling's name (``SELECT *, 1 AS x, 2 AS x`` — a
    # duplicate ``EXCEPT`` entry, a hard parse error). Materialize each sibling in the
    # inner subquery under a reserved, unique ``__giql_sibling_N`` name instead; the
    # outer star EXCEPTs those reserved names and the outer projection references them,
    # aliased back to the name they would carry un-clustered, at their written position —
    # so base columns are never stripped, the inner window references stay unambiguous,
    # and both the values *and the column order* match the un-clustered projection.
    # Materializing *every* sibling — not only aliased ones — is what pins column order:
    # a non-aliased sibling left to ride the star would always surface at the star's
    # position, silently reordering it ahead of any aliased sibling written before it
    # (#190; column-order pinning #192 round 3).
    outer_has_star = any(is_star_projection(item) for item in query.expressions)
    sibling_inner_names: dict[int, str] = {}
    if outer_has_star:
        for item in query.expressions:
            if is_star_projection(item):
                continue
            if isinstance(item, GIQLCluster):
                continue
            if isinstance(item, exp.Alias) and isinstance(item.this, GIQLCluster):
                continue
            sibling_inner_names[id(item)] = (
                f"{CLUSTER_SIBLING_PREFIX}{len(sibling_inner_names)}"
            )

    # Build CTE SELECT expressions (all original except CLUSTER, plus
    # __giql_is_new_cluster)
    cte_expressions = []
    for expression in query.expressions:
        # Skip CLUSTER expressions
        if isinstance(expression, GIQLCluster):
            continue
        elif isinstance(expression, exp.Alias) and isinstance(
            expression.this, GIQLCluster
        ):
            continue
        elif id(expression) in sibling_inner_names:
            # A star's sibling: materialize it under its reserved inner name so the outer
            # star can EXCEPT it (without stripping a same-named base column) and the
            # outer projection can re-project it at the user's written slot. Materialize
            # the aliased item's inner value, or the bare column / unnamed expression
            # itself when it has no alias.
            value = expression.this if isinstance(expression, exp.Alias) else expression
            cte_expressions.append(
                exp.alias_(
                    value.copy(),
                    sibling_inner_names[id(expression)],
                    quoted=False,
                )
            )
        else:
            cte_expressions.append(expression)

    # Ensure required columns for window functions are included
    required_cols = {chrom_col, start_col, end_col}
    if stranded:
        required_cols.add(strand_col)

    # The predicate is evaluated inside the __giql_lag_calc CTE, so every column
    # it references (current-row columns and PREV() arguments alike) must
    # be projected into that CTE. Folding them into required_cols makes the
    # scope dependency explicit and keeps the columns available even when a
    # later operator wraps this query in a further subquery.
    if predicate_expr is not None:
        required_cols |= {col.name for col in predicate_expr.find_all(exp.Column)}

    # Check if required columns are already in the select list
    selected_cols = set()
    for expr in cte_expressions:
        if is_star_projection(expr):
            # A bare ``*`` or a qualified ``rel.*`` covers every column. A qualified
            # ``rel.*`` parses as a Column wrapping a Star, so it must be caught here,
            # before the plain-Column branch below — whose ``.name`` for a star is the
            # literal ``'*'`` (never a required-column name), so it would leave chrom/
            # start/end seen as missing and re-inject them as chrom_1/start_1/end_1
            # duplicates (#185).
            selected_cols = required_cols  # Assume all are covered
            break
        elif isinstance(expr, exp.Column):
            selected_cols.add(expr.name)
        elif isinstance(expr, exp.Alias):
            # Don't count aliases as the source column
            pass

    # Add missing required columns
    # Sort the residual so the injected-column order is deterministic across runs
    # (set-difference iteration order is PYTHONHASHSEED-dependent; the legacy port
    # left it nondeterministic — see review #144 A2).
    for col in sorted(required_cols - selected_cols):
        cte_expressions.append(exp.column(col, quoted=True))

    # Add __giql_is_new_cluster calculation. The reserved __giql_ prefix keeps the
    # synthesized flag from colliding with a user column of the same name (#161).
    cte_expressions.append(exp.alias_(case_expr, CLUSTER_FLAG_COL, quoted=False))

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

    # Create outer query with SUM over __giql_is_new_cluster
    sum_window = exp.Window(
        this=exp.Sum(this=exp.column(CLUSTER_FLAG_COL)),
        partition_by=partition_cols,
        order=exp.Order(expressions=order_by),
    )

    # When the projection carries a star, every sibling was materialized in the inner
    # __giql_lag_calc subquery under a reserved ``__giql_sibling_N`` name (the CTE loop
    # above). The outer star EXCEPTs those reserved names so it does not surface them,
    # and each sibling is re-projected at its written position as a reference to its
    # reserved column, aliased back to the output name it would carry un-clustered.
    # Referencing the materialized column rather than re-evaluating the original
    # expression keeps it correct even when that expression's inputs were EXCEPTed from
    # the star inside the subquery, and — because the reserved name cannot collide —
    # keeps it correct for any name. Re-projecting *every* sibling (not only aliased
    # ones) pins each to the user's column order (#190; column-order pinning #192).
    #
    # Without a star nothing re-surfaces the item, so it is always kept as-is.
    star_hidden_names = list(sibling_inner_names.values())

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
        elif is_star_projection(expression):
            # A star projection (bare ``*`` or qualified ``rel.*``) expands over the
            # inner __giql_lag_calc subquery. It is normalized to a *bare* outer star
            # (dropping any dangling ``rel.`` qualifier — an executability fix), the
            # synthesized __giql_is_new_cluster flag is EXCEPTed from it when
            # hide_reserved is set (the default; MERGE opts out — see
            # expand_cluster_query), and each aliased sibling re-projected below is
            # EXCEPTed too so the star does not double-surface it (#190). See
            # _outer_star_over_lag_calc for how the concerns are gated.
            rewritten = _outer_star_over_lag_calc(
                expression, CLUSTER_FLAG_COL, hide_reserved, star_hidden_names
            )
            new_expressions.append(expression if rewritten is None else rewritten)
        elif not outer_has_star:
            # A non-star, non-CLUSTER item with no sibling star: keep it unchanged.
            new_expressions.append(expression)
        elif id(expression) in sibling_inner_names:
            # A sibling of a star: re-project the reserved column the CTE materialized it
            # under, aliased back to the output name it would carry un-clustered, at its
            # written position — so it surfaces exactly once (its reserved name is
            # EXCEPTed from the star above), collision-free for any name, in the
            # user's column order rather than jumping to the star's position. The output
            # name is the user's alias, a bare column's own name, or — for an unnamed
            # expression, which carries no output name — its rendered text, so a
            # positional consumer sees it in the written slot (#190; ordering #192).
            out_name = expression.output_name or expression.sql()
            new_expressions.append(
                exp.alias_(
                    exp.column(sibling_inner_names[id(expression)], quoted=True),
                    out_name,
                    quoted=True,
                )
            )

    # Build new query
    new_query = exp.Select()
    new_query.select(*new_expressions, copy=False)

    # Wrap CTE in subquery and set as FROM clause. The alias carries the reserved
    # __giql_ prefix for namespace consistency with the other synthesized names
    # (#161); a derived-table alias never actually collides with a user relation
    # (it lives in the enclosing query's scope), so this is hygiene, not a fix.
    subquery = exp.Subquery(
        this=cte_select,
        alias=exp.TableAlias(this=exp.Identifier(this="__giql_lag_calc")),
    )
    new_query.from_(subquery, copy=False)

    # Copy ORDER BY from original to outer query
    if query.args.get("order"):
        new_query.order_by(*query.args["order"].expressions, copy=False)

    return new_query


def _outer_star_over_lag_calc(
    expression: exp.Expression,
    reserved: str,
    hide_reserved: bool,
    extra_hidden: list[str] | None = None,
) -> exp.Star | None:
    """Rewrite a star projection for the outer query over ``__giql_lag_calc``.

    Returns a *bare* :class:`~sqlglot.expressions.Star` for a bare ``*`` or a qualified
    ``rel.*`` outer projection, or ``None`` for any non-star item. (The caller only ever
    invokes this on the star branch, so the ``None`` return is defensive — non-star
    items are handled by the caller, not passed here.) Three concerns are handled, on
    separate gates:

    **Normalization to a bare star — always.** The CLUSTER rewrite's outer query selects
    from the single renamed ``__giql_lag_calc`` subquery, so the user's ``rel.``
    qualifier no longer resolves there — an un-rewritten ``rel.*`` would dangle
    (``Binder Error: Referenced table "rel" not found``). CLUSTER runs over a single
    relation (joins are rejected upstream by
    :func:`giql.expanders._genomic.genomic_columns`), so every column ``rel.*`` would
    have named is exactly the subquery's projection, making ``rel.*`` equivalent to
    ``*``. The qualifier is therefore *always* dropped to a bare star. This is an
    **executability** fix (without it the query does not bind), so it must NOT be gated
    on the cosmetic *hide_reserved* flag (#185). Any ``EXCEPT`` the user's own star
    carried is already applied *inside* ``__giql_lag_calc`` (the inner subquery keeps
    the original projection node), so the outer bare star never needs to carry it.

    **Flag-hiding — only when *hide_reserved*.** The subquery's ``*`` also carries the
    synthesized ``reserved`` flag (``__giql_is_new_cluster``), so a bare outer ``*``
    would surface it — and, under a preserved DISTINCT, dedup over it. When
    *hide_reserved* is set (the default; MERGE opts out — its explicit projection never
    surfaces the flag), the bare star is emitted as ``* EXCEPT (reserved)`` to hide it
    (#184).

    The hiding is applied *eagerly*, at the same projection level as the outer star —
    NOT via an outer-wrapping statement finalizer like NEAREST's
    (:func:`giql.expanders.nearest._make_reserved_column_finalizer`). It must be:
    ``transplant`` re-attaches the user's ``DISTINCT`` onto this outer query, so a
    finalizer that wrapped it in ``SELECT * EXCEPT (...) FROM (<this>)`` would dedup
    over the flag *before* stripping it (splitting rows that should collapse). NEAREST
    can wrap because its reserved columns surface through an arbitrary enclosing SELECT
    it does not own; CLUSTER owns its outer star at build time and must excise the flag
    in place. A future consolidation (#187) must preserve this level constraint.

    **Aliased-sibling hiding — when *extra_hidden* is given.** When the projection pairs
    the star with aliased sibling items (``SELECT *, expr AS name, CLUSTER(...)``), each
    sibling was materialized in ``__giql_lag_calc`` under a reserved ``__giql_sibling_N``
    name and would be surfaced a second time by the outer star. The caller re-projects
    each such sibling (aliased back to the user's name) at its own position and passes
    the reserved names as *extra_hidden* so they are EXCEPTed from this star too — the
    star then carries only the base columns, and the sibling surfaces exactly once, in
    the user's column order. The reserved names are used precisely so this EXCEPT can
    never strip a same-named base column the star should keep (#190).

    ``* EXCEPT`` is the portable exclusion form GIQL already relies on (the coordinate
    canonicalizer and the NEAREST fallback finalizer emit it too), so sqlglot renders
    the target spelling (``EXCLUDE`` for DuckDB) without a capability branch here.
    """
    if not is_star_projection(expression):
        return None
    hidden = []
    if hide_reserved:
        hidden.append(exp.column(reserved, quoted=True))
    if extra_hidden:
        hidden.extend(exp.column(name, quoted=True) for name in extra_hidden)
    if hidden:
        return exp.Star(except_=hidden)
    return exp.Star()


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
