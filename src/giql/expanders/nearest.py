"""The NEAREST operator expander (epic #137, issue #142).

NEAREST is the first operator whose expansion is genuinely capability-driven.
The portable form is a correlated ``LATERAL`` subquery: each outer row drives a
``SELECT ... FROM <target> WHERE <chrom prefilter> ORDER BY ABS(distance)
LIMIT k`` whose reference endpoints are outer-table columns. DuckDB and the
generic target plan that directly (``supports_lateral == True``).

Apache DataFusion has no correlated-``LATERAL`` physical plan
(``supports_lateral == False``). For it the same k-nearest / ``max_distance`` /
``stranded`` / ``signed`` semantics are reproduced with a **decorrelated
window-function fallback**: the target is cross-joined against the outer
relation, each candidate is ranked with
``ROW_NUMBER() OVER (PARTITION BY <reference key> ORDER BY ABS(distance))``, and
the surrounding ``CROSS JOIN LATERAL`` is rewritten into a plain join that
re-associates the top-``k`` ranked candidates back to every outer row sharing
that reference key. Ranking depends only on the reference value, so ranking once
per distinct reference value and re-joining is set-equivalent to the per-row
LATERAL form — deterministic when the ``(start, end)`` tiebreaker distinguishes
candidates tied at the k-th distance (verified by the cross-target result
oracle).

A literal-reference (standalone) NEAREST is already an uncorrelated subquery, so
every target — DataFusion included — uses the LATERAL/standalone form unchanged;
only the *correlated* shape needs the fallback.

The expander assembles its distance/passthrough/filter SQL as string fragments
(the ``_nearest_*`` helpers below and
:func:`giql.expanders._distance.generate_distance_case`, shared in spirit with
DISTANCE, #140), then parses the assembled fragments into AST so the emitted SQL
is reserialized by the active target's serializer.
"""

from __future__ import annotations

from sqlglot import exp
from sqlglot import parse_one

from giql.canonical import decanonical_end
from giql.canonical import decanonical_start
from giql.dialect import GIQLDialect
from giql.expander import EXPAND_ALIAS_PREFIX
from giql.expander import ExpansionContext
from giql.expander import register
from giql.expanders._distance import generate_distance_case
from giql.expanders._params import coerce_bool_param
from giql.expressions import GIQLNearest
from giql.range_parser import RangeParser
from giql.resolver import META_KEY
from giql.resolver import OperatorResolution
from giql.resolver import ResolvedInterval
from giql.resolver import ResolvedRef
from giql.table import Table
from giql.targets import Capabilities
from giql.targets import GenericTarget

#: Reserved column names the window-function fallback synthesizes inside its
#: ranked subquery. They are derived from the expander's reserved
#: ``EXPAND_ALIAS_PREFIX`` (``__giql_x_``) — rather than hardcoded — so they stay
#: clear of user identifiers and track the prefix if it ever changes, mirroring
#: the other reserved internal prefixes.
_RANK_COL = f"{EXPAND_ALIAS_PREFIX}rn"
_REF_KEY_PREFIX = f"{EXPAND_ALIAS_PREFIX}rk_"


def _nearest_output_encoding(
    expression: GIQLNearest, target_ref: ResolvedRef
) -> Table | None:
    """Return the target's declared encoding for NEAREST's row passthrough.

    ``CanonicalizeCoordinates`` (pass 2) records the original
    :class:`~giql.table.Table` on the resolution when it wraps a non-canonical
    target in a ``__giql_canon_*`` CTE (blanking the slot's own ``table``). For
    an unwrapped target — a canonical registered table, or any target when the
    pass did not run — the slot's own ``table`` carries the (identity) encoding.

    :param expression:
        GIQLNearest expression node
    :param target_ref:
        The resolved target reference (post pass 2)
    :return:
        The target's declared :class:`~giql.table.Table`, or ``None``
    """
    resolution = expression.meta.get(META_KEY)
    if isinstance(resolution, OperatorResolution):
        preserved = resolution.output_tables.get("this")
        if preserved is not None:
            return preserved
    return target_ref.table


def _nearest_passthrough(
    table_name: str,
    target_start: str,
    target_end: str,
    output_table: Table | None,
    capabilities: Capabilities,
) -> str:
    """Project the target's full row, de-canonicalizing the interval columns.

    NEAREST passes the whole target row through (``SELECT {table_name}.*``)
    alongside the synthesized, encoding-invariant ``distance`` column. When the
    target's declared encoding is canonical 0-based half-open the row passes
    through as a plain ``{table_name}.*`` — the byte-identical identity fast
    path. When it is non-canonical the interval columns, canonical inside the
    ``__giql_canon_*`` CTE the target was rewritten to, are de-canonicalized
    back into that encoding so the passed-through interval matches the target's
    own convention.

    The emit strategy is chosen from *capabilities*, following the precedent of
    :func:`giql.expanders.disjoin._disjoin_passthrough` (issue #145) — the same
    two emit forms:

    * ``{table_name}.* REPLACE (...)`` when ``supports_star_replace`` holds —
      substitutes start/end in place;
    * the portable ``{table_name}.* EXCEPT (start, end), <start>, <end>`` form
      otherwise (the generic baseline / DataFusion family). Row-equivalent but
      not column-order-equivalent, and not DuckDB-runnable.

    :param table_name:
        The relation the row is selected from (the canon CTE name when wrapped,
        else the registered table name) — also the column qualifier.
    :param target_start:
        Physical start column name
    :param target_end:
        Physical end column name
    :param output_table:
        The target's declared :class:`~giql.table.Table`, or ``None``
    :param capabilities:
        The active target's :class:`~giql.targets.Capabilities`.
    :return:
        The passthrough projection fragment
    """
    if output_table is None or (
        output_table.coordinate_system == "0based"
        and output_table.interval_type == "half_open"
    ):
        return f"{table_name}.*"
    pt_start = decanonical_start(f'{table_name}."{target_start}"', output_table)
    pt_end = decanonical_end(f'{table_name}."{target_end}"', output_table)
    if capabilities.supports_star_replace:
        return (
            f"{table_name}.* REPLACE "
            f'({pt_start} AS "{target_start}", {pt_end} AS "{target_end}")'
        )
    # Portable form for engines without ``* REPLACE`` (generic / DataFusion):
    # drop the interval columns from the star and re-project them recomputed.
    return (
        f'{table_name}.* EXCEPT ("{target_start}", "{target_end}"), '
        f'{pt_start} AS "{target_start}", {pt_end} AS "{target_end}"'
    )


def _raise_nearest_reference_error(
    expression: GIQLNearest,
    resolution: OperatorResolution | None,
) -> None:
    """Raise the historical diagnostic for an unresolved NEAREST reference.

    ResolveOperatorRefs (pass 1) defers a reference slot it cannot resolve; this
    re-raises the historical pre-pass error verbatim. The implicit-outer failures
    rely on the :class:`~giql.resolver.SlotDeferral` the pass records (the
    ancestor walk that distinguished them now lives in the pass); a literal-range
    parse failure is reproduced by re-parsing.

    :param expression:
        GIQLNearest expression node
    :param resolution:
        The attached resolution metadata, if any
    :raises ValueError:
        Always — with the matching historical message
    """
    reference = expression.args.get("reference")

    if reference is None:
        # An absent reference is a correlated (implicit-outer) placement that the
        # resolver could not tie to a registered outer table; consult the
        # recorded deferral for the specific historical message.
        deferral = (
            resolution.deferral("reference") if resolution is not None else None
        )
        if deferral is not None and deferral.reason == "implicit_outer_unregistered":
            raise ValueError(
                f"Outer table '{deferral.detail}' not found in tables. "
                "Please specify reference parameter explicitly."
            )
        raise ValueError(
            "Could not find outer table in LATERAL join context. "
            "Please specify reference parameter explicitly."
        )

    # An explicit reference that deferred is a literal range that failed to
    # parse (column references always resolve). Re-parse to surface the
    # original parse error in the historical message.
    reference_sql = reference.sql(dialect=GIQLDialect)
    range_str = reference_sql.strip("'\"")
    try:
        RangeParser.parse(range_str).to_zero_based_half_open()
    except Exception as e:
        raise ValueError(
            f"Could not parse reference genomic range: {range_str}. Error: {e}"
        )
    raise ValueError(f"Could not parse reference genomic range: {range_str}.")


def _nearest_params(
    expression: GIQLNearest,
) -> tuple[int, int | None, bool, bool]:
    """Unpack the (k, max_distance, stranded, signed) parameters of a NEAREST."""
    k = expression.args.get("k")
    k_value = int(str(k)) if k else 1

    max_distance = expression.args.get("max_distance")
    max_dist_value = int(str(max_distance)) if max_distance else None

    is_stranded = coerce_bool_param(expression.args.get("stranded"))
    is_signed = coerce_bool_param(expression.args.get("signed"))
    return k_value, max_dist_value, is_stranded, is_signed


def _distance_and_filters(
    expression: GIQLNearest,
    table_name: str,
    target_ref: ResolvedRef,
    ref: ResolvedInterval,
    capabilities: Capabilities,
    ref_fragments: tuple[str, str, str, str | None] | None = None,
) -> tuple[str, str, list[str], str]:
    """Build the shared distance SQL, the qualified target columns, and WHERE.

    Returns ``(distance_expr, abs_distance_expr, where_clauses, passthrough)`` —
    the fragments common to the LATERAL/standalone form and the decorrelated
    fallback. Distance math, the chromosome pre-filter, the optional strand match,
    and the optional ``max_distance`` filter all reproduce the legacy
    ``giqlnearest_sql`` emitter exactly. Each form derives its deterministic
    ORDER BY tiebreaker from the target columns itself.

    ``capabilities`` is the active target's :class:`~giql.targets.Capabilities`,
    forwarded to :func:`_nearest_passthrough` to choose the target's
    de-canonicalization emit form (``* REPLACE`` vs the portable ``* EXCEPT``);
    both call sites pass ``ctx.capabilities``.

    ``ref_fragments`` optionally overrides the reference ``(chrom, start, end,
    strand)`` SQL fragments. The LATERAL form consumes the resolution's
    outer-qualified fragments verbatim; the fallback passes fragments pointing at
    its renamed, pre-projected reference relation so the cross-joined columns
    carry names distinct from the target's (DataFusion's planner cannot resolve a
    window ordering over a join with duplicate column names).
    """
    target_chrom, target_start, target_end = target_ref.cols
    _k_value, max_dist_value, is_stranded, is_signed = _nearest_params(expression)

    output_table = _nearest_output_encoding(expression, target_ref)
    passthrough = _nearest_passthrough(
        table_name, target_start, target_end, output_table, capabilities
    )

    if ref_fragments is not None:
        ref_chrom, ref_start, ref_end, ref_strand_frag = ref_fragments
    else:
        ref_chrom, ref_start, ref_end, ref_strand_frag = (
            ref.chrom,
            ref.start,
            ref.end,
            ref.strand,
        )

    ref_strand = None
    target_strand = None
    if is_stranded:
        ref_strand = ref_strand_frag
        if output_table and output_table.strand_col:
            target_strand = f'{table_name}."{output_table.strand_col}"'

    target_chrom_expr = f'{table_name}."{target_chrom}"'
    target_start_expr = f'{table_name}."{target_start}"'
    target_end_expr = f'{table_name}."{target_end}"'

    distance_expr = generate_distance_case(
        ref_chrom,
        ref_start,
        ref_end,
        ref_strand,
        target_chrom_expr,
        target_start_expr,
        target_end_expr,
        target_strand,
        stranded=is_stranded,
        signed=is_signed,
    )
    abs_distance_expr = f"ABS({distance_expr})"

    where_clauses = [f"{ref_chrom} = {target_chrom_expr}"]
    if is_stranded and ref_strand and target_strand:
        where_clauses.append(f"{ref_strand} = {target_strand}")
    if max_dist_value is not None:
        where_clauses.append(f"({abs_distance_expr}) <= {max_dist_value}")

    return distance_expr, abs_distance_expr, where_clauses, passthrough


def _lateral_form(
    expression: GIQLNearest,
    ctx: ExpansionContext,
    table_name: str,
    target_ref: ResolvedRef,
    ref: ResolvedInterval,
) -> exp.Expression:
    """The portable LATERAL/standalone subquery.

    Builds a two-level subquery: an inner ``SELECT <passthrough>, <distance> AS
    distance FROM <target> WHERE ...`` that materializes the distance, wrapped by
    an outer ``SELECT * FROM (<inner>) AS x ORDER BY ABS(x.distance), x.<start>,
    x.<end> LIMIT k`` that orders on the *precomputed* ``distance`` column. For a
    correlated placement the parent ``LATERAL`` correlates it to the outer row;
    for a standalone (literal-reference) placement it stands alone.

    Splitting the distance computation (inner) from the ordering (outer) is
    load-bearing for cross-engine support:

    * DuckDB's correlated-``LATERAL`` binder will not resolve a SELECT-list alias
      named ``distance`` from inside an ``ORDER BY`` that also projects
      ``<target>.*``, so the order key must reference a *materialized* column
      (``x.distance``) from the wrapping level rather than an alias in the same
      SELECT.
    * DataFusion's planner, given the distance ``CASE`` re-emitted inline in the
      ``ORDER BY`` over the chromosome-equality prefiltered scan, rewrites the
      filtered ``chrom`` to a self-comparison in one copy of the CASE but not the
      other and trips ``SanityCheckPlan``; ordering on the materialized column
      avoids re-deriving the key.

    A ``(start, end)`` tiebreaker follows ``ABS(distance)`` so rows tied at the
    k-th distance order deterministically — and identically across engines and
    against the decorrelated fallback's ranking — *when ``(start, end)``
    distinguishes the tied candidates*. Two target rows sharing both distance and
    ``(start, end)`` remain order-ambiguous (no key here breaks that residual
    tie); the two forms are set-equivalent up to such coordinate-duplicate ties
    (#142 A5).
    """
    k_value, *_ = _nearest_params(expression)
    (
        distance_expr,
        _abs_distance_expr,
        where_clauses,
        passthrough,
    ) = _distance_and_filters(
        expression, table_name, target_ref, ref, ctx.capabilities
    )
    where_sql = " AND ".join(where_clauses)
    # The wrapping level reads the inner row's *bare* column names (the passthrough
    # projected ``<target>.*``), so the tiebreaker qualifies them by the wrapper
    # alias, not the original ``table_name."col"``.
    _chrom, target_start_col, target_end_col = target_ref.cols
    wrapper = ctx.alias()
    inner = (
        f"SELECT {passthrough}, {distance_expr} AS distance "
        f"FROM {table_name} WHERE {where_sql}"
    )
    sql = (
        f"(SELECT * FROM ({inner}) AS {wrapper} "
        f'ORDER BY ABS({wrapper}."distance"), '
        f'{wrapper}."{target_start_col}", {wrapper}."{target_end_col}" '
        f"LIMIT {k_value})"
    )
    return parse_one(sql, dialect=GIQLDialect)


def _outer_relation(ref: ResolvedInterval) -> tuple[str, str]:
    """Return ``(physical_relation, alias)`` for the correlated reference table.

    The reference endpoints are alias-qualified fragments (``a."chrom"``). The
    alias is the outer table's correlation name in the query; the physical
    relation comes from the reference's backing :class:`~giql.table.Table`. Both
    are needed to re-introduce the outer relation inside the decorrelated
    subquery the fallback builds.

    The ``else alias`` branch (``ref.table is None``) only fires for a reference
    whose backing table the resolver could not attach. For a *correlated*
    NEAREST — the only shape that reaches the fallback — the reference is an
    outer-row column, so the resolver always attaches its table and this branch
    is not reached in practice. Falling back to ``alias`` (yielding a
    cosmetically redundant ``FROM <alias> AS <alias>``, where ``relation ==
    alias``) keeps the emitted SQL valid rather than emitting an empty relation
    name should that invariant ever not hold; the assert pins the expectation.
    """
    parsed = parse_one(ref.chrom, dialect=GIQLDialect)
    alias = parsed.table if isinstance(parsed, exp.Column) else ""
    if ref.table is not None:
        relation = ref.table.name
    else:
        # Defensive fallback only: a correlated reference always carries a
        # resolved backing table, so relation == alias here would be a redundant
        # self-alias rather than a real two-name relation.
        assert alias, (
            "correlated NEAREST fallback expected the reference to carry either a "
            "backing table or an alias-qualified column"
        )
        relation = alias
    return relation, alias


def _fallback_form(
    expression: GIQLNearest,
    ctx: ExpansionContext,
    table_name: str,
    target_ref: ResolvedRef,
    ref: ResolvedInterval,
) -> exp.Expression:
    """The decorrelated window-function fallback for non-LATERAL targets.

    Rewrites the surrounding ``<outer> AS a CROSS JOIN LATERAL (nearest) AS b``
    into ``<outer> AS a JOIN (<ranked subquery>) AS b ON <ref-key match> AND
    b.<rn> <= k``. The ranked subquery cross-joins the target against the outer
    relation and ranks candidates per distinct reference key with
    ``ROW_NUMBER()``; the join re-associates the top-k back to every outer row
    sharing that key, reproducing the per-row LATERAL semantics. Swaps the parent
    ``LATERAL`` for the decorrelated subquery in place and returns the (now
    detached) NEAREST node, so the pass's own ``node.replace`` is a no-op.

    The no-op return relies on NEAREST having no nestable inner GIQL operator: a
    detached node carrying a still-pending descendant would strand that
    descendant's later ``node.replace``. NEAREST's only operands are a registered
    target table and an interval reference, neither of which is an expandable
    operator, so nothing pending hangs off the node this detaches.
    """
    lateral = expression.parent
    # Internal invariants the surrounding-AST rewrite depends on. The fallback
    # only runs for a correlated NEAREST, whose pass-1 placement is always a CROSS
    # JOIN LATERAL under a JOIN; a violation is an internal pipeline bug, not user
    # error, so fail loudly with a clear message rather than dereferencing None.
    # (The LATERAL *alias* is NOT an internal invariant — it is optional user
    # input — so it is synthesized below rather than asserted; see B3.)
    assert isinstance(lateral, exp.Lateral), (
        "correlated NEAREST fallback expected its parent to be a LATERAL, got "
        f"{type(lateral).__name__}"
    )
    join = lateral.parent
    assert isinstance(join, exp.Join), (
        "correlated NEAREST fallback expected the LATERAL to sit under a JOIN, got "
        f"{type(join).__name__}"
    )
    # The decorrelated join references the LATERAL's alias on both sides (its ON
    # clause and the replacement subquery's name). A correlated NEAREST written
    # *without* a LATERAL alias is legitimate user input that transpiles fine on
    # lateral-capable engines, so the fallback must not require one: synthesize a
    # collision-safe alias from the run's sequence when the LATERAL carries none.
    # This closes the DuckDB/DataFusion behavior gap (and, unlike an ``assert``,
    # survives ``python -O``, which would otherwise strip the guard and leave a
    # ``NoneType`` deref). The synthesized name is internal to the rewritten join.
    lateral_alias = lateral.args.get("alias")
    if lateral_alias is None or not lateral_alias.name:
        alias = ctx.alias()
        alias_node = exp.TableAlias(this=exp.to_identifier(alias))
    else:
        alias = lateral_alias.name
        alias_node = lateral_alias.copy()

    relation, outer_alias = _outer_relation(ref)
    k_value, _max, is_stranded, _signed = _nearest_params(expression)
    # Bare target column names: the candidate subquery exposes the target row via
    # ``target.*``, so its tiebreaker columns are referenced by name, not by the
    # ``table_name."col"`` qualifier the distance math uses.
    _target_chrom, target_start_col, target_end_col = target_ref.cols

    # Pre-project the outer relation's reference columns under fresh, non-target
    # names into a renamed derived relation. Cross-joining *this* (rather than the
    # raw outer table) keeps every reference column distinct from the target's
    # columns: DataFusion's planner cannot resolve a window ordering over a join
    # whose two sides share column names (e.g. both expose ``start`` / ``end``).
    #
    # The reference key identifies one distinct reference interval, which the
    # ranking partitions by and the join re-associates on. Position
    # (chrom/start/end) alone identifies it in the unstranded case; in stranded
    # mode strand joins the key too, because two outer rows at the same position
    # but opposite strands must each get their own strand-filtered nearest. The
    # ref relation is de-duplicated on the key with DISTINCT so ranking happens
    # once per distinct reference and the join fans the top-k back out to every
    # outer row sharing it — exactly the per-row LATERAL semantics, even when the
    # outer table holds duplicate reference rows.
    # Mint the synthetic relation aliases from the run's collision-safe sequence
    # (rather than hardcoding ``__giql_x_ref`` / ``__giql_x_cand``) so two NEAREST
    # operators in one query never reuse the same derived-relation name. The
    # reserved *column* names below stay derived from EXPAND_ALIAS_PREFIX.
    ref_relation_alias = ctx.alias()
    candidate = ctx.alias()
    strand_name = f"{_REF_KEY_PREFIX}strand"
    stranded_key = is_stranded and ref.strand is not None

    key_names = [f"{_REF_KEY_PREFIX}chrom", f"{_REF_KEY_PREFIX}start",
                 f"{_REF_KEY_PREFIX}end"]
    source_frags = [ref.chrom, ref.start, ref.end]
    if stranded_key:
        key_names.append(strand_name)
        source_frags.append(ref.strand)

    ref_projection = ", ".join(
        f'{frag} AS "{name}"' for name, frag in zip(key_names, source_frags)
    )
    ref_relation = (
        f"(SELECT DISTINCT {ref_projection} FROM {relation} AS {outer_alias})"
        f" AS {ref_relation_alias}"
    )

    # Reference fragments now point at the renamed relation's safe columns.
    renamed = [f'{ref_relation_alias}."{name}"' for name in key_names]
    renamed_strand = (
        f'{ref_relation_alias}."{strand_name}"' if stranded_key else None
    )
    ref_fragments = (renamed[0], renamed[1], renamed[2], renamed_strand)

    (
        distance_expr,
        _abs_distance_expr,
        where_clauses,
        passthrough,
    ) = _distance_and_filters(
        expression,
        table_name,
        target_ref,
        ref,
        ctx.capabilities,
        ref_fragments=ref_fragments,
    )

    # Surface the reference-key columns so the rewritten join can match each
    # ranked candidate back to its outer row(s). Ranking depends only on these
    # values, so partitioning by them and re-joining is identical to the per-row
    # LATERAL form even when outer rows share a reference value.
    key_cols = list(zip(key_names, renamed))
    key_projection = ", ".join(f'{frag} AS "{name}"' for name, frag in key_cols)
    where_sql = " AND ".join(where_clauses)

    # Compute the candidate set (cross join + distance + reference keys) in an
    # inner subquery, then add ROW_NUMBER() in the enclosing one. Keeping the
    # join and the window in *separate* query levels is load-bearing on
    # DataFusion: fused into one level its optimizer mis-derives the window's sort
    # order from the chromosome-equality prefilter and trips ``SanityCheckPlan``.
    inner = (
        f"SELECT {passthrough}, {distance_expr} AS distance, {key_projection} "
        f"FROM {table_name} CROSS JOIN {ref_relation} "
        f"WHERE {where_sql}"
    )
    partition = ", ".join(f'{candidate}."{name}"' for name, _ in key_cols)
    # A ``(start, end)`` tiebreaker follows ``ABS(distance)`` so rows tied at the
    # k-th distance rank identically here and in the LATERAL form when
    # ``(start, end)`` distinguishes them, making the two set-equivalent up to ties
    # among candidates sharing both distance and coordinates (no engine-dependent
    # tie ordering otherwise).
    ranked = (
        f"(SELECT {candidate}.*, "
        f"ROW_NUMBER() OVER (PARTITION BY {partition} "
        f"ORDER BY ABS({candidate}.distance), "
        f'{candidate}."{target_start_col}", {candidate}."{target_end_col}") '
        f'AS "{_RANK_COL}" '
        f"FROM ({inner}) AS {candidate})"
    )
    ranked_subquery = parse_one(ranked, dialect=GIQLDialect)

    # Match each ranked candidate back to the *outer* relation by its reference
    # value (the original outer-qualified fragments, e.g. ``a."chrom"``), not the
    # renamed inner columns which exist only inside the subquery.
    on_parts = [
        f'{alias}."{name}" = {src}' for name, src in zip(key_names, source_frags)
    ]
    on_parts.append(f'{alias}."{_RANK_COL}" <= {k_value}')
    on_sql = " AND ".join(on_parts)
    on_expr = parse_one(on_sql, dialect=GIQLDialect)

    subquery = exp.Subquery(this=ranked_subquery.this, alias=alias_node)

    # Convert ``<outer> CROSS JOIN LATERAL (nearest) AS b`` into
    # ``<outer> JOIN (ranked) AS b ON <ref-key match> AND b.rn <= k``. Swap the
    # whole LATERAL out for the decorrelated subquery and drop the CROSS kind so
    # the ON clause attaches as a plain (inner) join.
    lateral.replace(subquery)
    join.set("kind", None)
    join.set("side", None)
    join.set("on", on_expr)

    # The LATERAL (and the NEAREST node within it) is now detached; returning the
    # node unchanged makes the pass's ``node.replace`` a no-op.
    return expression


def expand_nearest(node: exp.Expression, ctx: ExpansionContext) -> exp.Expression:
    """Expand a NEAREST node to LATERAL or the decorrelated window-function form.

    Selects on ``ctx.capabilities.supports_lateral`` and whether the node is
    correlated (its parent is a ``LATERAL``). Lateral-capable targets and every
    standalone (literal-reference) placement get the portable LATERAL/standalone
    subquery; a correlated NEAREST on a target without LATERAL support gets the
    decorrelated window-function fallback.
    """
    assert isinstance(node, GIQLNearest)
    resolution = ctx.resolution

    target_ref = resolution.slot("this") if resolution is not None else None
    if not isinstance(target_ref, ResolvedRef):
        # An unresolved target means it is not a registered table; raise the
        # historical diagnostic (verbatim from the removed giqlnearest_sql).
        target = node.this
        if isinstance(target, exp.Table):
            target_name = target.name
        elif isinstance(target, exp.Column):
            target_name = target.table if target.table else str(target.this)
        else:
            target_name = str(target)
        raise ValueError(
            f"Target table '{target_name}' not found in tables. "
            "Register the table before transpiling."
        )
    table_name = target_ref.name

    ref = resolution.slot("reference")
    if not isinstance(ref, ResolvedInterval):
        _raise_nearest_reference_error(node, resolution)

    # A literal-range reference is uncorrelated even under CROSS JOIN LATERAL: its
    # endpoints are constants, not outer-row columns, so the subquery stands alone
    # and every target — DataFusion included — takes the LATERAL/standalone form.
    # Only a genuinely correlated reference (a column / implicit-outer endpoint)
    # needs the decorrelated window-function fallback on a lateral-incapable
    # target. Gating on parentage alone would mis-route a literal range into
    # ``_fallback_form``, which dereferences a non-existent outer relation.
    correlated = isinstance(node.parent, exp.Lateral) and ref.kind != "literal_range"
    if correlated and not ctx.capabilities.supports_lateral:
        return _fallback_form(node, ctx, table_name, target_ref, ref)
    return _lateral_form(node, ctx, table_name, target_ref, ref)


# The generic registration covers every target through the registry's fallback
# chain; the expander branches on ctx.capabilities.supports_lateral internally,
# so no per-target override is needed.
register(GenericTarget, GIQLNearest)(expand_nearest)
