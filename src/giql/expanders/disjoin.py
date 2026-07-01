"""The DISJOIN operator expander (epic #137, wave 3 — issue #143).

DISJOIN splits each target interval at every reference breakpoint strictly
interior to it, passing the full target row through and appending the
sub-interval as ``disjoin_chrom`` / ``disjoin_start`` / ``disjoin_end`` in the
target table's declared coordinate system. A coverage filter drops sub-intervals
overlapping no reference interval; an omitted reference defaults to the target
set.

This module is the AST-expansion replacement for the legacy
``giqldisjoin_sql`` generator emitter. It assembles the same WITH-CTE
subquery, parses it back into a sqlglot :class:`~sqlglot.expressions.Expression`,
and returns that node so the active target's serializer renders it — dissolving
the emit-time string special-case left by epic #114.

Three behaviours are capability-driven / corrected here relative to the legacy
emitter:

* **star-REPLACE portability (#143).** The full-row passthrough de-canonicalizes
  the interval columns back into the target's declared encoding. On a target that
  supports ``SELECT * REPLACE (...)`` (``supports_star_replace`` — DuckDB) the
  passthrough emits ``t.* REPLACE (...)``. On a target without it (the generic
  baseline and the DataFusion family) it emits the portable
  ``t.* EXCEPT (...)`` form plus the two recomputed interval columns. This
  portable form is **not** SQL-92 — it requires ``SELECT * EXCEPT (...)``
  support, which the DataFusion family provides but a strict SQL-92 engine (and
  DuckDB, which prefers ``* REPLACE``) does not — so the generic non-canonical
  passthrough is **not DuckDB-runnable** and runs only on an ``* EXCEPT``-capable
  engine.
* **passthrough column order diverges across targets (issue #143).** The two
  passthrough forms are *row-equivalent* but **not column-order-equivalent**.
  ``* REPLACE`` substitutes the start/end columns *in place*, so the row keeps
  its original column order. ``* EXCEPT (start, end)`` drops those two columns
  and the recomputed ones are *re-appended* at the end of the projection, so on
  the generic / DataFusion path the start/end columns move to the right of any
  passthrough column. A ``SELECT *`` over a non-canonical DISJOIN is therefore
  column-order-divergent between DuckDB and the generic target; only the row
  *set* is identical. Select **explicit** columns (in a chosen order) when the
  output column order must agree across targets.
* **duplicate-column fix (#153).** Every UNION branch of the ``__giql_dj_cuts``
  CTE aliases all four projected columns (``kc`` / ``ks`` / ``ke`` / ``pos``), so
  the ``t."end"`` column and the end-cut expression never collide on an output
  name when ``end`` de-canonicalizes to the bare physical column (the default
  0-based half-open identity case). DuckDB tolerated the prior unaliased
  duplicate; DataFusion rejected it as a non-unique projection name.

Input canonicalization is owned by ``CanonicalizeCoordinates`` (pass 2, #122):
every non-canonical interval-bearing operand is rewritten to a canonical
``__giql_canon_*`` CTE before this pass runs, so the expander consumes
already-canonical 0-based half-open columns and applies no input-canonicalization
arithmetic. The output round-trip back to the target's declared encoding stays
here, driven by the original encoding pass 2 preserves on the resolution.
"""

from __future__ import annotations

from sqlglot import exp
from sqlglot import parse_one

from giql.canonical import decanonical_end
from giql.canonical import decanonical_start
from giql.constants import DJ_PREFIX
from giql.dialect import GIQLDialect
from giql.expander import ExpansionContext
from giql.expander import register
from giql.expressions import GIQLDisjoin
from giql.resolver import OperatorResolution
from giql.resolver import ResolvedRef
from giql.table import Table
from giql.targets import GenericTarget


@register(GenericTarget, GIQLDisjoin)
def expand_disjoin(node: GIQLDisjoin, ctx: ExpansionContext) -> exp.Expression:
    """Expand a DISJOIN node to its WITH-CTE subquery AST for the active target.

    Registered for :class:`~giql.targets.GenericTarget`, so it is the portable
    fallback every target resolves to through the registry's generic chain. The
    one capability branch — ``ctx.capabilities.supports_star_replace`` — selects
    the passthrough projection form, so a single expander covers both the
    ``* REPLACE`` (DuckDB) and ``* EXCEPT`` (DataFusion / generic) targets.

    Parameters
    ----------
    node : GIQLDisjoin
        The :class:`~giql.expressions.GIQLDisjoin` node being expanded.
    ctx : ExpansionContext
        The expansion context carrying the node's pass-1/pass-2 resolution
        metadata, the active target and its capabilities, and the registered
        tables.

    Returns
    -------
    exp.Expression
        The parsed AST of the parenthesized WITH-CTE subquery that replaces the
        DISJOIN node.
    """
    sql = _build_disjoin_sql(node, ctx)
    return parse_one(sql, dialect=GIQLDialect)


def _build_disjoin_sql(node: GIQLDisjoin, ctx: ExpansionContext) -> str:
    """Assemble the DISJOIN WITH-CTE subquery SQL string for the active target.

    Mirrors the legacy ``giqldisjoin_sql`` emitter's shape one named
    ``__giql_dj_*`` CTE at a time, with the #153 alias fix and the
    capability-driven passthrough applied.
    """
    resolution = ctx.resolution
    target_ref, ref, ref_from = _disjoin_resolution(node, resolution)
    target_name = target_ref.name
    target_chrom, target_start, target_end = target_ref.cols

    ref_chrom, ref_start, ref_end = ref.cols
    is_self_reference = ref.coverage_skippable

    # The target's *declared* encoding, which disjoin_* output and the
    # passed-through interval round-trip back into. Pass 2 preserves it on the
    # resolution when it wraps a non-canonical target (the slot's own Table is
    # then None); a canonical target keeps the (identity) encoding on its slot.
    output_table = _disjoin_output_encoding(resolution, target_ref)

    # Post-pass every operand is canonical 0-based half-open, so the physical
    # columns are consumed verbatim with no input-canonicalization arithmetic.
    t_chrom = f't."{target_chrom}"'
    t_start = f't."{target_start}"'
    t_end = f't."{target_end}"'

    # Reference endpoints: unqualified for the breakpoint CTE, qualified by 'r'
    # for the coverage EXISTS filter.
    bp_start = f'"{ref_start}"'
    bp_end = f'"{ref_end}"'
    r_start = f'r."{ref_start}"'
    r_end = f'r."{ref_end}"'

    # disjoin_start / disjoin_end are emitted in the target's declared coordinate
    # system so an output row carries one convention; the cut math stays canonical
    # internally.
    out_start = decanonical_start("s.seg_start", output_table)
    out_end = decanonical_end("s.seg_end", output_table)
    passthrough = _disjoin_passthrough(ctx, target_start, target_end, output_table)

    # Build the WITH clause one named fragment per __giql_dj_* CTE so each block
    # reads on its own. The `seg_end > seg_start` guard in the final WHERE is
    # belt-and-suspenders: UNION already dedupes cut positions, so LEAD cannot
    # produce a zero-length segment unless it becomes UNION ALL.
    ref_cte = f"__giql_dj_ref AS (SELECT * FROM {ref_from})"
    tgt_cte = f"__giql_dj_tgt AS (SELECT * FROM {target_name})"
    bp_cte = (
        "__giql_dj_bp AS ("
        f'SELECT "{ref_chrom}" AS chrom, {bp_start} AS pos FROM __giql_dj_ref '
        "UNION "
        f'SELECT "{ref_chrom}" AS chrom, {bp_end} AS pos FROM __giql_dj_ref)'
    )
    # #153: alias all four columns in EVERY UNION branch. The output names already
    # come from branch 1, so this is behaviour-preserving on DuckDB and makes each
    # branch's projection internally unique for strict engines (DataFusion), where
    # an unaliased t."end" would otherwise collide with the end-cut t."end".
    cuts_cte = (
        "__giql_dj_cuts AS ("
        f'SELECT t."{target_chrom}" AS kc, t."{target_start}" AS ks, '
        f't."{target_end}" AS ke, {t_start} AS pos FROM __giql_dj_tgt AS t '
        "UNION "
        f'SELECT t."{target_chrom}" AS kc, t."{target_start}" AS ks, '
        f't."{target_end}" AS ke, {t_end} AS pos FROM __giql_dj_tgt AS t '
        "UNION "
        f'SELECT t."{target_chrom}" AS kc, t."{target_start}" AS ks, '
        f't."{target_end}" AS ke, bp.pos AS pos '
        "FROM __giql_dj_tgt AS t JOIN __giql_dj_bp AS bp "
        f"ON bp.chrom = {t_chrom} AND bp.pos > {t_start} "
        f"AND bp.pos < {t_end})"
    )
    segs_cte = (
        "__giql_dj_segs AS ("
        "SELECT kc, ks, ke, pos AS seg_start, "
        "LEAD(pos) OVER (PARTITION BY kc, ks, ke ORDER BY pos) AS seg_end "
        "FROM __giql_dj_cuts)"
    )
    # In self-reference mode the coverage EXISTS is provably always true: every
    # emitted segment lies inside its parent target row, and that row is itself a
    # member of the reference set. Skip the clause so the planner does not waste
    # work on a no-op semi-join. The __giql_dj_ref CTE itself stays live because
    # __giql_dj_bp still draws breakpoints from it.
    where_clauses = ["s.seg_end IS NOT NULL", "s.seg_end > s.seg_start"]
    if not is_self_reference:
        where_clauses.append(
            f'EXISTS (SELECT 1 FROM __giql_dj_ref AS r WHERE r."{ref_chrom}" = s.kc '
            f"AND {r_start} <= s.seg_start AND {r_end} > s.seg_start)"
        )
    where_sql = " AND ".join(where_clauses)
    final_select = (
        f"SELECT {passthrough}, s.kc AS disjoin_chrom, "
        f"{out_start} AS disjoin_start, "
        f"{out_end} AS disjoin_end FROM __giql_dj_tgt AS t "
        f'JOIN __giql_dj_segs AS s ON t."{target_chrom}" = s.kc '
        f'AND t."{target_start}" = s.ks AND t."{target_end}" = s.ke '
        f"WHERE {where_sql}"
    )
    return (
        f"(WITH {ref_cte}, {tgt_cte}, {bp_cte}, "
        f"{cuts_cte}, {segs_cte} {final_select})"
    )


def _disjoin_passthrough(
    ctx: ExpansionContext,
    target_start: str,
    target_end: str,
    output_table: Table | None,
) -> str:
    """Project the target's full row, de-canonicalizing the interval columns.

    When the target's declared encoding is canonical 0-based half-open the row
    passes through as a plain ``t.*`` — the identity fast path, portable on every
    engine. When it is non-canonical the interval columns (canonical inside
    ``__giql_dj_tgt``) are de-canonicalized back into that encoding:

    * ``t.* REPLACE (...)`` on a target that supports it (``supports_star_replace``
      — DuckDB), substituting start/end **in place**;
    * the portable ``t.* EXCEPT (start, end), <start>, <end>`` form otherwise
      (the generic baseline / DataFusion family), which every ``* EXCEPT``-capable
      engine plans. This form is **not** SQL-92 and is **not DuckDB-runnable**: it
      requires ``SELECT * EXCEPT`` support.

    Column-order divergence (issue #143)
    ------------------------------------
    The two forms produce row-equivalent but **not column-order-equivalent**
    output. ``* REPLACE`` keeps the original column order; ``* EXCEPT`` drops
    start/end and **re-appends** the recomputed columns at the end of the
    projection, so they move to the right of any passthrough column on the
    generic / DataFusion path. A ``SELECT *`` over a non-canonical DISJOIN is
    therefore order-divergent across targets — only the row set agrees. Select
    explicit columns when the cross-target column order must match.
    """
    if output_table is None or (
        output_table.coordinate_system == "0based"
        and output_table.interval_type == "half_open"
    ):
        return "t.*"
    pt_start = decanonical_start(f't."{target_start}"', output_table)
    pt_end = decanonical_end(f't."{target_end}"', output_table)
    if ctx.capabilities.supports_star_replace:
        return (
            f't.* REPLACE ({pt_start} AS "{target_start}", {pt_end} AS "{target_end}")'
        )
    # Portable substitution: drop the two interval columns from the star and
    # re-project them recomputed. EXCEPT removes them from the row, the trailing
    # projections add them back in the target's encoding under their own names.
    return (
        f't.* EXCEPT ("{target_start}", "{target_end}"), '
        f'{pt_start} AS "{target_start}", {pt_end} AS "{target_end}"'
    )


def _disjoin_output_encoding(
    resolution: OperatorResolution | None, target_ref: ResolvedRef
) -> Table | None:
    """Return the target's declared encoding for DISJOIN's output round-trip.

    ``CanonicalizeCoordinates`` (pass 2) records the original
    :class:`~giql.table.Table` on the resolution when it wraps a non-canonical
    target (blanking the slot's own ``table``). For an unwrapped target — a
    canonical registered table, or any target when the pass did not run — the
    slot's own ``table`` carries the (identity) encoding.
    """
    if isinstance(resolution, OperatorResolution):
        preserved = resolution.output_tables.get("this")
        if preserved is not None:
            return preserved
    return target_ref.table


def _disjoin_resolution(
    node: GIQLDisjoin, resolution: OperatorResolution | None
) -> tuple[ResolvedRef, ResolvedRef, str]:
    """Unpack the DISJOIN resolution attached by ResolveOperatorRefs (pass 1).

    Returns ``(target_ref, ref, ref_from)`` where ``ref_from`` is the text
    following ``FROM`` inside the ``__giql_dj_ref`` CTE. A subquery reference
    carries no name, so it is rendered from the AST node as an aliased derived
    table; registered tables and CTEs are selected from by name.

    The resolver pass deliberately leaves unresolvable slots unresolved; for
    those, and for a resolved name using the reserved ``__giql_dj_`` prefix, this
    re-raises DISJOIN's historical diagnostics verbatim.
    """
    target_ref = (
        resolution.slot("this")
        if isinstance(resolution, OperatorResolution)
        else None
    )

    # An unresolved target means it is not a registered table.
    if target_ref is None:
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

    # The __giql_dj_ prefix names the operator's internal CTEs; a target table
    # using it would collide with them.
    if target_ref.name.startswith(DJ_PREFIX):
        raise ValueError(
            f"DISJOIN target {target_ref.name!r} uses the reserved "
            "'__giql_dj_' prefix, which names the operator's internal "
            "CTEs. Rename the table."
        )

    reference = node.args.get("reference")
    ref = resolution.slot("reference")
    if ref is not None:
        if ref.kind == "subquery":
            # Serializing the subquery with sqlglot's stock generator is safe
            # because nested GIQL operators expand deepest-first: any GIQL
            # construct inside this operand is already
            # rewritten to standard AST by the time this outer DISJOIN expands,
            # so the stock generator only ever round-trips plain SQL here.
            ref_sql = reference.sql(dialect=GIQLDialect)
            return target_ref, ref, f"{ref_sql} AS __giql_dj_rs"
        return target_ref, ref, ref.name

    # Unresolved reference: re-classify it and raise the matching historical
    # diagnostic. An omitted reference always resolves, so reference is non-None.
    if not isinstance(reference, (exp.Table, exp.Column, exp.Identifier)):
        raise ValueError(
            "DISJOIN reference must be a table name, a CTE, or a "
            f"(SELECT ...) subquery; got {type(reference).__name__}: "
            f"{reference}"
        )
    ref_name = reference.name
    if ref_name.startswith(DJ_PREFIX):
        raise ValueError(
            f"DISJOIN reference {ref_name!r} uses the reserved "
            "'__giql_dj_' prefix, which names the operator's internal "
            "CTEs. Rename the reference relation."
        )
    raise ValueError(
        f"DISJOIN reference {ref_name!r} is neither a registered table "
        "nor a CTE defined in this query. Register the table or define "
        "the CTE before transpiling."
    )
