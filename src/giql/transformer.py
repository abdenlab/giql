"""Query transformers for GIQL operations.

This module contains the DuckDB IEJoin pre-pass transformer, which rewrites a
column-to-column INTERSECTS join into a per-chromosome dynamic-SQL pattern DuckDB
plans through its range-join family. Every other target — and every shape the
IEJoin path declines — leaves the INTERSECTS in place for the pass-3
``(GenericTarget, Intersects)`` expander, which renders the naive overlap
predicate (a plain ``ON`` condition the engine plans as a range join). The generic
binned equi-join transformer was dropped in favor of that naive predicate (#167).
CLUSTER and MERGE were relocated to the operator-expander registry
(``giql.expanders.cluster`` / ``giql.expanders.merge``) in epic #137 (#144).
"""

from dataclasses import dataclass
from typing import ClassVar
from typing import Literal
from uuid import uuid4

from sqlglot import exp

from giql.canonical import canonical_end
from giql.canonical import canonical_start
from giql.constants import DEFAULT_CHROM_COL
from giql.constants import DEFAULT_END_COL
from giql.constants import DEFAULT_START_COL
from giql.expressions import Intersects
from giql.table import Table
from giql.table import Tables


class _UnqualifiedProjectionError(Exception):
    """Signal an unrecoverable projection-resolution failure in the DuckDB IEJoin path.

    Raised by :class:`IntersectsDuckDBIEJoinTransformer` when a query engages
    the dialect path but contains a SELECT-list shape the dialect cannot
    translate (bare ``*``, unqualified column, expression-form projection,
    unknown table qualifier, ``a.* AS x``, window aggregate, ``FILTER`` clause,
    arithmetic over an aggregate, or a scalar subquery in the projection list),
    or an extras predicate that references an unqualified or unknown-aliased
    column. Caught and re-raised as :class:`ValueError` at the dialect's public
    boundary.
    """


# These predicates support :class:`IntersectsDuckDBIEJoinTransformer`; they live
# at module scope rather than on the class so a future target-specific INTERSECTS
# override can reuse them without importing the transformer. Style §2's strict
# "functions after classes" ordering is relaxed here under the topic-interleaved
# exception.


def _is_column_intersects(node: Intersects) -> bool:
    """Return True if *node* is column-to-column :class:`Intersects`.

    Both operands must be table-qualified column references.
    """
    return (
        isinstance(node.this, exp.Column)
        and bool(node.this.table)
        and isinstance(node.expression, exp.Column)
        and bool(node.expression.table)
    )


def _find_column_intersects_in(expr: exp.Expression) -> Intersects | None:
    """Return the first column-to-column :class:`Intersects` in *expr*, or ``None``."""
    return next(
        (n for n in expr.find_all(Intersects) if _is_column_intersects(n)),
        None,
    )


def _normalize_alias(name: str, *, quoted: bool = False) -> str:
    """Return *name* case-folded unless it was a quoted identifier.

    DuckDB (and standard SQL) treat unquoted identifiers case-insensitively.
    Comparing aliases via raw Python equality rejects valid mixed-case queries
    that the naive-predicate plan accepts; normalize via :py:meth:`str.casefold` so the
    dialect matches DuckDB's identifier semantics. Quoted identifiers stay
    case-sensitive per the SQL standard.
    """
    return name if quoted else name.casefold()


def _has_outer_join_intersects(query: exp.Select) -> bool:
    """Return True if any outer JOIN has a column-to-column :class:`Intersects`.

    Handles two surfaces the user can write the INTERSECTS under an outer
    join: directly in the join's ``ON`` clause, or in the top-level
    ``WHERE`` while the join keeps its ``LEFT``/``RIGHT``/``FULL`` side
    modifier (``LEFT JOIN ... ON TRUE WHERE a.interval INTERSECTS
    b.interval``). The dialect's inner-join IEJoin rewrite would silently
    drop the outer-join's unmatched-row guarantee in the second shape,
    so we fall back whenever an outer-join is present and an
    INTERSECTS sits anywhere it could influence the join result.
    """
    if not isinstance(query, exp.Select):
        return False
    outer_joins = [
        join for join in query.args.get("joins") or [] if join.args.get("side")
    ]
    if not outer_joins:
        return False
    for join in outer_joins:
        on = join.args.get("on")
        if on is not None and _find_column_intersects_in(on):
            return True
    where = query.args.get("where")
    if where and _find_column_intersects_in(where.this):
        return True
    return False


def _strip_intersects(
    expr: exp.Expression | None, intersects: Intersects
) -> exp.Expression | None:
    """Return the parts of an AND tree that aren't the given :class:`Intersects`."""
    if expr is None or expr is intersects:
        return None
    if isinstance(expr, exp.And):
        if expr.this is intersects:
            return expr.expression
        if expr.expression is intersects:
            return expr.this
        left = _strip_intersects(expr.this, intersects)
        right = _strip_intersects(expr.expression, intersects)
        if left is None:
            return right
        if right is None:
            return left
        return exp.And(this=left, expression=right)
    return expr


def _count_column_intersects(query: exp.Select) -> int:
    """Count column-to-column :class:`Intersects` predicates anywhere in *query*.

    Walks the full AST so the dialect's "exactly one INTERSECTS" gate cannot
    be circumvented by an INTERSECTS hidden inside a subquery (e.g. inside
    a SELECT-list scalar subquery or inside ``EXISTS (SELECT ...)``). The
    SELECT-list / modifier-clause subquery gates in ``transform_to_sql``
    already reject those shapes, but this widened count is defense-in-depth
    against future gate relaxations.
    """
    if not isinstance(query, exp.Select):
        return 0
    return sum(1 for n in query.find_all(Intersects) if _is_column_intersects(n))


@dataclass(frozen=True)
class _IEJoinSides:
    """Resolved two-sided join shape for the DuckDB IEJoin dialect path.

    Bundles the left/right user aliases (normalized via :func:`_normalize_alias`)
    with the un-aliased :class:`~sqlglot.expressions.Table` AST nodes and a
    ``kind`` discriminator identifying the join variant. The hardcoded inner
    aliases ``a`` and ``b`` used inside each per-chromosome subquery are
    exposed as :data:`left_inner_alias` and :data:`right_inner_alias`.

    The :attr:`is_left_only_join` property is True for SEMI and ANTI shapes,
    which restrict the outer SELECT to left-side projections only.
    """

    left_user_alias: str
    right_user_alias: str
    left_table: exp.Table
    right_table: exp.Table
    kind: Literal["INNER", "SEMI", "ANTI"]

    left_inner_alias: ClassVar[str] = "a"
    right_inner_alias: ClassVar[str] = "b"

    @property
    def is_left_only_join(self) -> bool:
        """True when the join's output only references the left side (SEMI / ANTI)."""
        return self.kind in ("SEMI", "ANTI")

    @classmethod
    def from_intersects(
        cls,
        from_table: exp.Table,
        join: exp.Join,
        intersects: Intersects,
    ) -> "_IEJoinSides | None":
        """Resolve sides, run the FROM-side orientation swap, normalize aliases.

        Returns ``None`` when the join's table operand isn't a base
        :class:`~sqlglot.expressions.Table`, when neither side of the
        :class:`Intersects` predicate references the FROM table's alias, or
        when the join's ``kind`` is outside the supported set
        (``INNER`` / ``CROSS`` / ``SEMI`` / ``ANTI``; ``CROSS`` folds to
        ``INNER`` because per-chromosome partitioning makes them equivalent
        inside the inner subquery).
        """
        join_table = join.this
        if not isinstance(join_table, exp.Table):
            return None
        from_alias = _normalize_alias(from_table.alias_or_name)
        left_alias = _normalize_alias(intersects.this.table)
        right_alias = _normalize_alias(intersects.expression.table)
        if left_alias == from_alias:
            l_user_alias = from_alias
            l_table = from_table
            r_user_alias = right_alias
            r_table = join_table
        elif right_alias == from_alias:
            # Inner subquery template hardcodes "a" to the FROM side; swap so
            # that contract holds regardless of which operand the user put
            # first in INTERSECTS.
            l_user_alias = from_alias
            l_table = from_table
            r_user_alias = left_alias
            r_table = join_table
        else:
            return None
        raw_kind = (join.args.get("kind") or "").upper()
        if raw_kind in ("", "CROSS"):
            kind: Literal["INNER", "SEMI", "ANTI"] = "INNER"
        elif raw_kind in ("SEMI", "ANTI"):
            kind = raw_kind  # type: ignore[assignment]
        else:
            return None
        return cls(
            left_user_alias=l_user_alias,
            right_user_alias=r_user_alias,
            left_table=l_table,
            right_table=r_table,
            kind=kind,
        )


class IntersectsDuckDBIEJoinTransformer:
    """Transform column-to-column INTERSECTS joins into a DuckDB IEJoin pattern.

    Emits a two-step dynamic SQL output:

    1. ``SET VARIABLE __giql_iejoin_<token> = COALESCE(<dynamic-sql>, <empty-schema>);``
       where ``<dynamic-sql>`` is a ``string_agg`` over the chromosome
       partition that emits one ``SELECT ... FROM (... WHERE chrom = '<c>')
       a JOIN (... WHERE chrom = '<c>') b ON <range-overlap>`` per
       chromosome, joined by ``UNION ALL``.
    2. ``SELECT <projections> FROM query(getvariable('__giql_iejoin_<token>'))``

    Each per-chromosome subquery uses pure inequality predicates which DuckDB
    plans through its range-join family (``IE_JOIN`` or ``PIECEWISE_MERGE_JOIN``).
    The dialect supports INNER, SEMI, and ANTI variants — the per-chromosome
    join keyword switches accordingly, and SEMI / ANTI restrict the outer
    SELECT to left-side projections only.

    See :doc:`/transpilation/performance` for the full join-shape support
    matrix, fallback enumeration, and hard-error projection shapes.

    Notes:

    - ``a.*`` / ``b.*`` star projections expand to the genomic columns
      declared by the corresponding :class:`~giql.table.Table` config
      (chrom / start / end, plus strand when set). This narrows the result
      schema relative to the generic naive-predicate plan, which delegates star expansion
      to DuckDB and projects every base-table column. Users with additional
      non-genomic columns should list them explicitly.
    - ``SEMI JOIN`` and ``ANTI JOIN`` outputs are left-side only. Any
      ``b.col`` / ``b.*`` reference in the outer SELECT, aggregate
      arguments, GROUP BY, HAVING, or ORDER BY raises
      :class:`_UnqualifiedProjectionError` (wrapped to :class:`ValueError`
      at the public boundary).
    - ``ANTI JOIN`` partitions on the left table's distinct chromosomes
      only (not the chromosome INTERSECT used by INNER / SEMI) so that
      left rows on chromosomes absent from the right table are preserved.
    """

    def __init__(self, tables: Tables):
        """Initialize transformer.

        :param tables:
            Table configurations for column mapping.
        """
        self.tables = tables

    def transform_to_sql(self, query: exp.Expression) -> str | None:
        """Return DuckDB IEJoin SQL for *query* if applicable, else None.

        The dialect rewrites the *whole* query as
        ``SET VARIABLE __giql_iejoin_<token> = ...; SELECT ... FROM
        query(getvariable(...))``, so any top-level clause :meth:`_build_sql`
        does not read would be silently dropped. The check below is a
        whitelist: engage the dialect path only for a query that is
        exactly ``SELECT <qualified projections> FROM <registered base
        table> [JOIN <registered base table>] {ON|WHERE} <one
        column-to-column INTERSECTS>`` with no other top-level clause.
        Everything else returns ``None`` so the generic naive-predicate plan handles it.
        """
        if not isinstance(query, exp.Select):
            return None

        # Top-level WITH (CTEs) still falls back to the naive-predicate plan.
        # ORDER BY / LIMIT / OFFSET / DISTINCT / GROUP BY / HAVING ride
        # on the outer SELECT wrapper. ``DISTINCT ON (...)`` is the one
        # DISTINCT variant we don't try to translate — it would interact
        # with the outer alias rewrite in ways the naive-predicate plan handles
        # for free.
        if query.args.get("with_"):
            return None
        distinct_node = query.args.get("distinct")
        if distinct_node is not None and distinct_node.args.get("on"):
            return None

        if _has_outer_join_intersects(query):
            return None

        if _count_column_intersects(query) != 1:
            return None

        # The wrapper-relation rewriter is not scope-aware: any subquery in
        # GROUP BY / HAVING / ORDER BY (or inside a SELECT-list aggregate
        # argument) would have its inner column refs substituted with
        # wrapper-relation aliases that don't exist in the subquery's own
        # scope. Reject those shapes here. ``EXISTS (SELECT ...)`` parses
        # as ``Exists(Select(...))`` without an ``exp.Subquery`` wrapper,
        # so we also check ``exp.Select`` directly. Bare scalar subqueries
        # in the SELECT list (``SELECT (SELECT SUM(x) FROM ...)``) get a
        # targeted error from :meth:`_resolve_projections` rather than a
        # silent fallback — only subqueries hidden inside aggregate args
        # fall back here.
        for clause_key in ("group", "having", "order"):
            modifier = query.args.get(clause_key)
            if modifier is None:
                continue
            if (
                modifier.find(exp.Subquery) is not None
                or modifier.find(exp.Select) is not None
            ):
                return None
        for sel in query.expressions:
            target = sel.this if isinstance(sel, exp.Alias) else sel
            while isinstance(target, exp.Paren):
                target = target.this
            if isinstance(target, exp.AggFunc) and (
                target.find(exp.Subquery) is not None
                or target.find(exp.Select) is not None
            ):
                return None

        # The supported shape connects exactly two tables via one INTERSECTS,
        # so the top level has exactly one join (either an explicit JOIN ...
        # ON or an implicit comma-join in WHERE). A third comma/CROSS-joined
        # table would otherwise be silently dropped by _build_sql.
        joins = query.args.get("joins") or []
        if len(joins) != 1:
            return None
        the_join = joins[0]

        # NATURAL needs schema introspection we don't have at transpile time
        # (Table configs only expose chrom/start/end/strand; user columns
        # like name/score are invisible). The naive-predicate plan defers to
        # DuckDB, which does see the schema, so falling back IS the optimal plan.
        if the_join.args.get("method"):
            return None

        # USING admits only the single-column form whose column matches both
        # tables' ``chrom_col`` (the per-chromosome partition is exactly the
        # equi-join that ``USING(chrom)`` requests). The chrom-equivalence
        # check itself happens in ``_build_sql`` after ``l_table`` /
        # ``r_table`` are resolved. Multi-column USING and USING(non-chrom)
        # fall back; inline support for those is a documented follow-up.
        using_cols = the_join.args.get("using") or []
        if using_cols and len(using_cols) != 1:
            return None

        # INNER (the default), CROSS (functionally equivalent to INNER once
        # the chromosome partition is in place), SEMI, and ANTI all engage.
        # Everything else falls back so the naive-predicate plan can preserve its
        # semantics (e.g. RIGHT / FULL outer joins, which the naive overlap
        # predicate expresses as a plain outer ``ON`` condition).
        kind_str = the_join.args.get("kind")
        if kind_str and kind_str.upper() not in ("INNER", "CROSS", "SEMI", "ANTI"):
            return None

        intersects, the_join = self._find_target_join(query)
        if intersects is None or the_join is None:
            return None

        from_table = query.args["from_"].this
        if not isinstance(from_table, exp.Table):
            return None

        # Reject non-base-table operands. A GIQL table-function operand such
        # as ``DISJOIN(genes)`` parses as an ``exp.Table`` with an empty
        # ``name`` and would be interpolated as an empty identifier into
        # broken SQL. (A CTE-named operand is already handled by the
        # ``with_`` guard above; an unregistered base table is fine — the
        # dialect uses default columns just like the naive-predicate path does.)
        if not from_table.name:
            return None
        join_table = the_join.this
        if not isinstance(join_table, exp.Table) or not join_table.name:
            return None

        sides = _IEJoinSides.from_intersects(from_table, the_join, intersects)
        if sides is None:
            return None

        if sides.left_user_alias == sides.right_user_alias:
            # Same alias on both sides of the INTERSECTS — the inner subquery
            # would alias the same underlying relation as both "a" and "b",
            # which the naive-predicate plan renders correctly in place (no
            # inner subquery, so no ambiguous double-alias).
            return None

        # Compare fully-qualified identifiers (catalog/db/name) so two
        # distinct same-named tables in different schemas are not treated
        # as a self-join. Same qualified identifier still falls back so the
        # naive-predicate plan handles the legitimate self-join shape.
        if self._qualified_table_ident(sides.left_table) == self._qualified_table_ident(
            sides.right_table
        ):
            return None

        # Extra predicates are inlined into each per-chromosome subquery's
        # join ON (see :meth:`_build_sql`). Soft-fallback residuals (subquery,
        # aggregate, window, or an INTERSECTS still wrapped in OR/NOT/Paren)
        # cause ``_build_sql`` to return None. User-mistake residuals
        # (unqualified or unknown-aliased columns, or columns with a
        # catalog/schema qualifier) raise :class:`_UnqualifiedProjectionError`,
        # which we wrap to :class:`ValueError` here so the public surface is
        # uniform.
        try:
            return self._build_sql(query, intersects, sides, the_join)
        except _UnqualifiedProjectionError as exc:
            raise ValueError(str(exc)) from exc

    @staticmethod
    def _sql_escape(s: str) -> str:
        """Return *s* with single quotes doubled for safe SQL string literal use."""
        return s.replace("'", "''")

    @staticmethod
    def _sql_quote_ident(name: str) -> str:
        """Render *name* as a DuckDB-quoted identifier.

        Delegates to sqlglot's identifier emitter so embedded double quotes
        are doubled correctly and DuckDB-specific quoting rules (case,
        reserved words) are honored. Used at every identifier-interpolation
        site in the dialect's dynamic SQL.
        """
        return exp.to_identifier(name, quoted=True).sql(dialect="duckdb")

    @staticmethod
    def _qualified_table_ident(table: exp.Table) -> str:
        """Render *table*'s catalog/db/name as a qualified DuckDB identifier.

        Drops the user-supplied alias (``peaks AS a``) so the rendered string
        is suitable both as a bare ``FROM`` operand (where the caller adds
        its own alias) and as a subquery source. Catalog and schema qualifiers
        survive — sqlglot's :meth:`~sqlglot.expressions.Expression.sql` emitter
        handles identifier quoting per DuckDB's rules (reserved words and
        special characters get quoted; safe unquoted identifiers stay
        unquoted).
        """
        no_alias = table.copy()
        no_alias.set("alias", None)
        return no_alias.sql(dialect="duckdb")

    @staticmethod
    def _classify_extras(
        extras: list[exp.Expression],
    ) -> Literal["inline", "fallback"]:
        """Return ``"fallback"`` if any extra needs the naive-predicate plan.

        Soft-fallback shapes: an INTERSECTS still embedded in the residual
        (because :func:`_strip_intersects` couldn't peel an OR/NOT/Paren
        wrapper), a subquery, an aggregate function, or a window function.
        Returns ``"inline"`` when every extra is safe to inline into the
        per-chromosome subquery's join ON clause.
        """
        for predicate in extras:
            if predicate.find(Intersects) is not None:
                return "fallback"  # INTERSECTS wrapped in OR/NOT/Paren
            if predicate.find(exp.Subquery) is not None:
                return "fallback"
            if predicate.find(exp.Select) is not None:
                return "fallback"
            if predicate.find(exp.AggFunc) is not None:
                return "fallback"
            if predicate.find(exp.Window) is not None:
                return "fallback"
        return "inline"

    @staticmethod
    def _validate_extra_qualifiers(
        extras: list[exp.Expression],
        sides: "_IEJoinSides",
    ) -> None:
        """Raise for unqualified, unknown-aliased, or schema-qualified columns.

        Surfaces user mistakes at transpile time rather than deferring to a
        downstream binder error. Three failure modes raise:

        - **Unqualified column:** no ``.table`` qualifier.
        - **Unknown alias:** ``.table`` doesn't match either side of the join.
        - **Catalog/schema-qualified column:** ``mycat.myschema.a.score``
          would silently mis-bind once the alias-rewriter remaps ``a`` to
          the inner subquery's hardcoded ``a`` alias (the ``mycat.myschema``
          qualifier survives the rewrite and addresses a different relation
          than intended).
        """
        valid_aliases = {sides.left_user_alias, sides.right_user_alias}
        for predicate in extras:
            for col in predicate.find_all(exp.Column):
                col_table = _normalize_alias(col.table) if col.table else ""
                if not col_table:
                    raise _UnqualifiedProjectionError(
                        f"dialect='duckdb' cannot inline extra predicate "
                        f"{predicate.sql()!r}: column {col.sql()!r} must be "
                        f"qualified with {sides.left_user_alias!r} or "
                        f"{sides.right_user_alias!r}."
                    )
                if col_table not in valid_aliases:
                    raise _UnqualifiedProjectionError(
                        f"dialect='duckdb' cannot inline extra predicate "
                        f"{predicate.sql()!r}: column references unknown "
                        f"table qualifier {col.table!r}; expected "
                        f"{sides.left_user_alias!r} or "
                        f"{sides.right_user_alias!r}."
                    )
                if col.args.get("db") or col.args.get("catalog"):
                    raise _UnqualifiedProjectionError(
                        f"dialect='duckdb' cannot inline extra predicate "
                        f"{predicate.sql()!r}: column {col.sql()!r} carries "
                        f"a catalog/schema qualifier; use the alias-only "
                        f"form ({sides.left_user_alias!r}.<col> or "
                        f"{sides.right_user_alias!r}.<col>)."
                    )

    @staticmethod
    def _rewrite_refs_for_per_chrom_subquery(
        node: exp.Expression,
        sides: "_IEJoinSides",
    ) -> str:
        """Render a residual predicate against the inner subquery's ``a``/``b`` scope.

        The user's AST references their original aliases (``peaks a`` /
        ``genes b`` or ``peaks p`` / ``genes g``), but the per-chromosome
        inner subquery is emitted with hardcoded aliases ``a`` and ``b``
        regardless of what the user wrote. Rewrite each column reference
        accordingly before emitting the residual SQL into the inner join's
        ON clause. Alias matching is case-folded so mixed-case identifiers
        match the DuckDB-equivalent contract.

        Assumes residuals have already been validated to reference only
        ``sides.left_user_alias`` or ``sides.right_user_alias``; columns
        without a table qualifier or pointing at an unknown alias must
        have been rejected upstream.
        """

        def replace(n: exp.Expression) -> exp.Expression:
            if not isinstance(n, exp.Column):
                return n
            col_table = _normalize_alias(n.table) if n.table else ""
            if col_table == sides.left_user_alias:
                new = n.copy()
                new.set("table", exp.to_identifier(sides.left_inner_alias))
                return new
            if col_table == sides.right_user_alias:
                new = n.copy()
                new.set("table", exp.to_identifier(sides.right_inner_alias))
                return new
            return n

        rewritten = node.transform(replace, copy=True)
        return rewritten.sql(dialect="duckdb")

    @staticmethod
    def _rewrite_refs_for_wrapper_relation(
        node: exp.Expression,
        projection_map: dict[tuple[str, str], str],
        sides: "_IEJoinSides",
        clause_label: str,
    ) -> str:
        """Rewrite table-qualified column refs to the wrapper relation's inner alias.

        The dialect's outer ``SELECT`` reads from ``query(getvariable(...)) AS
        <alias>``. The wrapper relation's columns are whatever the inner
        per-chromosome subqueries projected — i.e., the ``__giql_p<n>`` names
        on the left-hand side of each ``AS`` in the inner SELECT. Any user
        ``ORDER BY`` / ``GROUP BY`` / ``HAVING`` (and any aggregate-function
        argument) written against the original table aliases (``a.start``,
        ``b.score``) must be rewritten to those inner alias names so DuckDB
        can resolve them against the wrapper relation. Alias matching is
        case-folded for DuckDB-equivalent semantics.

        An upstream pre-allocation step ensures every referenced column has
        an entry in ``projection_map``; missing keys raise
        :class:`_UnqualifiedProjectionError` as a guard against future scope
        relaxations that would otherwise produce silent miscompiles.

        Scope caveat: this walker is not scope-aware — it rewrites every
        :class:`~sqlglot.expressions.Column` whose ``.table`` matches either
        side's alias anywhere inside ``node``, including correlated or
        non-correlated subqueries that may rebind those aliases to a
        different table. The dispatcher gates already reject any subquery
        in GROUP BY / HAVING / ORDER BY / SELECT list so the unsafe shape
        can't reach here today. A future relaxation of those gates must
        add a scope-aware walker here.
        """
        valid_aliases = (sides.left_user_alias, sides.right_user_alias)

        def replace(n: exp.Expression) -> exp.Expression:
            if not isinstance(n, exp.Column):
                return n
            col_table = _normalize_alias(n.table) if n.table else ""
            if col_table not in valid_aliases:
                return n
            key = (col_table, n.name)
            if key not in projection_map:
                raise _UnqualifiedProjectionError(
                    f"dialect='duckdb' cannot resolve {n.sql()!r} in "
                    f"{clause_label}: ensure the column is referenced "
                    f"with a table qualifier matching the join's left "
                    f"or right alias."
                )
            return exp.column(projection_map[key])

        rewritten = node.transform(replace, copy=True)
        return rewritten.sql(dialect="duckdb")

    def _extract_extra_predicates(
        self,
        query: exp.Select,
        intersects: Intersects,
    ) -> list[exp.Expression]:
        """Return the non-INTERSECTS predicate residuals from joins and WHERE.

        Walks the AND tree under each JOIN ON and the WHERE, strips the
        target INTERSECTS node, and returns the residual predicates in
        document order. An empty list means the query has no extras
        beyond the target INTERSECTS.

        Note: :func:`_strip_intersects` only descends ``exp.And``. Predicates
        whose AND tree wraps the INTERSECTS in ``exp.Or`` / ``exp.Not`` /
        ``exp.Paren`` surface here with the INTERSECTS still embedded;
        :meth:`_classify_extras` then routes them to the naive-predicate plan, while
        :meth:`_validate_extra_qualifiers` enforces qualifier rules for the
        remaining inlinable residuals.
        """
        residuals: list[exp.Expression] = []
        for join in query.args.get("joins") or []:
            on = join.args.get("on")
            if on is None:
                continue
            residual = _strip_intersects(on, intersects)
            if residual is not None:
                residuals.append(residual)
        where = query.args.get("where")
        if where:
            residual = _strip_intersects(where.this, intersects)
            if residual is not None:
                residuals.append(residual)
        return residuals

    def _find_target_join(
        self, query: exp.Select
    ) -> tuple[Intersects | None, exp.Join | None]:
        """Locate the single INTERSECTS join target.

        Returns ``(intersects_node, join)`` where ``join`` is the
        :class:`~sqlglot.expressions.Join` AST node carrying the target
        INTERSECTS — either directly via its ``ON`` clause, or via an
        implicit cross-join whose ``INTERSECTS`` predicate lives in the
        top-level ``WHERE``. Alias matching on the WHERE-branch is
        case-folded for DuckDB-equivalent semantics.
        """
        for join in query.args.get("joins") or []:
            on = join.args.get("on")
            if on:
                intersects = _find_column_intersects_in(on)
                if intersects:
                    if not isinstance(join.this, exp.Table):
                        return None, None
                    return intersects, join

        where = query.args.get("where")
        if where:
            intersects = _find_column_intersects_in(where.this)
            if intersects:
                from_table = query.args["from_"].this
                if not isinstance(from_table, exp.Table):
                    return None, None
                from_alias = _normalize_alias(from_table.alias_or_name)
                left_alias = _normalize_alias(intersects.this.table)
                right_alias = _normalize_alias(intersects.expression.table)
                target_alias = right_alias if left_alias == from_alias else left_alias
                for join in query.args.get("joins") or []:
                    if isinstance(join.this, exp.Table):
                        if _normalize_alias(join.this.alias_or_name) == target_alias:
                            return intersects, join

        return None, None

    def _build_sql(
        self,
        query: exp.Select,
        intersects: Intersects,
        sides: "_IEJoinSides",
        the_join: exp.Join,
    ) -> str | None:
        """Render the multi-statement DuckDB IEJoin SQL string for *query*.

        Returns ``None`` when a soft-fallback condition fires (the residual
        of the join carries a shape the naive-predicate plan can handle but the
        dialect cannot inline, or a single-column USING references a column
        that is not both tables' ``chrom_col``). Raises
        :class:`_UnqualifiedProjectionError` when a user-mistake condition
        fires; callers translating to the public surface wrap the raise to
        :class:`ValueError`.
        """
        extras = self._extract_extra_predicates(query, intersects)
        if extras and self._classify_extras(extras) == "fallback":
            return None
        if extras:
            self._validate_extra_qualifiers(extras, sides)

        # Tables registry is keyed by the bare table name; an unregistered
        # table uses default column names like the naive-predicate plan does.
        l_table = self.tables.get(sides.left_table.name)
        r_table = self.tables.get(sides.right_table.name)
        l_chrom = l_table.chrom_col if l_table else DEFAULT_CHROM_COL
        l_start = l_table.start_col if l_table else DEFAULT_START_COL
        l_end = l_table.end_col if l_table else DEFAULT_END_COL
        r_chrom = r_table.chrom_col if r_table else DEFAULT_CHROM_COL
        r_start = r_table.start_col if r_table else DEFAULT_START_COL
        r_end = r_table.end_col if r_table else DEFAULT_END_COL

        # USING(<col>) admission check (the gate already restricted to
        # single-column USING). The dialect's per-chromosome partition IS
        # the equi-join that USING(<chrom_col>) requests; if the USING
        # column doesn't match both tables' chrom_col (or the chrom_cols
        # diverge), fall back so the naive-predicate plan preserves the
        # USING equi-join for the engine to resolve.
        using_cols = the_join.args.get("using") or []
        if using_cols:
            using_name = _normalize_alias(using_cols[0].name)
            if (
                _normalize_alias(l_chrom) != _normalize_alias(r_chrom)
                or _normalize_alias(l_chrom) != using_name
            ):
                return None

        (
            inner_projections,
            outer_projections,
            projection_map,
        ) = self._resolve_projections(
            query,
            sides,
            l_table,
            r_table,
            (l_chrom, l_start, l_end),
            (r_chrom, r_start, r_end),
        )

        # Per-call random token (full uuid4 hex = 128 bits) so the SET
        # VARIABLE name is collision-resistant even across many transpile()
        # calls interleaved in one DuckDB session. DuckDB session variables
        # are global session state, so token collision would silently
        # rebind the wrapper query to a different intersection.
        token = uuid4().hex
        var_name = f"__giql_iejoin_{token}"

        # Canonicalize endpoints to 0-based half-open and use strict
        # operators. Mixing inclusive operators with raw closed-closed
        # coordinates would silently report non-overlapping touching
        # intervals as overlapping.
        q = self._sql_quote_ident
        l_start_expr = canonical_start(f"a.{q(l_start)}", l_table)
        l_end_expr = canonical_end(f"a.{q(l_end)}", l_table)
        r_start_expr = canonical_start(f"b.{q(r_start)}", r_table)
        r_end_expr = canonical_end(f"b.{q(r_end)}", r_table)

        # Render each table operand qualified (catalog/db/name) without the
        # user's alias so the empty-schema fallback can supply its own
        # ``a``/``b`` alias without producing ``FROM peaks AS a a`` shapes.
        l_table_ident = self._qualified_table_ident(sides.left_table)
        r_table_ident = self._qualified_table_ident(sides.right_table)

        # Join keyword switches by ``sides.kind``: SEMI / ANTI engage
        # DuckDB's IE_JOIN with the corresponding semi/anti semantics; INNER
        # is the default. (CROSS folded to INNER inside ``_IEJoinSides``.)
        join_keyword = {
            "INNER": "JOIN",
            "SEMI": "SEMI JOIN",
            "ANTI": "ANTI JOIN",
        }[sides.kind]

        inner_select_list = ", ".join(inner_projections)
        per_chrom_select_prefix = f"SELECT {inner_select_list} "
        from_left = f"FROM (SELECT * FROM {l_table_ident} WHERE {q(l_chrom)} = "
        from_right_open = (
            f") a {join_keyword} (SELECT * FROM {r_table_ident} WHERE {q(r_chrom)} = "
        )
        on_clause_predicates = [
            f"{l_start_expr} < {r_end_expr}",
            f"{l_end_expr} > {r_start_expr}",
        ]
        for residual in extras:
            on_clause_predicates.append(
                self._rewrite_refs_for_per_chrom_subquery(residual, sides)
            )
        on_clause = ") b ON " + " AND ".join(on_clause_predicates)

        # The dynamic SQL builder, expressed as a SQL string-concat that
        # ``string_agg`` aggregates per-chromosome. Single quotes are
        # doubled because the result is itself an SQL string literal that
        # contains SQL.
        chrom_literal = "'''' || replace(chrom, '''', '''''') || ''''"

        esc = self._sql_escape
        per_chrom_sql_expr = (
            f"'{esc(per_chrom_select_prefix)}{esc(from_left)}' "
            f"|| {chrom_literal} "
            f"|| '{esc(from_right_open)}' "
            f"|| {chrom_literal} "
            f"|| '{esc(on_clause)}'"
        )

        # Empty-schema fallback when no chromosomes pass the partition
        # filter. Project from the source tables under WHERE FALSE so
        # DuckDB resolves every output column to its real declared type.
        # SEMI / ANTI outputs reference only the left side, so we drop
        # the right operand from the empty-schema FROM clause.
        empty_select_list = ", ".join(inner_projections)
        if sides.is_left_only_join:
            empty_schema = (
                f"SELECT {empty_select_list} FROM {l_table_ident} a WHERE FALSE"
            )
        else:
            empty_schema = (
                f"SELECT {empty_select_list} "
                f"FROM {l_table_ident} a, {r_table_ident} b WHERE FALSE"
            )

        # Partition source: INTERSECT for INNER / SEMI (a chromosome only
        # on one side produces no matches anyway), left-distinct for ANTI
        # (a chromosome present only on the left must still emit its rows
        # under ANTI semantics).
        if sides.kind == "ANTI":
            chrom_partition_subquery = (
                f"SELECT DISTINCT {q(l_chrom)} AS chrom FROM {l_table_ident}"
            )
        else:
            chrom_partition_subquery = (
                f"SELECT DISTINCT {q(l_chrom)} AS chrom FROM {l_table_ident} "
                f"INTERSECT SELECT DISTINCT {q(r_chrom)} AS chrom "
                f"FROM {r_table_ident}"
            )

        set_var_stmt = (
            f"SET VARIABLE {var_name} = COALESCE((\n"
            f"  SELECT string_agg(\n"
            f"    {per_chrom_sql_expr},\n"
            f"    ' UNION ALL '\n"
            f"  )\n"
            f"  FROM ({chrom_partition_subquery})\n"
            f"), '{esc(empty_schema)}')"
        )

        # Constant wrapper alias — uniqueness is provided by the
        # token-bearing session variable; the outer SELECT's relation
        # alias is never user-visible and doesn't need to vary per call.
        wrapper_alias = "__giql_iejoin_wrapper"
        outer_projection = ", ".join(outer_projections)
        distinct_kw = "DISTINCT " if query.args.get("distinct") else ""
        outer_select_parts = [
            f"SELECT {distinct_kw}{outer_projection} "
            f"FROM query(getvariable('{var_name}')) AS {wrapper_alias}"
        ]

        # The rewriter translates any ``a.col`` / ``b.col`` references in
        # the GROUP BY / HAVING / ORDER BY clauses to the wrapper
        # relation's inner-alias column names (``__giql_p<n>``), since
        # the user's original aliases are not in scope in the outer SELECT.
        group_node = query.args.get("group")
        if group_node is not None:
            outer_select_parts.append(
                self._rewrite_refs_for_wrapper_relation(
                    group_node, projection_map, sides, "GROUP BY"
                )
            )

        having_node = query.args.get("having")
        if having_node is not None:
            outer_select_parts.append(
                self._rewrite_refs_for_wrapper_relation(
                    having_node, projection_map, sides, "HAVING"
                )
            )

        order_node = query.args.get("order")
        if order_node is not None:
            outer_select_parts.append(
                self._rewrite_refs_for_wrapper_relation(
                    order_node, projection_map, sides, "ORDER BY"
                )
            )

        limit_node = query.args.get("limit")
        if limit_node is not None:
            outer_select_parts.append(limit_node.sql(dialect="duckdb"))

        offset_node = query.args.get("offset")
        if offset_node is not None:
            outer_select_parts.append(offset_node.sql(dialect="duckdb"))

        outer_select = " ".join(outer_select_parts)

        return set_var_stmt + ";\n" + outer_select

    def _resolve_projections(
        self,
        query: exp.Select,
        sides: "_IEJoinSides",
        l_table: Table | None,
        r_table: Table | None,
        l_cols: tuple[str, str, str],
        r_cols: tuple[str, str, str],
    ) -> tuple[list[str], list[str], dict[tuple[str, str], str]]:
        """Resolve the SELECT list (and any modifier-referenced columns).

        Returns ``(inner_projections, outer_projections, projection_map)``:

        * ``inner_projections`` is the per-chromosome subquery's SELECT
          list — every column the query touches in any clause, projected
          once each under a unique ``__giql_p<n>`` alias.
        * ``outer_projections`` is the outer SELECT's projection list,
          rendered using the inner aliases (qualified columns, aggregate
          calls, ``a.*``/``b.*`` expansions including ``strand_col`` when
          the corresponding :class:`~giql.table.Table` config declares one).
        * ``projection_map`` binds each ``(normalized_alias, column_name)``
          referenced anywhere in the query to its inner alias name so the
          wrapper-relation rewriter and aggregate-argument rewriting can
          translate user-written references to the wrapper relation.

        Raises :class:`_UnqualifiedProjectionError` for SELECT-list shapes
        the dialect cannot translate (bare ``*``, unqualified column,
        unknown table qualifier, ``a.* AS x``, window aggregate, ``FILTER``
        clause, scalar subquery, arithmetic over an aggregate, aggregate
        over an unqualified column), and (when ``sides.is_left_only_join``)
        any reference to the right side of the join.
        """
        l_chrom, l_start, l_end = l_cols
        r_chrom, r_start, r_end = r_cols
        q = self._sql_quote_ident
        l_alias = sides.left_user_alias
        r_alias = sides.right_user_alias

        def star_columns(
            table: Table | None, defaults: tuple[str, str, str]
        ) -> tuple[str, ...]:
            """Expand ``a.*`` / ``b.*`` to the genomic columns declared by
            the :class:`~giql.table.Table` config (chrom / start / end, plus
            ``strand_col`` when set). Without this, ``strand_col`` would
            silently disappear under ``dialect='duckdb'`` even when the
            user has registered it.
            """
            chrom, start, end = defaults
            if table is not None and table.strand_col:
                return (chrom, start, end, table.strand_col)
            return (chrom, start, end)

        inner_projections: list[str] = []
        projection_map: dict[tuple[str, str], str] = {}

        def reject_right_side_if_left_only(tbl_normalized: str, source_sql: str) -> None:
            """Raise when a SEMI / ANTI join projects from the right side."""
            if sides.is_left_only_join and tbl_normalized == r_alias:
                raise _UnqualifiedProjectionError(
                    f"dialect='duckdb' with kind={sides.kind!r} only "
                    f"projects left-side columns; got {source_sql!r} "
                    f"which references the right side "
                    f"({sides.right_user_alias!r})."
                )

        def allocate(tbl: str, col: str) -> str:
            """Return the inner alias for ``(tbl, col)``, allocating if new.

            Allocation order is deterministic per query but depends on
            (a) the fixed pre-walk order in :meth:`_resolve_projections` —
            modifier columns first (``("group", "having", "order")``),
            then aggregate-argument columns, then the SELECT list — and
            (b) sqlglot's ``find_all`` traversal order within each
            sub-clause. Tests therefore should NOT couple to a specific
            ``__giql_p<n>`` ↔ source-column mapping (e.g.
            ``assert projection_map[(a, start)] == "__giql_p3"``).
            """
            tbl = _normalize_alias(tbl)
            key = (tbl, col)
            existing = projection_map.get(key)
            if existing is not None:
                return existing
            if tbl == l_alias:
                side = sides.left_inner_alias
            elif tbl == r_alias:
                side = sides.right_inner_alias
            else:
                raise _UnqualifiedProjectionError(
                    f"Unknown table qualifier {tbl!r} in projection; "
                    f"expected {l_alias!r} or {r_alias!r}."
                )
            alias = f"__giql_p{len(inner_projections)}"
            inner_projections.append(f"{side}.{q(col)} AS {alias}")
            projection_map[key] = alias
            return alias

        def render_aggregate(agg_node: exp.Expression, user_alias: str | None) -> None:
            """Validate aggregate-argument qualifiers and emit the outer projection."""
            for col in agg_node.find_all(exp.Column):
                col_table = _normalize_alias(col.table) if col.table else ""
                if col_table not in (l_alias, r_alias):
                    raise _UnqualifiedProjectionError(
                        "dialect='duckdb' requires qualified aggregate "
                        f"arguments; got {agg_node.sql()!r}. Use a "
                        "qualified column reference like "
                        f"'{type(agg_node).__name__.upper()}(a.col)' or "
                        f"'{type(agg_node).__name__.upper()}(b.col)'."
                    )
                reject_right_side_if_left_only(col_table, agg_node.sql())
            rewritten_agg = self._rewrite_refs_for_wrapper_relation(
                agg_node, projection_map, sides, "aggregate argument"
            )
            outer_name = (
                user_alias
                if user_alias is not None
                else f"__giql_agg_{len(outer_projections)}"
            )
            outer_projections.append(f"{rewritten_agg} AS {q(outer_name)}")

        # Pre-allocate inner aliases for every column referenced in
        # GROUP BY / HAVING / ORDER BY so the wrapper relation exposes
        # them even if the user didn't put them in the SELECT list.
        for clause_key in ("group", "having", "order"):
            modifier_node = query.args.get(clause_key)
            if modifier_node is None:
                continue
            for col in modifier_node.find_all(exp.Column):
                col_table = _normalize_alias(col.table) if col.table else ""
                if col_table in (l_alias, r_alias):
                    reject_right_side_if_left_only(col_table, col.sql())
                    allocate(col_table, col.name)

        # Pre-allocate inner aliases for columns referenced inside any
        # aggregate function in the SELECT list (we'll rewrite the
        # aggregate's argument to use the inner alias when rendering the
        # outer projection below). Peel any ``exp.Paren`` wrapper first
        # so a paren-wrapped aggregate like ``(SUM(a.score)) AS s`` is
        # treated identically to its bare form.
        for sel in query.expressions:
            target = sel.this if isinstance(sel, exp.Alias) else sel
            while isinstance(target, exp.Paren):
                target = target.this
            if isinstance(target, exp.AggFunc):
                for col in target.find_all(exp.Column):
                    col_table = _normalize_alias(col.table) if col.table else ""
                    if col_table in (l_alias, r_alias):
                        reject_right_side_if_left_only(col_table, col.sql())
                        allocate(col_table, col.name)

        outer_projections: list[str] = []

        for sel in query.expressions:
            target = sel.this if isinstance(sel, exp.Alias) else sel
            user_alias = sel.alias if isinstance(sel, exp.Alias) else None
            # Peel paren wrappers at the top of dispatch so paren-wrapped
            # aggregates like ``(SUM(a.score)) AS s`` route through the
            # AggFunc branch below, and paren-wrapped diagnostic shapes
            # (``(SUM(a.score) OVER (...))``) hit the targeted error
            # branches below rather than the generic catch-all.
            while isinstance(target, exp.Paren):
                target = target.this

            if isinstance(target, exp.Star) and not isinstance(sel, exp.Column):
                raise _UnqualifiedProjectionError(
                    "dialect='duckdb' requires qualified projections; bare '*' "
                    "is not supported. Use 'a.*', 'b.*', or qualified columns."
                )

            if isinstance(target, exp.Column) and isinstance(target.this, exp.Star):
                # Reject ``a.* AS x`` — most SQL engines reject
                # star-with-alias outright, and the dialect has no
                # reasonable interpretation (apply ``x`` as a prefix?
                # as a single composite column name?).
                if user_alias is not None:
                    raise _UnqualifiedProjectionError(
                        "dialect='duckdb' does not support star projections "
                        f"with a user alias; got {sel.sql()!r}. Either drop "
                        "the alias (use ``a.*``) or list the columns "
                        "explicitly with per-column aliases."
                    )
                tbl = _normalize_alias(target.table) if target.table else ""
                if tbl == l_alias:
                    cols = star_columns(l_table, (l_chrom, l_start, l_end))
                elif tbl == r_alias:
                    cols = star_columns(r_table, (r_chrom, r_start, r_end))
                else:
                    raise _UnqualifiedProjectionError(
                        f"Unknown table qualifier {target.table!r} in projection; "
                        f"expected {l_alias!r} or {r_alias!r}."
                    )
                reject_right_side_if_left_only(tbl, sel.sql())
                # Expand to the genomic columns the Table config knows
                # about; arbitrary additional columns require explicit
                # listing since schema introspection isn't available at
                # transpile time.
                for col in cols:
                    alias = allocate(tbl, col)
                    outer_projections.append(f"{alias} AS {q(col)}")
                continue

            if isinstance(target, exp.AggFunc):
                render_aggregate(target, user_alias)
                continue

            if isinstance(target, exp.Column) and target.table:
                tbl = _normalize_alias(target.table)
                reject_right_side_if_left_only(tbl, target.sql())
                col = target.name
                alias = allocate(tbl, col)
                outer_name = user_alias if user_alias is not None else col
                outer_projections.append(f"{alias} AS {q(outer_name)}")
                continue

            # Diagnostic branches for shapes the dialect deliberately
            # rejects. Specialize the error message when we can recognize
            # the unsupported shape so the user isn't told "use a
            # qualified column" when they actually wrote an aggregate
            # variant the dispatcher couldn't classify. Direct-isinstance
            # branches catch bare shapes; the ``target.find(...)`` branches
            # below catch nested forms (``CAST(SUM(a.x) OVER (...) AS
            # DOUBLE)`` etc.) before they fall through to the generic
            # AggFunc-found catch-all.
            if isinstance(target, exp.Window):
                raise _UnqualifiedProjectionError(
                    "dialect='duckdb' does not support window aggregates "
                    f"in the SELECT list; got {target.sql()!r}. Either "
                    "drop the OVER (...) clause, or omit dialect='duckdb' "
                    "to use the generic naive-predicate plan."
                )
            if isinstance(target, exp.Filter):
                raise _UnqualifiedProjectionError(
                    "dialect='duckdb' does not support FILTER (WHERE ...) "
                    f"clauses on aggregates; got {target.sql()!r}. Either "
                    "rewrite the FILTER as a WHERE predicate alongside "
                    "the INTERSECTS, or omit dialect='duckdb'."
                )
            if isinstance(target, exp.Subquery):
                # A scalar subquery in the SELECT list (e.g.
                # ``(SELECT SUM(x) FROM other_table)``) is shaped like an
                # aggregate but doesn't fit the dialect's projection rules.
                raise _UnqualifiedProjectionError(
                    "dialect='duckdb' does not support scalar subqueries in "
                    f"the SELECT list; got {target.sql()!r}. Either rewrite "
                    "without the subquery, or omit dialect='duckdb' to use "
                    "the generic naive-predicate plan."
                )
            if target.find(exp.Window) is not None:
                raise _UnqualifiedProjectionError(
                    "dialect='duckdb' does not support window aggregates "
                    f"in the SELECT list; got {target.sql()!r}. Either "
                    "drop the OVER (...) clause, or omit dialect='duckdb' "
                    "to use the generic naive-predicate plan."
                )
            if target.find(exp.Filter) is not None:
                raise _UnqualifiedProjectionError(
                    "dialect='duckdb' does not support FILTER (WHERE ...) "
                    f"clauses on aggregates; got {target.sql()!r}. Either "
                    "rewrite the FILTER as a WHERE predicate alongside "
                    "the INTERSECTS, or omit dialect='duckdb'."
                )
            if target.find(exp.AggFunc) is not None:
                raise _UnqualifiedProjectionError(
                    "dialect='duckdb' does not support aggregates inside "
                    f"arithmetic or function expressions; got {target.sql()!r}. "
                    "Project the underlying aggregate as its own column and "
                    "do the arithmetic in a wrapping query, or omit "
                    "dialect='duckdb'."
                )
            raise _UnqualifiedProjectionError(
                "dialect='duckdb' requires qualified projections; got "
                f"{target.sql()!r}. Use a qualified column reference like "
                "'a.col' or 'b.col [AS x]'."
            )

        if not outer_projections:
            raise _UnqualifiedProjectionError(
                "dialect='duckdb' requires at least one qualified projection."
            )

        # An aggregate without any source-column argument (e.g.
        # ``COUNT(*)``) leaves ``inner_projections`` empty in queries that
        # have no other column references. The per-chromosome subquery
        # still needs at least one column to project, so emit a constant
        # placeholder that DuckDB optimizes away. (Without this the
        # generated SQL ``SELECT  FROM (...) a JOIN (...) b ON ...`` would
        # fail to parse.)
        if not inner_projections:
            inner_projections.append("1 AS __giql_placeholder")

        return inner_projections, outer_projections, projection_map
