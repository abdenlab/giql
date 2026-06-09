from dataclasses import dataclass

from sqlglot import exp
from sqlglot.generator import Generator

from giql.canonical import canonical_end
from giql.canonical import canonical_start
from giql.canonical import decanonical_end
from giql.canonical import decanonical_start
from giql.constants import CANONICAL_DEFAULT_COLS
from giql.constants import DEFAULT_CHROM_COL
from giql.constants import DEFAULT_END_COL
from giql.constants import DEFAULT_START_COL
from giql.constants import DEFAULT_STRAND_COL
from giql.expressions import Contains
from giql.expressions import GIQLDisjoin
from giql.expressions import GIQLDistance
from giql.expressions import GIQLNearest
from giql.expressions import Intersects
from giql.expressions import SpatialSetPredicate
from giql.expressions import Within
from giql.range_parser import ParsedRange
from giql.range_parser import RangeParser
from giql.table import Table
from giql.table import Tables

# Reserved identifier prefix for DISJOIN's internal CTE aliases. User-supplied
# target/reference names or enclosing CTE names sharing this prefix would
# collide with the generator's emitted SQL.
_DISJOIN_RESERVED_PREFIX = "__giql_dj_"

# Internal CTE aliases emitted by giqldisjoin_sql. Kept as named constants so
# a typo in the emitted SQL produces a Python NameError rather than a
# downstream engine "unknown CTE" error.
#
# - ``_DISJOIN_CTE_REF``: filtered/canonicalized reference rows.
# - ``_DISJOIN_CTE_TGT``: target rows brought into scope.
# - ``_DISJOIN_CTE_BP``: breakpoint positions ``(chrom, pos)``.
# - ``_DISJOIN_CTE_CUTS``: cut points per target row.
# - ``_DISJOIN_CTE_SEGS``: contiguous segments via ``LEAD()``.
# - ``_DISJOIN_CTE_RS``: alias for an inlined subquery reference.
_DISJOIN_CTE_REF = "__giql_dj_ref"
_DISJOIN_CTE_TGT = "__giql_dj_tgt"
_DISJOIN_CTE_BP = "__giql_dj_bp"
_DISJOIN_CTE_CUTS = "__giql_dj_cuts"
_DISJOIN_CTE_SEGS = "__giql_dj_segs"
_DISJOIN_CTE_RS = "__giql_dj_rs"


def _quote_relation(name: str) -> str:
    """Quote a relation identifier for emission into SQL.

    Escapes embedded double quotes by doubling them, matching the SQL
    standard for delimited identifiers.
    """
    return '"' + name.replace('"', '""') + '"'


class _NotABareName(Exception):
    """Internal marker raised by :func:`_extract_bare_name`.

    Callers catch this to attach operator/role-specific wording while
    bare-name validation stays centralized.
    """

    def __init__(self, node: exp.Expression) -> None:
        super().__init__()
        self.node = node


def _is_qualified_name(node: exp.Expression) -> bool:
    """Return ``True`` when ``node`` is a db/schema-qualified Table or a
    table-qualified Column. Used by callers to phrase ``is qualified``
    error messages after :func:`_extract_bare_name` raises.
    """
    if isinstance(node, exp.Table):
        return bool(node.db or node.catalog)
    if isinstance(node, exp.Column):
        return bool(node.table)
    return False


def _bare_name_or_none(node: exp.Expression | None) -> str | None:
    """Return :func:`_extract_bare_name`\\ 's result, or ``None`` on failure.

    Convenience wrapper for callers (e.g. ``_check_reserved_collisions``)
    that only care whether the node *is* a bare name, not why it isn't.
    """
    if node is None:
        return None
    try:
        return _extract_bare_name(node)
    except _NotABareName:
        return None


def _extract_bare_name(node: exp.Expression) -> str:
    """Return the bare identifier for a DISJOIN/NEAREST target or DISJOIN
    reference node.

    Quoted identifiers are normalized via sqlglot's unquoted ``.name`` so
    ``DISJOIN("X")`` matches a same-named enclosing ``WITH "X" AS (...)``
    whose alias is also collected unquoted.

    :param node:
        Accepts :class:`sqlglot.exp.Table`, unqualified
        :class:`sqlglot.exp.Column`, and bare :class:`sqlglot.exp.Identifier`.
    :return:
        The bare identifier.
    :raises _NotABareName:
        If the node is a qualified column, db/schema-qualified table, or
        any other shape that does not resolve to a single bare name.
        Callers catch this marker to phrase operator/role-specific
        wording.
    """
    if isinstance(node, exp.Table) and not (node.db or node.catalog):
        return node.name
    if isinstance(node, exp.Column) and not node.table:
        return node.name
    if isinstance(node, exp.Identifier):
        return node.name
    raise _NotABareName(node)


def _enclosing_cte_names(expression: exp.Expression) -> set[str]:
    """Collect CTE names visible from enclosing ``WITH`` clauses.

    Walks the ancestor chain and gathers the aliases of every
    :class:`exp.With`-clause CTE that is visible at ``expression``'s
    position per SQL scoping rules: a CTE is visible from the containing
    ``SELECT`` body, and from later siblings inside the same ``WITH``.
    From inside a CTE body only earlier siblings are visible (forward
    references are forbidden), and the CTE itself is visible only under
    ``WITH RECURSIVE``.

    ``WITH`` attaches to ``SELECT``, set operations (``UNION`` /
    ``INTERSECT`` / ``EXCEPT``), and DML/DDL forms (``INSERT`` /
    ``UPDATE`` / ``DELETE`` / ``MERGE`` / ``CREATE`` / ``SUBQUERY``);
    reading the ``with_`` arg from every ancestor covers all of them.
    Names are de-duplicated; innermost shadowing is not modeled —
    callers only test membership.

    :param expression:
        Any sqlglot expression node. The walk starts at this node and
        proceeds via ``.parent``.
    :return:
        Unordered set of CTE aliases visible at this node per SQL
        scoping. Empty when no enclosing ``WITH`` exists.
    """
    names: set[str] = set()
    last_cte: exp.CTE | None = None
    node: exp.Expression | None = expression
    while node is not None:
        if isinstance(node, exp.CTE):
            last_cte = node
        with_clause = node.args.get("with_")
        # Every WITH-bearing sqlglot node we care about stores the clause
        # under ``with_``; surface a schema rename loudly rather than
        # silently dropping visible CTE names. ``raise`` survives ``-O``,
        # which strips ``assert``.
        if with_clause is not None and not isinstance(with_clause, exp.With):
            raise TypeError(
                f"Expected exp.With under 'with_' arg, got {type(with_clause).__name__}"
            )
        if with_clause is not None:
            ctes = list(with_clause.expressions)
            if last_cte in ctes:
                recursive = bool(with_clause.args.get("recursive"))
                for cte in ctes:
                    if cte is last_cte:
                        if recursive and cte.alias:
                            names.add(cte.alias)
                        break
                    if cte.alias:
                        names.add(cte.alias)
            else:
                # Ascended from the SELECT body (or another non-CTE child);
                # all sibling CTEs in this WITH are visible.
                for cte in ctes:
                    if cte.alias:
                        names.add(cte.alias)
        node = node.parent
    return names


@dataclass(frozen=True, slots=True, eq=False)
class _ResolvedRelation:
    """A bare-name reference resolved against registered tables or enclosing
    CTEs.

    :ivar name:
        The unquoted identifier. Empty string for an inlined subquery
        reference (its :attr:`_ResolvedReference.subquery` field carries
        the AST node).
    :ivar cols:
        ``(chrom_col, start_col, end_col)``. Canonical defaults
        (:data:`giql.constants.CANONICAL_DEFAULT_COLS`) for a CTE or
        subquery relation; configured column names for a registered table.
    :ivar table:
        The registered :class:`Table` when the name maps to one, else
        ``None``. ``None`` signals a CTE or subquery relation, which by
        contract exposes canonical 0-based half-open columns.
    """

    name: str
    cols: tuple[str, str, str]
    table: Table | None


@dataclass(frozen=True, slots=True, eq=False)
class _ResolvedReference:
    """A resolved DISJOIN reference.

    :ivar relation:
        The relation-level fields (name, cols, table). For a subquery
        reference the relation's ``name`` is empty and its ``table`` is
        ``None`` — :attr:`subquery` is what the emitter uses to build the
        ``FROM`` clause.
    :ivar subquery:
        The inlined ``(SELECT ...)`` AST node when the user passed one as
        ``reference := (SELECT ...)``; otherwise ``None``.
    :ivar is_self_reference:
        ``True`` only when both target and reference provably resolve to
        the same registered :class:`Table`. The ``EXISTS`` coverage filter
        in ``giqldisjoin_sql`` is then a tautology and is omitted. For a
        CTE target the engine may inline the body twice with divergent
        results, so the optimization is unsafe and this flag stays
        ``False``; subquery references and CTE-shadowed names are also
        treated as distinct because row-set identity to the target is not
        provable from the AST.
    """

    relation: _ResolvedRelation
    subquery: exp.Subquery | None
    is_self_reference: bool


class BaseGIQLGenerator(Generator):
    """Base generator that outputs standard SQL.

    Works with any SQL database that supports:

    - Basic comparison operators (<, >, =, AND, OR)
    - String literals
    - Numeric comparisons

    This generator uses only SQL-92 compatible constructs, ensuring
    compatibility with virtually all SQL databases.
    """

    # Most databases support LATERAL joins (PostgreSQL 9.3+, DuckDB 0.7.0+)
    # SQLite does not support LATERAL, so it overrides this to False
    SUPPORTS_LATERAL = True

    @staticmethod
    def _extract_bool_param(param_expr: exp.Expression | None) -> bool:
        """Extract boolean value from a parameter expression.

        Handles exp.Boolean, exp.Literal, and string representations.
        """
        if not param_expr:
            return False
        elif isinstance(param_expr, exp.Boolean):
            return param_expr.this
        else:
            return str(param_expr).upper() in ("TRUE", "1", "YES")

    def __init__(self, tables: Tables | None = None, **kwargs):
        super().__init__(**kwargs)
        self.tables = tables or Tables()
        self._current_table = None  # Track current table for column resolution
        self._alias_to_table = {}  # Map aliases to table names

    def select_sql(self, expression: exp.Select) -> str:
        """Override SELECT generation to track table context and aliases."""
        # Build alias-to-table mapping
        self._alias_to_table = {}

        # Extract from FROM clause
        if expression.args.get("from_"):
            from_clause = expression.args["from_"]
            if isinstance(from_clause.this, exp.Table):
                table_name = from_clause.this.name
                self._current_table = table_name
                # Check if table has an alias
                if from_clause.this.alias:
                    self._alias_to_table[from_clause.this.alias] = table_name
                else:
                    # No alias, table referenced by name
                    self._alias_to_table[table_name] = table_name

        # Extract from JOINs
        if expression.args.get("joins"):
            for join in expression.args["joins"]:
                if isinstance(join.this, exp.Table):
                    table_name = join.this.name
                    # Check if table has an alias
                    if join.this.alias:
                        self._alias_to_table[join.this.alias] = table_name
                    else:
                        self._alias_to_table[table_name] = table_name

        # Call parent implementation
        return super().select_sql(expression)

    def intersects_sql(self, expression: Intersects) -> str:
        """Generate standard SQL for INTERSECTS.

        :param expression:
            INTERSECTS expression node
        :return:
            SQL predicate string
        """
        return self._generate_spatial_op(expression, "intersects")

    def contains_sql(self, expression: Contains) -> str:
        """Generate standard SQL for CONTAINS.

        :param expression:
            CONTAINS expression node
        :return:
            SQL predicate string
        """
        return self._generate_spatial_op(expression, "contains")

    def within_sql(self, expression: Within) -> str:
        """Generate standard SQL for WITHIN.

        :param expression:
            WITHIN expression node
        :return:
            SQL predicate string
        """
        return self._generate_spatial_op(expression, "within")

    def spatialsetpredicate_sql(self, expression: SpatialSetPredicate) -> str:
        """Generate SQL for spatial set predicates (ANY/ALL).

        :param expression:
            SpatialSetPredicate expression node
        :return:
            SQL predicate string
        """
        return self._generate_spatial_set(expression)

    def giqlnearest_sql(self, expression: GIQLNearest) -> str:
        """Generate SQL for NEAREST function.

        Detects mode (standalone vs correlated) and generates appropriate SQL:
        - Standalone: Direct query with ORDER BY + LIMIT
        - Correlated (LATERAL): Subquery for k-nearest neighbors

        :param expression:
            GIQLNearest expression node.
        :return:
            SQL string for NEAREST operation.
        :raises ValueError:
            If the target node is not a bare name (e.g. ``db.genes`` or a
            subquery); if a bare target resolves to neither a registered
            table nor an enclosing CTE; if a bare target name resolves to
            both a registered table and an enclosing CTE (ambiguous); if
            the target name resolves only to an enclosing CTE (CTE
            targets are unsupported for NEAREST — register the relation
            as a table); or if NEAREST runs in correlated mode against a
            dialect lacking ``LATERAL`` support (SQLite).
        """
        # Detect mode
        mode = self._detect_nearest_mode(expression)

        # Resolve target table
        target = self._resolve_target(
            expression, accept_cte=False, operator_label="NEAREST"
        )
        target_name = target.name
        target_chrom, target_start, target_end = target.cols
        target_table = target.table
        # Quoted form used in emitted SQL (FROM clause and column prefixes).
        # Quoting preserves case on engines that fold unquoted identifiers
        # and escapes any embedded double quotes per the SQL standard.
        target_ref = _quote_relation(target_name)

        # Resolve reference
        ref_chrom, ref_start, ref_end = self._resolve_nearest_reference(expression, mode)

        # Extract parameters
        k = expression.args.get("k")
        k_value = int(str(k)) if k else 1  # Default k=1

        max_distance = expression.args.get("max_distance")
        max_dist_value = int(str(max_distance)) if max_distance else None

        is_stranded = self._extract_bool_param(expression.args.get("stranded"))
        is_signed = self._extract_bool_param(expression.args.get("signed"))

        # Resolve strand columns if stranded mode
        ref_strand = None
        target_strand = None
        if is_stranded:
            # Get strand column for reference
            reference = expression.args.get("reference")
            if reference:
                reference_sql = self.sql(reference)

                # Check if reference is a literal string or column reference
                if reference_sql.startswith("'") or reference_sql.startswith('"'):
                    # Literal reference - parse for strand
                    range_str = reference_sql.strip("'\"")
                    from giql.range_parser import RangeParser

                    parsed_range = RangeParser.parse(range_str).to_zero_based_half_open()
                    if parsed_range.strand:
                        ref_strand = f"'{parsed_range.strand}'"
                else:
                    # Column reference - get strand column
                    ref_cols = self._get_column_refs(
                        reference_sql, None, include_strand=True
                    )
                    if len(ref_cols) == 4:
                        ref_strand = ref_cols[3]
            else:
                # Implicit reference in correlated mode - get strand from outer table
                outer_table = self._find_outer_table_in_lateral_join(expression)
                if outer_table and self.tables:
                    actual_table = self._alias_to_table.get(outer_table, outer_table)
                    table = self.tables.get(actual_table)
                    if table and table.strand_col:
                        ref_strand = f'{outer_table}."{table.strand_col}"'

            # Get strand column for target table
            if target_table and target_table.strand_col:
                target_strand = f'{target_ref}."{target_table.strand_col}"'

        # Distance math below assumes 0-based half-open.
        target_start_expr = canonical_start(
            f'{target_ref}."{target_start}"', target_table
        )
        target_end_expr = canonical_end(
            f'{target_ref}."{target_end}"', target_table
        )

        # Build distance calculation using CASE expression
        # For NEAREST: ORDER BY absolute distance, but RETURN signed distance
        distance_expr = self._generate_distance_case(
            ref_chrom,
            ref_start,
            ref_end,
            ref_strand,
            f'{target_ref}."{target_chrom}"',
            target_start_expr,
            target_end_expr,
            target_strand,
            stranded=is_stranded,
            signed=is_signed,
        )

        # Use absolute distance for ordering and filtering
        abs_distance_expr = f"ABS({distance_expr})"

        # Build WHERE clauses
        where_clauses = [
            f'{ref_chrom} = {target_ref}."{target_chrom}"'  # Chromosome pre-filter
        ]

        # Add strand matching for stranded mode
        if is_stranded and ref_strand and target_strand:
            where_clauses.append(f"{ref_strand} = {target_strand}")

        if max_dist_value is not None:
            where_clauses.append(f"({abs_distance_expr}) <= {max_dist_value}")

        where_sql = " AND ".join(where_clauses)

        # Generate SQL based on mode
        if mode == "standalone":
            # Standalone mode: direct ORDER BY + LIMIT
            # Return signed distance, but order by absolute distance
            sql = f"""(
                SELECT {target_ref}.*, {distance_expr} AS distance
                FROM {target_ref}
                WHERE {where_sql}
                ORDER BY {abs_distance_expr}
                LIMIT {k_value}
            )"""
        else:
            # Correlated mode: requires LATERAL join support
            if not self.SUPPORTS_LATERAL:
                raise ValueError(
                    "NEAREST in correlated mode (CROSS JOIN LATERAL) is not supported "
                    "in SQLite. SQLite does not support LATERAL joins. "
                    "\n\nAlternatives:"
                    "\n1. Use standalone mode: SELECT * FROM NEAREST(table, "
                    "reference='chr1:100-200', k=3)"
                    "\n2. Use DuckDB for queries requiring LATERAL joins"
                    "\n3. Manually write equivalent window function query"
                )

            # LATERAL mode: subquery for k-nearest neighbors
            # Return signed distance, but order by absolute distance
            sql = f"""(
                SELECT {target_ref}.*, {distance_expr} AS distance
                FROM {target_ref}
                WHERE {where_sql}
                ORDER BY {abs_distance_expr}
                LIMIT {k_value}
            )"""

        return sql.strip()

    def giqldisjoin_sql(self, expression: GIQLDisjoin) -> str:
        """Generate SQL for the DISJOIN table function.

        DISJOIN splits each target interval at every reference breakpoint
        strictly interior to it. The full target row passes through unchanged;
        the sub-interval is appended as ``disjoin_chrom`` / ``disjoin_start`` /
        ``disjoin_end`` in the target table's coordinate system. A coverage
        filter drops sub-intervals overlapping no reference interval. When no
        ``reference`` is given it defaults to the target set.

        The target may be a registered table or an enclosing-WITH CTE; a CTE
        target is assumed canonical (0-based half-open), and its appended
        ``disjoin_*`` columns are emitted in that convention.

        :param expression:
            GIQLDisjoin expression node.
        :return:
            SQL string (a parenthesized WITH-CTE subquery) for the DISJOIN
            table.
        :raises ValueError:
            If the target or reference name resolves to neither a registered
            table nor an enclosing CTE; if the target or reference name uses
            the reserved ``__giql_dj_`` prefix; or if any enclosing CTE alias
            uses the reserved ``__giql_dj_`` prefix.
        """
        # Reject any user-visible identifier that would shadow DISJOIN's
        # internal __giql_dj_* namespace. Run before resolution so a reserved
        # name reports as such even when no table/CTE of that name exists.
        self._check_reserved_collisions(expression)
        target = self._resolve_target(
            expression, accept_cte=True, operator_label="DISJOIN"
        )
        reference = self._resolve_disjoin_reference(expression, target)

        target_name = target.name
        target_chrom, target_start, target_end = target.cols
        target_table = target.table

        ref_chrom, ref_start, ref_end = reference.relation.cols
        ref_table = reference.relation.table

        # Endpoint qualification scheme: target endpoints are qualified by the
        # ``__giql_dj_tgt`` alias (``t``); reference endpoints appear
        # unqualified in the breakpoint CTE and qualified by ``r`` in the
        # coverage ``EXISTS`` filter.
        t_chrom = f't."{target_chrom}"'
        t_start = canonical_start(f't."{target_start}"', target_table)
        t_end = canonical_end(f't."{target_end}"', target_table)
        bp_start = canonical_start(f'"{ref_start}"', ref_table)
        bp_end = canonical_end(f'"{ref_end}"', ref_table)
        r_start = canonical_start(f'r."{ref_start}"', ref_table)
        r_end = canonical_end(f'r."{ref_end}"', ref_table)

        # disjoin_start / disjoin_end are emitted in the target table's
        # coordinate system so an output row carries one convention; the cut
        # math above stays canonical internally.
        out_start = decanonical_start("s.seg_start", target_table)
        out_end = decanonical_end("s.seg_end", target_table)

        # Build the FROM clause for the reference CTE. A subquery reference is
        # aliased with the reserved ``__giql_dj_rs`` so downstream column refs
        # resolve unambiguously; bare references quote the relation name.
        if reference.subquery is not None:
            ref_from = f"{self.sql(reference.subquery)} AS {_DISJOIN_CTE_RS}"
        else:
            ref_from = _quote_relation(reference.relation.name)

        ref_cte = f"{_DISJOIN_CTE_REF} AS (SELECT * FROM {ref_from})"
        tgt_cte = f"{_DISJOIN_CTE_TGT} AS (SELECT * FROM {_quote_relation(target_name)})"
        bp_cte = (
            f"{_DISJOIN_CTE_BP} AS ("
            f'SELECT "{ref_chrom}" AS chrom, {bp_start} AS pos '
            f"FROM {_DISJOIN_CTE_REF} "
            "UNION "
            f'SELECT "{ref_chrom}" AS chrom, {bp_end} AS pos '
            f"FROM {_DISJOIN_CTE_REF})"
        )
        cuts_cte = (
            f"{_DISJOIN_CTE_CUTS} AS ("
            f'SELECT t."{target_chrom}" AS kc, t."{target_start}" AS ks, '
            f't."{target_end}" AS ke, {t_start} AS pos '
            f"FROM {_DISJOIN_CTE_TGT} AS t "
            "UNION "
            f'SELECT t."{target_chrom}", t."{target_start}", t."{target_end}", '
            f"{t_end} FROM {_DISJOIN_CTE_TGT} AS t "
            "UNION "
            f'SELECT t."{target_chrom}", t."{target_start}", t."{target_end}", '
            f"bp.pos FROM {_DISJOIN_CTE_TGT} AS t JOIN {_DISJOIN_CTE_BP} AS bp "
            f"ON bp.chrom = {t_chrom} AND bp.pos > {t_start} "
            f"AND bp.pos < {t_end})"
        )
        segs_cte = (
            f"{_DISJOIN_CTE_SEGS} AS ("
            "SELECT kc, ks, ke, pos AS seg_start, "
            "LEAD(pos) OVER (PARTITION BY kc, ks, ke ORDER BY pos) AS seg_end "
            f"FROM {_DISJOIN_CTE_CUTS})"
        )
        # The ``seg_end > seg_start`` guard is belt-and-suspenders: UNION
        # already dedupes cut positions, so LEAD cannot produce a zero-length
        # segment unless it becomes UNION ALL.
        where_clauses = ["s.seg_end IS NOT NULL", "s.seg_end > s.seg_start"]
        if not reference.is_self_reference:
            # See ``_ResolvedReference.is_self_reference`` for when the
            # coverage EXISTS may be skipped. giql does not emit dialect-
            # specific MATERIALIZED hints, so we keep EXISTS rather than rely
            # on engine-specific materialization guarantees for CTE bodies.
            where_clauses.append(
                f"EXISTS (SELECT 1 FROM {_DISJOIN_CTE_REF} AS r "
                f'WHERE r."{ref_chrom}" = s.kc '
                f"AND {r_start} <= s.seg_start AND {r_end} > s.seg_start)"
            )
        where_sql = " AND ".join(where_clauses)
        final_select = (
            f"SELECT t.*, s.kc AS disjoin_chrom, {out_start} AS disjoin_start, "
            f"{out_end} AS disjoin_end FROM {_DISJOIN_CTE_TGT} AS t "
            f'JOIN {_DISJOIN_CTE_SEGS} AS s ON t."{target_chrom}" = s.kc '
            f'AND t."{target_start}" = s.ks AND t."{target_end}" = s.ke '
            f"WHERE {where_sql}"
        )
        return (
            f"(WITH {ref_cte}, {tgt_cte}, {bp_cte}, "
            f"{cuts_cte}, {segs_cte} {final_select})"
        )

    def giqldistance_sql(self, expression: GIQLDistance) -> str:
        """Generate SQL CASE expression for DISTANCE function.

        :param expression:
            GIQLDistance expression node
        :return:
            SQL CASE expression string calculating genomic distance
        """
        # Extract the two interval arguments
        interval_a = expression.this
        interval_b = expression.args.get("expression")

        stranded = self._extract_bool_param(expression.args.get("stranded"))
        signed = self._extract_bool_param(expression.args.get("signed"))

        # Get SQL representations
        interval_a_sql = self.sql(interval_a)
        interval_b_sql = self.sql(interval_b)

        # Check if we're dealing with column-to-column or column-to-literal
        if "." in interval_a_sql and not interval_a_sql.startswith("'"):
            # Column reference for interval_a
            if stranded:
                chrom_a, start_a, end_a, strand_a = self._get_column_refs(
                    interval_a_sql, None, include_strand=True
                )
            else:
                chrom_a, start_a, end_a = self._get_column_refs(interval_a_sql, None)
                strand_a = None
        else:
            # Literal range - not implemented yet for interval_a
            raise ValueError("Literal range as first argument not yet supported")

        if "." in interval_b_sql and not interval_b_sql.startswith("'"):
            # Column reference for interval_b
            if stranded:
                chrom_b, start_b, end_b, strand_b = self._get_column_refs(
                    interval_b_sql, None, include_strand=True
                )
            else:
                chrom_b, start_b, end_b = self._get_column_refs(interval_b_sql, None)
                strand_b = None
        else:
            # Literal range - not implemented yet
            raise ValueError("Literal range as second argument not yet supported")

        # Distance math below assumes 0-based half-open.
        table_a = self._resolve_table(interval_a_sql)
        table_b = self._resolve_table(interval_b_sql)
        start_a = canonical_start(start_a, table_a)
        end_a = canonical_end(end_a, table_a)
        start_b = canonical_start(start_b, table_b)
        end_b = canonical_end(end_b, table_b)

        # Generate CASE expression
        return self._generate_distance_case(
            chrom_a,
            start_a,
            end_a,
            strand_a,
            chrom_b,
            start_b,
            end_b,
            strand_b,
            stranded=stranded,
            signed=signed,
        )

    def _generate_distance_case(
        self,
        chrom_a: str,
        start_a: str,
        end_a: str,
        strand_a: str | None,
        chrom_b: str,
        start_b: str,
        end_b: str,
        strand_b: str | None,
        stranded: bool = False,
        signed: bool = False,
    ) -> str:
        """Generate SQL CASE expression for distance calculation.

        :param chrom_a:
            Chromosome column for interval A
        :param start_a:
            Start column for interval A
        :param end_a:
            End column for interval A
        :param strand_a:
            Strand column for interval A (None if not stranded)
        :param chrom_b:
            Chromosome column for interval B
        :param start_b:
            Start column for interval B
        :param end_b:
            End column for interval B
        :param strand_b:
            Strand column for interval B (None if not stranded)
        :param stranded:
            Whether to use strand-aware distance calculation
        :param signed:
            Whether to return signed distance (negative for upstream, positive for
            downstream)
        :return:
            SQL CASE expression
        """
        if not stranded or strand_a is None or strand_b is None:
            # Basic distance calculation without strand awareness
            if signed:
                # Signed distance: negative for upstream (B before A),
                # positive for downstream (B after A)
                return (
                    f"CASE WHEN {chrom_a} != {chrom_b} THEN NULL "
                    f"WHEN {start_a} < {end_b} AND {end_a} > {start_b} THEN 0 "
                    f"WHEN {end_a} <= {start_b} THEN ({start_b} - {end_a}) "
                    f"ELSE -({start_a} - {end_b}) END"
                )
            # Unsigned (absolute) distance
            return (
                f"CASE WHEN {chrom_a} != {chrom_b} THEN NULL "
                f"WHEN {start_a} < {end_b} AND {end_a} > {start_b} THEN 0 "
                f"WHEN {end_a} <= {start_b} THEN ({start_b} - {end_a}) "
                f"ELSE ({start_a} - {end_b}) END"
            )

        # Stranded distance calculation
        # Return NULL if either strand is '.', '?', or NULL
        # Calculate distance and multiply by -1 if first interval is on '-' strand
        if signed:
            # Stranded + signed: apply strand flip AND directional sign
            return (
                f"CASE WHEN {chrom_a} != {chrom_b} THEN NULL "
                f"WHEN {strand_a} IS NULL OR {strand_b} IS NULL THEN NULL "
                f"WHEN {strand_a} = '.' OR {strand_a} = '?' THEN NULL "
                f"WHEN {strand_b} = '.' OR {strand_b} = '?' THEN NULL "
                f"WHEN {start_a} < {end_b} AND {end_a} > {start_b} THEN 0 "
                f"WHEN {end_a} <= {start_b} THEN "
                f"CASE WHEN {strand_a} = '-' THEN -({start_b} - {end_a}) "
                f"ELSE ({start_b} - {end_a}) END "
                f"ELSE "
                f"CASE WHEN {strand_a} = '-' THEN ({start_a} - {end_b}) "
                f"ELSE -({start_a} - {end_b}) END END"
            )
        # Stranded but not signed: apply strand flip only
        return (
            f"CASE WHEN {chrom_a} != {chrom_b} THEN NULL "
            f"WHEN {strand_a} IS NULL OR {strand_b} IS NULL THEN NULL "
            f"WHEN {strand_a} = '.' OR {strand_a} = '?' THEN NULL "
            f"WHEN {strand_b} = '.' OR {strand_b} = '?' THEN NULL "
            f"WHEN {start_a} < {end_b} AND {end_a} > {start_b} THEN 0 "
            f"WHEN {end_a} <= {start_b} THEN "
            f"CASE WHEN {strand_a} = '-' THEN -({start_b} - {end_a}) "
            f"ELSE ({start_b} - {end_a}) END "
            f"ELSE "
            f"CASE WHEN {strand_a} = '-' THEN -({start_a} - {end_b}) "
            f"ELSE ({start_a} - {end_b}) END END"
        )

    def _generate_spatial_op(self, expression: exp.Binary, op_type: str) -> str:
        """Generate SQL for a spatial operation.

        :param expression:
            AST node (Intersects, Contains, or Within)
        :param op_type:
            'intersects', 'contains', or 'within'
        :return:
            SQL predicate string
        """
        left = self.sql(expression, "this")
        right_raw = self.sql(expression, "expression")

        # Check if right side is a column reference or a literal range string
        if "." in right_raw and not right_raw.startswith("'"):
            # Column-to-column join (e.g., a.interval INTERSECTS b.interval)
            return self._generate_column_join(left, right_raw, op_type)
        else:
            # Literal range string (e.g., interval INTERSECTS 'chr1:1000-2000')
            try:
                range_str = right_raw.strip("'\"")
                parsed_range = RangeParser.parse(range_str).to_zero_based_half_open()
                return self._generate_range_predicate(left, parsed_range, op_type)
            except Exception as e:
                raise ValueError(
                    f"Could not parse genomic range: {right_raw}. Error: {e}"
                )

    def _generate_range_predicate(
        self,
        column_ref: str,
        parsed_range: ParsedRange,
        op_type: str,
    ) -> str:
        """Generate SQL predicate for a range operation.

        :param column_ref:
            Column reference (e.g., 'v.interval' or 'interval')
        :param parsed_range:
            Parsed genomic range
        :param op_type:
            'intersects', 'contains', or 'within'
        :return:
            SQL predicate string
        """
        # Get column references
        chrom_col, raw_start_col, raw_end_col = self._get_column_refs(
            column_ref, self._current_table
        )
        table = self._resolve_table(column_ref, self._current_table)
        start_col = canonical_start(raw_start_col, table)
        end_col = canonical_end(raw_end_col, table)

        chrom = parsed_range.chromosome
        start = parsed_range.start
        end = parsed_range.end

        if op_type == "intersects":
            # Ranges overlap if: start1 < end2 AND end1 > start2
            return (
                f"({chrom_col} = '{chrom}' "
                f"AND {start_col} < {end} "
                f"AND {end_col} > {start})"
            )

        elif op_type == "contains":
            # Point query: start1 <= point < end1
            if end == start + 1:
                return (
                    f"({chrom_col} = '{chrom}' "
                    f"AND {start_col} <= {start} "
                    f"AND {end_col} > {start})"
                )
            # Range query: start1 <= start2 AND end1 >= end2
            else:
                return (
                    f"({chrom_col} = '{chrom}' "
                    f"AND {start_col} <= {start} "
                    f"AND {end_col} >= {end})"
                )

        elif op_type == "within":
            # Left within right: start1 >= start2 AND end1 <= end2
            return (
                f"({chrom_col} = '{chrom}' "
                f"AND {start_col} >= {start} "
                f"AND {end_col} <= {end})"
            )

        raise ValueError(f"Unknown operation: {op_type}")

    def _generate_column_join(self, left_col: str, right_col: str, op_type: str) -> str:
        """Generate SQL for column-to-column spatial joins.

        :param left_col:
            Left column reference (e.g., 'a.interval')
        :param right_col:
            Right column reference (e.g., 'b.interval')
        :param op_type:
            'intersects', 'contains', or 'within'
        :return:
            SQL predicate string
        """
        # Get column references for both sides
        # Pass None to let _get_column_refs extract and resolve table from column ref
        l_chrom, raw_l_start, raw_l_end = self._get_column_refs(left_col, None)
        r_chrom, raw_r_start, raw_r_end = self._get_column_refs(right_col, None)
        l_table = self._resolve_table(left_col)
        r_table = self._resolve_table(right_col)
        l_start = canonical_start(raw_l_start, l_table)
        l_end = canonical_end(raw_l_end, l_table)
        r_start = canonical_start(raw_r_start, r_table)
        r_end = canonical_end(raw_r_end, r_table)

        if op_type == "intersects":
            # Ranges overlap if: chrom1 = chrom2 AND start1 < end2 AND end1 > start2
            return (
                f"({l_chrom} = {r_chrom} "
                f"AND {l_start} < {r_end} "
                f"AND {l_end} > {r_start})"
            )

        elif op_type == "contains":
            # Left contains right: chrom1 = chrom2 AND start1 <= start2 AND end1 >= end2
            return (
                f"({l_chrom} = {r_chrom} "
                f"AND {l_start} <= {r_start} "
                f"AND {l_end} >= {r_end})"
            )

        elif op_type == "within":
            # Left within right: chrom1 = chrom2 AND start1 >= start2 AND end1 <= end2
            return (
                f"({l_chrom} = {r_chrom} "
                f"AND {l_start} >= {r_start} "
                f"AND {l_end} <= {r_end})"
            )

        raise ValueError(f"Unknown operation: {op_type}")

    def _generate_spatial_set(self, expression: SpatialSetPredicate) -> str:
        """Generate SQL for spatial set predicates (ANY/ALL).

        Examples:
            column INTERSECTS ANY(...) -> (condition1 OR condition2 OR ...)
            column INTERSECTS ALL(...) -> (condition1 AND condition2 AND ...)

        :param expression:
            SpatialSetPredicate expression node
        :return:
            SQL predicate string
        """
        column_ref = self.sql(expression, "this")
        operator = expression.args["operator"]
        quantifier = expression.args["quantifier"]
        ranges = expression.args["ranges"]

        # Parse all ranges
        parsed_ranges = []
        for range_expr in ranges:
            range_str = self.sql(range_expr).strip("'\"")
            parsed_range = RangeParser.parse(range_str).to_zero_based_half_open()
            parsed_ranges.append(parsed_range)

        op_type = operator.lower()

        # Generate conditions for each range
        conditions = []
        for parsed_range in parsed_ranges:
            condition = self._generate_range_predicate(column_ref, parsed_range, op_type)
            conditions.append(condition)

        # Combine with AND (for ALL) or OR (for ANY)
        combinator = " OR " if quantifier.upper() == "ANY" else " AND "
        return "(" + combinator.join(conditions) + ")"

    def _detect_nearest_mode(
        self, expression: GIQLNearest, parent_expression: exp.Expression | None = None
    ) -> str:
        """Detect whether NEAREST is in standalone or correlated mode.

        :param expression:
            GIQLNearest expression node
        :param parent_expression:
            Parent AST node (optional, used to detect LATERAL context)
        :return:
            "standalone" or "correlated"
        """
        # Check if reference parameter is explicitly provided
        reference = expression.args.get("reference")

        if reference:
            # Explicit reference means standalone mode
            return "standalone"

        # No explicit reference - check for LATERAL context
        # In correlated mode, NEAREST appears in a LATERAL join context
        # For now, default to correlated mode if no reference specified
        # (validation will catch missing reference errors later)
        return "correlated"

    def _find_outer_table_in_lateral_join(self, expression: GIQLNearest) -> str | None:
        """Find the outer table name in a LATERAL join context.

        Walks up the AST to find the JOIN clause and extracts the outer table
        that the LATERAL subquery is correlated with.

        :param expression:
            GIQLNearest expression node
        :return:
            Table name or alias of the outer table, or None if not found
        """
        # Walk up the AST to find the JOIN
        current = expression
        while current:
            parent = current.parent
            if not parent:
                break

            # Check if parent is a Lateral expression
            if isinstance(parent, exp.Lateral):
                # Continue up to find the Join
                current = parent
                continue

            # Check if parent is a Join
            if isinstance(parent, exp.Join):
                # The outer table is in the parent Select's FROM clause
                # or in previous joins
                select = parent.parent
                if isinstance(select, exp.Select):
                    # Get the FROM clause
                    from_expr = select.args.get("from_")
                    if from_expr:
                        # Extract table from FROM
                        table_expr = from_expr.this
                        if isinstance(table_expr, exp.Table):
                            # Return alias if it exists, otherwise table name
                            return table_expr.alias or table_expr.name
                        elif isinstance(table_expr, exp.Alias):
                            return table_expr.alias
                break

            current = parent

        return None

    def _resolve_nearest_reference(
        self, expression: GIQLNearest, mode: str
    ) -> tuple[str, str, str] | tuple[str, str, str, str]:
        """Resolve the reference position for NEAREST queries.

        :param expression:
            GIQLNearest expression node
        :param mode:
            "standalone" or "correlated"
        :return:
            Tuple of (chromosome, start, end) or (chromosome, start, end, strand)
            Returns SQL expressions (column refs for correlated, literals for standalone)
            Endpoints are canonicalized to 0-based half-open: literal references via
            :meth:`RangeParser.to_zero_based_half_open`, column references via
            :func:`giql.canonical.canonical_start` /
            :func:`giql.canonical.canonical_end`.
        :raises ValueError:
            If reference is missing in standalone mode or invalid format
        """
        reference = expression.args.get("reference")

        if mode == "standalone":
            if not reference:
                raise ValueError(
                    "NEAREST in standalone mode requires explicit reference parameter"
                )

            # Get SQL representation of reference
            reference_sql = self.sql(reference)

            # Check if it's a literal range string
            if reference_sql.startswith("'") or reference_sql.startswith('"'):
                # Parse literal genomic range
                range_str = reference_sql.strip("'\"")
                try:
                    parsed_range = RangeParser.parse(range_str).to_zero_based_half_open()
                    # Return as SQL literals
                    return (
                        f"'{parsed_range.chromosome}'",
                        str(parsed_range.start),
                        str(parsed_range.end),
                    )
                except Exception as e:
                    raise ValueError(
                        f"Could not parse reference genomic range: "
                        f"{range_str}. Error: {e}"
                    )
            else:
                # Column reference - resolve and canonicalize
                chrom, start, end = self._get_column_refs(reference_sql, None)
                ref_table = self._resolve_table(reference_sql)
                return (
                    chrom,
                    canonical_start(start, ref_table),
                    canonical_end(end, ref_table),
                )

        else:  # correlated mode
            if reference:
                # Explicit reference in correlated mode (e.g., peaks.interval)
                reference_sql = self.sql(reference)
                chrom, start, end = self._get_column_refs(reference_sql, None)
                ref_table = self._resolve_table(reference_sql)
                return (
                    chrom,
                    canonical_start(start, ref_table),
                    canonical_end(end, ref_table),
                )
            else:
                # Implicit reference - resolve from outer table in LATERAL join
                outer_table = self._find_outer_table_in_lateral_join(expression)
                if not outer_table:
                    raise ValueError(
                        "Could not find outer table in LATERAL join context. "
                        "Please specify reference parameter explicitly."
                    )

                # Look up the table to find the genomic column
                # Check if outer_table is an alias
                actual_table = self._alias_to_table.get(outer_table, outer_table)
                table = self.tables.get(actual_table)

                if not table:
                    raise ValueError(
                        f"Outer table '{outer_table}' not found in tables. "
                        "Please specify reference parameter explicitly."
                    )

                # Build column references using the outer table and genomic column
                reference_sql = f"{outer_table}.{table.genomic_col}"
                chrom, start, end = self._get_column_refs(reference_sql, None)
                return (
                    chrom,
                    canonical_start(start, table),
                    canonical_end(end, table),
                )

    def _resolve_cte_or_registered(
        self, name: str, expression: exp.Expression
    ) -> _ResolvedRelation | None:
        """Resolve ``name`` against enclosing CTEs first, then registered tables.

        A CTE shadowing a registered table of the same name yields the
        canonical default columns and ``None`` for the table config; a
        registered table contributes its configured column names and config.

        :param name:
            Bare identifier extracted from a target or reference node.
        :param expression:
            The operator expression; the ancestor walk for CTE collection
            starts from this node.
        :return:
            A :class:`_ResolvedRelation` on a hit, or ``None`` when the name
            matches neither an enclosing CTE nor a registered table.
        """
        if name in _enclosing_cte_names(expression):
            return _ResolvedRelation(name=name, cols=CANONICAL_DEFAULT_COLS, table=None)
        table = self.tables.get(name)
        if table is not None:
            return _ResolvedRelation(
                name=name,
                cols=(table.chrom_col, table.start_col, table.end_col),
                table=table,
            )
        return None

    def _resolve_target(
        self,
        expression: GIQLDisjoin | GIQLNearest,
        *,
        accept_cte: bool,
        operator_label: str,
    ) -> _ResolvedRelation:
        """Resolve a DISJOIN or NEAREST target to a :class:`_ResolvedRelation`.

        When ``accept_cte`` is true (DISJOIN), a name defined by an enclosing
        ``WITH`` shadows a same-named registered table (matching SQL scoping)
        and is assumed to expose canonical 0-based half-open ``chrom`` /
        ``start`` / ``end`` columns — the CTE contract is *unvalidated*: the
        resolver does not inspect the CTE body. When ``accept_cte`` is false
        (NEAREST), a CTE-only match is rejected, and a name matching *both*
        an enclosing CTE and a registered table is rejected as ambiguous —
        the emitted ``FROM "x"`` would be shadowed by the CTE under SQL
        lexical scoping while column references resolve against the
        registered table, so the two interpretations would silently diverge.

        :param expression:
            A :class:`GIQLDisjoin` or :class:`GIQLNearest` expression node.
        :param accept_cte:
            ``True`` for DISJOIN (CTE shadow allowed); ``False`` for NEAREST.
        :param operator_label:
            ``"DISJOIN"`` or ``"NEAREST"`` for operator-specific error
            wording.
        :return:
            The :class:`_ResolvedRelation`. ``table`` is ``None`` for a CTE
            target (DISJOIN only); never ``None`` for NEAREST.
        :raises ValueError:
            If the target node is not a bare name (qualified, subquery,
            etc.); if a bare name resolves to neither a registered table
            nor an enclosing CTE; for NEAREST, if a bare name matches an
            enclosing CTE (rejected) or matches both a registered table and
            an enclosing CTE (ambiguous).
        """
        target_node = expression.this
        try:
            name = _extract_bare_name(target_node)
        except _NotABareName:
            if _is_qualified_name(target_node):
                raise ValueError(
                    f"Target {target_node.sql()!r} is qualified; "
                    "pass a bare table or CTE name."
                ) from None
            raise ValueError(
                f"Target must be a bare table or CTE name; "
                f"got {type(target_node).__name__}: {target_node.sql()!r}."
            ) from None

        if accept_cte:
            resolved = self._resolve_cte_or_registered(name, expression)
            if resolved is None:
                raise ValueError(
                    f"{operator_label} target {name!r} is neither a "
                    "registered table nor a CTE defined in an enclosing "
                    "WITH clause. Pass the relation via "
                    "`transpile(..., tables=[...])` or define the CTE in an "
                    "enclosing WITH clause before transpiling."
                )
            return resolved

        # NEAREST policy: registered tables only.
        table = self.tables.get(name)
        is_cte = name in _enclosing_cte_names(expression)
        if table is not None and is_cte:
            raise ValueError(
                f"{operator_label} target {name!r} resolves to both a "
                "registered table and an enclosing CTE. The engine's SQL "
                "scoping would route NEAREST's rows from the CTE while "
                "column references resolve against the registered table — "
                "rename one of them to disambiguate."
            )
        if table is not None:
            return _ResolvedRelation(
                name=name,
                cols=(table.chrom_col, table.start_col, table.end_col),
                table=table,
            )
        if is_cte:
            raise ValueError(
                f"{operator_label} target {name!r} matches an enclosing CTE. "
                f"{operator_label} does not accept CTE targets; register the "
                "relation as a table via `transpile(..., tables=[...])`, or "
                "rephrase the query against a registered table."
            )
        raise ValueError(
            f"Target table {name!r} not found in tables. "
            "Register the table before transpiling."
        )

    def _resolve_disjoin_reference(
        self,
        expression: GIQLDisjoin,
        target: _ResolvedRelation,
    ) -> _ResolvedReference:
        """Resolve the reference relation for a DISJOIN query.

        The reference contributes the breakpoints at which target intervals
        are cut. When ``reference`` is omitted it defaults to the target set.

        A bare reference name is resolved against the query's CTEs first: a
        name defined by an enclosing ``WITH`` clause shadows a registered
        table of the same name (matching SQL scoping) and is assumed to
        expose canonical 0-based half-open ``chrom`` / ``start`` / ``end``
        columns. A name with no such CTE but a registered table contributes
        that table's column names and coordinate system. A name that is
        neither is an error. A ``(SELECT ...)`` subquery is inlined verbatim
        and assumed to project the canonical columns.

        :param expression:
            GIQLDisjoin expression node.
        :param target:
            The already-resolved target relation (used to default the
            reference and to determine whether the coverage EXISTS may be
            omitted).
        :return:
            A :class:`_ResolvedReference`. For a subquery reference,
            ``subquery`` carries the AST node and ``relation.name`` is
            empty. ``is_self_reference`` is ``True`` only when both sides
            provably resolve to the same registered :class:`Table` — i.e.
            an omitted reference against a registered target, or a bare
            name that maps to the same Table instance.
        :raises ValueError:
            If the reference is not a bare table/CTE name or a
            ``(SELECT ...)`` subquery, or if a bare name resolves to neither
            a registered table nor a CTE.
        """
        reference = expression.args.get("reference")

        # No reference: default to the target set (single-mode).
        if reference is None:
            return _ResolvedReference(
                relation=target,
                subquery=None,
                is_self_reference=target.table is not None,
            )

        # Subquery reference: inline as an aliased derived table. The subquery
        # may textually duplicate the target but proving equivalence is out
        # of scope, so treat it as a distinct relation.
        if isinstance(reference, exp.Subquery):
            return _ResolvedReference(
                relation=_ResolvedRelation(
                    name="", cols=CANONICAL_DEFAULT_COLS, table=None
                ),
                subquery=reference,
                is_self_reference=False,
            )

        if not isinstance(reference, (exp.Table, exp.Column, exp.Identifier)):
            raise ValueError(
                "DISJOIN reference must be a table name, a CTE, or a "
                f"(SELECT ...) subquery; got {type(reference).__name__}: "
                f"{reference.sql() if hasattr(reference, 'sql') else reference!r}"
            )

        try:
            ref_name = _extract_bare_name(reference)
        except _NotABareName:
            if _is_qualified_name(reference):
                raise ValueError(
                    f"Reference {reference.sql()!r} is qualified; "
                    "pass a bare table or CTE name, or a (SELECT ...) "
                    "subquery."
                ) from None
            raise ValueError(
                "Reference must be a bare table or CTE name, or a "
                f"(SELECT ...) subquery; got {type(reference).__name__}: "
                f"{reference.sql()!r}."
            ) from None

        resolved = self._resolve_cte_or_registered(ref_name, expression)
        if resolved is None:
            raise ValueError(
                f"DISJOIN reference {ref_name!r} is neither a registered "
                "table nor a CTE defined in an enclosing WITH clause. Pass "
                "the relation via `transpile(..., tables=[...])` or define "
                "the CTE in an enclosing WITH clause before transpiling. "
                "You may also inline the rows via `reference := (SELECT ...)`."
            )
        is_self_reference = (
            resolved.table is not None
            and resolved.name == target.name
            and resolved.table is target.table
        )
        return _ResolvedReference(
            relation=resolved,
            subquery=None,
            is_self_reference=is_self_reference,
        )

    def _check_reserved_collisions(self, expression: GIQLDisjoin) -> None:
        """Reject identifiers that would shadow DISJOIN's internal
        ``__giql_dj_*`` CTE namespace.

        Three sources are checked: every CTE alias visible in an enclosing
        ``WITH``, the target's bare name, and the reference's bare name
        (subquery references and non-bare reference shapes are skipped —
        their own resolution paths handle them). Surfacing the collision
        here keeps the diagnostic at transpile time; without it the
        duplicate-CTE name reaches the engine as an opaque parse error.

        :raises ValueError:
            If any of the three sources begins with ``__giql_dj_``.
        """
        prefix = _DISJOIN_RESERVED_PREFIX
        for cte_name in _enclosing_cte_names(expression):
            if cte_name.startswith(prefix):
                raise ValueError(
                    f"Enclosing CTE {cte_name!r} uses the reserved "
                    f"{prefix!r} prefix, which names DISJOIN's internal "
                    "CTEs — this CTE collides with DISJOIN's internal "
                    "namespace and would shadow its emitted SQL. Rename "
                    "the CTE before invoking DISJOIN."
                )
        target_name = _bare_name_or_none(expression.this)
        if target_name and target_name.startswith(prefix):
            raise ValueError(
                f"DISJOIN target {target_name!r} uses the reserved "
                f"{prefix!r} prefix, which names the operator's internal "
                "CTEs. Rename the relation."
            )
        ref_node = expression.args.get("reference")
        ref_name = _bare_name_or_none(ref_node) if ref_node is not None else None
        if ref_name and ref_name.startswith(prefix):
            raise ValueError(
                f"DISJOIN reference {ref_name!r} uses the reserved "
                f"{prefix!r} prefix, which names the operator's internal "
                "CTEs. Rename the relation."
            )

    def _get_column_refs(
        self,
        column_ref: str,
        table_name: str | None = None,
        include_strand: bool = False,
    ) -> tuple[str, str, str] | tuple[str, str, str, str]:
        """Get physical column names for genomic data.

        :param column_ref:
            Logical column reference (e.g., 'v.interval' or 'interval')
        :param table_name:
            Table name to look up schema (optional, overrides extraction from column_ref)
        :param include_strand:
            If True, return 4-tuple with strand column; otherwise return 3-tuple
        :return:
            Tuple of (chromosome_col, start_col, end_col) or
            (chromosome_col, start_col, end_col, strand_col) if include_strand=True
        """
        # Default column names
        chrom_col = DEFAULT_CHROM_COL
        start_col = DEFAULT_START_COL
        end_col = DEFAULT_END_COL
        strand_col = DEFAULT_STRAND_COL

        # Alias is kept verbatim for output formatting; table name is resolved
        # separately to look up the Table config (alias != name in joins).
        table_alias = column_ref.rsplit(".", 1)[0] if "." in column_ref else None
        table_name = self._resolve_table_name(column_ref, table_name)

        # Try to get custom column names from table config
        if table_name and self.tables:
            table = self.tables.get(table_name)
            if table:
                chrom_col = table.chrom_col
                start_col = table.start_col
                end_col = table.end_col
                if table.strand_col:
                    strand_col = table.strand_col

        # Format with table alias if present
        if table_alias:
            base_cols = (
                f'{table_alias}."{chrom_col}"',
                f'{table_alias}."{start_col}"',
                f'{table_alias}."{end_col}"',
            )
            if include_strand:
                return base_cols + (f'{table_alias}."{strand_col}"',)
            return base_cols
        else:
            base_cols = (
                f'"{chrom_col}"',
                f'"{start_col}"',
                f'"{end_col}"',
            )
            if include_strand:
                return base_cols + (f'"{strand_col}"',)
            return base_cols

    def _resolve_table_name(self, column_ref: str, table_name: str | None) -> str | None:
        """Resolve the underlying table name for a column reference.

        Precedence: explicit ``table_name`` (caller-provided context) > alias
        map (JOIN-side resolution via ``self._alias_to_table``) >
        ``self._current_table`` (FROM-clause fallback for unmapped dotted
        aliases). Undotted refs return ``None`` because their caller is
        expected to pass ``_current_table`` as the explicit ``table_name``
        when relevant.

        :param column_ref:
            Column reference, possibly aliased (e.g. ``a.interval``)
        :param table_name:
            Explicit table name; takes precedence if non-empty
        :return:
            ``table_name`` if non-empty; otherwise the alias mapping from
            ``self._alias_to_table`` (with ``self._current_table`` as fallback)
            if ``column_ref`` is dotted; otherwise None
        """
        if table_name:
            return table_name
        if "." in column_ref:
            alias = column_ref.rsplit(".", 1)[0]
            return self._alias_to_table.get(alias, self._current_table)
        return None

    def _resolve_table(
        self, column_ref: str, table_name: str | None = None
    ) -> Table | None:
        """Resolve the Table config that backs a column reference.

        :param column_ref:
            Column reference, possibly aliased (e.g. ``a.interval``)
        :param table_name:
            Explicit table name; if omitted, derived from ``column_ref`` alias
            or falls back to ``self._current_table``
        :return:
            The Table config if registered, otherwise None
        """
        table_name = self._resolve_table_name(column_ref, table_name)
        if table_name and self.tables:
            return self.tables.get(table_name)
        return None
