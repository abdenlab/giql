from sqlglot import exp
from sqlglot.generator import Generator

from giql.canonical import decanonical_end
from giql.canonical import decanonical_start
from giql.expressions import Contains
from giql.expressions import GIQLDisjoin
from giql.expressions import GIQLNearest
from giql.expressions import Intersects
from giql.expressions import SpatialSetPredicate
from giql.expressions import Within
from giql.range_parser import ParsedRange
from giql.range_parser import RangeParser
from giql.resolver import META_KEY
from giql.resolver import OperatorResolution
from giql.resolver import ResolvedColumn
from giql.resolver import ResolvedInterval
from giql.resolver import ResolvedRef
from giql.table import Table
from giql.table import Tables


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

    def __init__(self, tables: Tables | None = None, **kwargs):
        super().__init__(**kwargs)
        self.tables = tables or Tables()

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
            GIQLNearest expression node
        :return:
            SQL string for NEAREST operation
        """
        # Detect mode
        mode = self._detect_nearest_mode(expression)

        # Unpack the resolution metadata attached by ResolveOperatorRefs (pass 1).
        resolution = self._nearest_resolution(expression)

        # Target (already a registered-table ResolvedRef from the pass). An
        # unresolved target means it is not a registered table; raise the
        # historical diagnostic.
        target_ref = resolution.slot("this") if resolution is not None else None
        if not isinstance(target_ref, ResolvedRef):
            target = expression.this
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
        target_chrom, target_start, target_end = target_ref.cols

        # The target's *declared* encoding, which the passed-through target row
        # (SELECT {table_name}.*) must round-trip back into. CanonicalizeCoordinates
        # (pass 2) preserves it on the resolution when it wraps a non-canonical
        # target in a __giql_canon_* CTE (the slot's own Table is then None); a
        # canonical target is left unwrapped and its slot Table carries the
        # (identity) encoding. The synthesized `distance` column is encoding-
        # invariant (a count of bases) and needs no round-trip.
        output_table = self._nearest_output_encoding(expression, target_ref)
        passthrough = self._nearest_passthrough(
            table_name, target_start, target_end, output_table
        )

        # Reference interval (a ResolvedInterval from the pass). An unresolved
        # reference re-raises the generator's historical diagnostic. Input
        # canonicalization is owned by CanonicalizeCoordinates (pass 2, issue
        # #123): a literal range is already canonical, and a column / implicit-
        # outer reference's endpoints are canonicalized in place by the pass, so
        # the emitter consumes the fragments verbatim with no canonicalization.
        ref = resolution.slot("reference")
        if not isinstance(ref, ResolvedInterval):
            self._raise_nearest_reference_error(expression, mode, resolution)
        ref_chrom, ref_start, ref_end = ref.chrom, ref.start, ref.end

        # Extract parameters
        k = expression.args.get("k")
        k_value = int(str(k)) if k else 1  # Default k=1

        max_distance = expression.args.get("max_distance")
        max_dist_value = int(str(max_distance)) if max_distance else None

        is_stranded = self._extract_bool_param(expression.args.get("stranded"))
        is_signed = self._extract_bool_param(expression.args.get("signed"))

        # Resolve strand columns if stranded mode. The reference strand is
        # carried on the resolved interval (a literal's strand, an explicit
        # column's strand, or the outer table's strand for an implicit
        # reference — already gated to preserve the historical divergence).
        ref_strand = None
        target_strand = None
        if is_stranded:
            ref_strand = ref.strand
            # When pass 2 wraps a non-canonical target its slot Table is blanked,
            # so the strand column name comes from the *declared* encoding the
            # pass preserved (output_table). The canon CTE's SELECT * REPLACE
            # passes the strand column through unchanged under its physical name,
            # so the qualifier stays the relation NEAREST selects from.
            if output_table and output_table.strand_col:
                target_strand = f'{table_name}."{output_table.strand_col}"'

        # Distance math below assumes 0-based half-open. Input canonicalization
        # is owned by CanonicalizeCoordinates (pass 2, issue #123): a
        # non-canonical target is rewritten to a canonical __giql_canon_* CTE
        # before generation (table_name then names the CTE), so the target
        # endpoints are consumed verbatim with no in-emitter canonicalization. The
        # output round-trip of the passed-through target row stays here (see the
        # SELECT projection below).
        target_start_expr = f'{table_name}."{target_start}"'
        target_end_expr = f'{table_name}."{target_end}"'

        # Build distance calculation using CASE expression
        # For NEAREST: ORDER BY absolute distance, but RETURN signed distance
        distance_expr = self._generate_distance_case(
            ref_chrom,
            ref_start,
            ref_end,
            ref_strand,
            f'{table_name}."{target_chrom}"',
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
            f'{ref_chrom} = {table_name}."{target_chrom}"'  # Chromosome pre-filter
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
                SELECT {passthrough}, {distance_expr} AS distance
                FROM {table_name}
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
                SELECT {passthrough}, {distance_expr} AS distance
                FROM {table_name}
                WHERE {where_sql}
                ORDER BY {abs_distance_expr}
                LIMIT {k_value}
            )"""

        return sql.strip()

    def _nearest_output_encoding(
        self, expression: GIQLNearest, target_ref: ResolvedRef
    ) -> Table | None:
        """Return the target's declared encoding for NEAREST's row passthrough.

        ``CanonicalizeCoordinates`` (pass 2) records the original
        :class:`~giql.table.Table` on the resolution when it wraps a non-canonical
        target in a ``__giql_canon_*`` CTE (blanking the slot's own ``table``).
        For an unwrapped target — a canonical registered table, or any target when
        the pass did not run — the slot's own ``table`` carries the (identity)
        encoding.

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
        self,
        table_name: str,
        target_start: str,
        target_end: str,
        output_table: Table | None,
    ) -> str:
        """Project the target's full row, de-canonicalizing the interval columns.

        NEAREST passes the whole target row through (``SELECT {table_name}.*``)
        alongside the synthesized, encoding-invariant ``distance`` column. When the
        target's declared encoding is canonical 0-based half-open the row passes
        through as a plain ``{table_name}.*`` — the byte-identical identity fast
        path. When it is non-canonical the interval columns, canonical inside the
        ``__giql_canon_*`` CTE the target was rewritten to, are de-canonicalized
        back into that encoding via a star ``REPLACE`` so the passed-through
        interval matches the target's own convention. (Only non-canonical targets
        are wrapped, so the ``REPLACE`` appears only where a canonical CTE already
        shapes the SQL.)

        :param table_name:
            The relation the row is selected from (the canon CTE name when wrapped,
            else the registered table name) — also the column qualifier.
        :param target_start:
            Physical start column name
        :param target_end:
            Physical end column name
        :param output_table:
            The target's declared :class:`~giql.table.Table`, or ``None``
        :return:
            The passthrough projection fragment (``{table_name}.*`` or a star
            ``REPLACE``)
        """
        if output_table is None or (
            output_table.coordinate_system == "0based"
            and output_table.interval_type == "half_open"
        ):
            return f"{table_name}.*"
        pt_start = decanonical_start(f'{table_name}."{target_start}"', output_table)
        pt_end = decanonical_end(f'{table_name}."{target_end}"', output_table)
        return (
            f"{table_name}.* REPLACE "
            f'({pt_start} AS "{target_start}", {pt_end} AS "{target_end}")'
        )

    def giqldisjoin_sql(self, expression: GIQLDisjoin) -> str:
        """Generate SQL for the DISJOIN table function.

        DISJOIN splits each target interval at every reference breakpoint
        strictly interior to it. The full target row passes through unchanged;
        the sub-interval is appended as ``disjoin_chrom`` / ``disjoin_start`` /
        ``disjoin_end`` in the target table's coordinate system. A coverage
        filter drops sub-intervals overlapping no reference interval. When no
        ``reference`` is given it defaults to the target set.

        Input canonicalization is owned by ``CanonicalizeCoordinates`` (pass 2,
        issue #122): every non-canonical interval-bearing operand is rewritten to
        a canonical ``__giql_canon_*`` CTE before generation, so this emitter
        consumes already-canonical 0-based half-open columns and applies no
        in-emitter canonicalization arithmetic. The output round-trip back to the
        target's declared encoding stays here — the ``disjoin_*`` columns and the
        passed-through interval are synthesized at generation time and cannot be
        reached by a pass-2 outermost-projection rewrite — driven by the original
        encoding the pass preserves on the resolution.

        :param expression:
            GIQLDisjoin expression node
        :return:
            SQL string (a parenthesized WITH-CTE subquery) for the DISJOIN table
        """
        # Unpack the resolution metadata attached by ResolveOperatorRefs (pass 1)
        # and rewritten by CanonicalizeCoordinates (pass 2).
        target_ref, ref, ref_from = self._disjoin_resolution(expression)
        target_name = target_ref.name
        target_chrom, target_start, target_end = target_ref.cols
        ref_chrom, ref_start, ref_end = ref.cols
        is_self_reference = ref.coverage_skippable

        # The target's *declared* encoding, which disjoin_* output and the
        # passed-through interval must round-trip back into. Pass 2 preserves it
        # on the resolution when it wraps a non-canonical target (the slot's own
        # Table is then None); a canonical target is left unwrapped and its slot
        # Table carries the (identity) encoding.
        output_table = self._disjoin_output_encoding(expression, target_ref)

        # Post-pass every operand is canonical 0-based half-open (a registered
        # table is either identity-encoded or rewritten to a canonical CTE), so
        # the physical columns are consumed verbatim with no canonicalization.
        t_chrom = f't."{target_chrom}"'
        t_start = f't."{target_start}"'
        t_end = f't."{target_end}"'

        # Reference endpoints: unqualified for the breakpoint CTE, qualified by
        # 'r' for the coverage EXISTS filter.
        bp_start = f'"{ref_start}"'
        bp_end = f'"{ref_end}"'
        r_start = f'r."{ref_start}"'
        r_end = f'r."{ref_end}"'

        # disjoin_start / disjoin_end are emitted in the target's declared
        # coordinate system so an output row carries one convention; the cut math
        # stays canonical internally.
        out_start = decanonical_start("s.seg_start", output_table)
        out_end = decanonical_end("s.seg_end", output_table)
        passthrough = self._disjoin_passthrough(target_start, target_end, output_table)

        # Build the WITH clause one named fragment per __giql_dj_* CTE so each
        # block reads on its own. The `seg_end > seg_start` guard in the final
        # WHERE is belt-and-suspenders: UNION already dedupes cut positions, so
        # LEAD cannot produce a zero-length segment unless it becomes UNION ALL.
        ref_cte = f"__giql_dj_ref AS (SELECT * FROM {ref_from})"
        tgt_cte = f"__giql_dj_tgt AS (SELECT * FROM {target_name})"
        bp_cte = (
            "__giql_dj_bp AS ("
            f'SELECT "{ref_chrom}" AS chrom, {bp_start} AS pos FROM __giql_dj_ref '
            "UNION "
            f'SELECT "{ref_chrom}" AS chrom, {bp_end} AS pos FROM __giql_dj_ref)'
        )
        cuts_cte = (
            "__giql_dj_cuts AS ("
            f'SELECT t."{target_chrom}" AS kc, t."{target_start}" AS ks, '
            f't."{target_end}" AS ke, {t_start} AS pos FROM __giql_dj_tgt AS t '
            "UNION "
            f'SELECT t."{target_chrom}", t."{target_start}", t."{target_end}", '
            f"{t_end} FROM __giql_dj_tgt AS t "
            "UNION "
            f'SELECT t."{target_chrom}", t."{target_start}", t."{target_end}", '
            "bp.pos FROM __giql_dj_tgt AS t JOIN __giql_dj_bp AS bp "
            f"ON bp.chrom = {t_chrom} AND bp.pos > {t_start} "
            f"AND bp.pos < {t_end})"
        )
        segs_cte = (
            "__giql_dj_segs AS ("
            "SELECT kc, ks, ke, pos AS seg_start, "
            "LEAD(pos) OVER (PARTITION BY kc, ks, ke ORDER BY pos) AS seg_end "
            "FROM __giql_dj_cuts)"
        )
        # In self-reference mode the coverage EXISTS is provably always true:
        # every emitted segment lies inside its parent target row, and that
        # row is itself a member of the reference set. Skip the clause so the
        # planner does not waste work on a no-op semi-join. The __giql_dj_ref
        # CTE itself stays live because __giql_dj_bp still draws breakpoints
        # from it.
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

    def _disjoin_output_encoding(
        self, expression: GIQLDisjoin, target_ref: ResolvedRef
    ) -> Table | None:
        """Return the target's declared encoding for DISJOIN's output round-trip.

        ``CanonicalizeCoordinates`` (pass 2) records the original
        :class:`~giql.table.Table` on the resolution when it wraps a non-canonical
        target (blanking the slot's own ``table``). For an unwrapped target — a
        canonical registered table, or any target when the pass did not run — the
        slot's own ``table`` carries the (identity) encoding.

        :param expression:
            GIQLDisjoin expression node
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

    def _disjoin_passthrough(
        self, target_start: str, target_end: str, output_table: Table | None
    ) -> str:
        """Project the target's full row, de-canonicalizing the interval columns.

        When the target's declared encoding is canonical 0-based half-open the
        row passes through as a plain ``t.*`` — the byte-identical identity fast
        path. When it is non-canonical the interval columns, canonical inside
        ``__giql_dj_tgt``, are de-canonicalized back into that encoding via a star
        ``REPLACE`` so the passed-through interval matches the target's own
        convention. (Only non-canonical targets are wrapped, so the ``REPLACE``
        appears only where a canonical CTE already shapes the SQL.)

        :param target_start:
            Physical start column name
        :param target_end:
            Physical end column name
        :param output_table:
            The target's declared :class:`~giql.table.Table`, or ``None``
        :return:
            The passthrough projection fragment (``t.*`` or a star ``REPLACE``)
        """
        if output_table is None or (
            output_table.coordinate_system == "0based"
            and output_table.interval_type == "half_open"
        ):
            return "t.*"
        pt_start = decanonical_start(f't."{target_start}"', output_table)
        pt_end = decanonical_end(f't."{target_end}"', output_table)
        return (
            f't.* REPLACE ({pt_start} AS "{target_start}", {pt_end} AS "{target_end}")'
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

        .. note::

           KEEP IN SYNC: this method and the AST builder in
           ``giql.expanders.distance`` (``expand_distance``) produce the *same*
           distance CASE by two routes. DISTANCE itself migrated to the expander
           (epic #137, issue #140); this method survives only because NEAREST
           still calls it for its ORDER BY / filter math. Once NEAREST migrates,
           delete this method and retire the duplication. Until then, any change
           to the distance math here must be mirrored in the expander (and vice
           versa); the parity test in tests/test_distance_udf.py guards drift.

        Distances follow bedtools ``closest -d`` semantics: overlapping
        intervals report ``0``, book-ended (adjacent) intervals where
        ``A.end == B.start`` in half-open coordinates report ``1``, and a raw
        half-open gap of N bases reports ``N + 1``. The ``+ 1`` is applied to
        the absolute gap magnitude before any directional sign, so a downstream
        book-ended pair reports ``+1`` and an upstream one ``-1`` in signed mode.

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
                    f"WHEN {end_a} <= {start_b} THEN ({start_b} - {end_a} + 1) "
                    f"ELSE -({start_a} - {end_b} + 1) END"
                )
            # Unsigned (absolute) distance
            return (
                f"CASE WHEN {chrom_a} != {chrom_b} THEN NULL "
                f"WHEN {start_a} < {end_b} AND {end_a} > {start_b} THEN 0 "
                f"WHEN {end_a} <= {start_b} THEN ({start_b} - {end_a} + 1) "
                f"ELSE ({start_a} - {end_b} + 1) END"
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
                f"CASE WHEN {strand_a} = '-' THEN -({start_b} - {end_a} + 1) "
                f"ELSE ({start_b} - {end_a} + 1) END "
                f"ELSE "
                f"CASE WHEN {strand_a} = '-' THEN ({start_a} - {end_b} + 1) "
                f"ELSE -({start_a} - {end_b} + 1) END END"
            )
        # Stranded but not signed: apply strand flip only
        return (
            f"CASE WHEN {chrom_a} != {chrom_b} THEN NULL "
            f"WHEN {strand_a} IS NULL OR {strand_b} IS NULL THEN NULL "
            f"WHEN {strand_a} = '.' OR {strand_a} = '?' THEN NULL "
            f"WHEN {strand_b} = '.' OR {strand_b} = '?' THEN NULL "
            f"WHEN {start_a} < {end_b} AND {end_a} > {start_b} THEN 0 "
            f"WHEN {end_a} <= {start_b} THEN "
            f"CASE WHEN {strand_a} = '-' THEN -({start_b} - {end_a} + 1) "
            f"ELSE ({start_b} - {end_a} + 1) END "
            f"ELSE "
            f"CASE WHEN {strand_a} = '-' THEN -({start_a} - {end_b} + 1) "
            f"ELSE ({start_a} - {end_b} + 1) END END"
        )

    def _predicate_operand(self, expression: exp.Expression, arg: str) -> ResolvedColumn:
        """Return the :class:`ResolvedColumn` for a spatial predicate operand.

        Reads the column resolution attached to *expression* by the
        ``ResolveOperatorRefs`` pass (pass 1). The emitter consumes only the
        resolved metadata; all name/column resolution lives in the pass.

        :param expression:
            The spatial predicate node carrying the resolution metadata.
        :param arg:
            The operand slot key (``"this"`` or ``"expression"``).
        :return:
            The resolved column metadata.
        """
        resolution = expression.meta.get(META_KEY)
        if isinstance(resolution, OperatorResolution):
            resolved = resolution.column(arg)
            if resolved is not None:
                return resolved

        raise ValueError(
            f"Spatial predicate operand {arg!r} was not resolved; run the "
            "ResolveOperatorRefs pass (transpile pipeline) before generation."
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
        right_raw = self.sql(expression, "expression")

        # Check if right side is a column reference or a literal range string
        if "." in right_raw and not right_raw.startswith("'"):
            # Column-to-column join (e.g., a.interval INTERSECTS b.interval)
            left = self._predicate_operand(expression, "this")
            right = self._predicate_operand(expression, "expression")
            return self._generate_column_join(left, right, op_type)
        else:
            # Literal range string (e.g., interval INTERSECTS 'chr1:1000-2000')
            try:
                range_str = right_raw.strip("'\"")
                parsed_range = RangeParser.parse(range_str).to_zero_based_half_open()
                left = self._predicate_operand(expression, "this")
                return self._generate_range_predicate(left, parsed_range, op_type)
            except Exception as e:
                raise ValueError(
                    f"Could not parse genomic range: {right_raw}. Error: {e}"
                )

    def _generate_range_predicate(
        self,
        column: ResolvedColumn,
        parsed_range: ParsedRange,
        op_type: str,
    ) -> str:
        """Generate SQL predicate for a range operation.

        :param column:
            Resolved column operand (physical chrom/start/end fragments plus the
            backing :class:`~giql.table.Table` config for canonicalization).
        :param parsed_range:
            Parsed genomic range
        :param op_type:
            'intersects', 'contains', or 'within'
        :return:
            SQL predicate string
        """
        # The alias-qualified column fragments come pre-resolved on the
        # ResolvedColumn, already canonicalized to 0-based half-open by
        # CanonicalizeCoordinates (pass 2, issue #123). The predicate returns a
        # boolean, which is encoding-invariant, so no output de-canonicalization
        # is needed.
        chrom_col = column.chrom
        start_col = column.start
        end_col = column.end

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

    def _generate_column_join(
        self, left: ResolvedColumn, right: ResolvedColumn, op_type: str
    ) -> str:
        """Generate SQL for column-to-column spatial joins.

        :param left:
            Resolved left column operand (e.g., for 'a.interval').
        :param right:
            Resolved right column operand (e.g., for 'b.interval').
        :param op_type:
            'intersects', 'contains', or 'within'
        :return:
            SQL predicate string
        """
        # The alias-qualified chrom/start/end fragments come pre-resolved on the
        # ResolvedColumns, already canonicalized to 0-based half-open by
        # CanonicalizeCoordinates (pass 2, issue #123). The predicate returns a
        # boolean (encoding-invariant), so no output de-canonicalization is needed.
        l_chrom = left.chrom
        r_chrom = right.chrom
        l_start = left.start
        l_end = left.end
        r_start = right.start
        r_end = right.end

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
        operator = expression.args["operator"]
        quantifier = expression.args["quantifier"]
        ranges = expression.args["ranges"]

        # Resolve the (single) left column operand once; every range condition
        # compares against the same column. The set predicate's ranges are
        # always literals, so only this operand needs resolution.
        column = self._predicate_operand(expression, "this")

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
            condition = self._generate_range_predicate(column, parsed_range, op_type)
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

    def _nearest_resolution(self, expression: GIQLNearest) -> OperatorResolution | None:
        """Return the NEAREST resolution attached by ResolveOperatorRefs (pass 1).

        The transpile pipeline attaches an
        :class:`~giql.resolver.OperatorResolution` before generation, and it
        survives the generator's defensive tree copy. The emitter reads only the
        attached metadata; resolution lives entirely in the pass.

        :param expression:
            GIQLNearest expression node
        :return:
            The attached :class:`~giql.resolver.OperatorResolution`, or ``None``
            if resolution did not produce one.
        """
        resolution = expression.meta.get(META_KEY)
        return resolution if isinstance(resolution, OperatorResolution) else None

    def _raise_nearest_reference_error(
        self,
        expression: GIQLNearest,
        mode: str,
        resolution: OperatorResolution | None,
    ) -> None:
        """Raise the historical diagnostic for an unresolved NEAREST reference.

        ResolveOperatorRefs (pass 1) defers a reference slot it cannot resolve;
        this re-raises the generator's pre-pass error verbatim. The implicit-
        outer failures rely on the :class:`~giql.resolver.SlotDeferral` the pass
        records (the ancestor walk that distinguished them now lives in the
        pass); a literal-range parse failure is reproduced by re-parsing.

        :param expression:
            GIQLNearest expression node
        :param mode:
            "standalone" or "correlated"
        :param resolution:
            The attached resolution metadata, if any
        :raises ValueError:
            Always — with the matching historical message
        """
        reference = expression.args.get("reference")

        if reference is None:
            # Implicit-outer reference: standalone mode (unreachable in practice,
            # since an absent reference is classified as correlated) keeps the
            # historical message; otherwise consult the recorded deferral.
            if mode == "standalone":
                raise ValueError(
                    "NEAREST in standalone mode requires explicit reference parameter"
                )
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
        reference_sql = self.sql(reference)
        range_str = reference_sql.strip("'\"")
        try:
            RangeParser.parse(range_str).to_zero_based_half_open()
        except Exception as e:
            raise ValueError(
                f"Could not parse reference genomic range: {range_str}. Error: {e}"
            )
        raise ValueError(f"Could not parse reference genomic range: {range_str}.")

    def _disjoin_resolution(
        self, expression: GIQLDisjoin
    ) -> tuple[ResolvedRef, ResolvedRef, str]:
        """Unpack the DISJOIN resolution attached by ResolveOperatorRefs (pass 1).

        Reads the :class:`~giql.resolver.OperatorResolution` from
        ``expression.meta`` and returns ``(target_ref, ref, ref_from)`` where
        ``ref_from`` is the text following ``FROM`` inside the ``__giql_dj_ref``
        CTE. A subquery reference carries no name, so it is rendered from the
        AST node as an aliased derived table; registered tables and CTEs are
        selected from by name.

        The resolver pass deliberately leaves unresolvable slots unresolved
        (unregistered target; unsupported reference node type; reference name
        using the reserved ``__giql_dj_`` prefix or matching neither a
        registered table nor an in-query CTE). For those, and for a target name
        using the reserved prefix (which the resolver does resolve), this
        re-raises the generator's historical diagnostics verbatim.

        :param expression:
            GIQLDisjoin expression node
        :return:
            Tuple of ``(target_ref, ref, ref_from)``
        :raises ValueError:
            If the target or reference slot is unresolved, or a resolved name
            uses the reserved ``__giql_dj_`` prefix.
        """
        resolution = expression.meta.get(META_KEY)
        target_ref = (
            resolution.slot("this")
            if isinstance(resolution, OperatorResolution)
            else None
        )

        # An unresolved target means it is not a registered table.
        if target_ref is None:
            target = expression.this
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

        # The __giql_dj_ prefix names the operator's internal CTEs; a target
        # table using it would collide with them.
        if target_ref.name.startswith("__giql_dj_"):
            raise ValueError(
                f"DISJOIN target {target_ref.name!r} uses the reserved "
                "'__giql_dj_' prefix, which names the operator's internal "
                "CTEs. Rename the table."
            )

        reference = expression.args.get("reference")
        ref = resolution.slot("reference")
        if ref is not None:
            if ref.kind == "subquery":
                return target_ref, ref, f"{self.sql(reference)} AS __giql_dj_rs"
            return target_ref, ref, ref.name

        # Unresolved reference: re-classify it and raise the matching
        # historical diagnostic. An omitted reference always resolves (to the
        # target set), so reference is non-None here.
        if not isinstance(reference, (exp.Table, exp.Column, exp.Identifier)):
            raise ValueError(
                "DISJOIN reference must be a table name, a CTE, or a "
                f"(SELECT ...) subquery; got {type(reference).__name__}: "
                f"{reference}"
            )
        ref_name = reference.name
        if ref_name.startswith("__giql_dj_"):
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

    @staticmethod
    def _extract_bool_param(param_expr: exp.Expression | None) -> bool:
        """Extract boolean value from a parameter expression.

        Handles :class:`exp.Boolean`, :class:`exp.Literal`, and string
        representations.
        """
        if not param_expr:
            return False
        elif isinstance(param_expr, exp.Boolean):
            return param_expr.this
        else:
            return str(param_expr).upper() in ("TRUE", "1", "YES")
