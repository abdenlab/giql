from sqlglot import exp
from sqlglot.generator import Generator

from giql.canonical import decanonical_end
from giql.canonical import decanonical_start
from giql.dialect import GIQLDialect
from giql.expressions import GIQLNearest
from giql.range_parser import RangeParser
from giql.resolver import META_KEY
from giql.resolver import OperatorResolution
from giql.resolver import ResolvedColumn
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

    def __init__(self, tables: Tables | None = None, **kwargs):
        super().__init__(**kwargs)
        self.tables = tables or Tables()

    # INTERSECTS / CONTAINS / WITHIN and the ANY/ALL set predicates are migrated
    # to the ExpandOperators registry (#141): they expand to standard boolean AST
    # in ``giql.expanders.intersects`` before generation, so the generator no
    # longer carries ``intersects_sql`` / ``contains_sql`` / ``within_sql`` /
    # ``spatialsetpredicate_sql`` emitters or their ``_generate_spatial_*`` /
    # ``_predicate_operand`` helpers.

    @staticmethod
    def _nearest_output_encoding(
        expression: GIQLNearest, target_ref: ResolvedRef
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

    @staticmethod
    def _nearest_passthrough(
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
        # TODO(#142): this emits an unconditional ``* REPLACE`` (DuckDB-only).
        # When DataFusion gains correlated LATERAL, adopt the capability branch the
        # DISJOIN expander uses (``giql.expanders.disjoin._disjoin_passthrough``):
        # ``* REPLACE`` where ``supports_star_replace`` holds, the portable
        # ``* EXCEPT`` form otherwise, so a non-canonical NEAREST passthrough runs
        # on the DataFusion family too.
        return (
            f"{table_name}.* REPLACE "
            f'({pt_start} AS "{target_start}", {pt_end} AS "{target_end}")'
        )

    @staticmethod
    def _generate_distance_case(
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

    @staticmethod
    def _detect_nearest_mode(
        expression: GIQLNearest, parent_expression: exp.Expression | None = None
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

    @staticmethod
    def _raise_nearest_reference_error(
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
        reference_sql = reference.sql(dialect=GIQLDialect)
        range_str = reference_sql.strip("'\"")
        try:
            RangeParser.parse(range_str).to_zero_based_half_open()
        except Exception as e:
            raise ValueError(
                f"Could not parse reference genomic range: {range_str}. Error: {e}"
            )
        raise ValueError(f"Could not parse reference genomic range: {range_str}.")

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
