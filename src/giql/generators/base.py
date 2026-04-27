from sqlglot import exp
from sqlglot.generator import Generator

from giql.constants import DEFAULT_CHROM_COL
from giql.constants import DEFAULT_END_COL
from giql.constants import DEFAULT_START_COL
from giql.constants import DEFAULT_STRAND_COL
from giql.expressions import Contains
from giql.expressions import GIQLDistance
from giql.expressions import GIQLNearest
from giql.expressions import Intersects
from giql.expressions import SpatialSetPredicate
from giql.expressions import Within
from giql.range_parser import ParsedRange
from giql.range_parser import RangeParser
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
            GIQLNearest expression node
        :return:
            SQL string for NEAREST operation
        """
        # Detect mode
        mode = self._detect_nearest_mode(expression)

        # Resolve target table
        table_name, (target_chrom, target_start, target_end) = (
            self._resolve_target_table(expression)
        )

        # Resolve reference
        ref_chrom, ref_start, ref_end = self._resolve_nearest_reference(expression, mode)

        # Extract parameters
        k = expression.args.get("k")
        k_value = int(str(k)) if k else 1  # Default k=1

        max_distance = expression.args.get("max_distance")
        max_dist_value = int(str(max_distance)) if max_distance else None

        is_stranded = self._extract_bool_param(expression.args.get("stranded"))
        is_signed = self._extract_bool_param(expression.args.get("signed"))

        target_table = self.tables.get(table_name) if self.tables else None

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
                target_strand = f'{table_name}."{target_table.strand_col}"'

        # Distance math below assumes 0-based half-open.
        target_start_expr = self._canonical_start(
            f'{table_name}."{target_start}"', target_table
        )
        target_end_expr = self._canonical_end(
            f'{table_name}."{target_end}"', target_table
        )

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
                SELECT {table_name}.*, {distance_expr} AS distance
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
                SELECT {table_name}.*, {distance_expr} AS distance
                FROM {table_name}
                WHERE {where_sql}
                ORDER BY {abs_distance_expr}
                LIMIT {k_value}
            )"""

        return sql.strip()

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
        start_a = self._canonical_start(start_a, table_a)
        end_a = self._canonical_end(end_a, table_a)
        start_b = self._canonical_start(start_b, table_b)
        end_b = self._canonical_end(end_b, table_b)

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
        start_col = self._canonical_start(raw_start_col, table)
        end_col = self._canonical_end(raw_end_col, table)

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
        l_start = self._canonical_start(raw_l_start, l_table)
        l_end = self._canonical_end(raw_l_end, l_table)
        r_start = self._canonical_start(raw_r_start, r_table)
        r_end = self._canonical_end(raw_r_end, r_table)

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

    def _find_outer_table_in_lateral_join(
        self, expression: GIQLNearest
    ) -> str | None:
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
            ``RangeParser.to_zero_based_half_open``, column references via
            ``_canonical_start`` / ``_canonical_end``.
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
                return self._canonical_endpoints(chrom, start, end, ref_table)

        else:  # correlated mode
            if reference:
                # Explicit reference in correlated mode (e.g., peaks.interval)
                reference_sql = self.sql(reference)
                chrom, start, end = self._get_column_refs(reference_sql, None)
                ref_table = self._resolve_table(reference_sql)
                return self._canonical_endpoints(chrom, start, end, ref_table)
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
                return self._canonical_endpoints(chrom, start, end, table)

    def _resolve_target_table(
        self, expression: GIQLNearest
    ) -> tuple[str, tuple[str, str, str]]:
        """Resolve the target table name and its genomic column references.

        :param expression:
            GIQLNearest expression node
        :return:
            Tuple of (table_name, (chromosome_col, start_col, end_col))
        :raises ValueError:
            If target table is not found or doesn't have genomic columns
        """
        # Extract target table from 'this' argument
        target = expression.this

        if isinstance(target, exp.Table):
            table_name = target.name
        elif isinstance(target, exp.Column):
            # If it's a column reference, extract table name
            table_name = target.table if target.table else str(target.this)
        else:
            # Try to extract as string
            table_name = str(target)

        table = self.tables.get(table_name)
        if not table:
            raise ValueError(
                f"Target table '{table_name}' not found in tables. "
                "Register the table before transpiling."
            )

        # Get physical column names from table config
        return table_name, (table.chrom_col, table.start_col, table.end_col)

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

    @staticmethod
    def _canonical_start(raw_start: str, table: Table | None) -> str:
        """Wrap a raw start column expression to yield canonical 0-based half-open start.

        Conversion by ``coordinate_system``:

        - 0based: ``start`` (identity)
        - 1based: ``start - 1``
        """
        if table is None or table.coordinate_system == "0based":
            return raw_start
        return f"({raw_start} - 1)"

    @staticmethod
    def _canonical_end(raw_end: str, table: Table | None) -> str:
        """Wrap a raw end column expression to yield canonical 0-based half-open end.

        Conversion by ``(coordinate_system, interval_type)``:

        - 0based / half_open: ``end`` (identity)
        - 0based / closed:    ``end + 1``
        - 1based / half_open: ``end - 1``
        - 1based / closed:    ``end`` (identity)
        """
        if table is None:
            return raw_end
        key = (table.coordinate_system, table.interval_type)
        if key == ("0based", "closed"):
            return f"({raw_end} + 1)"
        if key == ("1based", "half_open"):
            return f"({raw_end} - 1)"
        return raw_end  # 0based/half_open and 1based/closed: identity

    @staticmethod
    def _canonical_endpoints(
        chrom: str,
        raw_start: str,
        raw_end: str,
        table: Table | None,
    ) -> tuple[str, str, str]:
        """Canonicalize a ``(chrom, start, end)`` triple to 0-based half-open.

        The chromosome expression passes through unchanged; start and end are
        wrapped via ``_canonical_start`` / ``_canonical_end``.
        """
        return (
            chrom,
            BaseGIQLGenerator._canonical_start(raw_start, table),
            BaseGIQLGenerator._canonical_end(raw_end, table),
        )

    def _resolve_table_name(
        self, column_ref: str, table_name: str | None
    ) -> str | None:
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
