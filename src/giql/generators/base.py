"""Base generator that outputs standard SQL.

Works with any SQL database that supports:
- Basic comparison operators (<, >, =, AND, OR)
- String literals
- Numeric comparisons

This generator uses only SQL-92 compatible constructs, ensuring compatibility
with virtually all SQL databases.
"""

from typing import Optional

from sqlglot import exp
from sqlglot.generator import Generator

from giql.constants import DEFAULT_CHROM_COL
from giql.constants import DEFAULT_END_COL
from giql.constants import DEFAULT_START_COL
from giql.constants import DEFAULT_STRAND_COL
from giql.expressions import Contains
from giql.expressions import GIQLDistance
from giql.expressions import Intersects
from giql.expressions import SpatialSetPredicate
from giql.expressions import Within
from giql.range_parser import ParsedRange
from giql.range_parser import RangeParser
from giql.schema import SchemaInfo


class BaseGIQLGenerator(Generator):
    """Base generator for standard SQL output.

    This generator uses only SQL-92 compatible constructs,
    ensuring compatibility with virtually all SQL databases.
    """

    def __init__(self, schema_info: Optional[SchemaInfo] = None, **kwargs):
        super().__init__(**kwargs)
        self.schema_info = schema_info or SchemaInfo()
        self._current_table = None  # Track current table for column resolution
        self._alias_to_table = {}  # Map aliases to table names

    def select_sql(self, expression: exp.Select) -> str:
        """Override SELECT generation to track table context and aliases."""
        # Build alias-to-table mapping
        self._alias_to_table = {}

        # Extract from FROM clause
        if expression.args.get("from"):
            from_clause = expression.args["from"]
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

        # Extract stranded parameter
        stranded_expr = expression.args.get("stranded")
        stranded = False
        if stranded_expr:
            if isinstance(stranded_expr, exp.Boolean):
                stranded = stranded_expr.this
            elif isinstance(stranded_expr, exp.Literal):
                stranded = str(stranded_expr.this).upper() == "TRUE"
            else:
                stranded = str(stranded_expr).upper() in ("TRUE", "1", "YES")

        # Get SQL representations
        interval_a_sql = self.sql(interval_a)
        interval_b_sql = self.sql(interval_b)

        # Check if we're dealing with column-to-column or column-to-literal
        if "." in interval_a_sql and not interval_a_sql.startswith("'"):
            # Column reference for interval_a
            if stranded:
                chrom_a, start_a, end_a, strand_a = self._get_column_refs(interval_a_sql, None, include_strand=True)
            else:
                chrom_a, start_a, end_a = self._get_column_refs(interval_a_sql, None)
                strand_a = None
        else:
            # Literal range - not implemented yet for interval_a
            raise ValueError("Literal range as first argument not yet supported")

        if "." in interval_b_sql and not interval_b_sql.startswith("'"):
            # Column reference for interval_b
            if stranded:
                chrom_b, start_b, end_b, strand_b = self._get_column_refs(interval_b_sql, None, include_strand=True)
            else:
                chrom_b, start_b, end_b = self._get_column_refs(interval_b_sql, None)
                strand_b = None
        else:
            # Literal range - not implemented yet
            raise ValueError("Literal range as second argument not yet supported")

        # Generate CASE expression
        return self._generate_distance_case(
            chrom_a, start_a, end_a, strand_a,
            chrom_b, start_b, end_b, strand_b,
            stranded=stranded
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
    ) -> str:
        """Generate SQL CASE expression for distance calculation.

        :param chrom_a: Chromosome column for interval A
        :param start_a: Start column for interval A
        :param end_a: End column for interval A
        :param strand_a: Strand column for interval A (None if not stranded)
        :param chrom_b: Chromosome column for interval B
        :param start_b: Start column for interval B
        :param end_b: End column for interval B
        :param strand_b: Strand column for interval B (None if not stranded)
        :param stranded: Whether to use strand-aware distance calculation
        :return: SQL CASE expression
        """
        if not stranded or strand_a is None or strand_b is None:
            # Basic distance calculation without strand awareness
            return f"""CASE WHEN {chrom_a} != {chrom_b} THEN NULL WHEN {start_a} < {end_b} AND {end_a} > {start_b} THEN 0 WHEN {end_a} <= {start_b} THEN ({start_b} - {end_a}) ELSE ({start_a} - {end_b}) END"""

        # Stranded distance calculation
        # Return NULL if either strand is '.', '?', or NULL
        # Calculate distance and multiply by -1 if first interval is on '-' strand
        return f"""CASE WHEN {chrom_a} != {chrom_b} THEN NULL WHEN {strand_a} IS NULL OR {strand_b} IS NULL THEN NULL WHEN {strand_a} = '.' OR {strand_a} = '?' THEN NULL WHEN {strand_b} = '.' OR {strand_b} = '?' THEN NULL WHEN {start_a} < {end_b} AND {end_a} > {start_b} THEN 0 WHEN {end_a} <= {start_b} THEN CASE WHEN {strand_a} = '-' THEN -({start_b} - {end_a}) ELSE ({start_b} - {end_a}) END ELSE CASE WHEN {strand_a} = '-' THEN -({start_a} - {end_b}) ELSE ({start_a} - {end_b}) END END"""

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
            # Column-to-column join (e.g., a.position INTERSECTS b.position)
            return self._generate_column_join(left, right_raw, op_type)
        else:
            # Literal range string (e.g., position INTERSECTS 'chr1:1000-2000')
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
            Column reference (e.g., 'v.position' or 'position')
        :param parsed_range:
            Parsed genomic range
        :param op_type:
            'intersects', 'contains', or 'within'
        :return:
            SQL predicate string
        """
        # Get column references
        chrom_col, start_col, end_col = self._get_column_refs(
            column_ref, self._current_table
        )

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
            Left column reference (e.g., 'a.position')
        :param right_col:
            Right column reference (e.g., 'b.position')
        :param op_type:
            'intersects', 'contains', or 'within'
        :return:
            SQL predicate string
        """
        # Get column references for both sides
        # Pass None to let _get_column_refs extract and resolve table from column ref
        l_chrom, l_start, l_end = self._get_column_refs(left_col, None)
        r_chrom, r_start, r_end = self._get_column_refs(right_col, None)

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

    def _get_column_refs(
        self, column_ref: str, table_name: str | None = None, include_strand: bool = False
    ) -> tuple[str, str, str] | tuple[str, str, str, str]:
        """Get physical column names for genomic data.

        :param column_ref:
            Logical column reference (e.g., 'v.position' or 'position')
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

        # Extract table alias/name from column reference if present
        table_alias = None
        if "." in column_ref:
            table_alias, _ = column_ref.rsplit(".", 1)
            # If no explicit table_name provided, resolve alias to table name
            if not table_name:
                # Look up actual table name from alias
                table_name = self._alias_to_table.get(table_alias, self._current_table)

        # Try to get custom column names from schema
        if table_name and self.schema_info:
            table_schema = self.schema_info.get_table(table_name)
            if table_schema:
                # Find the genomic column
                for col_info in table_schema.columns.values():
                    if col_info.is_genomic:
                        if col_info.chrom_col:
                            chrom_col = col_info.chrom_col
                        if col_info.start_col:
                            start_col = col_info.start_col
                        if col_info.end_col:
                            end_col = col_info.end_col
                        if col_info.strand_col:
                            strand_col = col_info.strand_col
                        break

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
