"""
Base generator that outputs standard SQL.

Works with any SQL database that supports:
- Basic comparison operators (<, >, =, AND, OR)
- String literals
- Numeric comparisons
"""

from typing import Optional

from sqlglot import exp
from sqlglot.generator import Generator

from ..expressions import Contains
from ..expressions import Intersects
from ..expressions import SpatialSetPredicate
from ..expressions import Within
from ..range_parser import ParsedRange
from ..range_parser import RangeParser
from ..schema import SchemaInfo


class BaseGIQLGenerator(Generator):
    """
    Base generator for standard SQL output.

    This generator uses only SQL-92 compatible constructs,
    ensuring compatibility with virtually all SQL databases.
    """

    def __init__(self, schema_info: Optional[SchemaInfo] = None, **kwargs):
        super().__init__(**kwargs)
        self.schema_info = schema_info or SchemaInfo()


    def intersects_sql(self, expression: Intersects) -> str:
        """Generate standard SQL for INTERSECTS."""
        return self._generate_spatial_op(expression, "intersects")

    def contains_sql(self, expression: Contains) -> str:
        """Generate standard SQL for CONTAINS."""
        return self._generate_spatial_op(expression, "contains")

    def within_sql(self, expression: Within) -> str:
        """Generate standard SQL for WITHIN."""
        return self._generate_spatial_op(expression, "within")

    def spatialsetpredicate_sql(self, expression: SpatialSetPredicate) -> str:
        """Generate SQL for spatial set predicates (ANY/ALL)."""
        return self._generate_spatial_set(expression)

    def _generate_spatial_op(self, expression: exp.Binary, op_type: str) -> str:
        """
        Generate SQL for a spatial operation.

        Args:
            expression: AST node (Intersects, Contains, or Within)
            op_type: 'intersects', 'contains', or 'within'
        """
        left = self.sql(expression, "this")
        right_raw = self.sql(expression, "expression")

        # Parse the genomic range
        try:
            range_str = right_raw.strip("'\"")
            parsed_range = RangeParser.parse(range_str).to_zero_based_half_open()
            return self._generate_range_predicate(left, parsed_range, op_type)
        except Exception as e:
            raise ValueError(f"Could not parse genomic range: {right_raw}. Error: {e}")

    def _generate_range_predicate(
        self,
        column_ref: str,
        parsed_range: ParsedRange,
        op_type: str,
    ) -> str:
        """
        Generate SQL predicate for a range operation.

        Args:
            column_ref: Column reference (e.g., 'v.position' or 'position')
            parsed_range: Parsed genomic range
            op_type: 'intersects', 'contains', or 'within'
        """
        # Get column references
        chrom_col, start_col, end_col = self._get_column_refs(column_ref)

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

    def _generate_spatial_set(self, expression: SpatialSetPredicate) -> str:
        """
        Generate SQL for spatial set predicates (ANY/ALL).

        Examples:
            column INTERSECTS ANY(...) -> (condition1 OR condition2 OR ...)
            column INTERSECTS ALL(...) -> (condition1 AND condition2 AND ...)
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

    def _get_column_refs(self, column_ref: str) -> tuple[str, str, str]:
        """
        Get physical column names for genomic data.

        Args:
            column_ref: Logical column reference (e.g., 'v.position')

        Returns:
            Tuple of (chromosome_col, start_col, end_col)
        """
        # TODO: Use schema_info to get actual column mappings
        # For now, use default column names
        if "." in column_ref:
            table_alias, _ = column_ref.rsplit(".", 1)
            return (
                f"{table_alias}.chromosome",
                f"{table_alias}.start_pos",
                f"{table_alias}.end_pos",
            )
        else:
            return "chromosome", "start_pos", "end_pos"
