"""Schema information for transpilation.

This module manages schema metadata for tables, including how genomic
ranges are physically stored in the database.
"""

from dataclasses import dataclass
from typing import Dict
from typing import Optional

from giql.range_parser import CoordinateSystem
from giql.range_parser import IntervalType


@dataclass
class ColumnInfo:
    """Information about a column."""

    name: str
    type: str
    is_genomic: bool = False
    # For genomic columns stored as separate fields
    chrom_col: Optional[str] = None
    start_col: Optional[str] = None
    end_col: Optional[str] = None
    strand_col: Optional[str] = None
    # Coordinate system configuration for genomic columns
    coordinate_system: CoordinateSystem = CoordinateSystem.ZERO_BASED
    interval_type: IntervalType = IntervalType.HALF_OPEN


@dataclass
class TableSchema:
    """Schema for a table."""

    name: str
    columns: Dict[str, ColumnInfo]


class SchemaInfo:
    """Manages schema information for transpilation.

    Tracks how genomic ranges are stored:
        - Separate columns (chromosome, start_pos, end_pos)
        - STRUCT types
        - Custom types
    """

    def __init__(self):
        self.tables: Dict[str, TableSchema] = {}

    def register_table(self, name: str, schema: TableSchema):
        """Register a table schema.

        :param name: Table name
        :param schema: TableSchema object
        """
        self.tables[name] = schema

    def get_table(self, name: str) -> Optional[TableSchema]:
        """Get table schema by name.

        :param name:
            Table name
        :return:
            TableSchema object or None if not found
        """
        return self.tables.get(name)

    def get_column_info(self, table: str, column: str) -> Optional[ColumnInfo]:
        """Get column information.

        :param table:
            Table name
        :param column:
            Column name
        :return:
            ColumnInfo object or None if not found
        """
        table_schema = self.get_table(table)
        if table_schema:
            return table_schema.columns.get(column)
        return None
