"""
Schema information for transpilation.
"""

from dataclasses import dataclass
from typing import Dict
from typing import Optional


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


@dataclass
class TableSchema:
    """Schema for a table."""

    name: str
    columns: Dict[str, ColumnInfo]


class SchemaInfo:
    """
    Manages schema information for transpilation.

    Tracks how genomic ranges are stored:
        - Separate columns (chromosome, start_pos, end_pos)
        - STRUCT types
        - Custom types
    """

    def __init__(self):
        self.tables: Dict[str, TableSchema] = {}

    def register_table(self, name: str, schema: TableSchema):
        """Register a table schema."""
        self.tables[name] = schema

    def get_table(self, name: str) -> Optional[TableSchema]:
        """Get table schema by name."""
        return self.tables.get(name)

    def get_column_info(self, table: str, column: str) -> Optional[ColumnInfo]:
        """Get column information."""
        table_schema = self.get_table(table)
        if table_schema:
            return table_schema.columns.get(column)
        return None
