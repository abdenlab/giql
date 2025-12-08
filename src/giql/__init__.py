"""GIQL - Genomic Interval Query Language.

A SQL dialect for genomic range queries with multi-database support.

This package provides:
    - GIQL dialect extending SQL with spatial operators
    - Query engine supporting multiple backends (DuckDB, SQLite)
    - Range parser for genomic coordinate strings
    - Schema management for genomic data
"""

from giql.constants import DEFAULT_CHROM_COL as DEFAULT_CHROM_COL
from giql.constants import DEFAULT_END_COL as DEFAULT_END_COL
from giql.constants import DEFAULT_GENOMIC_COL as DEFAULT_GENOMIC_COL
from giql.constants import DEFAULT_START_COL as DEFAULT_START_COL
from giql.constants import DEFAULT_STRAND_COL as DEFAULT_STRAND_COL
from giql.dialect import GIQLDialect as GIQLDialect
from giql.engine import DialectType as DialectType
from giql.engine import GIQLEngine as GIQLEngine
from giql.generators import BaseGIQLGenerator as BaseGIQLGenerator
from giql.generators import GIQLDuckDBGenerator as GIQLDuckDBGenerator
from giql.protocols import CursorLike as CursorLike
from giql.range_parser import CoordinateSystem as CoordinateSystem
from giql.range_parser import IntervalType as IntervalType
from giql.range_parser import ParsedRange as ParsedRange
from giql.range_parser import RangeParser as RangeParser
from giql.schema import ColumnInfo as ColumnInfo
from giql.schema import SchemaInfo as SchemaInfo
from giql.schema import TableSchema as TableSchema

__version__ = "0.1.0"


__all__ = [
    "BaseGIQLGenerator",
    "GIQLDuckDBGenerator",
    "GIQLEngine",
    "DialectType",
    "GIQLDialect",
    "RangeParser",
    "ParsedRange",
    "CoordinateSystem",
    "IntervalType",
    "SchemaInfo",
    "TableSchema",
    "ColumnInfo",
    "DEFAULT_CHROM_COL",
    "DEFAULT_START_COL",
    "DEFAULT_END_COL",
    "DEFAULT_STRAND_COL",
    "DEFAULT_GENOMIC_COL",
]
