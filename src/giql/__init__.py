"""GIQL - Genomic Interval Query Language.

A SQL dialect for genomic range queries with multi-database support.

This package provides:
    - GIQL dialect extending SQL with spatial operators
    - Query engine supporting multiple backends (DuckDB, SQLite)
    - Range parser for genomic coordinate strings
    - Schema management for genomic data
"""

__version__ = "0.1.0"

from giql.dialect import GIQLDialect as GIQLDialect
from giql.engine import DialectType as DialectType
from giql.engine import GIQLEngine as GIQLEngine
from giql.generators import BaseGIQLGenerator as BaseGIQLGenerator
from giql.generators import GIQLDuckDBGenerator as GIQLDuckDBGenerator
from giql.range_parser import CoordinateSystem as CoordinateSystem
from giql.range_parser import IntervalType as IntervalType
from giql.range_parser import ParsedRange as ParsedRange
from giql.range_parser import RangeParser as RangeParser
from giql.schema import ColumnInfo as ColumnInfo
from giql.schema import SchemaInfo as SchemaInfo
from giql.schema import TableSchema as TableSchema

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
]
