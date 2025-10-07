"""
GIQL - Genomic Interval Query Language

A SQL dialect for genomic range queries with multi-database support.
"""

__version__ = "0.1.0"

from giql.dialect import GIQLDialect
from giql.engine import DialectType
from giql.engine import GIQLEngine
from giql.range_parser import CoordinateSystem
from giql.range_parser import IntervalType
from giql.range_parser import ParsedRange
from giql.range_parser import RangeParser
from giql.schema import ColumnInfo
from giql.schema import SchemaInfo
from giql.schema import TableSchema

__all__ = [
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
