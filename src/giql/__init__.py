"""GIQL - Genomic Interval Query Language.

A SQL dialect for genomic range queries with multi-database support.

This package provides:
    - GIQL dialect extending SQL with spatial operators
    - Query engine supporting multiple backends (DuckDB, SQLite)
    - Range parser for genomic coordinate strings
    - Schema management for genomic data
"""

from giql.engine import GIQLEngine as GIQLEngine

__version__ = "0.1.0"


__all__ = [
    "GIQLEngine",
]
