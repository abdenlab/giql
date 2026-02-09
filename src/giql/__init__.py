"""GIQL - Genomic Interval Query Language.

A SQL dialect for genomic range queries.

This package provides:
    - GIQL dialect extending SQL with spatial operators (INTERSECTS, CONTAINS, WITHIN)
    - CLUSTER and MERGE operations for interval grouping
    - NEAREST operator for finding closest intervals
    - Range parser for genomic coordinate strings
    - Transpilation to standard SQL-92 compatible output
"""

from giql.table import Table
from giql.transpile import transpile

__version__ = "0.1.0"


__all__ = [
    "Table",
    "transpile",
]
