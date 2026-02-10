"""GIQL - Genomic Interval Query Language.

A SQL dialect for genomic range queries.
"""

from giql.table import Table
from giql.transpile import transpile

__version__ = "0.1.0"


__all__ = [
    "Table",
    "transpile",
]
