"""GIQL - Genomic Interval Query Language.

A SQL dialect for genomic range queries.
"""

from giql.constants import DEFAULT_BIN_SIZE
from giql.table import Table
from giql.transpile import transpile

__version__ = "0.1.0"


__all__ = [
    "DEFAULT_BIN_SIZE",
    "Table",
    "transpile",
]
