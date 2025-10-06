"""
DuckDB-specific generator with optimizations.
"""

from sqlglot.dialects.duckdb import DuckDB

from .base import BaseGIQLGenerator


class GIQLDuckDBGenerator(BaseGIQLGenerator, DuckDB.Generator):
    """
    DuckDB-specific optimizations.

    Can leverage:
    - Efficient list operations
    - STRUCT types
    - Columnar optimizations
    """

    def __init__(self, schema_info=None, **kwargs):
        BaseGIQLGenerator.__init__(self, schema_info=schema_info, **kwargs)
        DuckDB.Generator.__init__(self, **kwargs)
