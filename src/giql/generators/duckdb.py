from sqlglot.dialects.duckdb import DuckDB

from giql.generators.base import BaseGIQLGenerator


class GIQLDuckDBGenerator(BaseGIQLGenerator, DuckDB.Generator):
    """DuckDB-specific generator with optimizations."""
