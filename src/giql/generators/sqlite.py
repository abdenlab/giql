"""SQLite-specific generator.

This module provides SQLite-specific SQL generation for GIQL queries.
SQLite does not support LATERAL joins, so NEAREST uses window functions instead.
"""

from sqlglot.dialects.sqlite import SQLite

from giql.generators.base import BaseGIQLGenerator


class GIQLSQLiteGenerator(BaseGIQLGenerator, SQLite.Generator):
    """SQLite-specific SQL generator.

    Key differences from other dialects:
    - No LATERAL join support - uses window functions for NEAREST
    - Window functions available since SQLite 3.25.0 (2018-09-15)
    """

    # SQLite does not support LATERAL joins
    SUPPORTS_LATERAL = False

    def __init__(self, schema_info=None, **kwargs):
        super().__init__(schema_info=schema_info, **kwargs)
