from typing import Final

from sqlglot.dialects.sqlite import SQLite

from giql.generators.base import BaseGIQLGenerator


class GIQLSQLiteGenerator(BaseGIQLGenerator, SQLite.Generator):
    """SQLite-specific SQL generator.

    SQLite does not support LATERAL joins, so correlated NEAREST queries
    (without explicit reference) will raise an error. Use standalone mode
    with an explicit reference parameter instead.

    Example::

        -- This works (standalone mode with explicit reference):
        SELECT * FROM NEAREST(genes, reference='chr1:1000-2000', k=3)

        -- This fails (correlated mode requires LATERAL):
        SELECT * FROM peaks CROSS JOIN LATERAL NEAREST(genes, k=3)
    """

    SUPPORTS_LATERAL: Final = False
