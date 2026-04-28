"""Pytest fixtures for coordinate-space matrix integration tests.

These tests assert convention-invariance of GIQL DISTANCE / NEAREST against
itself via DuckDB. They do not invoke ``bedtools`` or ``pybedtools`` -- the
comparison reference is GIQL with the canonical 0-based half-open encoding.
"""

import pytest

from giql import transpile
from giql.table import Table

duckdb = pytest.importorskip("duckdb")

pytestmark = pytest.mark.integration

from tests.integration.bedtools.utils.duckdb_loader import (  # noqa: E402
    load_intervals,
)


@pytest.fixture(scope="function")
def duckdb_connection():
    """Provide a clean DuckDB connection for each test."""
    conn = duckdb.connect(":memory:")
    yield conn
    conn.close()


@pytest.fixture(scope="function")
def giql_query(duckdb_connection):
    """Provide a helper that loads encoded intervals, transpiles, and executes.

    The caller passes ``tables`` as a list of :class:`giql.table.Table`
    objects (with ``coordinate_system`` / ``interval_type`` already set) and
    ``**table_data`` mapping table name to a list of pre-encoded
    ``(chrom, stored_start, stored_end, name, score, strand)`` tuples.
    """

    def _run(query: str, *, tables: list[Table], **table_data):
        for name, rows in table_data.items():
            load_intervals(duckdb_connection, name, rows)
        sql = transpile(query, tables=tables)
        return duckdb_connection.execute(sql).fetchall()

    return _run
