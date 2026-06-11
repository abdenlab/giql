"""Pytest fixtures for CLUSTER/MERGE predicate integration tests.

These tests assert the functional behavior of the optional ``predicate :=``
argument by executing the generated SQL against DuckDB. They do not invoke
``bedtools`` or ``pybedtools`` -- the comparison reference is a partition
computed directly in the test, matching the DISJOIN coordinate-space approach.
"""

import pytest

from giql import transpile

duckdb = pytest.importorskip("duckdb")

pytestmark = pytest.mark.integration

from tests.integration.bedtools.utils.duckdb_loader import load_intervals  # noqa: E402


@pytest.fixture(scope="function")
def duckdb_connection():
    """Provide a clean DuckDB connection for each test."""
    conn = duckdb.connect(":memory:")
    yield conn
    conn.close()


@pytest.fixture(scope="function")
def giql_query(duckdb_connection):
    """Provide a helper that loads data, transpiles GIQL, and executes.

    Usage::

        rows = giql_query(
            "SELECT MERGE(interval, predicate := score = prev.score) FROM t",
            tables=["t"],
            t=[("chr1", 0, 10, "a", 5, "+"), ...],
        )
    """

    def _run(query: str, *, tables: list[str], **table_data):
        for name, intervals in table_data.items():
            load_intervals(duckdb_connection, name, intervals)
        sql = transpile(query, tables=tables)
        return duckdb_connection.execute(sql).fetchall()

    return _run
