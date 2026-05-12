"""Pytest fixtures for bedtools integration tests."""

import shutil

import pytest

from giql import transpile

duckdb = pytest.importorskip("duckdb")
pytest.importorskip("pybedtools")

if not shutil.which("bedtools"):
    pytest.skip(
        "bedtools binary not found in PATH",
        allow_module_level=True,
    )

from .utils.duckdb_loader import load_intervals  # noqa: E402


@pytest.fixture(scope="function")
def duckdb_connection():
    """Provide clean DuckDB connection for each test.

    Each test gets a fresh in-memory database with no shared state.
    """
    conn = duckdb.connect(":memory:")
    yield conn
    conn.close()


@pytest.fixture(scope="function")
def giql_query(duckdb_connection):
    """Provide a helper that loads data, transpiles GIQL, and executes.

    Usage::

        result = giql_query(
            "SELECT * FROM t WHERE interval INTERSECTS 'chr1:1-100'",
            tables=["t"],
            t=[GenomicInterval("chr1", 50, 150, "x", 0, "+")],
        )
    """

    def _run(query: str, *, tables: list[str], **table_data):
        for name, intervals in table_data.items():
            load_intervals(
                duckdb_connection,
                name,
                [i.to_tuple() for i in intervals],
            )
        sql = transpile(query, tables=tables)
        return duckdb_connection.execute(sql).fetchall()

    return _run
