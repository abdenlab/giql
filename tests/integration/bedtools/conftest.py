"""Pytest fixtures for bedtools integration tests."""

import shutil

import pytest

duckdb = pytest.importorskip("duckdb")
pytest.importorskip("pybedtools")

if not shutil.which("bedtools"):
    pytest.skip(
        "bedtools binary not found in PATH",
        allow_module_level=True,
    )


@pytest.fixture(scope="function")
def duckdb_connection():
    """Provide clean DuckDB connection for each test.

    Each test gets a fresh in-memory database with no shared state.
    """
    conn = duckdb.connect(":memory:")
    yield conn
    conn.close()
