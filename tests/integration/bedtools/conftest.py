"""Pytest fixtures for bedtools integration tests.

This module provides shared fixtures for:
- DuckDB connections
- Interval generators
"""

import pytest

from .utils.data_models import IntervalGeneratorConfig
from .utils.interval_generator import IntervalGenerator


@pytest.fixture(scope="function")
def duckdb_connection():
    """Provide clean DuckDB connection for each test.

    Yields:
        DuckDB connection to in-memory database

    Note:
        Each test gets a fresh database with no shared state.
        Connection is automatically closed after test.
    """
    try:
        import duckdb
    except ImportError:
        pytest.skip("DuckDB not installed. Install with: pip install duckdb")

    conn = duckdb.connect(":memory:")
    yield conn
    conn.close()


@pytest.fixture(scope="function")
def interval_generator():
    """Provide configured interval generator.

    Returns:
        IntervalGenerator with deterministic seed

    Note:
        Uses seed=42 for reproducible test data.
    """
    config = IntervalGeneratorConfig(seed=42)
    return IntervalGenerator(config)
