"""Pytest fixtures for GIQL tests."""

import pandas as pd
import pytest


@pytest.fixture(scope="session")
def to_df():
    """Fixture providing a helper to convert DuckDB results to DataFrames.

    Returns a function that materializes query results for testing.
    Session-scoped since it's a pure function with no state.

    Usage:
        result = to_df(conn.execute("SELECT ..."))
    """

    def _to_df(cursor):
        if cursor.description:
            columns = [desc[0] for desc in cursor.description]
            return pd.DataFrame(cursor.fetchall(), columns=columns)
        return pd.DataFrame()

    return _to_df
