"""Unit tests for DuckDB interval loading utility."""

import duckdb
import pytest

from .duckdb_loader import load_intervals

pytestmark = pytest.mark.integration


@pytest.fixture()
def conn():
    c = duckdb.connect(":memory:")
    yield c
    c.close()


def test_creates_table_with_correct_schema(conn):
    """
    GIVEN a DuckDB connection and interval tuples
    WHEN load_intervals() is called
    THEN table is created with columns: chrom, start, end, name, score, strand
    """
    load_intervals(conn, "test_table", [("chr1", 100, 200, "a1", 50, "+")])
    cols = conn.execute(
        "SELECT column_name FROM information_schema.columns "
        "WHERE table_name = 'test_table' ORDER BY ordinal_position"
    ).fetchall()
    col_names = [c[0] for c in cols]
    assert col_names == ["chrom", "start", "end", "name", "score", "strand"]


def test_inserts_all_rows(conn):
    """
    GIVEN multiple interval tuples
    WHEN load_intervals() is called and table is queried
    THEN all rows are present with correct values
    """
    intervals = [
        ("chr1", 100, 200, "a1", 50, "+"),
        ("chr2", 300, 400, "a2", 75, "-"),
    ]
    load_intervals(conn, "t", intervals)
    rows = conn.execute("SELECT * FROM t ORDER BY chrom").fetchall()
    assert len(rows) == 2
    assert rows[0] == ("chr1", 100, 200, "a1", 50, "+")
    assert rows[1] == ("chr2", 300, 400, "a2", 75, "-")


def test_null_handling(conn):
    """
    GIVEN tuples with None values for optional fields
    WHEN load_intervals() is called
    THEN NULL values stored correctly in DuckDB
    """
    load_intervals(conn, "t", [("chr1", 100, 200, None, None, None)])
    row = conn.execute("SELECT * FROM t").fetchone()
    assert row == ("chr1", 100, 200, None, None, None)


def test_multi_chromosome(conn):
    """
    GIVEN intervals across multiple chromosomes
    WHEN load_intervals() is called
    THEN all intervals inserted regardless of chromosome
    """
    intervals = [
        ("chr1", 100, 200, "a", 0, "+"),
        ("chr2", 100, 200, "b", 0, "+"),
        ("chrX", 100, 200, "c", 0, "+"),
    ]
    load_intervals(conn, "t", intervals)
    count = conn.execute("SELECT COUNT(*) FROM t").fetchone()[0]
    assert count == 3


def test_empty_dataset(conn):
    """
    GIVEN an empty list of intervals
    WHEN load_intervals() is called
    THEN DuckDB raises an error (executemany requires non-empty list)
    """
    with pytest.raises(duckdb.InvalidInputException):
        load_intervals(conn, "t", [])
