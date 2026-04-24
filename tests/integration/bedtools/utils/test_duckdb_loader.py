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


def test_load_intervals_should_create_table_with_default_schema(conn):
    """Test that load_intervals creates a table with the default GIQL schema.

    Given:
        A DuckDB connection and a single interval tuple.
    When:
        load_intervals is called with a target table name.
    Then:
        It should create a table with columns chrom, start, end, name, score, strand.
    """
    # Arrange, act, & assert
    load_intervals(conn, "test_table", [("chr1", 100, 200, "a1", 50, "+")])
    cols = conn.execute(
        "SELECT column_name FROM information_schema.columns "
        "WHERE table_name = 'test_table' ORDER BY ordinal_position"
    ).fetchall()
    col_names = [c[0] for c in cols]
    assert col_names == ["chrom", "start", "end", "name", "score", "strand"]


def test_load_intervals_should_insert_all_tuples(conn):
    """Test that load_intervals inserts every provided tuple.

    Given:
        A DuckDB connection and multiple interval tuples.
    When:
        load_intervals is called and the resulting table is queried.
    Then:
        It should persist each tuple with its exact field values.
    """
    # Arrange
    intervals = [
        ("chr1", 100, 200, "a1", 50, "+"),
        ("chr2", 300, 400, "a2", 75, "-"),
    ]

    # Act
    load_intervals(conn, "t", intervals)

    # Assert
    rows = conn.execute("SELECT * FROM t ORDER BY chrom").fetchall()
    assert len(rows) == 2
    assert rows[0] == ("chr1", 100, 200, "a1", 50, "+")
    assert rows[1] == ("chr2", 300, 400, "a2", 75, "-")


def test_load_intervals_should_store_nulls_when_optional_fields_are_none(conn):
    """Test that load_intervals preserves None values for optional fields.

    Given:
        A DuckDB connection and an interval tuple with None for name, score, and strand.
    When:
        load_intervals is called and the row is read back.
    Then:
        It should store the optional fields as SQL NULL values.
    """
    # Arrange, act, & assert
    load_intervals(conn, "t", [("chr1", 100, 200, None, None, None)])
    row = conn.execute("SELECT * FROM t").fetchone()
    assert row == ("chr1", 100, 200, None, None, None)


def test_load_intervals_should_insert_all_rows_when_intervals_span_multiple_chromosomes(conn):
    """Test that load_intervals loads intervals across different chromosomes.

    Given:
        A DuckDB connection and interval tuples referencing chr1, chr2, and chrX.
    When:
        load_intervals is called with the cross-chromosome dataset.
    Then:
        It should insert every row regardless of its chromosome label.
    """
    # Arrange
    intervals = [
        ("chr1", 100, 200, "a", 0, "+"),
        ("chr2", 100, 200, "b", 0, "+"),
        ("chrX", 100, 200, "c", 0, "+"),
    ]

    # Act
    load_intervals(conn, "t", intervals)

    # Assert
    count = conn.execute("SELECT COUNT(*) FROM t").fetchone()[0]
    assert count == 3


def test_load_intervals_should_create_empty_table_when_intervals_empty(conn):
    """Test that load_intervals accepts an empty interval list.

    Given:
        A DuckDB connection and an empty list of intervals.
    When:
        load_intervals is called with the empty list.
    Then:
        It should create the table with the default schema and zero rows.
    """
    # Arrange, act
    load_intervals(conn, "t", [])

    # Assert
    count = conn.execute("SELECT COUNT(*) FROM t").fetchone()[0]
    assert count == 0
