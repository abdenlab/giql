"""End-to-end execution tests for the DISJOIN operator.

Tests transpile DISJOIN queries and execute them against in-memory DuckDB and
SQLite to verify split correctness, the coverage filter, parent passthrough,
and degenerate-input handling on every supported in-process backend.
"""

import sqlite3

import duckdb
import pytest

from giql import transpile

_TABLE_COLUMNS = '(chrom VARCHAR, "start" INTEGER, "end" INTEGER, name VARCHAR)'


def _run_duckdb(query: str, tables: list[str], **table_data):
    """Transpile a DISJOIN query and execute it against in-memory DuckDB."""
    conn = duckdb.connect(":memory:")
    for name, rows in table_data.items():
        conn.execute(f"CREATE TABLE {name}{_TABLE_COLUMNS}")
        conn.executemany(f"INSERT INTO {name} VALUES (?, ?, ?, ?)", rows)
    result = conn.execute(transpile(query, tables=tables)).fetchall()
    conn.close()
    return result


def _run_sqlite(query: str, tables: list[str], **table_data):
    """Transpile a DISJOIN query and execute it against in-memory SQLite."""
    conn = sqlite3.connect(":memory:")
    for name, rows in table_data.items():
        conn.execute(f"CREATE TABLE {name}{_TABLE_COLUMNS}")
        conn.executemany(f"INSERT INTO {name} VALUES (?, ?, ?, ?)", rows)
    result = conn.execute(transpile(query, tables=tables)).fetchall()
    conn.close()
    return result


_ENGINES = {"duckdb": _run_duckdb, "sqlite": _run_sqlite}


@pytest.fixture(
    params=[
        "duckdb",
        pytest.param(
            "sqlite",
            marks=pytest.mark.skipif(
                sqlite3.sqlite_version_info < (3, 25),
                reason="DISJOIN emits LEAD(), which requires SQLite 3.25+",
            ),
        ),
    ]
)
def run(request):
    """Yield an engine-specific transpile-and-execute callable."""
    return _ENGINES[request.param]


class TestDisjoinExecution:
    """End-to-end DISJOIN correctness on every supported backend."""

    def test_disjoin_should_split_intervals_at_breakpoints(self, run):
        """Test that DISJOIN splits intervals at interior breakpoints.

        Given:
            Two overlapping target intervals A=[0,20) and B=[10,30).
        When:
            DISJOIN runs in self-mode.
        Then:
            It should split each interval at the other's interior breakpoint.
        """
        # Arrange & act
        result = run(
            "SELECT name, disjoin_start, disjoin_end FROM DISJOIN(features) "
            "ORDER BY name, disjoin_start",
            tables=["features"],
            features=[("chr1", 0, 20, "A"), ("chr1", 10, 30, "B")],
        )

        # Assert
        assert result == [
            ("A", 0, 10),
            ("A", 10, 20),
            ("B", 10, 20),
            ("B", 20, 30),
        ]

    def test_disjoin_should_yield_bioconductor_partition_when_distinct(self, run):
        """Test that DISTINCT sub-intervals yield the Bioconductor partition.

        Given:
            Overlapping target intervals in self-mode.
        When:
            Selecting DISTINCT sub-intervals from DISJOIN.
        Then:
            It should yield the globally non-overlapping partition.
        """
        # Arrange & act
        result = run(
            "SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(features) ORDER BY disjoin_start",
            tables=["features"],
            features=[("chr1", 0, 20, "A"), ("chr1", 10, 30, "B")],
        )

        # Assert
        assert result == [("chr1", 0, 10), ("chr1", 10, 20), ("chr1", 20, 30)]

    def test_disjoin_should_split_against_explicit_reference(self, run):
        """Test that DISJOIN cuts the target at an explicit reference.

        Given:
            A target interval and a separate reference set.
        When:
            DISJOIN runs with an explicit reference.
        Then:
            It should cut the target at the reference breakpoints.
        """
        # Arrange & act
        result = run(
            "SELECT disjoin_start, disjoin_end "
            "FROM DISJOIN(features, reference := refs) ORDER BY disjoin_start",
            tables=["features", "refs"],
            features=[("chr1", 0, 30, "T")],
            refs=[("chr1", 0, 10, "a"), ("chr1", 10, 30, "b")],
        )

        # Assert
        assert result == [(0, 10), (10, 30)]

    def test_disjoin_should_drop_pieces_overlapping_no_reference(self, run):
        """Test that sub-intervals overlapping no reference are dropped.

        Given:
            A target spanning a gap in the reference coverage.
        When:
            DISJOIN runs with an explicit reference.
        Then:
            It should drop sub-intervals overlapping no reference interval.
        """
        # Arrange & act
        result = run(
            "SELECT disjoin_start, disjoin_end "
            "FROM DISJOIN(features, reference := refs) ORDER BY disjoin_start",
            tables=["features", "refs"],
            features=[("chr1", 0, 30, "T")],
            refs=[("chr1", 0, 10, "a"), ("chr1", 20, 30, "b")],
        )

        # Assert
        assert result == [(0, 10), (20, 30)]

    def test_disjoin_should_yield_no_rows_when_target_is_a_point(self, run):
        """Test that a zero-length target produces no rows.

        Given:
            A zero-length target interval.
        When:
            DISJOIN runs.
        Then:
            It should produce no output rows.
        """
        # Arrange & act
        result = run(
            "SELECT * FROM DISJOIN(features)",
            tables=["features"],
            features=[("chr1", 5, 5, "P")],
        )

        # Assert
        assert result == []

    def test_disjoin_should_split_duplicate_targets_independently(self, run):
        """Test that duplicate target rows are split independently.

        Given:
            Two target rows with identical geometry.
        When:
            DISJOIN runs with a reference that cuts them.
        Then:
            It should split each duplicate row independently.
        """
        # Arrange & act
        result = run(
            "SELECT name, disjoin_start, disjoin_end "
            "FROM DISJOIN(features, reference := refs) ORDER BY name, disjoin_start",
            tables=["features", "refs"],
            features=[("chr1", 0, 10, "X"), ("chr1", 0, 10, "Y")],
            refs=[("chr1", 0, 5, "a"), ("chr1", 5, 10, "b")],
        )

        # Assert
        assert result == [("X", 0, 5), ("X", 5, 10), ("Y", 0, 5), ("Y", 5, 10)]

    def test_disjoin_should_not_cross_chromosome_boundaries(self, run):
        """Test that sub-intervals never cross a chromosome boundary.

        Given:
            Target intervals on different chromosomes.
        When:
            DISJOIN runs in self-mode.
        Then:
            It should never produce a sub-interval spanning a chromosome
            boundary.
        """
        # Arrange & act
        result = run(
            "SELECT disjoin_chrom, disjoin_start, disjoin_end FROM DISJOIN(features) "
            "ORDER BY disjoin_chrom, disjoin_start",
            tables=["features"],
            features=[("chr1", 0, 20, "A"), ("chr2", 5, 25, "B")],
        )

        # Assert
        assert result == [("chr1", 0, 20), ("chr2", 5, 25)]

    def test_disjoin_should_pass_through_non_interval_columns(self, run):
        """Test that non-interval columns pass through unchanged.

        Given:
            A target table carrying a non-interval column.
        When:
            DISJOIN runs.
        Then:
            It should carry the intact parent row on every output row
            alongside its sub-interval.
        """
        # Arrange & act
        result = run(
            'SELECT name, chrom, "start", "end", disjoin_start, disjoin_end '
            "FROM DISJOIN(features) ORDER BY disjoin_start, name",
            tables=["features"],
            features=[("chr1", 0, 20, "A"), ("chr1", 10, 30, "B")],
        )

        # Assert
        assert result == [
            ("A", "chr1", 0, 20, 0, 10),
            ("A", "chr1", 0, 20, 10, 20),
            ("B", "chr1", 10, 30, 10, 20),
            ("B", "chr1", 10, 30, 20, 30),
        ]

    def test_disjoin_should_accept_a_cte_as_reference(self, run):
        """Test that a CTE can be used as the reference.

        Given:
            A reference supplied as a CTE defined in the outer query.
        When:
            DISJOIN runs against that CTE.
        Then:
            It should split the target at the CTE's breakpoints.
        """
        # Arrange & act
        result = run(
            "WITH bins AS ("
            "  SELECT 'chr1' AS chrom, 0 AS \"start\", 10 AS \"end\" "
            "  UNION ALL SELECT 'chr1', 10, 20"
            ") "
            "SELECT disjoin_start, disjoin_end "
            "FROM DISJOIN(features, reference := bins) ORDER BY disjoin_start",
            tables=["features"],
            features=[("chr1", 0, 20, "T")],
        )

        # Assert
        assert result == [(0, 10), (10, 20)]

    def test_disjoin_should_not_cut_at_breakpoint_on_target_boundary(self, run):
        """Test that a breakpoint on a target boundary does not split it.

        Given:
            A reference breakpoint coinciding with a target boundary.
        When:
            DISJOIN runs.
        Then:
            It should not split the target at that boundary.
        """
        # Arrange & act
        result = run(
            "SELECT disjoin_start, disjoin_end "
            "FROM DISJOIN(features, reference := refs)",
            tables=["features", "refs"],
            features=[("chr1", 10, 20, "T")],
            refs=[("chr1", 10, 20, "r")],
        )

        # Assert
        assert result == [(10, 20)]
