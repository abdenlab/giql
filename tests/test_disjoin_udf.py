"""End-to-end execution tests for the DISJOIN operator on DuckDB.

Tests transpile DISJOIN queries and execute them against in-memory DuckDB
to verify split correctness, the coverage filter, parent passthrough, and
degenerate-input handling.
"""

import duckdb

from giql import transpile


def _run(query: str, tables: list[str], **table_data):
    """Transpile a DISJOIN query and execute it against in-memory DuckDB."""
    conn = duckdb.connect(":memory:")
    for name, rows in table_data.items():
        conn.execute(
            f'CREATE TABLE {name}(chrom VARCHAR, "start" INTEGER, '
            f'"end" INTEGER, name VARCHAR)'
        )
        conn.executemany(f"INSERT INTO {name} VALUES (?, ?, ?, ?)", rows)
    result = conn.execute(transpile(query, tables=tables)).fetchall()
    conn.close()
    return result


class TestDisjoinExecution:
    """End-to-end DISJOIN correctness on DuckDB."""

    def test_disjoin_should_split_intervals_at_breakpoints(self):
        """
        GIVEN two overlapping target intervals A=[0,20) and B=[10,30)
        WHEN DISJOIN runs in self-mode
        THEN each interval should be split at the other's interior breakpoint
        """
        result = _run(
            "SELECT name, disjoin_start, disjoin_end FROM DISJOIN(features) "
            "ORDER BY name, disjoin_start",
            tables=["features"],
            features=[("chr1", 0, 20, "A"), ("chr1", 10, 30, "B")],
        )

        assert result == [
            ("A", 0, 10),
            ("A", 10, 20),
            ("B", 10, 20),
            ("B", 20, 30),
        ]

    def test_disjoin_should_yield_bioconductor_partition_when_distinct(self):
        """
        GIVEN overlapping target intervals in self-mode
        WHEN selecting DISTINCT sub-intervals from DISJOIN
        THEN the result should be the globally non-overlapping partition
        """
        result = _run(
            "SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(features) ORDER BY disjoin_start",
            tables=["features"],
            features=[("chr1", 0, 20, "A"), ("chr1", 10, 30, "B")],
        )

        assert result == [("chr1", 0, 10), ("chr1", 10, 20), ("chr1", 20, 30)]

    def test_disjoin_should_split_against_explicit_reference(self):
        """
        GIVEN a target interval and a separate reference set
        WHEN DISJOIN runs with an explicit reference
        THEN the target should be cut at the reference breakpoints
        """
        result = _run(
            "SELECT disjoin_start, disjoin_end "
            "FROM DISJOIN(features, reference := refs) ORDER BY disjoin_start",
            tables=["features", "refs"],
            features=[("chr1", 0, 30, "T")],
            refs=[("chr1", 0, 10, "a"), ("chr1", 10, 30, "b")],
        )

        assert result == [(0, 10), (10, 30)]

    def test_disjoin_should_drop_pieces_overlapping_no_reference(self):
        """
        GIVEN a target spanning a gap in the reference coverage
        WHEN DISJOIN runs with an explicit reference
        THEN sub-intervals overlapping no reference interval should be dropped
        """
        result = _run(
            "SELECT disjoin_start, disjoin_end "
            "FROM DISJOIN(features, reference := refs) ORDER BY disjoin_start",
            tables=["features", "refs"],
            features=[("chr1", 0, 30, "T")],
            refs=[("chr1", 0, 10, "a"), ("chr1", 20, 30, "b")],
        )

        assert result == [(0, 10), (20, 30)]

    def test_disjoin_should_yield_no_rows_when_target_is_a_point(self):
        """
        GIVEN a zero-length target interval
        WHEN DISJOIN runs
        THEN it should produce no output rows
        """
        result = _run(
            "SELECT * FROM DISJOIN(features)",
            tables=["features"],
            features=[("chr1", 5, 5, "P")],
        )

        assert result == []

    def test_disjoin_should_split_duplicate_targets_independently(self):
        """
        GIVEN two target rows with identical geometry
        WHEN DISJOIN runs with a reference that cuts them
        THEN each duplicate row should be split independently
        """
        result = _run(
            "SELECT name, disjoin_start, disjoin_end "
            "FROM DISJOIN(features, reference := refs) ORDER BY name, disjoin_start",
            tables=["features", "refs"],
            features=[("chr1", 0, 10, "X"), ("chr1", 0, 10, "Y")],
            refs=[("chr1", 0, 5, "a"), ("chr1", 5, 10, "b")],
        )

        assert result == [("X", 0, 5), ("X", 5, 10), ("Y", 0, 5), ("Y", 5, 10)]

    def test_disjoin_should_not_cross_chromosome_boundaries(self):
        """
        GIVEN target intervals on different chromosomes
        WHEN DISJOIN runs in self-mode
        THEN sub-intervals should never span a chromosome boundary
        """
        result = _run(
            "SELECT disjoin_chrom, disjoin_start, disjoin_end FROM DISJOIN(features) "
            "ORDER BY disjoin_chrom, disjoin_start",
            tables=["features"],
            features=[("chr1", 0, 20, "A"), ("chr2", 5, 25, "B")],
        )

        assert result == [("chr1", 0, 20), ("chr2", 5, 25)]

    def test_disjoin_should_pass_through_non_interval_columns(self):
        """
        GIVEN a target table carrying a non-interval column
        WHEN DISJOIN runs
        THEN every output row should carry the intact parent row alongside
            its sub-interval
        """
        result = _run(
            'SELECT name, chrom, "start", "end", disjoin_start, disjoin_end '
            "FROM DISJOIN(features) ORDER BY disjoin_start, name",
            tables=["features"],
            features=[("chr1", 0, 20, "A"), ("chr1", 10, 30, "B")],
        )

        assert result == [
            ("A", "chr1", 0, 20, 0, 10),
            ("A", "chr1", 0, 20, 10, 20),
            ("B", "chr1", 10, 30, 10, 20),
            ("B", "chr1", 10, 30, 20, 30),
        ]

    def test_disjoin_should_accept_a_cte_as_reference(self):
        """
        GIVEN a reference supplied as a CTE defined in the outer query
        WHEN DISJOIN runs against that CTE
        THEN the target should be split at the CTE's breakpoints
        """
        result = _run(
            "WITH bins AS ("
            "  SELECT 'chr1' AS chrom, 0 AS \"start\", 10 AS \"end\" "
            "  UNION ALL SELECT 'chr1', 10, 20"
            ") "
            "SELECT disjoin_start, disjoin_end "
            "FROM DISJOIN(features, reference := bins) ORDER BY disjoin_start",
            tables=["features"],
            features=[("chr1", 0, 20, "T")],
        )

        assert result == [(0, 10), (10, 20)]

    def test_disjoin_should_not_cut_at_breakpoint_on_target_boundary(self):
        """
        GIVEN a reference breakpoint coinciding with a target boundary
        WHEN DISJOIN runs
        THEN the target should not be split at that boundary
        """
        result = _run(
            "SELECT disjoin_start, disjoin_end "
            "FROM DISJOIN(features, reference := refs)",
            tables=["features", "refs"],
            features=[("chr1", 10, 20, "T")],
            refs=[("chr1", 10, 20, "r")],
        )

        assert result == [(10, 20)]
