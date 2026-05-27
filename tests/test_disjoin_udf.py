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


def _create_and_insert(conn, name: str, rows, schema_ddl: str) -> None:
    """Create a table with ``schema_ddl`` and bulk-insert ``rows``."""
    conn.execute(f"CREATE TABLE {name}{schema_ddl}")
    if not rows:
        return
    placeholders = ", ".join("?" for _ in rows[0])
    conn.executemany(f"INSERT INTO {name} VALUES ({placeholders})", rows)


def _run_duckdb(
    query: str,
    tables: list,
    schemas: dict[str, str] | None = None,
    **table_data,
):
    """Transpile a DISJOIN query and execute it against in-memory DuckDB.

    ``schemas`` maps a table name to a CREATE-TABLE column list (including the
    enclosing parentheses) for tables whose physical layout deviates from the
    canonical ``(chrom, start, end, name)`` schema. Tables not present in
    ``schemas`` use the canonical schema.
    """
    schemas = schemas or {}
    conn = duckdb.connect(":memory:")
    for name, rows in table_data.items():
        _create_and_insert(conn, name, rows, schemas.get(name, _TABLE_COLUMNS))
    result = conn.execute(transpile(query, tables=tables)).fetchall()
    conn.close()
    return result


def _run_sqlite(
    query: str,
    tables: list,
    schemas: dict[str, str] | None = None,
    **table_data,
):
    """Transpile a DISJOIN query and execute it against in-memory SQLite.

    ``schemas`` behaves identically to :func:`_run_duckdb`.
    """
    schemas = schemas or {}
    conn = sqlite3.connect(":memory:")
    for name, rows in table_data.items():
        _create_and_insert(conn, name, rows, schemas.get(name, _TABLE_COLUMNS))
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
            '  SELECT \'chr1\' AS chrom, 0 AS "start", 10 AS "end" '
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


class TestDisjoinCteTargetSelfMode:
    """Execution tests for a DISJOIN whose target is an in-query CTE in self-mode."""

    def test_disjoin_with_cte_target_should_match_standalone_disjoin_rows(self, run):
        """Test that a CTE-target self-mode DISJOIN matches STANDALONE row-by-row.

        Given:
            A registered ``features`` table with two overlapping intervals.
        When:
            Running both ``SELECT * FROM DISJOIN(features)`` and a
            CTE-target variant that wraps ``features`` in an identity CTE.
        Then:
            Both queries should return identical disjoin segments.
        """
        # Arrange
        rows = [("chr1", 0, 20, "A"), ("chr1", 10, 30, "B")]
        projection = (
            "SELECT disjoin_chrom, disjoin_start, disjoin_end {body} "
            "ORDER BY disjoin_chrom, disjoin_start, disjoin_end"
        )

        # Act
        standalone = run(
            projection.format(body="FROM DISJOIN(features)"),
            tables=["features"],
            features=rows,
        )
        cte = run(
            'WITH x AS (SELECT chrom, "start", "end" FROM features) '
            + projection.format(body="FROM DISJOIN(x)"),
            tables=["features"],
            features=rows,
        )

        # Assert
        assert cte == standalone

    def test_disjoin_with_cte_target_should_resolve_columns_independent_of_projection_order(
        self, run
    ):
        """Test that the CTE's column ordering does not affect resolution.

        Given:
            A ``features`` table with two overlapping intervals.
        When:
            Running a self-mode DISJOIN against a CTE that reverses the
            column projection order.
        Then:
            It should yield the canonical Bioconductor partition.
        """
        # Arrange & act
        result = run(
            'WITH x AS (SELECT name, "end", "start", chrom FROM features) '
            "SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(x) ORDER BY disjoin_start",
            tables=["features"],
            features=[("chr1", 0, 20, "A"), ("chr1", 10, 30, "B")],
        )

        # Assert
        assert result == [("chr1", 0, 10), ("chr1", 10, 20), ("chr1", 20, 30)]

    def test_disjoin_with_cte_target_should_pass_through_extra_cte_columns(self, run):
        """Test that extra non-genomic CTE columns reach the output.

        Given:
            A CTE that projects an extra ``nl`` column computed from
            ``name``.
        When:
            Running a self-mode DISJOIN over the CTE.
        Then:
            The parent passthrough should carry the extra column.
        """
        # Arrange & act
        result = run(
            'WITH x AS (SELECT chrom, "start", "end", name, '
            "length(name) AS nl FROM features) "
            "SELECT disjoin_chrom, disjoin_start, disjoin_end, name, nl "
            "FROM DISJOIN(x) ORDER BY name, disjoin_start",
            tables=["features"],
            features=[("chr1", 0, 20, "A"), ("chr1", 10, 30, "BB")],
        )

        # Assert
        assert result == [
            ("chr1", 0, 10, "A", 1),
            ("chr1", 10, 20, "A", 1),
            ("chr1", 10, 20, "BB", 2),
            ("chr1", 20, 30, "BB", 2),
        ]

    def test_disjoin_with_values_only_cte_target_should_succeed_without_registered_table(
        self,
    ):
        """Test the headline use case: a VALUES-only CTE target.

        Given:
            No registered table; a CTE built from a ``VALUES`` clause.
        When:
            Running a self-mode DISJOIN over the CTE.
        Then:
            The query should execute successfully against DuckDB and yield
            the expected disjoint segments.
        """
        # Arrange & act
        result = _run_duckdb(
            'WITH x(chrom, "start", "end") AS ('
            "VALUES ('chr1', 0, 100), ('chr1', 50, 150)) "
            "SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(x) ORDER BY disjoin_start",
            tables=[],
        )

        # Assert
        assert result == [
            ("chr1", 0, 50),
            ("chr1", 50, 100),
            ("chr1", 100, 150),
        ]

    def test_disjoin_with_empty_cte_target_should_return_no_rows(self, run):
        """Test that an empty CTE target yields zero rows without error.

        Given:
            A non-empty registered ``features`` table and a CTE filtered to
            no rows.
        When:
            Running a self-mode DISJOIN over the empty CTE.
        Then:
            The result should be an empty list.
        """
        # Arrange & act
        result = run(
            "WITH x AS (SELECT * FROM features WHERE 1 = 0) SELECT * FROM DISJOIN(x)",
            tables=["features"],
            features=[("chr1", 0, 20, "A"), ("chr1", 10, 30, "B")],
        )

        # Assert
        assert result == []

    def test_disjoin_with_null_endpoint_cte_target_should_return_no_rows(self):
        """Test that NULL endpoints in the CTE propagate as zero output rows.

        Given:
            A CTE projecting ``NULL`` for the ``start`` column.
        When:
            Running a self-mode DISJOIN against DuckDB.
        Then:
            The result should be an empty list — NULL endpoints fail the
            ``s.seg_end IS NOT NULL`` and ``s.seg_end > s.seg_start`` guards.
        """
        # Arrange & act
        result = _run_duckdb(
            'WITH x AS (SELECT chrom, CAST(NULL AS INTEGER) AS "start", '
            '"end", name FROM features) SELECT * FROM DISJOIN(x)',
            tables=["features"],
            features=[("chr1", 0, 20, "A"), ("chr1", 10, 30, "B")],
        )

        # Assert
        assert result == []


class TestDisjoinCteTargetReferenceMode:
    """Execution tests for a CTE-target DISJOIN with an explicit reference."""

    def test_disjoin_with_filtering_cte_target_and_registered_reference(self, run):
        """Test that a filtering CTE target composes with a registered reference.

        Given:
            ``features`` with rows ``A`` and ``B``, ``mask`` covering only
            parts of ``A``; a CTE projects only ``A``.
        When:
            Running a DISJOIN with the CTE target and the registered mask.
        Then:
            Only the mask-covered piece of ``A`` should survive.
        """
        # Arrange & act
        result = run(
            "WITH filtered AS (SELECT * FROM features WHERE name = 'A') "
            "SELECT name, disjoin_start, disjoin_end "
            "FROM DISJOIN(filtered, reference := mask) "
            "ORDER BY disjoin_start",
            tables=["features", "mask"],
            features=[("chr1", 0, 30, "A"), ("chr1", 60, 90, "B")],
            mask=[("chr1", 0, 20, "m1"), ("chr1", 70, 80, "m2")],
        )

        # Assert
        assert result == [("A", 0, 20)]

    def test_disjoin_with_two_distinct_cte_target_and_reference(self, run):
        """Test that a CTE target plus a distinct CTE reference work together.

        Given:
            Two registered tables wrapped in distinct CTEs as target and
            reference.
        When:
            Running a DISJOIN with both CTE inputs.
        Then:
            The output should reflect ``target ∩ reference-coverage`` cuts.
        """
        # Arrange & act
        result = run(
            "WITH tgt_cte AS (SELECT * FROM features), "
            "ref_cte AS (SELECT * FROM regions) "
            "SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(tgt_cte, reference := ref_cte) "
            "ORDER BY disjoin_chrom, disjoin_start",
            tables=["features", "regions"],
            features=[("chr1", 0, 30, "A"), ("chr2", 0, 30, "X")],
            regions=[("chr1", 10, 25, "r"), ("chr2", 0, 15, "s")],
        )

        # Assert
        assert result == [("chr1", 10, 25), ("chr2", 0, 15)]

    def test_disjoin_with_cte_target_and_subquery_reference(self, run):
        """Test that a CTE target accepts a subquery reference.

        Given:
            A CTE target and a subquery reference; mask covers two pieces of
            the target.
        When:
            Running the DISJOIN.
        Then:
            Only the covered pieces should survive.
        """
        # Arrange & act
        result = run(
            "WITH x AS (SELECT * FROM features) "
            "SELECT disjoin_start, disjoin_end "
            "FROM DISJOIN(x, reference := (SELECT * FROM mask)) "
            "ORDER BY disjoin_start",
            tables=["features", "mask"],
            features=[("chr1", 0, 30, "A")],
            mask=[("chr1", 0, 20, "m1"), ("chr1", 25, 30, "m2")],
        )

        # Assert
        assert result == [(0, 20), (25, 30)]

    def test_disjoin_with_same_cte_for_target_and_reference_keeps_exists(self, run):
        """Test that a single CTE used on both sides yields the self-mode partition.

        Given:
            A single CTE used as both target and reference; output equals the
            self-mode DISJOIN even though EXISTS is emitted (conservative
            duplicate handling).
        When:
            Running the DISJOIN.
        Then:
            The rows should match a self-mode DISJOIN over the same CTE.
        """
        # Arrange
        rows = [("chr1", 0, 30, "A"), ("chr1", 20, 50, "B")]
        projection = (
            "SELECT name, disjoin_start, disjoin_end {body} ORDER BY name, disjoin_start"
        )

        # Act
        same_cte = run(
            "WITH x AS (SELECT * FROM features) "
            + projection.format(body="FROM DISJOIN(x, reference := x)"),
            tables=["features"],
            features=rows,
        )
        self_mode = run(
            "WITH x AS (SELECT * FROM features) "
            + projection.format(body="FROM DISJOIN(x)"),
            tables=["features"],
            features=rows,
        )

        # Assert
        assert same_cte == self_mode

    def test_disjoin_with_cte_target_and_custom_column_registered_reference(self):
        """Test that a CTE target works with a custom-column registered reference.

        Given:
            A canonical ``features`` table and a custom-column ``mask`` table
            registered with ``Table`` config.
        When:
            Running a DISJOIN with the CTE target and the registered mask.
        Then:
            The mask's coverage filter should drop uncovered pieces.
        """
        # Arrange
        from giql.table import Table

        # Act
        result = _run_duckdb(
            "WITH x AS (SELECT * FROM features) "
            "SELECT disjoin_start, disjoin_end "
            "FROM DISJOIN(x, reference := mask) "
            "ORDER BY disjoin_start",
            tables=[
                "features",
                Table("mask", chrom_col="seqid", start_col="lo", end_col="hi"),
            ],
            schemas={
                "mask": ("(seqid VARCHAR, lo INTEGER, hi INTEGER, label VARCHAR)"),
            },
            features=[("chr1", 0, 30, "A")],
            mask=[("chr1", 0, 20, "m1")],
        )

        # Assert
        assert result == [(0, 20)]

    def test_disjoin_with_cte_target_and_one_based_registered_reference(self):
        """Test that a CTE target works with a 1-based-closed registered reference.

        Given:
            A canonical ``features`` table and a 1-based-closed ``mask``
            table registered with ``Table`` config.
        When:
            Running a DISJOIN with the CTE target and the registered mask.
        Then:
            The output should be in canonical 0-based half-open coordinates.
        """
        # Arrange
        from giql.table import Table

        # Act
        result = _run_duckdb(
            "WITH x AS (SELECT * FROM features) "
            "SELECT disjoin_start, disjoin_end "
            "FROM DISJOIN(x, reference := mask) "
            "ORDER BY disjoin_start",
            tables=[
                "features",
                Table("mask", coordinate_system="1based", interval_type="closed"),
            ],
            features=[("chr1", 0, 30, "A")],
            mask=[("chr1", 1, 20, "m1")],
        )

        # Assert
        assert result == [(0, 20)]

    def test_disjoin_with_cte_target_and_empty_registered_reference(self, run):
        """Test that an empty registered reference drops every CTE-target row.

        Given:
            A non-empty CTE target and an empty registered mask.
        When:
            Running a DISJOIN with the CTE target and the empty mask.
        Then:
            The result should be an empty list — every segment fails the
            coverage filter.
        """
        # Arrange & act
        result = run(
            "WITH x AS (SELECT * FROM features) "
            "SELECT * FROM DISJOIN(x, reference := empty_mask)",
            tables=["features", "empty_mask"],
            features=[("chr1", 0, 30, "A")],
            empty_mask=[],
        )

        # Assert
        assert result == []

    def test_disjoin_with_filtering_cte_target_and_registered_reference_keeps_coverage(
        self, run
    ):
        """Test that CTE row-set filtering composes with reference-mode coverage.

        Given:
            A CTE filters ``features`` to row ``B`` only; ``mask`` has two
            ranges; only one overlaps row ``B``.
        When:
            Running the DISJOIN.
        Then:
            Only the mask-covered piece of ``B`` should survive.
        """
        # Arrange & act
        result = run(
            "WITH b_only AS (SELECT * FROM features WHERE name = 'B') "
            "SELECT name, disjoin_start, disjoin_end "
            "FROM DISJOIN(b_only, reference := mask) "
            "ORDER BY disjoin_start",
            tables=["features", "mask"],
            features=[("chr1", 0, 40, "A"), ("chr1", 50, 90, "B")],
            mask=[("chr1", 0, 15, "m1"), ("chr1", 60, 80, "m2")],
        )

        # Assert
        assert result == [("B", 60, 80)]


class TestDisjoinCteTargetEncoding:
    """Execution tests for CTE targets crossing custom columns or encodings."""

    def test_disjoin_with_custom_column_alias_in_cte_matches_canonical_output(self):
        """Test that aliasing custom columns inside the CTE matches canonical DISJOIN.

        Given:
            A custom-column ``features`` table; the CTE projects canonical
            ``chrom`` / ``start`` / ``end`` names.
        When:
            Running a self-mode DISJOIN.
        Then:
            The result should match a canonical-column equivalent.
        """
        # Arrange
        from giql.table import Table

        # Act
        result = _run_duckdb(
            'WITH x AS (SELECT seqid AS chrom, lo AS "start", hi AS "end" '
            "FROM features) "
            "SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(x) ORDER BY disjoin_chrom, disjoin_start",
            tables=[
                Table("features", chrom_col="seqid", start_col="lo", end_col="hi"),
            ],
            schemas={
                "features": ("(seqid VARCHAR, lo INTEGER, hi INTEGER, label VARCHAR)"),
            },
            features=[("chr1", 0, 100, "a"), ("chr1", 60, 180, "b")],
        )

        # Assert
        assert result == [("chr1", 0, 60), ("chr1", 60, 100), ("chr1", 100, 180)]

    def test_disjoin_with_user_canonicalized_one_based_cte_matches_zero_based_run(self):
        """Test that user-supplied 1-based canonicalization matches a 0-based table.

        Given:
            A 1-based-closed ``features_1based`` registered table and a 0-based
            canonical ``features_0based`` registered table holding the same
            logical intervals.
        When:
            The CTE applies ``start - 1`` canonicalization to ``features_1based``
            and DISJOIN runs on it; separately DISJOIN runs on
            ``features_0based``.
        Then:
            Both result sets should be identical.
        """
        # Arrange
        from giql.table import Table

        projection = (
            "SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end {body} "
            "ORDER BY disjoin_chrom, disjoin_start"
        )

        # Act
        cte_run = _run_duckdb(
            'WITH x AS (SELECT chrom, "start" - 1 AS "start", "end" AS "end" '
            "FROM features_1based) " + projection.format(body="FROM DISJOIN(x)"),
            tables=[
                Table(
                    "features_1based",
                    coordinate_system="1based",
                    interval_type="closed",
                ),
                "features_0based",
            ],
            features_1based=[("chr1", 1, 100, "a"), ("chr1", 61, 180, "b")],
            features_0based=[("chr1", 0, 100, "a"), ("chr1", 60, 180, "b")],
        )
        zero_run = _run_duckdb(
            projection.format(body="FROM DISJOIN(features_0based)"),
            tables=[
                Table(
                    "features_1based",
                    coordinate_system="1based",
                    interval_type="closed",
                ),
                "features_0based",
            ],
            features_1based=[("chr1", 1, 100, "a"), ("chr1", 61, 180, "b")],
            features_0based=[("chr1", 0, 100, "a"), ("chr1", 60, 180, "b")],
        )

        # Assert
        assert cte_run == zero_run

    def test_disjoin_with_uncanonicalized_one_based_cte_is_off_by_one(self):
        """Test the documented contract: uncanonicalized 1-based CTE produces off-by-one rows.

        Given:
            A 1-based-closed ``features`` registered table; the CTE projects
            raw values without canonicalization.
        When:
            Running a self-mode DISJOIN.
        Then:
            Output is in raw 1-based values treated as 0-based half-open —
            this is the documented contract for non-canonical CTE input.
        """
        # Arrange
        from giql.table import Table

        # Act
        result = _run_duckdb(
            "WITH x AS (SELECT * FROM features_1based) "
            "SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(x) ORDER BY disjoin_start",
            tables=[
                Table(
                    "features_1based",
                    coordinate_system="1based",
                    interval_type="closed",
                ),
            ],
            features_1based=[("chr1", 1, 100, "a"), ("chr1", 61, 180, "b")],
        )

        # Assert
        assert result == [("chr1", 1, 61), ("chr1", 61, 100), ("chr1", 100, 180)]

    def test_disjoin_with_cte_shadow_ignores_registered_one_based_config(self):
        """Test that a CTE-shadow target ignores the registered table's coord system.

        Given:
            A registered 1-based ``features`` shadowed by a CTE; the CTE
            projects rows verbatim.
        When:
            Running a self-mode DISJOIN.
        Then:
            Output reflects raw values (no ``-1`` shift) — the registered
            ``Table`` config is ignored once shadowed.
        """
        # Arrange
        from giql.table import Table

        # Act
        result = _run_duckdb(
            "WITH features AS (SELECT * FROM features) "
            "SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(features) ORDER BY disjoin_start",
            tables=[
                Table(
                    "features",
                    coordinate_system="1based",
                    interval_type="closed",
                ),
            ],
            features=[("chr1", 1, 20, "A"), ("chr1", 11, 30, "B")],
        )

        # Assert — raw 1-based values, treated as 0-based, partition at
        # the interior breakpoint 11 (canonical), shared boundary 20.
        assert result == [("chr1", 1, 11), ("chr1", 11, 20), ("chr1", 20, 30)]

    def test_disjoin_with_cte_target_always_emits_canonical_zero_based_output(self):
        """Test that ``disjoin_*`` columns are emitted canonical regardless of source.

        Given:
            Two equivalent setups — a canonical ``features`` table and a
            custom-column ``features_custom`` table whose CTE aliases columns
            back to canonical.
        When:
            Running a self-mode DISJOIN on each.
        Then:
            ``disjoin_start`` / ``disjoin_end`` should be identical 0-based
            half-open values across both runs.
        """
        # Arrange
        from giql.table import Table

        projection = (
            "SELECT DISTINCT disjoin_start, disjoin_end {body} ORDER BY disjoin_start"
        )

        # Act
        canonical = _run_duckdb(
            "WITH x AS (SELECT * FROM features) "
            + projection.format(body="FROM DISJOIN(x)"),
            tables=["features"],
            features=[("chr1", 0, 20, "A"), ("chr1", 10, 30, "B")],
        )
        custom_aliased = _run_duckdb(
            'WITH x AS (SELECT seqid AS chrom, lo AS "start", hi AS "end" '
            "FROM features_custom) " + projection.format(body="FROM DISJOIN(x)"),
            tables=[
                "features",
                Table(
                    "features_custom",
                    chrom_col="seqid",
                    start_col="lo",
                    end_col="hi",
                ),
            ],
            schemas={
                "features_custom": (
                    "(seqid VARCHAR, lo INTEGER, hi INTEGER, label VARCHAR)"
                ),
            },
            features=[("chr1", 0, 20, "A"), ("chr1", 10, 30, "B")],
            features_custom=[("chr1", 0, 20, "A"), ("chr1", 10, 30, "B")],
        )

        # Assert
        assert canonical == custom_aliased == [(0, 10), (10, 20), (20, 30)]

    def test_disjoin_with_type_cast_in_cte_target_succeeds(self, run):
        """Test that explicit type casting inside the CTE does not break DISJOIN.

        Given:
            A CTE that casts ``start`` and ``end`` to ``BIGINT``.
        When:
            Running a self-mode DISJOIN.
        Then:
            The query should execute and yield the expected partition.
        """
        # Arrange & act
        result = run(
            'WITH x AS (SELECT chrom, CAST("start" AS BIGINT) AS "start", '
            'CAST("end" AS BIGINT) AS "end", name FROM features) '
            "SELECT DISTINCT disjoin_start, disjoin_end FROM DISJOIN(x) "
            "ORDER BY disjoin_start",
            tables=["features"],
            features=[("chr1", 0, 20, "A"), ("chr1", 10, 30, "B")],
        )

        # Assert
        assert result == [(0, 10), (10, 20), (20, 30)]


class TestDisjoinCteTargetFiltering:
    """Execution tests for CTE targets that filter, join, or aggregate input rows."""

    def test_disjoin_with_cte_where_filter_drops_unmatched_rows(self, run):
        """Test that a CTE WHERE filter removes rows before DISJOIN sees them.

        Given:
            A multi-chromosome ``features`` table.
        When:
            A CTE filters to one chromosome and DISJOIN runs over it.
        Then:
            Only rows from the filtered chromosome partition appear.
        """
        # Arrange & act
        result = run(
            "WITH x AS (SELECT * FROM features WHERE chrom = 'chr1') "
            "SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(x) ORDER BY disjoin_start",
            tables=["features"],
            features=[
                ("chr1", 0, 20, "A"),
                ("chr1", 10, 30, "B"),
                ("chr2", 0, 40, "C"),
            ],
        )

        # Assert
        assert result == [("chr1", 0, 10), ("chr1", 10, 20), ("chr1", 20, 30)]

    def test_disjoin_with_cte_interval_filter_drops_short_rows(self, run):
        """Test that a CTE filtering by interval length excludes short rows.

        Given:
            Three rows with various lengths; the CTE keeps only the long ones.
        When:
            DISJOIN runs over the CTE.
        Then:
            Only breakpoints from long rows contribute.
        """
        # Arrange & act
        result = run(
            'WITH x AS (SELECT * FROM features WHERE "end" - "start" > 50) '
            "SELECT DISTINCT disjoin_start, disjoin_end "
            "FROM DISJOIN(x) ORDER BY disjoin_start",
            tables=["features"],
            features=[
                ("chr1", 0, 10, "short"),
                ("chr1", 5, 80, "long1"),
                ("chr1", 40, 100, "long2"),
            ],
        )

        # Assert
        assert result == [(5, 40), (40, 80), (80, 100)]

    def test_disjoin_with_cte_joining_two_tables(self, run):
        """Test that DISJOIN works on a CTE that joins two interval tables.

        Given:
            A CTE that joins ``features`` and ``annotations`` on chromosome
            overlap.
        When:
            DISJOIN runs over the CTE.
        Then:
            Only overlap-joined rows contribute breakpoints.
        """
        # Arrange & act
        result = run(
            "WITH x AS ("
            'SELECT a.chrom, a."start", a."end", a.name '
            "FROM features AS a JOIN annotations AS b "
            'ON a.chrom = b.chrom AND a."start" < b."end" AND a."end" > b."start"'
            ") SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(x) ORDER BY disjoin_chrom, disjoin_start",
            tables=["features", "annotations"],
            features=[
                ("chr1", 0, 100, "f1"),
                ("chr1", 200, 300, "f2"),
                ("chr2", 0, 50, "f3"),
            ],
            annotations=[
                ("chr1", 10, 30, "a1"),
                ("chr1", 250, 260, "a2"),
                ("chr2", 100, 200, "a3"),
            ],
        )

        # Assert
        assert result == [("chr1", 0, 100), ("chr1", 200, 300)]

    def test_disjoin_with_cte_unioning_two_tables(self, run):
        """Test that DISJOIN works on a CTE built via UNION ALL.

        Given:
            Two interval tables combined via UNION ALL inside the CTE.
        When:
            DISJOIN runs over the CTE.
        Then:
            The output partitions the combined interval set.
        """
        # Arrange & act
        result = run(
            "WITH x AS ("
            'SELECT chrom, "start", "end", name FROM features_a '
            'UNION ALL SELECT chrom, "start", "end", name FROM features_b'
            ") SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(x) ORDER BY disjoin_start",
            tables=["features_a", "features_b"],
            features_a=[("chr1", 0, 20, "a1")],
            features_b=[("chr1", 10, 30, "b1")],
        )

        # Assert
        assert result == [("chr1", 0, 10), ("chr1", 10, 20), ("chr1", 20, 30)]

    def test_disjoin_with_cte_aggregating_rows(self, run):
        """Test that DISJOIN works on a CTE that aggregates rows per chromosome.

        Given:
            A CTE that groups rows by ``chrom`` and reports the min/max.
        When:
            DISJOIN runs over the aggregate CTE.
        Then:
            Each merged row passes through unchanged in self-mode (no
            interior breakpoints).
        """
        # Arrange & act
        result = run(
            "WITH x AS ("
            'SELECT chrom, MIN("start") AS "start", MAX("end") AS "end", '
            "'merged' AS name FROM features GROUP BY chrom"
            ") SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(x) ORDER BY disjoin_chrom, disjoin_start",
            tables=["features"],
            features=[
                ("chr1", 0, 20, "a"),
                ("chr1", 30, 50, "b"),
                ("chr2", 100, 150, "c"),
                ("chr2", 200, 250, "d"),
            ],
        )

        # Assert
        assert result == [("chr1", 0, 50), ("chr2", 100, 250)]

    def test_disjoin_with_cte_ordered_and_limited(self, run):
        """Test that ORDER BY + LIMIT inside the CTE restricts DISJOIN's input.

        Given:
            A CTE that orders by ``start`` and limits to the first 2 rows.
        When:
            DISJOIN runs over the limited CTE.
        Then:
            Only the first two intervals contribute breakpoints.
        """
        # Arrange & act
        result = run(
            'WITH x AS (SELECT * FROM features ORDER BY "start" LIMIT 2) '
            "SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(x) ORDER BY disjoin_start",
            tables=["features"],
            features=[
                ("chr1", 0, 20, "a"),
                ("chr1", 10, 30, "b"),
                ("chr1", 50, 80, "c"),
                ("chr1", 60, 90, "d"),
            ],
        )

        # Assert
        assert result == [("chr1", 0, 10), ("chr1", 10, 20), ("chr1", 20, 30)]

    def test_disjoin_called_twice_over_same_cte_produces_identical_partitions(self, run):
        """Test that two DISJOIN calls over the same CTE return identical rows.

        Given:
            A query that calls DISJOIN twice over the same CTE inside a
            UNION ALL with a tag column.
        When:
            Running the query.
        Then:
            Both tagged groups should hold the same partition rows.
        """
        # Arrange & act
        result = run(
            "WITH x AS (SELECT * FROM features) "
            "SELECT 'first' AS tag, disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(x) "
            "UNION ALL "
            "SELECT 'second', disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(x) "
            "ORDER BY tag, disjoin_start, disjoin_end",
            tables=["features"],
            features=[("chr1", 0, 20, "A"), ("chr1", 10, 30, "B")],
        )

        # Assert — same row set under each tag (parents repeat across the
        # two distinct sub-intervals their breakpoints produce).
        first = [tuple(row[1:]) for row in result if row[0] == "first"]
        second = [tuple(row[1:]) for row in result if row[0] == "second"]
        assert first == second
        assert sorted(set(first)) == [
            ("chr1", 0, 10),
            ("chr1", 10, 20),
            ("chr1", 20, 30),
        ]


_CANONICAL_FEATURES = [
    ("chr1", 0, 100, "a"),
    ("chr1", 60, 180, "b"),
    ("chr2", 10, 40, "c"),
]


class TestDisjoinCteComposition:
    """Execution tests for a CTE-target DISJOIN composed with downstream SQL."""

    def test_outer_where_filters_cte_target_disjoin_output(self, run):
        """Test that an outer WHERE filters segments derived from the CTE target.

        Given:
            A canonical ``features`` table wrapped in a CTE.
        When:
            An outer WHERE filters segments by length.
        Then:
            Only segments matching the predicate appear.
        """
        # Arrange & act
        result = run(
            'WITH x AS (SELECT chrom, "start", "end" FROM features) '
            "SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(x) WHERE disjoin_end - disjoin_start >= 50 "
            "ORDER BY disjoin_chrom, disjoin_start",
            tables=["features"],
            features=_CANONICAL_FEATURES,
        )

        # Assert
        assert result == [("chr1", 0, 60), ("chr1", 100, 180)]

    def test_outer_order_by_limit_applied_to_cte_target_disjoin(self, run):
        """Test that outer ORDER BY / LIMIT applies to CTE-target DISJOIN output.

        Given:
            A canonical CTE-target DISJOIN.
        When:
            Outer ORDER BY + LIMIT 1 selects the first segment.
        Then:
            Only the leftmost segment appears.
        """
        # Arrange & act
        result = run(
            'WITH x AS (SELECT chrom, "start", "end" FROM features) '
            "SELECT disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(x) "
            "ORDER BY disjoin_start, disjoin_chrom, disjoin_end LIMIT 1",
            tables=["features"],
            features=_CANONICAL_FEATURES,
        )

        # Assert
        assert result == [("chr1", 0, 60)]

    def test_outer_group_by_aggregates_cte_target_disjoin_output(self, run):
        """Test that GROUP BY aggregates segments by chromosome.

        Given:
            A canonical CTE-target DISJOIN.
        When:
            Outer GROUP BY counts segments per chromosome.
        Then:
            Per-chromosome counts reflect the DISJOIN partition.
        """
        # Arrange & act
        result = run(
            'WITH x AS (SELECT chrom, "start", "end" FROM features) '
            "SELECT disjoin_chrom, COUNT(*) AS n FROM DISJOIN(x) "
            "GROUP BY disjoin_chrom ORDER BY disjoin_chrom",
            tables=["features"],
            features=_CANONICAL_FEATURES,
        )

        # Assert
        assert result == [("chr1", 4), ("chr2", 1)]

    def test_outer_join_against_cte_target_disjoin_output(self, run):
        """Test that segments can be overlap-joined against another table.

        Given:
            A canonical CTE-target DISJOIN and an ``annotations`` table.
        When:
            Outer JOIN on chromosome + interval overlap.
        Then:
            Only segments overlapping an annotation appear, paired with the
            matching annotation name.
        """
        # Arrange & act
        result = run(
            'WITH x AS (SELECT chrom, "start", "end" FROM features) '
            "SELECT d.disjoin_chrom, d.disjoin_start, d.disjoin_end, a.name "
            "FROM DISJOIN(x) AS d JOIN annotations AS a "
            'ON a.chrom = d.disjoin_chrom AND a."start" < d.disjoin_end '
            'AND a."end" > d.disjoin_start '
            "ORDER BY d.disjoin_chrom, d.disjoin_start, a.name",
            tables=["features", "annotations"],
            features=_CANONICAL_FEATURES,
            annotations=[("chr1", 0, 30, "g1"), ("chr2", 0, 1000, "g2")],
        )

        # Assert
        assert result == [("chr1", 0, 60, "g1"), ("chr2", 10, 40, "g2")]

    def test_outer_scalar_aggregate_sums_cte_target_disjoin_lengths(self, run):
        """Test that an outer SUM aggregates segment lengths.

        Given:
            A canonical CTE-target DISJOIN.
        When:
            Outer SUM totals segment lengths.
        Then:
            The total equals the sum of segment widths.
        """
        # Arrange & act
        result = run(
            'WITH x AS (SELECT chrom, "start", "end" FROM features) '
            "SELECT SUM(disjoin_end - disjoin_start) AS total_bp FROM DISJOIN(x)",
            tables=["features"],
            features=_CANONICAL_FEATURES,
        )

        # Assert
        assert result == [(250,)]

    def test_cte_target_disjoin_nested_in_outer_cte(self, run):
        """Test that DISJOIN-on-CTE composes with a wrapping outer CTE.

        Given:
            A two-level CTE pipeline: inner CTE feeds DISJOIN, outer CTE
            re-projects the output.
        When:
            Selecting from the outer CTE.
        Then:
            Rows match the inner DISJOIN output.
        """
        # Arrange & act
        result = run(
            'WITH x AS (SELECT chrom, "start", "end" FROM features), '
            "y AS (SELECT * FROM DISJOIN(x)) "
            "SELECT disjoin_chrom, disjoin_start, disjoin_end FROM y "
            "ORDER BY disjoin_chrom, disjoin_start, disjoin_end",
            tables=["features"],
            features=_CANONICAL_FEATURES,
        )

        # Assert
        assert result == [
            ("chr1", 0, 60),
            ("chr1", 60, 100),
            ("chr1", 60, 100),
            ("chr1", 100, 180),
            ("chr2", 10, 40),
        ]

    def test_cte_target_disjoin_nested_in_derived_subquery(self, run):
        """Test that DISJOIN-on-CTE composes with a derived-subquery wrapper.

        Given:
            A derived subquery wrapping the CTE-target DISJOIN.
        When:
            Selecting from the derived subquery alias.
        Then:
            Rows match the inner DISJOIN output.
        """
        # Arrange & act
        result = run(
            'WITH x AS (SELECT chrom, "start", "end" FROM features) '
            "SELECT disjoin_chrom, disjoin_start, disjoin_end "
            "FROM (SELECT * FROM DISJOIN(x)) AS s "
            "ORDER BY disjoin_chrom, disjoin_start, disjoin_end",
            tables=["features"],
            features=_CANONICAL_FEATURES,
        )

        # Assert
        assert result == [
            ("chr1", 0, 60),
            ("chr1", 60, 100),
            ("chr1", 60, 100),
            ("chr1", 100, 180),
            ("chr2", 10, 40),
        ]

    def test_cte_target_disjoin_combined_in_union_arm(self, run):
        """Test that CTE-target DISJOIN combines with another arm via UNION ALL.

        Given:
            One UNION arm runs CTE-target DISJOIN; the other selects from
            ``annotations``.
        When:
            Combining via UNION ALL on chromosome.
        Then:
            All chromosomes from both arms appear.
        """
        # Arrange & act
        result = run(
            'WITH x AS (SELECT chrom, "start", "end" FROM features) '
            "SELECT disjoin_chrom AS c FROM DISJOIN(x) "
            "UNION ALL SELECT chrom AS c FROM annotations "
            "ORDER BY c",
            tables=["features", "annotations"],
            features=_CANONICAL_FEATURES,
            annotations=[("chr1", 0, 30, "g1"), ("chr2", 0, 1000, "g2")],
        )

        # Assert
        assert result == [
            ("chr1",),
            ("chr1",),
            ("chr1",),
            ("chr1",),
            ("chr1",),
            ("chr2",),
            ("chr2",),
        ]

    def test_cte_target_disjoin_with_window_function(self, run):
        """Test that a window function over DISJOIN segments works.

        Given:
            A canonical CTE-target DISJOIN.
        When:
            COUNT(*) OVER (PARTITION BY disjoin_chrom) is applied.
        Then:
            Each row carries the per-chromosome segment count.
        """
        # Arrange & act
        result = run(
            'WITH x AS (SELECT chrom, "start", "end" FROM features) '
            "SELECT disjoin_chrom, disjoin_start, disjoin_end, "
            "COUNT(*) OVER (PARTITION BY disjoin_chrom) AS chrom_n "
            "FROM DISJOIN(x) "
            "ORDER BY disjoin_chrom, disjoin_start, disjoin_end",
            tables=["features"],
            features=_CANONICAL_FEATURES,
        )

        # Assert
        assert result == [
            ("chr1", 0, 60, 4),
            ("chr1", 60, 100, 4),
            ("chr1", 60, 100, 4),
            ("chr1", 100, 180, 4),
            ("chr2", 10, 40, 1),
        ]

    def test_chained_disjoin_via_cte_indirection_is_idempotent_on_partition(self, run):
        """Test that DISJOIN over an already-disjoint CTE is idempotent.

        Given:
            An inner DISJOIN produces sub-intervals; an outer CTE re-projects
            them under canonical names; an outer DISJOIN runs on that CTE.
        When:
            Running the chained query.
        Then:
            The outer DISJOIN's output equals the inner partition row-by-row
            (idempotence on an already-disjoint set).
        """
        # Arrange & act
        result = run(
            "WITH first AS ("
            'SELECT disjoin_chrom AS chrom, disjoin_start AS "start", '
            'disjoin_end AS "end" FROM DISJOIN(features)'
            "), second AS (SELECT * FROM first) "
            "SELECT disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(second) "
            "ORDER BY disjoin_chrom, disjoin_start, disjoin_end",
            tables=["features"],
            features=[("chr1", 0, 20, "A"), ("chr1", 10, 30, "B")],
        )

        # Assert
        assert result == [
            ("chr1", 0, 10),
            ("chr1", 10, 20),
            ("chr1", 10, 20),
            ("chr1", 20, 30),
        ]

    def test_outer_group_by_having_filters_cte_target_disjoin_groups(self, run):
        """Test that HAVING filters group-by output of CTE-target DISJOIN.

        Given:
            A canonical CTE-target DISJOIN.
        When:
            Outer GROUP BY + HAVING COUNT(*) > 1.
        Then:
            Only chromosomes with multiple segments appear.
        """
        # Arrange & act
        result = run(
            'WITH x AS (SELECT chrom, "start", "end" FROM features) '
            "SELECT disjoin_chrom, COUNT(*) AS n FROM DISJOIN(x) "
            "GROUP BY disjoin_chrom HAVING COUNT(*) > 1 ORDER BY disjoin_chrom",
            tables=["features"],
            features=_CANONICAL_FEATURES,
        )

        # Assert
        assert result == [("chr1", 4)]
