"""Integration tests asserting DISJOIN is convention-invariant.

For a fixed set of canonical 0-based half-open intervals, encoding the same
logical intervals under any of the four ``(coordinate_system, interval_type)``
combinations -- on the target, the reference, or both -- must produce the same
logical partition. DISJOIN emits ``disjoin_start`` / ``disjoin_end`` in the
target table's coordinate system, so the expected rows are re-encoded under
that convention.
"""

import pytest

from .encodings import CONVENTIONS
from .encodings import encode
from .encodings import make_table

pytestmark = pytest.mark.integration


def _row(chrom: str, start: int, end: int, name: str) -> tuple:
    """Build a 6-tuple suitable for ``load_intervals`` with strand ``+``."""
    return (chrom, start, end, name, 100, "+")


class TestDisjoinCoordinateSpace:
    """Convention-invariance of DISJOIN end-to-end via DuckDB."""

    @pytest.mark.parametrize(
        ("coordinate_system", "interval_type"),
        CONVENTIONS,
        ids=lambda v: str(v),
    )
    def test_disjoin_yields_the_partition_in_the_target_encoding(
        self, giql_query, coordinate_system, interval_type
    ):
        """Test DISJOIN self-mode yields the partition in the target's encoding.

        Given:
            Two overlapping intervals (canonical [0, 20) and [10, 30))
            re-encoded under one convention.
        When:
            DISJOIN(features) runs in self-mode and DISTINCT sub-intervals
            are selected.
        Then:
            The sub-intervals should be the partition {[0, 10), [10, 20),
            [20, 30)} re-encoded under the target table's convention.
        """
        # Arrange
        s_a, e_a = encode(0, 20, coordinate_system, interval_type)
        s_b, e_b = encode(10, 30, coordinate_system, interval_type)
        tables = [make_table("features", coordinate_system, interval_type)]

        # Act
        result = giql_query(
            "SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(features) ORDER BY disjoin_start",
            tables=tables,
            features=[_row("chr1", s_a, e_a, "a"), _row("chr1", s_b, e_b, "b")],
        )

        # Assert
        expected = [
            ("chr1", *encode(0, 10, coordinate_system, interval_type)),
            ("chr1", *encode(10, 20, coordinate_system, interval_type)),
            ("chr1", *encode(20, 30, coordinate_system, interval_type)),
        ]
        assert result == expected

    @pytest.mark.parametrize(
        ("coordinate_system", "interval_type"),
        CONVENTIONS,
        ids=lambda v: str(v),
    )
    def test_disjoin_yields_target_encoded_subintervals_when_reference_differs(
        self, giql_query, coordinate_system, interval_type
    ):
        """Test DISJOIN emits target-encoded sub-intervals across mixed encodings.

        Given:
            A target interval (canonical [0, 30)) under one convention and a
            reference set (canonical [0, 10) and [10, 30)) under a fixed,
            different convention.
        When:
            DISJOIN(features, reference := refs) runs.
        Then:
            The sub-intervals {[0, 10), [10, 30)} should be re-encoded under
            the target table's convention regardless of either side's storage
            encoding.
        """
        # Arrange
        ref_system, ref_type = "1based", "closed"
        s_t, e_t = encode(0, 30, coordinate_system, interval_type)
        r1s, r1e = encode(0, 10, ref_system, ref_type)
        r2s, r2e = encode(10, 30, ref_system, ref_type)
        tables = [
            make_table("features", coordinate_system, interval_type),
            make_table("refs", ref_system, ref_type),
        ]

        # Act
        result = giql_query(
            "SELECT disjoin_start, disjoin_end "
            "FROM DISJOIN(features, reference := refs) ORDER BY disjoin_start",
            tables=tables,
            features=[_row("chr1", s_t, e_t, "t")],
            refs=[_row("chr1", r1s, r1e, "r1"), _row("chr1", r2s, r2e, "r2")],
        )

        # Assert
        expected = [
            encode(0, 10, coordinate_system, interval_type),
            encode(10, 30, coordinate_system, interval_type),
        ]
        assert result == expected

    @pytest.mark.parametrize(
        ("coordinate_system", "interval_type"),
        CONVENTIONS,
        ids=lambda v: str(v),
    )
    def test_disjoin_self_mode_and_explicit_self_reference_yield_identical_partitions(
        self, giql_query, coordinate_system, interval_type
    ):
        """Test omitted-reference and explicit-self-reference produce identical rows.

        Given:
            Two overlapping intervals (canonical [0, 20) and [10, 30))
            re-encoded under one convention.
        When:
            ``DISJOIN(features)`` and ``DISJOIN(features, reference := features)``
            are unioned (tagged) in one query.
        Then:
            Both tag groups should return the same partition
            {[0, 10), [10, 20), [20, 30)} re-encoded under the target
            convention — proving the two self-reference paths take the same
            optimized branch.
        """
        # Arrange
        s_a, e_a = encode(0, 20, coordinate_system, interval_type)
        s_b, e_b = encode(10, 30, coordinate_system, interval_type)
        tables = [make_table("features", coordinate_system, interval_type)]
        rows = [_row("chr1", s_a, e_a, "a"), _row("chr1", s_b, e_b, "b")]
        expected_partition = [
            ("chr1", *encode(0, 10, coordinate_system, interval_type)),
            ("chr1", *encode(10, 20, coordinate_system, interval_type)),
            ("chr1", *encode(20, 30, coordinate_system, interval_type)),
        ]

        # Act
        result = giql_query(
            "SELECT 'omitted' AS tag, disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(features) "
            "UNION ALL "
            "SELECT 'explicit' AS tag, disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(features, reference := features) "
            "ORDER BY tag, disjoin_start",
            tables=tables,
            features=rows,
        )

        # Assert
        omitted = sorted({row[1:] for row in result if row[0] == "omitted"})
        explicit = sorted({row[1:] for row in result if row[0] == "explicit"})
        assert omitted == expected_partition
        assert explicit == expected_partition

    def test_disjoin_drops_segments_outside_partial_reference_coverage(self, giql_query):
        """Test DISJOIN drops sub-intervals not covered by the reference set.

        Given:
            A target interval [0, 100) and a reference set covering only
            [0, 30) and [70, 100) (gap [30, 70)), canonical 0-based half-open.
        When:
            DISJOIN(features, reference := refs) runs.
        Then:
            The result should contain exactly {[0, 30), [70, 100)} and the
            uncovered [30, 70) segment should be dropped — proving the EXISTS
            coverage filter is preserved for distinct references.
        """
        # Arrange
        tables = [
            make_table("features", "0based", "half_open"),
            make_table("refs", "0based", "half_open"),
        ]

        # Act
        result = giql_query(
            "SELECT DISTINCT disjoin_start, disjoin_end "
            "FROM DISJOIN(features, reference := refs) ORDER BY disjoin_start",
            tables=tables,
            features=[_row("chr1", 0, 100, "t")],
            refs=[_row("chr1", 0, 30, "r1"), _row("chr1", 70, 100, "r2")],
        )

        # Assert
        assert result == [(0, 30), (70, 100)]

    def test_disjoin_drops_segments_outside_partial_reference_coverage_when_target_and_reference_use_different_encodings(
        self, giql_query
    ):
        """Test EXISTS canonicalization across mixed target/reference encodings.

        Given:
            A target interval canonical [0, 100) stored 0-based half-open and
            a partial-coverage reference {[0, 30), [70, 100)} stored 1-based
            closed.
        When:
            DISJOIN(features, reference := refs) runs.
        Then:
            The result should contain the same canonical surviving sub-intervals
            {[0, 30), [70, 100)} re-encoded under the target convention — the
            EXISTS body still canonicalizes the reference endpoints correctly.
        """
        # Arrange
        target_system, target_type = "0based", "half_open"
        ref_system, ref_type = "1based", "closed"
        s_t, e_t = encode(0, 100, target_system, target_type)
        r1s, r1e = encode(0, 30, ref_system, ref_type)
        r2s, r2e = encode(70, 100, ref_system, ref_type)
        tables = [
            make_table("features", target_system, target_type),
            make_table("refs", ref_system, ref_type),
        ]

        # Act
        result = giql_query(
            "SELECT DISTINCT disjoin_start, disjoin_end "
            "FROM DISJOIN(features, reference := refs) ORDER BY disjoin_start",
            tables=tables,
            features=[_row("chr1", s_t, e_t, "t")],
            refs=[_row("chr1", r1s, r1e, "r1"), _row("chr1", r2s, r2e, "r2")],
        )

        # Assert
        expected = [
            encode(0, 30, target_system, target_type),
            encode(70, 100, target_system, target_type),
        ]
        assert result == expected

    def test_disjoin_keeps_coverage_filter_when_cte_shadows_target_name(
        self, giql_query
    ):
        """Test CTE-shadowed self-reference still applies the coverage filter.

        Given:
            A registered ``features`` table with two intervals and an outer
            WITH clause defining a CTE ``features`` that is a strict subset
            (one row).
        When:
            ``WITH features AS (SELECT * FROM features WHERE name = 'a')
            SELECT ... FROM DISJOIN(features, reference := features)`` is
            unioned (tagged) with the registered-table self-mode partition.
        Then:
            The CTE-shadowed result should be just {[0, 50)} (one segment),
            strictly smaller than the registered-table partition
            {[0, 30), [30, 50), [50, 100)} — the CTE shadows both target and
            reference and EXISTS remains active.
        """
        # Arrange
        tables = [make_table("features", "0based", "half_open")]
        rows = [_row("chr1", 0, 50, "a"), _row("chr1", 30, 100, "b")]

        # Act
        result = giql_query(
            "SELECT 'shadowed' AS tag, disjoin_chrom, disjoin_start, disjoin_end "
            "FROM (WITH features AS (SELECT * FROM features WHERE name = 'a') "
            "SELECT * FROM DISJOIN(features, reference := features)) AS s "
            "UNION ALL "
            "SELECT 'unshadowed' AS tag, disjoin_chrom, disjoin_start, disjoin_end "
            "FROM DISJOIN(features) "
            "ORDER BY tag, disjoin_start",
            tables=tables,
            features=rows,
        )

        # Assert
        shadowed = sorted({row[1:] for row in result if row[0] == "shadowed"})
        unshadowed = sorted({row[1:] for row in result if row[0] == "unshadowed"})
        assert shadowed == [("chr1", 0, 50)]
        assert unshadowed == [("chr1", 0, 30), ("chr1", 30, 50), ("chr1", 50, 100)]

    def test_disjoin_keeps_coverage_filter_when_reference_is_partial_coverage_subquery(
        self, giql_query
    ):
        """Test a partial-coverage subquery reference still drops uncovered segments.

        Given:
            A target with one wide row [0, 100) and one narrow row [20, 50)
            on chr1; a subquery reference selecting only rows with
            ``"start" >= 20`` (which yields only the narrow row).
        When:
            ``DISJOIN(features, reference := (SELECT * FROM features WHERE
            "start" >= 20))`` runs.
        Then:
            Only the [20, 50) sub-interval should survive — proving subquery
            references retain the EXISTS coverage filter end-to-end.
        """
        # Arrange
        tables = [make_table("features", "0based", "half_open")]

        # Act
        result = giql_query(
            "SELECT DISTINCT disjoin_start, disjoin_end "
            "FROM DISJOIN(features, reference := ("
            'SELECT * FROM features WHERE "start" >= 20'
            ")) ORDER BY disjoin_start",
            tables=tables,
            features=[_row("chr1", 0, 100, "wide"), _row("chr1", 20, 50, "partial")],
        )

        # Assert
        assert result == [(20, 50)]
