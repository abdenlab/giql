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
    def test_disjoin_should_yield_the_partition_in_the_target_encoding(
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
    def test_disjoin_should_yield_target_encoded_subintervals_when_reference_differs(
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
