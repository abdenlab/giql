"""Integration tests asserting DISTANCE is convention-invariant.

For a fixed pair of canonical 0-based half-open intervals, encoding the same
logical pair under any of the four ``(coordinate_system, interval_type)``
combinations -- on either or both sides of the join -- must produce the same
numeric distance through the GIQL transpiler and DuckDB.
"""

import pytest

from .encodings import CONVENTIONS
from .encodings import encode
from .encodings import make_table

pytestmark = pytest.mark.integration


def _row(chrom: str, start: int, end: int, name: str) -> tuple:
    """Build a 6-tuple suitable for ``load_intervals`` with strand ``+``."""
    return (chrom, start, end, name, 100, "+")


class TestDistanceCoordinateSpace:
    """Convention-invariance of DISTANCE end-to-end via DuckDB."""

    @pytest.mark.parametrize(
        ("coordinate_system", "interval_type"),
        CONVENTIONS,
        ids=lambda v: str(v),
    )
    def test_distance_returns_canonical_gap_when_both_sides_share_convention(
        self, giql_query, coordinate_system, interval_type
    ):
        """Test DISTANCE returns 101 for a 100 bp gap in any same-convention encoding.

        Given:
            Two non-overlapping intervals (canonical [100, 200) and
            [300, 400)) re-encoded under the same convention on both sides.
        When:
            DISTANCE(a.interval, b.interval) is executed.
        Then:
            It should return 101 (the canonical 0-based half-open gap of 100
            plus the bedtools closest -d parity + 1), regardless of how the
            intervals are stored.
        """
        # Arrange
        s_a, e_a = encode(100, 200, coordinate_system, interval_type)
        s_b, e_b = encode(300, 400, coordinate_system, interval_type)
        tables = [
            make_table("intervals_a", coordinate_system, interval_type),
            make_table("intervals_b", coordinate_system, interval_type),
        ]

        # Act
        result = giql_query(
            """
            SELECT DISTANCE(a.interval, b.interval) AS dist
            FROM intervals_a a, intervals_b b
            """,
            tables=tables,
            intervals_a=[_row("chr1", s_a, e_a, "a1")],
            intervals_b=[_row("chr1", s_b, e_b, "b1")],
        )

        # Assert
        assert result == [(101,)]

    @pytest.mark.parametrize(
        ("coordinate_system", "interval_type"),
        CONVENTIONS,
        ids=lambda v: str(v),
    )
    def test_distance_returns_zero_when_intervals_overlap(
        self, giql_query, coordinate_system, interval_type
    ):
        """Test DISTANCE returns 0 for overlapping intervals in any encoding.

        Given:
            Canonical intervals [100, 300) and [200, 400) re-encoded under
            the same convention on both sides.
        When:
            DISTANCE is executed.
        Then:
            It should return 0 regardless of encoding.
        """
        # Arrange
        s_a, e_a = encode(100, 300, coordinate_system, interval_type)
        s_b, e_b = encode(200, 400, coordinate_system, interval_type)
        tables = [
            make_table("intervals_a", coordinate_system, interval_type),
            make_table("intervals_b", coordinate_system, interval_type),
        ]

        # Act
        result = giql_query(
            """
            SELECT DISTANCE(a.interval, b.interval) AS dist
            FROM intervals_a a, intervals_b b
            """,
            tables=tables,
            intervals_a=[_row("chr1", s_a, e_a, "a1")],
            intervals_b=[_row("chr1", s_b, e_b, "b1")],
        )

        # Assert
        assert result == [(0,)]

    @pytest.mark.parametrize(
        ("coordinate_system", "interval_type"),
        CONVENTIONS,
        ids=lambda v: str(v),
    )
    def test_distance_returns_one_for_adjacent_half_open_boundary(
        self, giql_query, coordinate_system, interval_type
    ):
        """Test DISTANCE returns 1 when intervals touch at the half-open boundary.

        Given:
            Canonical adjacent intervals [100, 200) and [200, 300)
            re-encoded under the same convention on both sides.
        When:
            DISTANCE is executed.
        Then:
            It should return 1 (bedtools closest -d parity for book-ended
            features) regardless of encoding.
        """
        # Arrange
        s_a, e_a = encode(100, 200, coordinate_system, interval_type)
        s_b, e_b = encode(200, 300, coordinate_system, interval_type)
        tables = [
            make_table("intervals_a", coordinate_system, interval_type),
            make_table("intervals_b", coordinate_system, interval_type),
        ]

        # Act
        result = giql_query(
            """
            SELECT DISTANCE(a.interval, b.interval) AS dist
            FROM intervals_a a, intervals_b b
            """,
            tables=tables,
            intervals_a=[_row("chr1", s_a, e_a, "a1")],
            intervals_b=[_row("chr1", s_b, e_b, "b1")],
        )

        # Assert
        assert result == [(1,)]

    @pytest.mark.parametrize(
        ("conv_a", "conv_b"),
        [(a, b) for a in CONVENTIONS for b in CONVENTIONS],
        ids=lambda v: str(v),
    )
    def test_distance_is_invariant_under_mixed_conventions(
        self, giql_query, conv_a, conv_b
    ):
        """Test DISTANCE returns 101 for any pair of (convention_a, convention_b).

        Given:
            Canonical intervals [100, 200) and [300, 400) encoded in
            convention_a on table A and convention_b on table B.
        When:
            DISTANCE is executed.
        Then:
            It should return 101 (the 100 bp half-open gap plus the bedtools
            parity + 1) for every (convention_a, convention_b) pair,
            demonstrating that each side is canonicalized independently.
        """
        # Arrange
        s_a, e_a = encode(100, 200, *conv_a)
        s_b, e_b = encode(300, 400, *conv_b)
        tables = [
            make_table("intervals_a", *conv_a),
            make_table("intervals_b", *conv_b),
        ]

        # Act
        result = giql_query(
            """
            SELECT DISTANCE(a.interval, b.interval) AS dist
            FROM intervals_a a, intervals_b b
            """,
            tables=tables,
            intervals_a=[_row("chr1", s_a, e_a, "a1")],
            intervals_b=[_row("chr1", s_b, e_b, "b1")],
        )

        # Assert
        assert result == [(101,)]

    @pytest.mark.parametrize(
        ("coordinate_system", "interval_type"),
        CONVENTIONS,
        ids=lambda v: str(v),
    )
    def test_distance_returns_null_for_cross_chromosome_pair(
        self, giql_query, coordinate_system, interval_type
    ):
        """Test DISTANCE returns NULL across chromosomes regardless of encoding.

        Given:
            Canonical [100, 200) on chr1 and [100, 200) on chr2 re-encoded
            under the same convention on both sides.
        When:
            DISTANCE is executed.
        Then:
            It should return NULL regardless of encoding.
        """
        # Arrange
        s_a, e_a = encode(100, 200, coordinate_system, interval_type)
        s_b, e_b = encode(100, 200, coordinate_system, interval_type)
        tables = [
            make_table("intervals_a", coordinate_system, interval_type),
            make_table("intervals_b", coordinate_system, interval_type),
        ]

        # Act
        result = giql_query(
            """
            SELECT DISTANCE(a.interval, b.interval) AS dist
            FROM intervals_a a, intervals_b b
            """,
            tables=tables,
            intervals_a=[_row("chr1", s_a, e_a, "a1")],
            intervals_b=[_row("chr2", s_b, e_b, "b1")],
        )

        # Assert
        assert result == [(None,)]

    @pytest.mark.parametrize(
        ("coordinate_system", "interval_type"),
        CONVENTIONS,
        ids=lambda v: str(v),
    )
    def test_signed_distance_returns_positive_when_b_is_downstream(
        self, giql_query, coordinate_system, interval_type
    ):
        """Test signed DISTANCE returns +101 when B is downstream in any encoding.

        Given:
            Canonical A=[100, 200) and B=[300, 400) re-encoded under the
            same convention on both sides, queried with ``signed := true``.
        When:
            DISTANCE runs.
        Then:
            It should return +101 (100 bp gap + bedtools parity 1) regardless
            of encoding.
        """
        # Arrange
        s_a, e_a = encode(100, 200, coordinate_system, interval_type)
        s_b, e_b = encode(300, 400, coordinate_system, interval_type)
        tables = [
            make_table("intervals_a", coordinate_system, interval_type),
            make_table("intervals_b", coordinate_system, interval_type),
        ]

        # Act
        result = giql_query(
            """
            SELECT DISTANCE(a.interval, b.interval, signed := true) AS dist
            FROM intervals_a a, intervals_b b
            """,
            tables=tables,
            intervals_a=[_row("chr1", s_a, e_a, "a1")],
            intervals_b=[_row("chr1", s_b, e_b, "b1")],
        )

        # Assert
        assert result == [(101,)]

    @pytest.mark.parametrize(
        ("coordinate_system", "interval_type"),
        CONVENTIONS,
        ids=lambda v: str(v),
    )
    def test_signed_distance_returns_negative_when_b_is_upstream(
        self, giql_query, coordinate_system, interval_type
    ):
        """Test signed DISTANCE returns -101 when B is upstream in any encoding.

        Given:
            Canonical A=[300, 400) and B=[100, 200) re-encoded under the
            same convention on both sides, queried with ``signed := true``.
        When:
            DISTANCE runs.
        Then:
            It should return -101 (100 bp gap + bedtools parity 1, negated)
            regardless of encoding.
        """
        # Arrange
        s_a, e_a = encode(300, 400, coordinate_system, interval_type)
        s_b, e_b = encode(100, 200, coordinate_system, interval_type)
        tables = [
            make_table("intervals_a", coordinate_system, interval_type),
            make_table("intervals_b", coordinate_system, interval_type),
        ]

        # Act
        result = giql_query(
            """
            SELECT DISTANCE(a.interval, b.interval, signed := true) AS dist
            FROM intervals_a a, intervals_b b
            """,
            tables=tables,
            intervals_a=[_row("chr1", s_a, e_a, "a1")],
            intervals_b=[_row("chr1", s_b, e_b, "b1")],
        )

        # Assert
        assert result == [(-101,)]

    @pytest.mark.parametrize(
        ("conv_a", "conv_b"),
        [(a, b) for a in CONVENTIONS for b in CONVENTIONS],
        ids=lambda v: str(v),
    )
    def test_distance_book_ended_is_one_under_mixed_conventions(
        self, giql_query, conv_a, conv_b
    ):
        """Test DISTANCE returns 1 for a book-ended pair for any convention pair.

        Given:
            Canonical book-ended intervals [100, 200) and [200, 300) encoded in
            convention_a on table A and convention_b on table B.
        When:
            DISTANCE is executed.
        Then:
            It should return 1 (bedtools parity for touching features) for every
            (convention_a, convention_b) pair, proving the +1 lands on the
            canonical gap regardless of how each side is stored.
        """
        # Arrange
        s_a, e_a = encode(100, 200, *conv_a)
        s_b, e_b = encode(200, 300, *conv_b)
        tables = [
            make_table("intervals_a", *conv_a),
            make_table("intervals_b", *conv_b),
        ]

        # Act
        result = giql_query(
            """
            SELECT DISTANCE(a.interval, b.interval) AS dist
            FROM intervals_a a, intervals_b b
            """,
            tables=tables,
            intervals_a=[_row("chr1", s_a, e_a, "a1")],
            intervals_b=[_row("chr1", s_b, e_b, "b1")],
        )

        # Assert
        assert result == [(1,)]

    @pytest.mark.parametrize(
        ("coordinate_system", "interval_type"),
        CONVENTIONS,
        ids=lambda v: str(v),
    )
    def test_distance_one_bp_gap_returns_two_in_any_encoding(
        self, giql_query, coordinate_system, interval_type
    ):
        """Test DISTANCE returns 2 for a 1 bp gap in any same-convention encoding.

        Given:
            Canonical intervals [100, 200) and [201, 300) (a 1 bp half-open
            gap) re-encoded under the same convention on both sides.
        When:
            DISTANCE is executed.
        Then:
            It should return 2 (half-open gap 1 + bedtools parity 1) regardless
            of encoding.
        """
        # Arrange
        s_a, e_a = encode(100, 200, coordinate_system, interval_type)
        s_b, e_b = encode(201, 300, coordinate_system, interval_type)
        tables = [
            make_table("intervals_a", coordinate_system, interval_type),
            make_table("intervals_b", coordinate_system, interval_type),
        ]

        # Act
        result = giql_query(
            """
            SELECT DISTANCE(a.interval, b.interval) AS dist
            FROM intervals_a a, intervals_b b
            """,
            tables=tables,
            intervals_a=[_row("chr1", s_a, e_a, "a1")],
            intervals_b=[_row("chr1", s_b, e_b, "b1")],
        )

        # Assert
        assert result == [(2,)]

    @pytest.mark.parametrize(
        ("conv_a", "conv_b"),
        [(a, b) for a in CONVENTIONS for b in CONVENTIONS],
        ids=lambda v: str(v),
    )
    def test_distance_one_bp_gap_is_two_under_mixed_conventions(
        self, giql_query, conv_a, conv_b
    ):
        """Test DISTANCE returns 2 for a 1 bp gap for any convention pair.

        Given:
            Canonical intervals [100, 200) and [201, 300) (a 1 bp half-open
            gap) encoded in convention_a on table A and convention_b on table B.
        When:
            DISTANCE is executed.
        Then:
            It should return 2 for every (convention_a, convention_b) pair, the
            tightest parity check across independent per-side canonicalization.
        """
        # Arrange
        s_a, e_a = encode(100, 200, *conv_a)
        s_b, e_b = encode(201, 300, *conv_b)
        tables = [
            make_table("intervals_a", *conv_a),
            make_table("intervals_b", *conv_b),
        ]

        # Act
        result = giql_query(
            """
            SELECT DISTANCE(a.interval, b.interval) AS dist
            FROM intervals_a a, intervals_b b
            """,
            tables=tables,
            intervals_a=[_row("chr1", s_a, e_a, "a1")],
            intervals_b=[_row("chr1", s_b, e_b, "b1")],
        )

        # Assert
        assert result == [(2,)]

    @pytest.mark.parametrize(
        ("coordinate_system", "interval_type"),
        CONVENTIONS,
        ids=lambda v: str(v),
    )
    def test_stranded_distance_book_ended_is_one_in_any_encoding(
        self, giql_query, coordinate_system, interval_type
    ):
        """Test stranded DISTANCE returns 1 for a book-ended pair in any encoding.

        Given:
            Canonical book-ended intervals [100, 200) and [200, 300), both on
            '+' strand, re-encoded under the same convention on both sides.
        When:
            DISTANCE with stranded := true is executed.
        Then:
            It should return 1 regardless of encoding, exercising the parity +1
            in the stranded branch under coordinate canonicalization.
        """
        # Arrange
        s_a, e_a = encode(100, 200, coordinate_system, interval_type)
        s_b, e_b = encode(200, 300, coordinate_system, interval_type)
        tables = [
            make_table("intervals_a", coordinate_system, interval_type),
            make_table("intervals_b", coordinate_system, interval_type),
        ]

        # Act
        result = giql_query(
            """
            SELECT DISTANCE(a.interval, b.interval, stranded := true) AS dist
            FROM intervals_a a, intervals_b b
            """,
            tables=tables,
            intervals_a=[_row("chr1", s_a, e_a, "a1")],
            intervals_b=[_row("chr1", s_b, e_b, "b1")],
        )

        # Assert
        assert result == [(1,)]
