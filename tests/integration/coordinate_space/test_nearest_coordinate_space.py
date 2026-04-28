"""Integration tests asserting NEAREST is convention-invariant.

For a fixed canonical configuration, encoding the same logical intervals
under any of the four ``(coordinate_system, interval_type)`` combinations
must select the same neighbor row(s) through GIQL transpilation and DuckDB.
"""

import pytest

from .encodings import CONVENTIONS
from .encodings import encode
from .encodings import make_table

pytestmark = pytest.mark.integration


def _row(chrom: str, start: int, end: int, name: str) -> tuple:
    """Build a 6-tuple suitable for ``load_intervals`` with strand ``+``."""
    return (chrom, start, end, name, 100, "+")


# Canonical layout reused by NC-001 / NC-005 to keep the matched-name baseline
# stable across same- and mixed-convention parametrizations.
#
#   A: [200, 300)
#   B candidates: b_far  [100, 150) -> gap 50
#                 b_near [310, 350) -> gap 10  <-- nearest
#                 b_mid  [500, 600) -> gap 200
_CANONICAL_A = [(200, 300, "a1")]
_CANONICAL_B = [
    (100, 150, "b_far"),
    (310, 350, "b_near"),
    (500, 600, "b_mid"),
]


def _encode_rows(rows: list[tuple], coordinate_system: str, interval_type: str):
    """Re-encode canonical (start, end, name) tuples for the given convention."""
    encoded = []
    for canonical_start, canonical_end, name in rows:
        s, e = encode(canonical_start, canonical_end, coordinate_system, interval_type)
        encoded.append(_row("chr1", s, e, name))
    return encoded


class TestNearestCoordinateSpace:
    """Convention-invariance of NEAREST end-to-end via DuckDB."""

    @pytest.mark.parametrize(
        ("coordinate_system", "interval_type"),
        CONVENTIONS,
        ids=lambda v: str(v),
    )
    def test_correlated_nearest_picks_canonical_neighbor_in_any_encoding(
        self, giql_query, coordinate_system, interval_type
    ):
        """Test correlated NEAREST k=1 selects the same B row in any same-convention encoding.

        Given:
            Canonical A=[200, 300) and three B candidates re-encoded under
            the same convention on both sides.
        When:
            ``NEAREST(intervals_b, reference := a.interval, k := 1)`` runs
            via LATERAL join.
        Then:
            It should return ``b_near`` regardless of encoding.
        """
        # Arrange
        tables = [
            make_table("intervals_a", coordinate_system, interval_type),
            make_table("intervals_b", coordinate_system, interval_type),
        ]
        intervals_a = _encode_rows(_CANONICAL_A, coordinate_system, interval_type)
        intervals_b = _encode_rows(_CANONICAL_B, coordinate_system, interval_type)

        # Act
        result = giql_query(
            """
            SELECT a.name, b.name
            FROM intervals_a a
            CROSS JOIN LATERAL NEAREST(
                intervals_b,
                reference := a.interval,
                k := 1
            ) b
            """,
            tables=tables,
            intervals_a=intervals_a,
            intervals_b=intervals_b,
        )

        # Assert
        assert result == [("a1", "b_near")]

    @pytest.mark.parametrize(
        ("coordinate_system", "interval_type"),
        CONVENTIONS,
        ids=lambda v: str(v),
    )
    def test_correlated_nearest_k3_returns_all_three_candidates_in_any_encoding(
        self, giql_query, coordinate_system, interval_type
    ):
        """Test correlated NEAREST k=3 returns all candidates in any encoding.

        Given:
            Canonical A=[200, 300) and three B candidates re-encoded under
            the same convention on both sides.
        When:
            NEAREST with ``k := 3`` runs via LATERAL join.
        Then:
            It should return all three B rows regardless of encoding.
        """
        # Arrange
        tables = [
            make_table("intervals_a", coordinate_system, interval_type),
            make_table("intervals_b", coordinate_system, interval_type),
        ]
        intervals_a = _encode_rows(_CANONICAL_A, coordinate_system, interval_type)
        intervals_b = _encode_rows(_CANONICAL_B, coordinate_system, interval_type)

        # Act
        result = giql_query(
            """
            SELECT b.name
            FROM intervals_a a
            CROSS JOIN LATERAL NEAREST(
                intervals_b,
                reference := a.interval,
                k := 3
            ) b
            """,
            tables=tables,
            intervals_a=intervals_a,
            intervals_b=intervals_b,
        )

        # Assert
        assert {row[0] for row in result} == {"b_near", "b_far", "b_mid"}

    @pytest.mark.parametrize(
        ("coordinate_system", "interval_type"),
        CONVENTIONS,
        ids=lambda v: str(v),
    )
    def test_correlated_nearest_max_distance_filters_by_canonical_gap(
        self, giql_query, coordinate_system, interval_type
    ):
        """Test correlated NEAREST honors max_distance against the canonical gap.

        Given:
            Canonical A=[200, 300) and the same three B candidates,
            re-encoded under the same convention on both sides; ``k := 3``
            and ``max_distance := 50`` is supplied.
        When:
            NEAREST runs via LATERAL join.
        Then:
            It should return only ``b_near`` (10 bp) and ``b_far`` (50 bp,
            inclusive boundary) regardless of encoding; ``b_mid`` (200 bp)
            should be excluded.
        """
        # Arrange
        tables = [
            make_table("intervals_a", coordinate_system, interval_type),
            make_table("intervals_b", coordinate_system, interval_type),
        ]
        intervals_a = _encode_rows(_CANONICAL_A, coordinate_system, interval_type)
        intervals_b = _encode_rows(_CANONICAL_B, coordinate_system, interval_type)

        # Act
        result = giql_query(
            """
            SELECT b.name
            FROM intervals_a a
            CROSS JOIN LATERAL NEAREST(
                intervals_b,
                reference := a.interval,
                k := 3,
                max_distance := 50
            ) b
            """,
            tables=tables,
            intervals_a=intervals_a,
            intervals_b=intervals_b,
        )

        # Assert
        assert sorted(row[0] for row in result) == ["b_far", "b_near"]

    @pytest.mark.parametrize(
        ("coordinate_system", "interval_type"),
        CONVENTIONS,
        ids=lambda v: str(v),
    )
    def test_standalone_nearest_with_literal_reference_resolves_across_encodings(
        self, giql_query, coordinate_system, interval_type
    ):
        """Test standalone NEAREST with a 0-based half-open literal works under any target encoding.

        Given:
            Three canonical B candidates re-encoded under the given
            convention; the literal reference ``'chr1:350-360'`` is parsed
            as 0-based half-open [350, 360).
        When:
            ``SELECT * FROM NEAREST(intervals_b, reference := ..., k := 2)``
            runs in standalone mode.
        Then:
            It should return ``b_near`` (10 bp) and ``b_far`` (200 bp)
            regardless of target encoding; ``b_mid`` (140 bp) is the
            second-nearest, so the expected pair is ``{b_near, b_mid}``.
        """
        # Arrange
        tables = [make_table("intervals_b", coordinate_system, interval_type)]
        intervals_b = _encode_rows(_CANONICAL_B, coordinate_system, interval_type)

        # Act
        result = giql_query(
            """
            SELECT name
            FROM NEAREST(
                intervals_b,
                reference := 'chr1:350-360',
                k := 2
            )
            """,
            tables=tables,
            intervals_b=intervals_b,
        )

        # Assert
        # Distances from [350, 360): b_near [310,350) -> 0 (touching),
        # b_mid [500,600) -> 140, b_far [100,150) -> 200.
        assert sorted(row[0] for row in result) == ["b_mid", "b_near"]

    @pytest.mark.parametrize(
        ("conv_a", "conv_b"),
        [(a, b) for a in CONVENTIONS for b in CONVENTIONS],
        ids=lambda v: str(v),
    )
    def test_correlated_nearest_is_invariant_under_mixed_conventions(
        self, giql_query, conv_a, conv_b
    ):
        """Test correlated NEAREST k=1 picks the same B row for any mixed-convention pair.

        Given:
            Canonical A=[200, 300) on table A encoded in convention_a, and
            three canonical B candidates on table B encoded in convention_b.
        When:
            NEAREST k=1 runs via LATERAL join.
        Then:
            It should return ``b_near`` for every (convention_a,
            convention_b) pair.
        """
        # Arrange
        tables = [
            make_table("intervals_a", *conv_a),
            make_table("intervals_b", *conv_b),
        ]
        intervals_a = _encode_rows(_CANONICAL_A, *conv_a)
        intervals_b = _encode_rows(_CANONICAL_B, *conv_b)

        # Act
        result = giql_query(
            """
            SELECT a.name, b.name
            FROM intervals_a a
            CROSS JOIN LATERAL NEAREST(
                intervals_b,
                reference := a.interval,
                k := 1
            ) b
            """,
            tables=tables,
            intervals_a=intervals_a,
            intervals_b=intervals_b,
        )

        # Assert
        assert result == [("a1", "b_near")]

    @pytest.mark.parametrize(
        ("coordinate_system", "interval_type"),
        CONVENTIONS,
        ids=lambda v: str(v),
    )
    def test_correlated_nearest_picks_adjacent_neighbor_when_touching_half_open_boundary(
        self, giql_query, coordinate_system, interval_type
    ):
        """Test correlated NEAREST picks the touching candidate (gap 0) in any encoding.

        Given:
            Canonical A=[100, 200) with two B candidates -- b_adjacent at
            canonical [200, 300) (touching half-open boundary, gap 0) and
            b_far at [500, 600) -- re-encoded under the same convention on
            both sides.
        When:
            NEAREST k=1 runs via LATERAL join.
        Then:
            It should return ``b_adjacent`` regardless of encoding.
        """
        # Arrange
        tables = [
            make_table("intervals_a", coordinate_system, interval_type),
            make_table("intervals_b", coordinate_system, interval_type),
        ]
        intervals_a = _encode_rows(
            [(100, 200, "a1")], coordinate_system, interval_type
        )
        intervals_b = _encode_rows(
            [(200, 300, "b_adjacent"), (500, 600, "b_far")],
            coordinate_system,
            interval_type,
        )

        # Act
        result = giql_query(
            """
            SELECT a.name, b.name
            FROM intervals_a a
            CROSS JOIN LATERAL NEAREST(
                intervals_b,
                reference := a.interval,
                k := 1
            ) b
            """,
            tables=tables,
            intervals_a=intervals_a,
            intervals_b=intervals_b,
        )

        # Assert
        assert result == [("a1", "b_adjacent")]
