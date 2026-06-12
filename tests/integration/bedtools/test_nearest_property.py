"""Property-based correctness tests for the NEAREST distance shift (#134).

These tests use Hypothesis to generate random reference intervals and
candidate sets, then verify three invariants of the bedtools-parity +1
applied to NEAREST: the reported distance equals the half-open gap plus
one (zero for overlaps), the k-nearest *selection* is unchanged by the
uniform monotonic shift, and ``max_distance`` filters on the shifted
distance.
"""

import pytest
from hypothesis import HealthCheck
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from giql import transpile

from .utils.data_models import GenomicInterval
from .utils.duckdb_loader import load_intervals

duckdb = pytest.importorskip("duckdb")


def _expected_distance(ref: GenomicInterval, cand: GenomicInterval) -> int:
    """Compute the bedtools-parity distance between two half-open intervals.

    Overlapping intervals are 0; otherwise the raw half-open gap plus one.
    """
    if ref.start < cand.end and ref.end > cand.start:
        return 0
    if cand.start >= ref.end:
        return cand.start - ref.end + 1
    return ref.start - cand.end + 1


@st.composite
def _candidates(draw, min_size=1, max_size=20):
    """Generate a list of same-chromosome candidate intervals with unique names."""
    n = draw(st.integers(min_value=min_size, max_value=max_size))
    candidates = []
    for i in range(n):
        start = draw(st.integers(min_value=0, max_value=120_000))
        length = draw(st.integers(min_value=1, max_value=5_000))
        candidates.append(
            GenomicInterval("chr1", start, start + length, f"c{i}", 100, "+")
        )
    return candidates


@st.composite
def _reference(draw):
    """Generate a single canonical (0-based half-open) reference interval."""
    start = draw(st.integers(min_value=0, max_value=100_000))
    length = draw(st.integers(min_value=1, max_value=5_000))
    return GenomicInterval("chr1", start, start + length, "ref", 100, "+")


def _run_nearest(reference, candidates, k, max_distance=None):
    """Run standalone NEAREST in DuckDB and return [(name, distance), ...]."""
    conn = duckdb.connect(":memory:")
    try:
        load_intervals(conn, "candidates", [c.to_tuple() for c in candidates])
        max_clause = "" if max_distance is None else f", max_distance := {max_distance}"
        sql = transpile(
            f"""
            SELECT name, distance
            FROM NEAREST(
                candidates,
                reference := 'chr1:{reference.start}-{reference.end}',
                k := {k}{max_clause}
            )
            """,
            tables=["candidates"],
        )
        return conn.execute(sql).fetchall()
    finally:
        conn.close()


@given(reference=_reference(), candidates=_candidates())
@settings(
    max_examples=50,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_nearest_reported_distance_equals_gap_plus_one(reference, candidates):
    """
    GIVEN a random reference interval and a random candidate set
    WHEN NEAREST returns every candidate with its distance
    THEN each reported distance equals the half-open gap + 1 (0 for overlaps)
    """
    by_name = {c.name: c for c in candidates}

    rows = _run_nearest(reference, candidates, k=len(candidates))

    for name, distance in rows:
        assert distance == _expected_distance(reference, by_name[name]), (
            f"{name}: reported {distance}, expected "
            f"{_expected_distance(reference, by_name[name])}"
        )


@given(data=st.data(), reference=_reference(), candidates=_candidates())
@settings(
    max_examples=50,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_nearest_k_selection_is_unchanged_by_the_shift(data, reference, candidates):
    """
    GIVEN a random reference, candidate set, and k
    WHEN NEAREST returns the k nearest candidates
    THEN the multiset of reported distances equals the k smallest distances by
        raw half-open gap, so the uniform +1 does not change which candidates
        are selected
    """
    k = data.draw(st.integers(min_value=1, max_value=len(candidates)))
    expected_sorted = sorted(_expected_distance(reference, c) for c in candidates)

    rows = _run_nearest(reference, candidates, k=k)

    reported_sorted = sorted(distance for _, distance in rows)
    assert reported_sorted == expected_sorted[:k]


@given(
    data=st.data(),
    reference=_reference(),
    candidates=_candidates(),
)
@settings(
    max_examples=50,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_nearest_max_distance_filters_on_the_shifted_distance(
    data, reference, candidates
):
    """
    GIVEN a random reference, candidate set, and threshold
    WHEN NEAREST is run with that max_distance over all candidates
    THEN the returned candidates are exactly those whose gap + 1 is <= the
        threshold, so no candidate at threshold + 1 ever survives
    """
    expected = {c.name: _expected_distance(reference, c) for c in candidates}
    threshold = data.draw(st.integers(min_value=0, max_value=max(expected.values()) + 5))

    rows = _run_nearest(reference, candidates, k=len(candidates), max_distance=threshold)

    returned = {name for name, _ in rows}
    survivors = {name for name, dist in expected.items() if dist <= threshold}
    assert returned == survivors
