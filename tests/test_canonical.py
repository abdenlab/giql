"""Unit tests for ``giql.canonical`` coordinate canonicalization helpers.

These tests pin down the externally observable contract of
:func:`giql.canonical.canonical_start` / :func:`giql.canonical.canonical_end`
and their inverses :func:`giql.canonical.decanonical_start` /
:func:`giql.canonical.decanonical_end` — pure functions that wrap a raw
column reference with the offset required to convert between a source
:class:`giql.table.Table`'s declared encoding and the canonical 0-based
half-open form the transpiler emits.
"""

import pytest

from giql.canonical import canonical_end
from giql.canonical import canonical_start
from giql.canonical import decanonical_end
from giql.canonical import decanonical_start
from giql.table import Table

hypothesis = pytest.importorskip("hypothesis")
from hypothesis import given  # noqa: E402
from hypothesis import settings  # noqa: E402
from hypothesis import strategies as st  # noqa: E402


def test_canonical_start_should_return_input_unchanged_when_table_is_none():
    """Test that a missing table leaves the start expression untouched.

    Given:
        A raw start expression and ``table=None``.
    When:
        ``canonical_start(raw, None)`` is called.
    Then:
        It should return the original raw expression unchanged.
    """
    # Arrange
    raw = "a.start"

    # Act
    result = canonical_start(raw, None)

    # Assert
    assert result == raw


def test_canonical_start_should_return_input_unchanged_when_table_is_zero_based():
    """Test that a 0-based table leaves the start expression untouched.

    Given:
        A :class:`giql.table.Table` with ``coordinate_system="0based"``.
    When:
        ``canonical_start(raw, table)`` is called.
    Then:
        It should return ``raw`` unchanged because the source already
        uses the canonical 0-based start convention.
    """
    # Arrange
    raw = "a.start"
    table = Table(name="features", coordinate_system="0based")

    # Act
    result = canonical_start(raw, table)

    # Assert
    assert result == raw


def test_canonical_start_should_subtract_one_when_table_is_one_based():
    """Test that a 1-based table subtracts one from the start expression.

    Given:
        A :class:`giql.table.Table` with ``coordinate_system="1based"``.
    When:
        ``canonical_start(raw, table)`` is called.
    Then:
        It should return ``f"({raw} - 1)"`` so the canonical 0-based
        half-open form is emitted downstream.
    """
    # Arrange
    raw = "a.start"
    table = Table(name="features", coordinate_system="1based")

    # Act
    result = canonical_start(raw, table)

    # Assert
    assert result == f"({raw} - 1)"


def test_canonical_end_should_return_input_unchanged_when_table_is_none():
    """Test that a missing table leaves the end expression untouched.

    Given:
        A raw end expression and ``table=None``.
    When:
        ``canonical_end(raw, None)`` is called.
    Then:
        It should return the original raw expression unchanged.
    """
    # Arrange
    raw = "a.end"

    # Act
    result = canonical_end(raw, None)

    # Assert
    assert result == raw


def test_canonical_end_should_return_input_unchanged_when_table_is_zero_based_half_open():
    """Test that 0-based half-open tables leave the end expression untouched.

    Given:
        A :class:`giql.table.Table` with ``coordinate_system="0based"``
        and ``interval_type="half_open"``.
    When:
        ``canonical_end(raw, table)`` is called.
    Then:
        It should return ``raw`` unchanged because the source already
        matches the canonical convention.
    """
    # Arrange
    raw = "a.end"
    table = Table(
        name="features",
        coordinate_system="0based",
        interval_type="half_open",
    )

    # Act
    result = canonical_end(raw, table)

    # Assert
    assert result == raw


def test_canonical_end_should_add_one_when_table_is_zero_based_closed():
    """Test that 0-based closed tables add one to the end expression.

    Given:
        A :class:`giql.table.Table` with ``coordinate_system="0based"``
        and ``interval_type="closed"``.
    When:
        ``canonical_end(raw, table)`` is called.
    Then:
        It should return ``f"({raw} + 1)"`` to convert closed to
        half-open.
    """
    # Arrange
    raw = "a.end"
    table = Table(
        name="features",
        coordinate_system="0based",
        interval_type="closed",
    )

    # Act
    result = canonical_end(raw, table)

    # Assert
    assert result == f"({raw} + 1)"


def test_canonical_end_should_subtract_one_when_table_is_one_based_half_open():
    """Test that 1-based half-open tables subtract one from the end expression.

    Given:
        A :class:`giql.table.Table` with ``coordinate_system="1based"``
        and ``interval_type="half_open"``.
    When:
        ``canonical_end(raw, table)`` is called.
    Then:
        It should return ``f"({raw} - 1)"`` to convert 1-based
        half-open to canonical 0-based half-open.
    """
    # Arrange
    raw = "a.end"
    table = Table(
        name="features",
        coordinate_system="1based",
        interval_type="half_open",
    )

    # Act
    result = canonical_end(raw, table)

    # Assert
    assert result == f"({raw} - 1)"


def test_canonical_end_should_return_input_unchanged_when_table_is_one_based_closed():
    """Test that 1-based closed tables leave the end expression untouched.

    Given:
        A :class:`giql.table.Table` with ``coordinate_system="1based"``
        and ``interval_type="closed"``.
    When:
        ``canonical_end(raw, table)`` is called.
    Then:
        It should return ``raw`` unchanged because a 1-based closed
        end equals a 0-based half-open end numerically.
    """
    # Arrange
    raw = "a.end"
    table = Table(
        name="features",
        coordinate_system="1based",
        interval_type="closed",
    )

    # Act
    result = canonical_end(raw, table)

    # Assert
    assert result == raw


@settings(max_examples=100, deadline=None)
@given(
    raw=st.text(
        alphabet=st.characters(min_codepoint=0x20, max_codepoint=0x7E),
        min_size=1,
        max_size=20,
    ),
    coord_system=st.sampled_from(["0based", "1based"]),
    interval_type=st.sampled_from(["half_open", "closed"]),
)
def test_canonical_end_should_return_one_of_three_documented_forms(
    raw, coord_system, interval_type
):
    """Test that canonical_end output matches one of three documented shapes.

    Given:
        A Hypothesis-generated raw column reference (ASCII-printable
        text of length 1..20) and a sampled
        ``(coordinate_system, interval_type)`` enum pair drawn from
        the four allowed combinations.
    When:
        ``canonical_end(raw, Table(...))`` is called.
    Then:
        It should return exactly one of ``raw``, ``f"({raw} + 1)"``,
        or ``f"({raw} - 1)"`` — no other shape is ever produced.
    """
    # Arrange
    table = Table(
        name="features",
        coordinate_system=coord_system,
        interval_type=interval_type,
    )
    allowed = {raw, f"({raw} + 1)", f"({raw} - 1)"}

    # Act
    result = canonical_end(raw, table)

    # Assert
    assert result in allowed


def test_decanonical_start_should_return_input_unchanged_when_table_is_none():
    """Test that a missing table leaves the canonical start untouched.

    Given:
        A canonical 0-based start expression and ``table=None``.
    When:
        ``decanonical_start(canon, None)`` is called.
    Then:
        It should return the canonical expression unchanged.
    """
    # Arrange
    canon = "a.start"

    # Act
    result = decanonical_start(canon, None)

    # Assert
    assert result == canon


def test_decanonical_start_should_return_input_unchanged_when_table_is_zero_based():
    """Test that a 0-based table leaves the canonical start untouched.

    Given:
        A :class:`giql.table.Table` with ``coordinate_system="0based"``.
    When:
        ``decanonical_start(canon, table)`` is called.
    Then:
        It should return ``canon`` unchanged because the canonical
        form already matches the source convention.
    """
    # Arrange
    canon = "a.start"
    table = Table(name="features", coordinate_system="0based")

    # Act
    result = decanonical_start(canon, table)

    # Assert
    assert result == canon


def test_decanonical_start_should_add_one_when_table_is_one_based():
    """Test that a 1-based table adds one to the canonical start.

    Given:
        A :class:`giql.table.Table` with ``coordinate_system="1based"``.
    When:
        ``decanonical_start(canon, table)`` is called.
    Then:
        It should return ``f"({canon} + 1)"`` so the value lands back
        in the table's 1-based encoding.
    """
    # Arrange
    canon = "(a.start - 1)"
    table = Table(name="features", coordinate_system="1based")

    # Act
    result = decanonical_start(canon, table)

    # Assert
    assert result == f"({canon} + 1)"


@settings(max_examples=100, deadline=None)
@given(
    x=st.integers(min_value=0, max_value=10**9),
    coord_system=st.sampled_from(["0based", "1based"]),
)
def test_decanonical_start_should_invert_canonical_start_numerically(x, coord_system):
    """Test that decanonical_start is the numeric inverse of canonical_start.

    Given:
        A Hypothesis-generated integer start coordinate and a sampled
        ``coordinate_system``.
    When:
        ``decanonical_start(canonical_start(str(x), table), table)`` is
        evaluated as a Python expression.
    Then:
        It should evaluate back to ``x`` — the two helpers are exact
        inverses across the full coordinate-system enum.
    """
    # Arrange
    table = Table(name="features", coordinate_system=coord_system)

    # Act
    composed = decanonical_start(canonical_start(str(x), table), table)

    # Assert
    assert eval(composed) == x


def test_decanonical_end_should_return_input_unchanged_when_table_is_none():
    """Test that a missing table leaves the canonical end untouched.

    Given:
        A canonical 0-based half-open end expression and ``table=None``.
    When:
        ``decanonical_end(canon, None)`` is called.
    Then:
        It should return the canonical expression unchanged.
    """
    # Arrange
    canon = "a.end"

    # Act
    result = decanonical_end(canon, None)

    # Assert
    assert result == canon


def test_decanonical_end_should_return_input_unchanged_when_table_is_zero_based_half_open():
    """Test that 0-based half-open tables leave the canonical end untouched.

    Given:
        A :class:`giql.table.Table` with ``coordinate_system="0based"``
        and ``interval_type="half_open"``.
    When:
        ``decanonical_end(canon, table)`` is called.
    Then:
        It should return ``canon`` unchanged because the canonical
        form already matches the source convention.
    """
    # Arrange
    canon = "a.end"
    table = Table(
        name="features",
        coordinate_system="0based",
        interval_type="half_open",
    )

    # Act
    result = decanonical_end(canon, table)

    # Assert
    assert result == canon


def test_decanonical_end_should_subtract_one_when_table_is_zero_based_closed():
    """Test that 0-based closed tables subtract one from the canonical end.

    Given:
        A :class:`giql.table.Table` with ``coordinate_system="0based"``
        and ``interval_type="closed"``.
    When:
        ``decanonical_end(canon, table)`` is called.
    Then:
        It should return ``f"({canon} - 1)"`` to convert canonical
        half-open back to the source's closed form.
    """
    # Arrange
    canon = "(a.end + 1)"
    table = Table(
        name="features",
        coordinate_system="0based",
        interval_type="closed",
    )

    # Act
    result = decanonical_end(canon, table)

    # Assert
    assert result == f"({canon} - 1)"


def test_decanonical_end_should_add_one_when_table_is_one_based_half_open():
    """Test that 1-based half-open tables add one to the canonical end.

    Given:
        A :class:`giql.table.Table` with ``coordinate_system="1based"``
        and ``interval_type="half_open"``.
    When:
        ``decanonical_end(canon, table)`` is called.
    Then:
        It should return ``f"({canon} + 1)"`` to convert canonical
        half-open back to the source's 1-based half-open form.
    """
    # Arrange
    canon = "(a.end - 1)"
    table = Table(
        name="features",
        coordinate_system="1based",
        interval_type="half_open",
    )

    # Act
    result = decanonical_end(canon, table)

    # Assert
    assert result == f"({canon} + 1)"


def test_decanonical_end_should_return_input_unchanged_when_table_is_one_based_closed():
    """Test that 1-based closed tables leave the canonical end untouched.

    Given:
        A :class:`giql.table.Table` with ``coordinate_system="1based"``
        and ``interval_type="closed"``.
    When:
        ``decanonical_end(canon, table)`` is called.
    Then:
        It should return ``canon`` unchanged because a 1-based closed
        end equals a canonical 0-based half-open end numerically.
    """
    # Arrange
    canon = "a.end"
    table = Table(
        name="features",
        coordinate_system="1based",
        interval_type="closed",
    )

    # Act
    result = decanonical_end(canon, table)

    # Assert
    assert result == canon


@settings(max_examples=100, deadline=None)
@given(
    x=st.integers(min_value=1, max_value=10**9),
    coord_system=st.sampled_from(["0based", "1based"]),
    interval_type=st.sampled_from(["half_open", "closed"]),
)
def test_decanonical_end_should_invert_canonical_end_numerically(
    x, coord_system, interval_type
):
    """Test that decanonical_end is the numeric inverse of canonical_end.

    Given:
        A Hypothesis-generated integer end coordinate and a sampled
        ``(coordinate_system, interval_type)`` pair drawn from the
        four allowed combinations.
    When:
        ``decanonical_end(canonical_end(str(x), table), table)`` is
        evaluated as a Python expression.
    Then:
        It should evaluate back to ``x`` — the two helpers are exact
        inverses across the full enum cross-product.
    """
    # Arrange
    table = Table(
        name="features",
        coordinate_system=coord_system,
        interval_type=interval_type,
    )

    # Act
    composed = decanonical_end(canonical_end(str(x), table), table)

    # Assert
    assert eval(composed) == x
