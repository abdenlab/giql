"""Encoding helpers for coordinate-space matrix integration tests.

A *canonical* interval in this module is always 0-based half-open. The
:func:`encode` helper translates a canonical ``(start, end)`` pair into the
integer values that should be stored under each of the four
``(coordinate_system, interval_type)`` combinations supported by GIQL.

The translation is the inverse of GIQL's coordinate canonicalization (see
:mod:`giql.canonical`): re-encoding a canonical interval and feeding the stored
values back through the canonicalization helpers must yield the original
canonical values.
"""

from giql.table import Table

CONVENTIONS: list[tuple[str, str]] = [
    ("0based", "half_open"),
    ("0based", "closed"),
    ("1based", "half_open"),
    ("1based", "closed"),
]


def encode(
    canonical_start: int,
    canonical_end: int,
    coordinate_system: str,
    interval_type: str,
) -> tuple[int, int]:
    """Translate a canonical 0-based half-open interval to a stored encoding.

    :param canonical_start:
        Canonical start (0-based half-open).
    :param canonical_end:
        Canonical end (0-based half-open).
    :param coordinate_system:
        Either ``"0based"`` or ``"1based"``.
    :param interval_type:
        Either ``"half_open"`` or ``"closed"``.
    :return:
        ``(stored_start, stored_end)`` for the given convention.
    :raises ValueError:
        If ``coordinate_system`` or ``interval_type`` is not a recognized
        value.
    """
    if coordinate_system == "0based":
        stored_start = canonical_start
    elif coordinate_system == "1based":
        stored_start = canonical_start + 1
    else:
        raise ValueError(f"Unknown coordinate_system: {coordinate_system!r}")

    key = (coordinate_system, interval_type)
    if key == ("0based", "half_open"):
        stored_end = canonical_end
    elif key == ("0based", "closed"):
        stored_end = canonical_end - 1
    elif key == ("1based", "half_open"):
        stored_end = canonical_end + 1
    elif key == ("1based", "closed"):
        stored_end = canonical_end
    else:
        raise ValueError(f"Unknown interval_type: {interval_type!r}")

    return stored_start, stored_end


def make_table(name: str, coordinate_system: str, interval_type: str) -> Table:
    """Build a :class:`giql.table.Table` for the given convention."""
    return Table(
        name,
        coordinate_system=coordinate_system,
        interval_type=interval_type,
    )
