"""Convert raw column expressions to and from canonical 0-based half-open coordinates.

Genomic intervals can be expressed in 0-based half-open, 0-based closed,
1-based half-open, or 1-based closed coordinates. The transpiler emits SQL
in 0-based half-open by canonicalizing every column endpoint at the source.
These helpers wrap a raw column reference with the offset needed for that
canonicalization based on the source table's declared coordinate system.
The ``decanonical_*`` helpers apply the inverse offset to emit an endpoint
back in the table's own coordinate system. (Literal ranges go through
:meth:`giql.range_parser.RangeParser.to_zero_based_half_open`.)
"""

from giql.table import Table


def canonical_start(raw_start: str, table: Table | None) -> str:
    """Wrap a raw start column expression to canonical 0-based half-open start.

    When ``table`` is ``None``, the raw expression is returned unchanged
    (treated as 0-based half-open).

    Conversion by :attr:`Table.coordinate_system`:

    - ``0based``: ``start`` (identity)
    - ``1based``: ``start - 1``
    """
    if table is None or table.coordinate_system == "0based":
        return raw_start
    return f"({raw_start} - 1)"


def canonical_end(raw_end: str, table: Table | None) -> str:
    """Wrap a raw end column expression to canonical 0-based half-open end.

    When ``table`` is ``None``, the raw expression is returned unchanged
    (treated as 0-based half-open).

    Conversion by (:attr:`Table.coordinate_system`, :attr:`Table.interval_type`):

    - ``0based`` / ``half_open``: ``end`` (identity)
    - ``0based`` / ``closed``:    ``end + 1``
    - ``1based`` / ``half_open``: ``end - 1``
    - ``1based`` / ``closed``:    ``end`` (identity)
    """
    if table is None:
        return raw_end
    key = (table.coordinate_system, table.interval_type)
    if key == ("0based", "closed"):
        return f"({raw_end} + 1)"
    if key == ("1based", "half_open"):
        return f"({raw_end} - 1)"
    return raw_end


def decanonical_start(canonical_expr: str, table: Table | None) -> str:
    """Convert a canonical 0-based half-open start to the table's encoding.

    Inverse of :func:`canonical_start`. When ``table`` is ``None``, the
    expression is returned unchanged (treated as 0-based half-open).

    Conversion by :attr:`Table.coordinate_system`:

    - ``0based``: ``start`` (identity)
    - ``1based``: ``start + 1``
    """
    if table is None or table.coordinate_system == "0based":
        return canonical_expr
    return f"({canonical_expr} + 1)"


def decanonical_end(canonical_expr: str, table: Table | None) -> str:
    """Convert a canonical 0-based half-open end to the table's encoding.

    Inverse of :func:`canonical_end`. When ``table`` is ``None``, the
    expression is returned unchanged (treated as 0-based half-open).

    Conversion by (:attr:`Table.coordinate_system`, :attr:`Table.interval_type`):

    - ``0based`` / ``half_open``: ``end`` (identity)
    - ``0based`` / ``closed``:    ``end - 1``
    - ``1based`` / ``half_open``: ``end + 1``
    - ``1based`` / ``closed``:    ``end`` (identity)
    """
    if table is None:
        return canonical_expr
    key = (table.coordinate_system, table.interval_type)
    if key == ("0based", "closed"):
        return f"({canonical_expr} - 1)"
    if key == ("1based", "half_open"):
        return f"({canonical_expr} + 1)"
    return canonical_expr
