"""Result comparison logic for GIQL vs bedtools outputs.

This module provides functions for:
- Comparing GIQL and bedtools results with appropriate tolerance
- Order-independent row sorting
- Epsilon-based float comparison
"""

from typing import Any
from typing import List
from typing import Tuple

from .data_models import ComparisonResult


def _sort_key(row: Tuple) -> Tuple:
    """Generate sort key for order-independent comparison."""
    return tuple("" if v is None else v for v in row)


def _values_match(val1: Any, val2: Any, epsilon: float = 1e-9) -> bool:
    """Compare two values with appropriate tolerance."""
    if val1 is None and val2 is None:
        return True
    if val1 is None or val2 is None:
        return False

    if isinstance(val1, float) or isinstance(val2, float):
        try:
            return abs(float(val1) - float(val2)) <= epsilon
        except (ValueError, TypeError):
            return False

    return val1 == val2


def compare_results(
    giql_rows: List[Tuple],
    bedtools_rows: List[Tuple],
    epsilon: float = 1e-9,
) -> ComparisonResult:
    """Compare GIQL and bedtools results with appropriate tolerance.

    Comparison rules:
    - Integer positions/counts: exact match required
    - Floating-point values: epsilon tolerance
    - Row ordering: order-independent (sorts both result sets)
    """
    giql_count = len(giql_rows)
    bedtools_count = len(bedtools_rows)

    giql_sorted = sorted(giql_rows, key=_sort_key)
    bedtools_sorted = sorted(bedtools_rows, key=_sort_key)

    differences = []

    if giql_count != bedtools_count:
        differences.append(
            f"Row count mismatch: GIQL has {giql_count} rows, "
            f"bedtools has {bedtools_count} rows"
        )

    max_rows = max(giql_count, bedtools_count)
    for i in range(max_rows):
        if i >= giql_count:
            differences.append(
                f"Row {i}: Missing in GIQL, present in bedtools: {bedtools_sorted[i]}"
            )
            continue
        if i >= bedtools_count:
            differences.append(
                f"Row {i}: Present in GIQL, missing in bedtools: {giql_sorted[i]}"
            )
            continue

        giql_row = giql_sorted[i]
        bedtools_row = bedtools_sorted[i]

        if len(giql_row) != len(bedtools_row):
            differences.append(
                f"Row {i}: Column count mismatch "
                f"(GIQL: {len(giql_row)} cols, "
                f"bedtools: {len(bedtools_row)} cols)"
            )
            continue

        for col_idx, (giql_val, bedtools_val) in enumerate(zip(giql_row, bedtools_row)):
            if not _values_match(giql_val, bedtools_val, epsilon):
                differences.append(
                    f"Row {i}, col {col_idx}: "
                    f"GIQL={giql_val!r} != bedtools={bedtools_val!r}"
                )

    match = len(differences) == 0

    return ComparisonResult(
        match=match,
        giql_row_count=giql_count,
        bedtools_row_count=bedtools_count,
        differences=differences,
        comparison_metadata={"epsilon": epsilon, "sorted": True},
    )
