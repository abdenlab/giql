"""Pybedtools wrapper for genomic interval operations.

This module provides functions for:
- Creating BedTool objects from interval data
- Executing bedtools operations via pybedtools
- Converting results to comparable formats
"""

from typing import List
from typing import Tuple

import pybedtools

# Strand flag constants for bedtools commands
STRAND_SAME = True  # Require same strand (pybedtools uses True for -s)
# Require opposite strands (pybedtools uses "opposite" for -S)
STRAND_OPPOSITE = "opposite"


class BedtoolsError(Exception):
    """Raised when bedtools operation fails."""

    pass


def create_bedtool(intervals: List[Tuple]) -> pybedtools.BedTool:
    """Create BedTool object from interval tuples.

    Args:
        intervals: List of tuples, each containing:
            - (chrom, start, end) for BED3 format
            - (chrom, start, end, name, score, strand) for BED6 format

    Returns:
        pybedtools.BedTool object

    Example:
        >>> intervals = [("chr1", 100, 200, "a1", 100, "+")]
        >>> bt = create_bedtool(intervals)
    """
    # Convert tuples to BED format strings
    bed_strings = []
    for interval in intervals:
        if len(interval) == 3:
            # BED3 format
            bed_strings.append(f"{interval[0]}\t{interval[1]}\t{interval[2]}")
        elif len(interval) >= 6:
            # BED6 format
            chrom, start, end, name, score, strand = interval[:6]
            # Handle None values
            name = name if name is not None else "."
            score = score if score is not None else 0
            strand = strand if strand is not None else "."
            bed_strings.append(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}")
        else:
            raise ValueError(f"Invalid interval format: {interval}")

    bed_string = "\n".join(bed_strings)
    return pybedtools.BedTool(bed_string, from_string=True)


def intersect(
    intervals_a: List[Tuple],
    intervals_b: List[Tuple],
    strand_mode: str | None = None,
) -> List[Tuple]:
    """Find overlapping intervals using bedtools intersect.

    Args:
        intervals_a: First set of intervals
        intervals_b: Second set of intervals
        strand_mode: Strand requirement ('same', 'opposite', or None for ignore)

    Returns:
        List of tuples matching intervals_a format

    Example:
        >>> a = [("chr1", 100, 200, "a1", 100, "+")]
        >>> b = [("chr1", 150, 250, "b1", 100, "+")]
        >>> result = intersect(a, b)
    """
    try:
        bt_a = create_bedtool(intervals_a)
        bt_b = create_bedtool(intervals_b)

        # Build kwargs for intersect
        # Use -u (unique) to return each interval from A only once
        # This matches GIQL's DISTINCT behavior
        kwargs = {"u": True}
        if strand_mode == "same":
            kwargs["s"] = True
        elif strand_mode == "opposite":
            kwargs["S"] = True

        # Perform intersection
        result = bt_a.intersect(bt_b, **kwargs)

        # Convert to tuples
        return bedtool_to_tuples(result)

    except Exception as e:
        raise BedtoolsError(f"Intersect operation failed: {e}")


def merge(intervals: List[Tuple], strand_mode: str | None = None) -> List[Tuple]:
    """Merge overlapping intervals using bedtools merge.

    Args:
        intervals: List of intervals to merge
        strand_mode: Strand requirement ('same' to merge per-strand, None for ignore)

    Returns:
        List of tuples in BED3 format (chrom, start, end)

    Example:
        >>> intervals = [
        ...     ("chr1", 100, 200, "a1", 100, "+"),
        ...     ("chr1", 180, 300, "a2", 100, "+"),
        ... ]
        >>> result = merge(intervals)
        >>> # Returns: [("chr1", 100, 300)]
    """
    try:
        bt = create_bedtool(intervals)

        # Sort before merging (required by bedtools merge)
        bt_sorted = bt.sort()

        # Build kwargs for merge
        kwargs = {}
        if strand_mode == "same":
            kwargs["s"] = True

        # Perform merge
        result = bt_sorted.merge(**kwargs)

        # Convert to tuples (merge returns BED3 format)
        return bedtool_to_tuples(result, format="bed3")

    except Exception as e:
        raise BedtoolsError(f"Merge operation failed: {e}")


def closest(
    intervals_a: List[Tuple],
    intervals_b: List[Tuple],
    strand_mode: str | None = None,
    k: int = 1,
) -> List[Tuple]:
    """Find closest intervals using bedtools closest.

    Args:
        intervals_a: Query intervals
        intervals_b: Database intervals to search
        strand_mode: Strand requirement ('same', 'opposite', or None for ignore)
        k: Number of closest intervals to report (default: 1)

    Returns:
        List of tuples with format: (a_fields..., b_fields..., distance)

    Example:
        >>> a = [("chr1", 100, 200, "a1", 100, "+")]
        >>> b = [("chr1", 300, 400, "b1", 100, "+")]
        >>> result = closest(a, b)
        >>> # Returns intervals from a and b with distance
    """
    try:
        bt_a = create_bedtool(intervals_a)
        bt_b = create_bedtool(intervals_b)

        # Sort inputs (required for -t flag)
        bt_a = bt_a.sort()
        bt_b = bt_b.sort()

        # Build kwargs for closest
        kwargs = {"d": True, "t": "first"}  # Report distance, break ties by taking first
        if k > 1:
            kwargs["k"] = k
        if strand_mode == "same":
            kwargs["s"] = True
        elif strand_mode == "opposite":
            kwargs["S"] = True

        # Perform closest
        result = bt_a.closest(bt_b, **kwargs)

        # Convert to tuples (closest returns concatenated fields + distance)
        return bedtool_to_tuples(result, format="closest")

    except Exception as e:
        raise BedtoolsError(f"Closest operation failed: {e}")


def bedtool_to_tuples(bedtool: pybedtools.BedTool, format: str = "bed6") -> List[Tuple]:
    """Convert BedTool object to list of tuples.

    Args:
        bedtool: pybedtools.BedTool object
        format: Expected format ('bed3', 'bed6', or 'closest')

    Returns:
        List of tuples matching the format

    Note:
        - bed3: (chrom, start, end)
        - bed6: (chrom, start, end, name, score, strand)
        - closest: (chrom_a, start_a, end_a, name_a, score_a, strand_a,
                    chrom_b, start_b, end_b, name_b, score_b, strand_b, distance)
    """
    rows = []

    for interval in bedtool:
        fields = interval.fields

        if format == "bed3":
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            rows.append((chrom, start, end))

        elif format == "bed6":
            if len(fields) < 6:
                # Pad with defaults if needed
                while len(fields) < 6:
                    if len(fields) == 3:
                        fields.append(".")  # name
                    elif len(fields) == 4:
                        fields.append("0")  # score
                    elif len(fields) == 5:
                        fields.append(".")  # strand

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3] if fields[3] != "." else None
            score = int(fields[4]) if fields[4] != "." else None
            strand = fields[5] if fields[5] != "." else None

            rows.append((chrom, start, end, name, score, strand))

        elif format == "closest":
            # Closest returns: a_fields + b_fields + distance
            # For BED6: 6 fields for a, 6 fields for b, 1 distance = 13 total
            if len(fields) >= 13:
                # Parse all fields as-is, converting appropriate ones to int
                row = []
                for i, field in enumerate(fields):
                    # Positions (1, 2, 7, 8) and distance (12) should be int
                    if i in (1, 2, 7, 8, 12):
                        row.append(int(field))
                    # Scores (4, 10) should be int if not "."
                    elif i in (4, 10):
                        row.append(int(field) if field != "." else None)
                    # Names (3, 9) and strands (5, 11) should be None if "."
                    elif i in (3, 5, 9, 11):
                        row.append(field if field != "." else None)
                    else:
                        row.append(field)
                rows.append(tuple(row))
            else:
                raise ValueError(
                    f"Unexpected number of fields for closest: {len(fields)}"
                )

        else:
            raise ValueError(f"Unsupported format: {format}")

    return rows


def add_strand_flag(kwargs: dict, strand_mode: str | None) -> dict:
    """Add strand flag to bedtools kwargs.

    Args:
        kwargs: Base kwargs dictionary
        strand_mode: Strand requirement ('same', 'opposite', or None for ignore)

    Returns:
        Updated kwargs dictionary with strand flag

    Example:
        >>> kwargs = add_strand_flag({}, "same")
        >>> # Returns: {"s": True}
    """
    updated_kwargs = kwargs.copy()

    if strand_mode == "same":
        updated_kwargs["s"] = True
    elif strand_mode == "opposite":
        updated_kwargs["S"] = True
    # None or other values = ignore strand (no flag added)

    return updated_kwargs
