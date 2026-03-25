"""Pybedtools wrapper for genomic interval operations.

This module provides functions for:
- Creating BedTool objects from interval data
- Executing bedtools operations via pybedtools
- Converting results to comparable formats
"""

from typing import List
from typing import Tuple

import pybedtools


class BedtoolsError(Exception):
    """Raised when bedtools operation fails."""

    pass


def create_bedtool(intervals: List[Tuple]) -> pybedtools.BedTool:
    """Create BedTool object from interval tuples.

    Args:
        intervals: List of tuples, each containing:
            - (chrom, start, end) for BED3 format
            - (chrom, start, end, name, score, strand) for BED6 format
    """
    bed_strings = []
    for interval in intervals:
        if len(interval) == 3:
            bed_strings.append(f"{interval[0]}\t{interval[1]}\t{interval[2]}")
        elif len(interval) >= 6:
            chrom, start, end, name, score, strand = interval[:6]
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
        strand_mode: Strand requirement ('same', 'opposite', or None)
    """
    try:
        bt_a = create_bedtool(intervals_a)
        bt_b = create_bedtool(intervals_b)

        kwargs = {"u": True}
        if strand_mode == "same":
            kwargs["s"] = True
        elif strand_mode == "opposite":
            kwargs["S"] = True

        result = bt_a.intersect(bt_b, **kwargs)
        return bedtool_to_tuples(result)

    except Exception as e:
        raise BedtoolsError(f"Intersect operation failed: {e}")


def merge(intervals: List[Tuple], strand_mode: str | None = None) -> List[Tuple]:
    """Merge overlapping intervals using bedtools merge.

    Args:
        intervals: List of intervals to merge
        strand_mode: 'same' to merge per-strand, None to ignore
    """
    try:
        bt = create_bedtool(intervals)
        bt_sorted = bt.sort()

        kwargs = {}
        if strand_mode == "same":
            kwargs["s"] = True

        result = bt_sorted.merge(**kwargs)
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
        strand_mode: Strand requirement ('same', 'opposite', or None)
        k: Number of closest intervals to report
    """
    try:
        bt_a = create_bedtool(intervals_a)
        bt_b = create_bedtool(intervals_b)

        bt_a = bt_a.sort()
        bt_b = bt_b.sort()

        kwargs = {"d": True, "t": "first"}
        if k > 1:
            kwargs["k"] = k
        if strand_mode == "same":
            kwargs["s"] = True
        elif strand_mode == "opposite":
            kwargs["S"] = True

        result = bt_a.closest(bt_b, **kwargs)
        return bedtool_to_tuples(result, format="closest")

    except Exception as e:
        raise BedtoolsError(f"Closest operation failed: {e}")


def bedtool_to_tuples(bedtool: pybedtools.BedTool, format: str = "bed6") -> List[Tuple]:
    """Convert BedTool object to list of tuples.

    Args:
        bedtool: pybedtools.BedTool object
        format: Expected format ('bed3', 'bed6', or 'closest')
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
            while len(fields) < 6:
                if len(fields) == 3:
                    fields.append(".")
                elif len(fields) == 4:
                    fields.append("0")
                elif len(fields) == 5:
                    fields.append(".")

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3] if fields[3] != "." else None
            score = int(fields[4]) if fields[4] != "." else None
            strand = fields[5] if fields[5] != "." else None

            rows.append((chrom, start, end, name, score, strand))

        elif format == "closest":
            if len(fields) >= 13:
                row = []
                for i, field in enumerate(fields):
                    if i in (1, 2, 7, 8, 12):
                        row.append(int(field))
                    elif i in (4, 10):
                        row.append(int(field) if field != "." else None)
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
