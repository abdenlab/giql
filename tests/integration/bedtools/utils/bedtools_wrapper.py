"""Pybedtools wrapper for genomic interval operations."""

import pybedtools


class BedtoolsError(Exception):
    """Raised when bedtools operation fails."""


def create_bedtool(intervals: list[tuple]) -> pybedtools.BedTool:
    """Create BedTool object from interval tuples."""
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
    intervals_a: list[tuple],
    intervals_b: list[tuple],
    strand_mode: str | None = None,
    *,
    loj: bool = False,
) -> list[tuple]:
    """Find overlapping intervals using bedtools intersect.

    When *loj* is True, use left outer join mode (-loj): every A
    interval appears in the output, paired with overlapping B
    intervals or with null-placeholder fields when there is no match.
    """
    try:
        bt_a = create_bedtool(intervals_a)
        bt_b = create_bedtool(intervals_b)

        if loj:
            kwargs = {"loj": True}
        else:
            kwargs = {"u": True}
        if strand_mode == "same":
            kwargs["s"] = True
        elif strand_mode == "opposite":
            kwargs["S"] = True

        result = bt_a.intersect(bt_b, **kwargs)

        if loj:
            return bedtool_to_tuples(result, bed_format="loj")
        return bedtool_to_tuples(result)

    except Exception as e:
        raise BedtoolsError(f"Intersect operation failed: {e}")


def merge(intervals: list[tuple], strand_mode: str | None = None) -> list[tuple]:
    """Merge overlapping intervals using bedtools merge."""
    try:
        bt = create_bedtool(intervals)
        bt_sorted = bt.sort()

        kwargs = {}
        if strand_mode == "same":
            kwargs["s"] = True

        result = bt_sorted.merge(**kwargs)
        return bedtool_to_tuples(result, bed_format="bed3")

    except Exception as e:
        raise BedtoolsError(f"Merge operation failed: {e}")


def closest(
    intervals_a: list[tuple],
    intervals_b: list[tuple],
    strand_mode: str | None = None,
    k: int = 1,
) -> list[tuple]:
    """Find closest intervals using bedtools closest."""
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
        return bedtool_to_tuples(result, bed_format="closest")

    except Exception as e:
        raise BedtoolsError(f"Closest operation failed: {e}")


def bedtool_to_tuples(
    bedtool: pybedtools.BedTool, bed_format: str = "bed6"
) -> list[tuple]:
    """Convert BedTool object to list of tuples.

    Args:
        bedtool: pybedtools.BedTool object
        bed_format: Expected format ('bed3', 'bed6', 'loj', or 'closest')

    LOJ format assumes BED6(A)+BED6(B) (12 fields):
        Fields 0-5: A interval
        Fields 6-11: B interval (all '.' / -1 when unmatched)

    Closest format assumes BED6+BED6+distance (13 fields):
        Fields 0-5: A interval (chrom, start, end, name, score, strand)
        Fields 6-11: B interval (chrom, start, end, name, score, strand)
        Field 12: distance (integer)
    """
    rows = []

    for interval in bedtool:
        fields = interval.fields

        if bed_format == "bed3":
            rows.append((fields[0], int(fields[1]), int(fields[2])))

        elif bed_format == "bed6":
            while len(fields) < 6:
                if len(fields) == 3:
                    fields.append(".")
                elif len(fields) == 4:
                    fields.append("0")
                elif len(fields) == 5:
                    fields.append(".")

            rows.append(
                (
                    fields[0],
                    int(fields[1]),
                    int(fields[2]),
                    fields[3] if fields[3] != "." else None,
                    int(fields[4]) if fields[4] != "." else None,
                    fields[5] if fields[5] != "." else None,
                )
            )

        elif bed_format == "loj":
            if len(fields) < 12:
                raise ValueError(f"Unexpected number of fields for loj: {len(fields)}")

            def _loj_field(val, as_int=False):
                if val == "." or val == "-1":
                    return None
                return int(val) if as_int else val

            rows.append(
                (
                    fields[0],
                    int(fields[1]),
                    int(fields[2]),
                    fields[3] if fields[3] != "." else None,
                    int(fields[4]) if fields[4] != "." else None,
                    fields[5] if fields[5] != "." else None,
                    _loj_field(fields[6]),
                    _loj_field(fields[7], as_int=True),
                    _loj_field(fields[8], as_int=True),
                    _loj_field(fields[9]),
                    _loj_field(fields[10], as_int=True),
                    _loj_field(fields[11]),
                )
            )

        elif bed_format == "closest":
            if len(fields) < 13:
                raise ValueError(
                    f"Unexpected number of fields for closest: {len(fields)}"
                )
            row = []
            for i, field_val in enumerate(fields):
                if i in (1, 2, 7, 8, 12):
                    row.append(int(field_val))
                elif i in (4, 10):
                    row.append(int(field_val) if field_val != "." else None)
                elif i in (3, 5, 9, 11):
                    row.append(field_val if field_val != "." else None)
                else:
                    row.append(field_val)
            rows.append(tuple(row))

        else:
            raise ValueError(f"Unsupported format: {bed_format}")

    return rows
