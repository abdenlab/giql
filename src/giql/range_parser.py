"""Parse genomic range strings into structured data.

Supported formats:
    - Simple: 'chr1:1000-2000'
    - Explicit half-open: 'chr1:[1000,2000)'
    - Explicit closed: 'chr1:[1001,2000]'
    - With strand: 'chr1:1000-2000:+'
    - Points: 'chr1:1500'
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum
from typing import Literal
from typing import Optional


class CoordinateSystem(Enum):
    """Coordinate system for genomic ranges."""

    ZERO_BASED = "0based"
    ONE_BASED = "1based"


class IntervalType(Enum):
    """Interval endpoint handling."""

    HALF_OPEN = "half_open"  # [start, end)
    CLOSED = "closed"  # [start, end]


@dataclass
class ParsedRange:
    """Structured representation of a genomic range."""

    chromosome: str
    start: int
    end: int
    interval_type: IntervalType
    strand: Optional[Literal["+", "-"]] = None

    def to_zero_based_half_open(self) -> ParsedRange:
        """Convert to canonical 0-based half-open representation.

        Conversions:
            - Closed [1000, 1999] -> Half-open [1000, 2000)

        :return: ParsedRange in 0-based half-open format
        """
        if self.interval_type == IntervalType.HALF_OPEN:
            return self

        # Closed to half-open: make end exclusive
        return ParsedRange(
            chromosome=self.chromosome,
            start=self.start,
            end=self.end + 1,
            interval_type=IntervalType.HALF_OPEN,
            strand=self.strand,
        )

    def length(self) -> int:
        """Calculate range length.

        :return: Length of the genomic range in base pairs
        """
        if self.interval_type == IntervalType.HALF_OPEN:
            return self.end - self.start
        else:
            return self.end - self.start + 1


class RangeParser:
    """Parse genomic range strings."""

    # chr1:1000-2000 or chr1:1000-2000:+
    SIMPLE_PATTERN = re.compile(
        r"^(?P<chr>[\w.]+):(?P<start>\d+)-(?P<end>\d+)(?::(?P<strand>[+-]))?$"
    )

    # chr1:[1000,2000) or chr1:[1000,2000]:+
    EXPLICIT_PATTERN = re.compile(
        r"^(?P<chr>[\w.]+):\[(?P<start>\d+),(?P<end>\d+)(?P<bracket>[\)\]])(?::(?P<strand>[+-]))?$"
    )

    # chr1:1500
    POINT_PATTERN = re.compile(r"^(?P<chr>[\w.]+):(?P<pos>\d+)$")

    @classmethod
    def parse(cls, range_str: str) -> ParsedRange:
        """Parse a genomic range string.

        :param range_str: String like 'chr1:1000-2000'
        :return: ParsedRange object
        :raises ValueError: If the string cannot be parsed
        """
        range_str = range_str.strip().strip("'\"")

        # Try point format
        match = cls.POINT_PATTERN.match(range_str)
        if match:
            return cls._parse_point(match)

        # Try explicit format
        match = cls.EXPLICIT_PATTERN.match(range_str)
        if match:
            return cls._parse_explicit(match)

        # Try simple format
        match = cls.SIMPLE_PATTERN.match(range_str)
        if match:
            return cls._parse_simple(match)

        raise ValueError(f"Invalid genomic range format: {range_str}")

    @classmethod
    def _parse_point(cls, match) -> ParsedRange:
        """Parse point format: chr1:1500 -> [1500, 1501).

        :param match: Regex match object
        :return: ParsedRange representing a single base position
        """
        chromosome = match.group("chr")
        position = int(match.group("pos"))

        return ParsedRange(
            chromosome=chromosome,
            start=position,
            end=position + 1,
            interval_type=IntervalType.HALF_OPEN,
            strand=None,
        )

    @classmethod
    def _parse_explicit(cls, match) -> ParsedRange:
        """Parse explicit format: chr1:[1000,2000).

        :param match: Regex match object
        :return: ParsedRange with explicit interval type
        :raises ValueError: If start >= end
        """
        chromosome = match.group("chr")
        start = int(match.group("start"))
        end = int(match.group("end"))
        bracket = match.group("bracket")
        strand = match.group("strand")

        if start >= end:
            raise ValueError(f"Start must be less than end: {start} >= {end}")

        interval_type = IntervalType.HALF_OPEN if bracket == ")" else IntervalType.CLOSED

        return ParsedRange(
            chromosome=chromosome,
            start=start,
            end=end,
            interval_type=interval_type,
            strand=strand,
        )

    @classmethod
    def _parse_simple(cls, match) -> ParsedRange:
        """Parse simple format: chr1:1000-2000.

        :param match: Regex match object
        :return: ParsedRange in half-open format
        :raises ValueError: If start >= end
        """
        chromosome = match.group("chr")
        start = int(match.group("start"))
        end = int(match.group("end"))
        strand = match.group("strand")

        if start >= end:
            raise ValueError(f"Start must be less than end: {start} >= {end}")

        return ParsedRange(
            chromosome=chromosome,
            start=start,
            end=end,
            interval_type=IntervalType.HALF_OPEN,
            strand=strand,
        )
