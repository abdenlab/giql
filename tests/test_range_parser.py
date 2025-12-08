import pytest

from giql.range_parser import IntervalType
from giql.range_parser import ParsedRange
from giql.range_parser import RangeParser


class TestRangeParser:
    def test_parse_simple_range(self):
        """
        GIVEN a simple range string
        WHEN parsing the range
        THEN should return a ParsedRange with correct values
        """
        result = RangeParser.parse("chr1:1000-2000")
        assert result.chromosome == "chr1"
        assert result.start == 1000
        assert result.end == 2000
        assert result.interval_type == IntervalType.HALF_OPEN
        assert result.strand is None

    def test_parse_explicit_half_open(self):
        """
        GIVEN an explicit half-open range string
        WHEN parsing the range
        THEN should return a ParsedRange with HALF_OPEN interval type
        """
        result = RangeParser.parse("chr1:[1000,2000)")
        assert result.interval_type == IntervalType.HALF_OPEN
        assert result.end == 2000

    def test_parse_explicit_closed(self):
        """
        GIVEN an explicit closed range string
        WHEN parsing the range
        THEN should return a ParsedRange with CLOSED interval type
        """
        result = RangeParser.parse("chr1:[1001,2000]")
        assert result.interval_type == IntervalType.CLOSED
        assert result.end == 2000

    def test_parse_with_strand(self):
        """
        GIVEN range strings with strand information
        WHEN parsing the ranges
        THEN should correctly parse the strand
        """
        result = RangeParser.parse("chr1:1000-2000:+")
        assert result.strand == "+"

        result = RangeParser.parse("chr1:1000-2000:-")
        assert result.strand == "-"

    def test_parse_point(self):
        """
        GIVEN a point range string
        WHEN parsing the range
        THEN should return a ParsedRange representing a single position
        """
        result = RangeParser.parse("chr1:1500")
        assert result.start == 1500
        assert result.end == 1501
        assert result.interval_type == IntervalType.HALF_OPEN

    def test_to_zero_based_half_open(self):
        """
        GIVEN a closed interval ParsedRange
        WHEN converting to zero-based half-open
        THEN should increment the end position
        """
        closed = ParsedRange("chr1", 1001, 2000, IntervalType.CLOSED)
        converted = closed.to_zero_based_half_open()
        assert converted.end == 2001
        assert converted.interval_type == IntervalType.HALF_OPEN

    def test_range_length(self):
        """
        GIVEN ParsedRange objects with different interval types
        WHEN calculating length
        THEN should return correct length for each type
        """
        half_open = ParsedRange("chr1", 1000, 2000, IntervalType.HALF_OPEN)
        assert half_open.length() == 1000

        closed = ParsedRange("chr1", 1000, 2000, IntervalType.CLOSED)
        assert closed.length() == 1001

    def test_invalid_range(self):
        """
        GIVEN invalid range strings
        WHEN parsing the ranges
        THEN should raise ValueError
        """
        with pytest.raises(ValueError):
            RangeParser.parse("invalid")

        with pytest.raises(ValueError):
            RangeParser.parse("chr1:2000-1000")

    def test_chromosome_formats(self):
        """
        GIVEN various chromosome naming conventions
        WHEN parsing the ranges
        THEN should correctly parse all formats
        """
        assert RangeParser.parse("chr1:100-200").chromosome == "chr1"
        assert RangeParser.parse("1:100-200").chromosome == "1"
        assert RangeParser.parse("chrX:100-200").chromosome == "chrX"
        assert RangeParser.parse("chrM:100-200").chromosome == "chrM"
