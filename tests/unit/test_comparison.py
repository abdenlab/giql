"""Unit tests for result comparison logic."""

from hypothesis import given
from hypothesis import strategies as st

from tests.integration.bedtools.utils.comparison import compare_results


class TestCompareResults:
    def test_exact_match(self):
        """
        GIVEN two identical lists of tuples
        WHEN compare_results() is called
        THEN returns match=True with no differences
        """
        rows = [("chr1", 100, 200), ("chr1", 300, 400)]
        result = compare_results(rows, rows)
        assert result.match is True
        assert result.differences == []

    def test_order_independent(self):
        """
        GIVEN same tuples in different order
        WHEN compare_results() is called
        THEN returns match=True
        """
        a = [("chr1", 300, 400), ("chr1", 100, 200)]
        b = [("chr1", 100, 200), ("chr1", 300, 400)]
        result = compare_results(a, b)
        assert result.match is True

    def test_row_count_mismatch(self):
        """
        GIVEN lists with different row counts
        WHEN compare_results() is called
        THEN returns match=False with row count difference
        """
        a = [("chr1", 100, 200)]
        b = [("chr1", 100, 200), ("chr1", 300, 400)]
        result = compare_results(a, b)
        assert result.match is False
        assert any("Row count" in d for d in result.differences)

    def test_integer_exact_match(self):
        """
        GIVEN rows with identical integer values
        WHEN compare_results() is called
        THEN returns match=True
        """
        a = [("chr1", 100, 200, 50)]
        b = [("chr1", 100, 200, 50)]
        result = compare_results(a, b)
        assert result.match is True

    def test_float_within_epsilon(self):
        """
        GIVEN rows with floats differing by less than epsilon
        WHEN compare_results() is called
        THEN returns match=True
        """
        a = [(1.0000000001,)]
        b = [(1.0,)]
        result = compare_results(a, b)
        assert result.match is True

    def test_float_beyond_epsilon(self):
        """
        GIVEN rows with floats differing by more than epsilon
        WHEN compare_results() is called
        THEN returns match=False
        """
        a = [(1.5,)]
        b = [(1.0,)]
        result = compare_results(a, b)
        assert result.match is False

    def test_custom_epsilon(self):
        """
        GIVEN rows with floats differing by 0.05
        WHEN compare_results() is called with epsilon=0.1
        THEN returns match=True
        """
        a = [(1.05,)]
        b = [(1.0,)]
        result = compare_results(a, b, epsilon=0.1)
        assert result.match is True

    def test_none_none_match(self):
        """
        GIVEN rows with None in the same positions
        WHEN compare_results() is called
        THEN returns match=True
        """
        a = [("chr1", None, 200)]
        b = [("chr1", None, 200)]
        result = compare_results(a, b)
        assert result.match is True

    def test_none_vs_value_mismatch(self):
        """
        GIVEN rows where one has None and other has a value
        WHEN compare_results() is called
        THEN returns match=False
        """
        a = [("chr1", None, 200)]
        b = [("chr1", 100, 200)]
        result = compare_results(a, b)
        assert result.match is False

    def test_column_count_mismatch(self):
        """
        GIVEN rows with different column counts
        WHEN compare_results() is called
        THEN returns match=False with column count difference
        """
        a = [("chr1", 100, 200)]
        b = [("chr1", 100)]
        result = compare_results(a, b)
        assert result.match is False
        assert any("Column count" in d for d in result.differences)

    def test_extra_giql_rows(self):
        """
        GIVEN GIQL has extra rows not in bedtools
        WHEN compare_results() is called
        THEN differences list the extra rows
        """
        a = [("chr1", 100, 200), ("chr1", 300, 400)]
        b = [("chr1", 100, 200)]
        result = compare_results(a, b)
        assert result.match is False
        assert any(
            "missing in bedtools" in d.lower() or "Present in GIQL" in d
            for d in result.differences
        )

    def test_extra_bedtools_rows(self):
        """
        GIVEN bedtools has extra rows not in GIQL
        WHEN compare_results() is called
        THEN differences list the missing rows
        """
        a = [("chr1", 100, 200)]
        b = [("chr1", 100, 200), ("chr1", 300, 400)]
        result = compare_results(a, b)
        assert result.match is False
        assert any("Missing in GIQL" in d for d in result.differences)

    def test_empty_comparison(self):
        """
        GIVEN both lists empty
        WHEN compare_results() is called
        THEN returns match=True with zero row counts
        """
        result = compare_results([], [])
        assert result.match is True
        assert result.giql_row_count == 0
        assert result.bedtools_row_count == 0

    def test_metadata_populated(self):
        """
        GIVEN any comparison
        WHEN compare_results() is called
        THEN comparison_metadata contains epsilon and sorted keys
        """
        result = compare_results([], [])
        assert "epsilon" in result.comparison_metadata
        assert "sorted" in result.comparison_metadata

    def test_row_counts_set(self):
        """
        GIVEN lists of different sizes
        WHEN compare_results() is called
        THEN giql_row_count and bedtools_row_count are set correctly
        """
        result = compare_results(
            [("a",), ("b",)],
            [("a",), ("b",), ("c",)],
        )
        assert result.giql_row_count == 2
        assert result.bedtools_row_count == 3

    def test_sorting_with_none_values(self):
        """
        GIVEN rows containing None values in different positions
        WHEN compare_results() is called
        THEN sorting handles None deterministically without errors
        """
        a = [("chr1", None, 200), ("chr1", 100, 200)]
        b = [("chr1", 100, 200), ("chr1", None, 200)]
        result = compare_results(a, b)
        assert result.match is True

    @given(
        rows=st.lists(
            st.tuples(
                st.sampled_from(["chr1", "chr2"]),
                st.integers(min_value=0, max_value=10000),
                st.integers(min_value=0, max_value=10000),
            ),
            min_size=0,
            max_size=20,
        )
    )
    def test_self_comparison_always_matches(self, rows):
        """
        GIVEN any list of tuples
        WHEN compare_results(rows, rows) is called
        THEN always returns match=True
        """
        result = compare_results(rows, rows)
        assert result.match is True
