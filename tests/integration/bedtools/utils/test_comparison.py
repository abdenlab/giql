"""Unit tests for result comparison logic."""

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from .comparison import compare_results

pytestmark = pytest.mark.integration


def test_compare_results_should_report_match_when_rows_identical():
    """Test that identical row lists compare as matching.

    Given:
        Two identical lists of tuples
    When:
        compare_results() is called
    Then:
        It should return match=True with no differences
    """
    # Arrange
    rows = [("chr1", 100, 200), ("chr1", 300, 400)]

    # Act
    result = compare_results(rows, rows)

    # Assert
    assert result.match is True
    assert result.differences == []


def test_compare_results_should_match_when_rows_in_different_order():
    """Test that row order does not affect match outcome.

    Given:
        Same tuples in different order
    When:
        compare_results() is called
    Then:
        It should return match=True
    """
    # Arrange
    a = [("chr1", 300, 400), ("chr1", 100, 200)]
    b = [("chr1", 100, 200), ("chr1", 300, 400)]

    # Act
    result = compare_results(a, b)

    # Assert
    assert result.match is True


def test_compare_results_should_report_mismatch_when_row_counts_differ():
    """Test that differing row counts produce a mismatch.

    Given:
        Lists with different row counts
    When:
        compare_results() is called
    Then:
        It should return match=False with a row count difference
    """
    # Arrange
    a = [("chr1", 100, 200)]
    b = [("chr1", 100, 200), ("chr1", 300, 400)]

    # Act
    result = compare_results(a, b)

    # Assert
    assert result.match is False
    assert any("Row count" in d for d in result.differences)


def test_compare_results_should_match_when_integer_values_identical():
    """Test that identical integer values compare as matching.

    Given:
        Rows with identical integer values
    When:
        compare_results() is called
    Then:
        It should return match=True
    """
    # Arrange
    a = [("chr1", 100, 200, 50)]
    b = [("chr1", 100, 200, 50)]

    # Act
    result = compare_results(a, b)

    # Assert
    assert result.match is True


def test_compare_results_should_match_when_floats_within_epsilon():
    """Test that floats within default epsilon compare as matching.

    Given:
        Rows with floats differing by less than epsilon
    When:
        compare_results() is called
    Then:
        It should return match=True
    """
    # Arrange
    a = [(1.0000000001,)]
    b = [(1.0,)]

    # Act
    result = compare_results(a, b)

    # Assert
    assert result.match is True


def test_compare_results_should_report_mismatch_when_floats_beyond_epsilon():
    """Test that floats beyond default epsilon produce a mismatch.

    Given:
        Rows with floats differing by more than epsilon
    When:
        compare_results() is called
    Then:
        It should return match=False
    """
    # Arrange
    a = [(1.5,)]
    b = [(1.0,)]

    # Act
    result = compare_results(a, b)

    # Assert
    assert result.match is False


def test_compare_results_should_match_when_custom_epsilon_tolerates_difference():
    """Test that a larger custom epsilon accommodates small float deltas.

    Given:
        Rows with floats differing by 0.05
    When:
        compare_results() is called with epsilon=0.1
    Then:
        It should return match=True
    """
    # Arrange
    a = [(1.05,)]
    b = [(1.0,)]

    # Act
    result = compare_results(a, b, epsilon=0.1)

    # Assert
    assert result.match is True


def test_compare_results_should_match_when_none_values_align():
    """Test that aligned None values compare as matching.

    Given:
        Rows with None in the same positions
    When:
        compare_results() is called
    Then:
        It should return match=True
    """
    # Arrange
    a = [("chr1", None, 200)]
    b = [("chr1", None, 200)]

    # Act
    result = compare_results(a, b)

    # Assert
    assert result.match is True


def test_compare_results_should_report_mismatch_when_none_vs_value():
    """Test that None paired with a concrete value produces a mismatch.

    Given:
        Rows where one has None and the other has a value
    When:
        compare_results() is called
    Then:
        It should return match=False
    """
    # Arrange
    a = [("chr1", None, 200)]
    b = [("chr1", 100, 200)]

    # Act
    result = compare_results(a, b)

    # Assert
    assert result.match is False


def test_compare_results_should_report_mismatch_when_column_counts_differ():
    """Test that differing column counts produce a mismatch.

    Given:
        Rows with different column counts
    When:
        compare_results() is called
    Then:
        It should return match=False with a column count difference
    """
    # Arrange
    a = [("chr1", 100, 200)]
    b = [("chr1", 100)]

    # Act
    result = compare_results(a, b)

    # Assert
    assert result.match is False
    assert any("Column count" in d for d in result.differences)


def test_compare_results_should_list_extra_rows_when_giql_has_more():
    """Test that extra GIQL rows are reported in differences.

    Given:
        GIQL has extra rows not in bedtools
    When:
        compare_results() is called
    Then:
        It should list the extra rows in differences
    """
    # Arrange
    a = [("chr1", 100, 200), ("chr1", 300, 400)]
    b = [("chr1", 100, 200)]

    # Act
    result = compare_results(a, b)

    # Assert
    assert result.match is False
    assert any(
        "missing in bedtools" in d.lower() or "Present in GIQL" in d
        for d in result.differences
    )


def test_compare_results_should_list_missing_rows_when_bedtools_has_more():
    """Test that extra bedtools rows are reported as missing in GIQL.

    Given:
        bedtools has extra rows not in GIQL
    When:
        compare_results() is called
    Then:
        It should list the missing rows in differences
    """
    # Arrange
    a = [("chr1", 100, 200)]
    b = [("chr1", 100, 200), ("chr1", 300, 400)]

    # Act
    result = compare_results(a, b)

    # Assert
    assert result.match is False
    assert any("Missing in GIQL" in d for d in result.differences)


def test_compare_results_should_match_when_both_empty():
    """Test that two empty lists compare as matching with zero counts.

    Given:
        Both lists empty
    When:
        compare_results() is called
    Then:
        It should return match=True with zero row counts
    """
    # Arrange
    # (no inputs to arrange beyond the empty lists passed below)

    # Act
    result = compare_results([], [])

    # Assert
    assert result.match is True
    assert result.giql_row_count == 0
    assert result.bedtools_row_count == 0


def test_compare_results_should_populate_metadata_with_epsilon_and_sorted():
    """Test that comparison metadata includes epsilon and sorted keys.

    Given:
        Any comparison
    When:
        compare_results() is called
    Then:
        It should populate comparison_metadata with epsilon and sorted keys
    """
    # Arrange
    # (no inputs to arrange beyond the empty lists passed below)

    # Act
    result = compare_results([], [])

    # Assert
    assert "epsilon" in result.comparison_metadata
    assert "sorted" in result.comparison_metadata


def test_compare_results_should_set_row_counts_when_sizes_differ():
    """Test that row counts are populated from the input list sizes.

    Given:
        Lists of different sizes
    When:
        compare_results() is called
    Then:
        It should set giql_row_count and bedtools_row_count correctly
    """
    # Arrange
    # (inputs are supplied inline in the Act step)

    # Act
    result = compare_results(
        [("a",), ("b",)],
        [("a",), ("b",), ("c",)],
    )

    # Assert
    assert result.giql_row_count == 2
    assert result.bedtools_row_count == 3


def test_compare_results_should_match_when_sorting_handles_none_values():
    """Test that sorting with None values completes without errors.

    Given:
        Rows containing None values in different positions
    When:
        compare_results() is called
    Then:
        It should handle None deterministically and return match=True
    """
    # Arrange
    a = [("chr1", None, 200), ("chr1", 100, 200)]
    b = [("chr1", 100, 200), ("chr1", None, 200)]

    # Act
    result = compare_results(a, b)

    # Assert
    assert result.match is True


@settings(max_examples=50)
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
def test_compare_results_should_always_match_when_comparing_rows_to_themselves(rows):
    """Test that self-comparison always yields a match.

    Given:
        Any list of tuples
    When:
        compare_results(rows, rows) is called
    Then:
        It should always return match=True
    """
    # Arrange
    # (rows supplied by Hypothesis)

    # Act
    result = compare_results(rows, rows)

    # Assert
    assert result.match is True
