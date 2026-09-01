"""Tests for the Table genomic-table configuration."""

import dataclasses

import pytest

from giql import Table


class TestTable:
    """Tests for Table."""

    def test___post_init___should_raise_when_coordinate_system_is_unknown(self):
        """Test that an unknown coordinate_system is rejected at construction.

        Given:
            A coordinate_system outside the documented ``0based`` / ``1based``
            pair.
        When:
            A Table is constructed with it.
        Then:
            It should raise ``ValueError`` naming the offending value.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError, match="coordinate_system"):
            Table("peaks", coordinate_system="2based")

    def test___post_init___should_raise_when_interval_type_is_unknown(self):
        """Test that an unknown interval_type is rejected at construction.

        Given:
            An interval_type outside the documented ``half_open`` / ``closed``
            pair.
        When:
            A Table is constructed with it.
        Then:
            It should raise ``ValueError`` naming the offending value.
        """
        # Arrange, act, & assert
        with pytest.raises(ValueError, match="interval_type"):
            Table("peaks", interval_type="inclusive")

    def test_table_should_raise_when_a_field_is_assigned_after_construction(self):
        """Test that a constructed Table cannot be mutated.

        Given:
            A validly constructed Table.
        When:
            One of its fields is assigned.
        Then:
            It should raise ``FrozenInstanceError``.

        The validation above runs only in ``__post_init__``, so a mutable Table
        would let a caller install a value the constructor rejects. That value is
        not inert: ``coordinate_system`` selects the coordinate translation the
        emitted SQL performs, so an unvalidated one produces a predicate matching
        no coordinate system at all, with no error raised anywhere.
        """
        # Arrange
        table = Table("peaks")

        # Act & assert
        with pytest.raises(dataclasses.FrozenInstanceError):
            table.coordinate_system = "2based"

    def test_table_should_build_a_variant_through_replace(self):
        """Test that dataclasses.replace re-validates the variant it builds.

        Given:
            A validly constructed Table.
        When:
            ``dataclasses.replace`` derives a variant from it.
        Then:
            A valid field value should produce a new Table, and an invalid one
            should raise -- ``replace`` re-runs ``__post_init__``, which is what
            makes it the supported way to vary a frozen config.
        """
        # Arrange
        table = Table("peaks")

        # Act
        variant = dataclasses.replace(table, coordinate_system="1based")

        # Assert
        assert variant.coordinate_system == "1based"
        assert table.coordinate_system == "0based"
        with pytest.raises(ValueError, match="coordinate_system"):
            dataclasses.replace(table, coordinate_system="2based")

    def test_table_should_be_hashable_and_value_equal(self):
        """Test that two Tables with identical fields are interchangeable.

        Given:
            Two separately constructed Tables carrying identical field values.
        When:
            They are compared and used as set members.
        Then:
            They should be equal and collapse to one set entry.

        A plain dataclass sets ``__hash__`` to None once it generates ``__eq__``,
        so Table was unhashable before it was frozen.
        """
        # Arrange
        first = Table("peaks", chrom_col="contig")
        second = Table("peaks", chrom_col="contig")

        # Act & assert
        assert first == second
        assert len({first, second}) == 1
