"""Tests for giql.table module."""

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from giql.table import Table
from giql.table import Tables


class TestTable:
    """Tests for the Table dataclass."""

    def test___init___should_use_default_values_when_only_name_provided(self):
        """Test Table uses default values when only `name` is provided.

        Given:
            Only the required arg `name`
        When:
            Table is instantiated
        Then:
            It should set all fields to their default values
        """
        # Arrange / Act
        table = Table(name="peaks")

        # Assert
        assert table.name == "peaks"
        assert table.genomic_col == "interval"
        assert table.chrom_col == "chrom"
        assert table.start_col == "start"
        assert table.end_col == "end"
        assert table.strand_col == "strand"
        assert table.coordinate_system == "0based"
        assert table.interval_type == "half_open"

    def test___init___should_reflect_custom_values_when_all_fields_provided(self):
        """Test Table reflects custom values when all fields are provided.

        Given:
            All fields provided with custom values
        When:
            Table is instantiated
        Then:
            It should populate all fields with the custom values
        """
        # Arrange / Act
        table = Table(
            name="variants",
            genomic_col="position",
            chrom_col="chr",
            start_col="pos_start",
            end_col="pos_end",
            strand_col="direction",
            coordinate_system="1based",
            interval_type="closed",
        )

        # Assert
        assert table.name == "variants"
        assert table.genomic_col == "position"
        assert table.chrom_col == "chr"
        assert table.start_col == "pos_start"
        assert table.end_col == "pos_end"
        assert table.strand_col == "direction"
        assert table.coordinate_system == "1based"
        assert table.interval_type == "closed"

    def test___init___should_allow_none_when_strand_col_is_none(self):
        """Test Table allows strand_col to be None.

        Given:
            strand_col=None
        When:
            Table is instantiated
        Then:
            It should set strand_col to None
        """
        # Arrange / Act
        table = Table(name="peaks", strand_col=None)

        # Assert
        assert table.strand_col is None

    def test___init___should_accept_1based_when_coordinate_system_is_1based(self):
        """Test Table accepts the 1based coordinate system.

        Given:
            coordinate_system="1based"
        When:
            Table is instantiated
        Then:
            It should set coordinate_system to "1based"
        """
        # Arrange / Act
        table = Table(name="peaks", coordinate_system="1based")

        # Assert
        assert table.coordinate_system == "1based"

    def test___init___should_accept_closed_when_interval_type_is_closed(self):
        """Test Table accepts the closed interval type.

        Given:
            interval_type="closed"
        When:
            Table is instantiated
        Then:
            It should set interval_type to "closed"
        """
        # Arrange / Act
        table = Table(name="peaks", interval_type="closed")

        # Assert
        assert table.interval_type == "closed"

    def test___init___should_raise_when_coordinate_system_invalid(self):
        """Test Table raises when coordinate_system is invalid.

        Given:
            coordinate_system="invalid"
        When:
            Table is instantiated
        Then:
            It should raise ValueError mentioning coordinate_system
        """
        # Arrange / Act / Assert
        with pytest.raises(ValueError, match="coordinate_system"):
            Table(name="peaks", coordinate_system="invalid")

    def test___init___should_raise_when_interval_type_invalid(self):
        """Test Table raises when interval_type is invalid.

        Given:
            interval_type="invalid"
        When:
            Table is instantiated
        Then:
            It should raise ValueError mentioning interval_type
        """
        # Arrange / Act / Assert
        with pytest.raises(ValueError, match="interval_type"):
            Table(name="peaks", interval_type="invalid")

    @given(
        coordinate_system=st.sampled_from(["0based", "1based"]),
        interval_type=st.sampled_from(["half_open", "closed"]),
    )
    @settings(max_examples=20)
    def test___init___should_not_raise_when_params_are_valid(
        self, coordinate_system, interval_type
    ):
        """Test Table never raises for any valid parameter combination.

        Given:
            Any Table with valid coordinate_system and interval_type
        When:
            Table is instantiated
        Then:
            It should not raise and all fields should be accessible
        """
        # Arrange / Act
        table = Table(
            name="test",
            coordinate_system=coordinate_system,
            interval_type=interval_type,
        )

        # Assert
        assert table.coordinate_system == coordinate_system
        assert table.interval_type == interval_type


class TestTables:
    """Tests for the Tables container class."""

    def test_get_should_return_none_when_name_absent(self):
        """Test get returns None for an unregistered name.

        Given:
            A fresh Tables instance
        When:
            get is called with an unregistered name
        Then:
            It should return None
        """
        # Arrange
        tables = Tables()

        # Act / Assert
        assert tables.get("unknown") is None

    def test_get_should_return_table_when_name_registered(self):
        """Test get returns the Table for a registered name.

        Given:
            A Tables instance with one registered table
        When:
            get is called with the registered name
        Then:
            It should return the registered Table object
        """
        # Arrange
        tables = Tables()
        table = Table(name="peaks")
        tables.register("peaks", table)

        # Act / Assert
        assert tables.get("peaks") is table

    def test_register_should_store_all_tables_when_called_multiple_times(self):
        """Test register stores every table when called with distinct names.

        Given:
            A Tables instance with one registered table
        When:
            register is called with a new name and Table
        Then:
            It should make both tables retrievable via get
        """
        # Arrange
        tables = Tables()
        peaks = Table(name="peaks")
        variants = Table(name="variants")
        tables.register("peaks", peaks)
        tables.register("variants", variants)

        # Act / Assert
        assert tables.get("peaks") is peaks
        assert tables.get("variants") is variants

    def test_register_should_overwrite_when_name_already_registered(self):
        """Test register overwrites an existing entry with the same name.

        Given:
            A Tables instance with a registered table
        When:
            register is called with the same name and a different Table
        Then:
            It should make get return the new Table
        """
        # Arrange
        tables = Tables()
        old_table = Table(name="peaks")
        new_table = Table(name="peaks", chrom_col="chr")
        tables.register("peaks", old_table)
        tables.register("peaks", new_table)

        # Act / Assert
        assert tables.get("peaks") is new_table

    def test___contains___should_return_true_when_name_registered(self):
        """Test __contains__ returns True for a registered name.

        Given:
            A Tables instance with registered tables
        When:
            the in operator is used with a registered name
        Then:
            It should return True
        """
        # Arrange
        tables = Tables()
        tables.register("peaks", Table(name="peaks"))

        # Act / Assert
        assert "peaks" in tables

    def test___contains___should_return_false_when_name_absent(self):
        """Test __contains__ returns False for an unregistered name.

        Given:
            A Tables instance with registered tables
        When:
            the in operator is used with an unregistered name
        Then:
            It should return False
        """
        # Arrange
        tables = Tables()
        tables.register("peaks", Table(name="peaks"))

        # Act / Assert
        assert "unknown" not in tables

    def test___iter___should_yield_all_registered(self):
        """Test __iter__ yields all registered Table objects.

        Given:
            A Tables instance with registered tables
        When:
            iterated with a for loop
        Then:
            It should yield all registered Table objects
        """
        # Arrange
        tables = Tables()
        peaks = Table(name="peaks")
        variants = Table(name="variants")
        tables.register("peaks", peaks)
        tables.register("variants", variants)

        # Act
        result = []
        for table in tables:
            result.append(table)

        # Assert
        assert len(result) == 2
        assert peaks in result
        assert variants in result

    def test___iter___should_yield_nothing_when_empty(self):
        """Test __iter__ yields nothing when no tables are registered.

        Given:
            A fresh Tables instance with no tables
        When:
            iterated with a for loop
        Then:
            It should yield nothing
        """
        # Arrange
        tables = Tables()

        # Act
        result = []
        for table in tables:
            result.append(table)

        # Assert
        assert result == []
