"""Tests for giql.table module."""

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from giql.table import Table
from giql.table import Tables


class TestTable:
    """Tests for the Table dataclass."""

    def test_default_values(self):
        """
        GIVEN only the required arg `name`
        WHEN Table is instantiated
        THEN all fields have their default values.
        """
        table = Table(name="peaks")

        assert table.name == "peaks"
        assert table.genomic_col == "interval"
        assert table.chrom_col == "chrom"
        assert table.start_col == "start"
        assert table.end_col == "end"
        assert table.strand_col == "strand"
        assert table.coordinate_system == "0based"
        assert table.interval_type == "half_open"

    def test_all_custom_values(self):
        """
        GIVEN all fields provided with custom values
        WHEN Table is instantiated
        THEN all fields reflect the custom values.
        """
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

        assert table.name == "variants"
        assert table.genomic_col == "position"
        assert table.chrom_col == "chr"
        assert table.start_col == "pos_start"
        assert table.end_col == "pos_end"
        assert table.strand_col == "direction"
        assert table.coordinate_system == "1based"
        assert table.interval_type == "closed"

    def test_strand_col_none(self):
        """
        GIVEN strand_col=None
        WHEN Table is instantiated
        THEN strand_col is None.
        """
        table = Table(name="peaks", strand_col=None)

        assert table.strand_col is None

    def test_coordinate_system_1based(self):
        """
        GIVEN coordinate_system="1based"
        WHEN Table is instantiated
        THEN coordinate_system is "1based".
        """
        table = Table(name="peaks", coordinate_system="1based")

        assert table.coordinate_system == "1based"

    def test_interval_type_closed(self):
        """
        GIVEN interval_type="closed"
        WHEN Table is instantiated
        THEN interval_type is "closed".
        """
        table = Table(name="peaks", interval_type="closed")

        assert table.interval_type == "closed"

    def test_invalid_coordinate_system(self):
        """
        GIVEN coordinate_system="invalid"
        WHEN Table is instantiated
        THEN raises ValueError with message about valid options.
        """
        with pytest.raises(ValueError, match="coordinate_system"):
            Table(name="peaks", coordinate_system="invalid")

    def test_invalid_interval_type(self):
        """
        GIVEN interval_type="invalid"
        WHEN Table is instantiated
        THEN raises ValueError with message about valid options.
        """
        with pytest.raises(ValueError, match="interval_type"):
            Table(name="peaks", interval_type="invalid")

    @given(
        coordinate_system=st.sampled_from(["0based", "1based"]),
        interval_type=st.sampled_from(["half_open", "closed"]),
    )
    @settings(max_examples=20)
    def test_valid_params_never_raise(self, coordinate_system, interval_type):
        """
        GIVEN any Table with valid coordinate_system and interval_type
        WHEN Table is instantiated
        THEN no exception is raised and all fields are accessible.
        """
        table = Table(
            name="test",
            coordinate_system=coordinate_system,
            interval_type=interval_type,
        )

        assert table.coordinate_system == coordinate_system
        assert table.interval_type == interval_type


class TestTables:
    """Tests for the Tables container class."""

    def test_get_missing_key(self):
        """
        GIVEN a fresh Tables instance
        WHEN get is called with an unregistered name
        THEN returns None.
        """
        tables = Tables()

        assert tables.get("unknown") is None

    def test_get_existing_key(self):
        """
        GIVEN a Tables instance with one registered table
        WHEN get is called with the registered name
        THEN returns the Table object.
        """
        tables = Tables()
        table = Table(name="peaks")
        tables.register("peaks", table)

        assert tables.get("peaks") is table

    def test_register_multiple_tables(self):
        """
        GIVEN a Tables instance with one registered table
        WHEN register is called with a new name and Table
        THEN both tables are retrievable via get.
        """
        tables = Tables()
        peaks = Table(name="peaks")
        variants = Table(name="variants")
        tables.register("peaks", peaks)
        tables.register("variants", variants)

        assert tables.get("peaks") is peaks
        assert tables.get("variants") is variants

    def test_register_overwrites(self):
        """
        GIVEN a Tables instance with a registered table
        WHEN register is called with the same name and a different Table
        THEN get returns the new Table (overwrite).
        """
        tables = Tables()
        old_table = Table(name="peaks")
        new_table = Table(name="peaks", chrom_col="chr")
        tables.register("peaks", old_table)
        tables.register("peaks", new_table)

        assert tables.get("peaks") is new_table

    def test_contains(self):
        """
        GIVEN a Tables instance with registered tables
        WHEN the in operator is used
        THEN returns True for registered names, False for others.
        """
        tables = Tables()
        tables.register("peaks", Table(name="peaks"))

        assert "peaks" in tables
        assert "unknown" not in tables

    def test_iter(self):
        """
        GIVEN a Tables instance with registered tables
        WHEN iterated with a for loop
        THEN yields all registered Table objects.
        """
        tables = Tables()
        peaks = Table(name="peaks")
        variants = Table(name="variants")
        tables.register("peaks", peaks)
        tables.register("variants", variants)

        result = []
        for table in tables:
            result.append(table)

        assert len(result) == 2
        assert peaks in result
        assert variants in result

    def test_iter_empty(self):
        """
        GIVEN a fresh Tables instance with no tables
        WHEN iterated with a for loop
        THEN yields nothing (empty iteration).
        """
        tables = Tables()

        result = []
        for table in tables:
            result.append(table)

        assert result == []
