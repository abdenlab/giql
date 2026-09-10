"""Runner-internals tests for the polars-bio oracle engine (#77).

These exercise :func:`run_polars_bio` directly with trivial ``SELECT ... FROM t``
SQL (no ``transpile``) so any failure localizes to the runner: explicit dtype
mapping (including the 32-bit kind), empty-table registration, the
overwrite-on-re-register isolation the runner relies on, and output-shape parity
with the vanilla DataFusion runner.
"""

import pytest

pytest.importorskip("polars_bio")
pytest.importorskip("polars")
pytest.importorskip("datafusion")
pytest.importorskip("pyarrow")

from tests.integration._oracle import run_datafusion  # noqa: E402
from tests.integration._oracle import run_polars_bio  # noqa: E402

pytestmark = pytest.mark.integration

_INTERVAL_COLUMNS = (("chrom", "utf8"), ("start", "int64"), ("end", "int64"))
_INT32_COLUMNS = (("chrom", "utf8"), ("start", "int32"), ("end", "int32"))


class TestRunPolarsBio:
    """`run_polars_bio` schema handling and session isolation."""

    def test_run_polars_bio_should_return_zero_rows_when_table_empty(self):
        """Test run_polars_bio registers an empty table with its declared schema.

        Given:
            An empty interval table.
        When:
            run_polars_bio selects from it.
        Then:
            It should return zero rows rather than failing on schema inference.
        """
        # Arrange / Act
        rows = run_polars_bio("SELECT chrom FROM t", {"t": []}, _INTERVAL_COLUMNS)

        # Assert
        assert rows == []

    def test_run_polars_bio_should_round_trip_int32_columns(self):
        """Test the 32-bit column kind is honored at the engine.

        Given:
            An interval table declared with ``int32`` coordinate columns.
        When:
            run_polars_bio selects the coordinates and their Arrow type name.
        Then:
            It should return the rows unchanged and report a 32-bit type.
        """
        # Arrange
        data = {"t": [("chr1", 5, 10)]}

        # Act
        rows = run_polars_bio(
            'SELECT chrom, start, "end", arrow_typeof(start) AS ty FROM t',
            data,
            _INT32_COLUMNS,
        )

        # Assert
        assert rows == [("chr1", 5, 10, "Int32")]

    def test_run_polars_bio_should_overwrite_when_table_reregistered(self):
        """Test that re-registering a table name replaces the earlier contents.

        Given:
            A table registered with one row, then registered again under the
            same name with a different row.
        When:
            run_polars_bio selects from it after the second registration.
        Then:
            It should see only the second contents, the isolation guarantee the
            process-global session runner depends on.
        """
        # Arrange
        run_polars_bio("SELECT chrom FROM t", {"t": [("chrA", 1, 2)]}, _INTERVAL_COLUMNS)

        # Act
        rows = run_polars_bio(
            "SELECT chrom FROM t", {"t": [("chrB", 3, 4)]}, _INTERVAL_COLUMNS
        )

        # Assert
        assert rows == [("chrB",)]

    def test_run_polars_bio_should_match_run_datafusion_on_trivial_select(self):
        """Test both DataFusion-family runners produce identical normalized output.

        Given:
            The same table data and trivial SELECT for both engines.
        When:
            run_datafusion and run_polars_bio each execute it.
        Then:
            Their normalized results should be identical, the parity the oracle
            relies on for its differential comparison.
        """
        # Arrange
        data = {"t": [("chr1", 1, 2), ("chr2", 3, 4)]}
        sql = 'SELECT chrom, start, "end" FROM t'

        # Act
        vanilla = run_datafusion(sql, data, _INTERVAL_COLUMNS)
        bio = run_polars_bio(sql, data, _INTERVAL_COLUMNS)

        # Assert
        assert vanilla == bio
