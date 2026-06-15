"""Runner-internals tests for the cross-target oracle (#139, T4).

These exercise :func:`run_duckdb` / :func:`run_datafusion` directly with trivial
``SELECT ... FROM t`` SQL (no ``transpile``) so any failure localizes to the
runner -- schema construction, type mapping, reserved-word quoting, the
empty-table handling that differs between engines, custom columns, and
cross-runner output-shape parity.
"""

import pytest

pytest.importorskip("duckdb")
pytest.importorskip("datafusion")
pytest.importorskip("pyarrow")

from tests.integration._oracle import arrow_schema  # noqa: E402
from tests.integration._oracle import run_datafusion  # noqa: E402
from tests.integration._oracle import run_duckdb  # noqa: E402

pytestmark = pytest.mark.integration

_INTERVAL_COLUMNS = (("chrom", "utf8"), ("start", "int64"), ("end", "int64"))


class TestRunDuckDB:
    """`run_duckdb` schema, types, quoting, and empty handling."""

    def test_run_duckdb_returns_normalized_rows(self):
        """Test run_duckdb registers a table and returns normalized rows.

        Given:
            A two-row interval table and a plain SELECT.
        When:
            run_duckdb executes the SQL.
        Then:
            It should return both rows as sorted plain-tuple values.
        """
        # Arrange / Act
        rows = run_duckdb(
            'SELECT chrom, start, "end" FROM t',
            {"t": [("chr2", 5, 6), ("chr1", 1, 2)]},
            _INTERVAL_COLUMNS,
        )

        # Assert
        assert rows == [("chr1", 1, 2), ("chr2", 5, 6)]

    def test_run_duckdb_maps_utf8_and_int64_types(self):
        """Test run_duckdb maps utf8 to VARCHAR and int64 to BIGINT.

        Given:
            A row whose chrom is a string and coordinates are ints.
        When:
            run_duckdb executes a SELECT.
        Then:
            The returned values should preserve str and int Python types.
        """
        # Arrange / Act
        rows = run_duckdb(
            "SELECT chrom, start FROM t", {"t": [("chr1", 1, 2)]}, _INTERVAL_COLUMNS
        )

        # Assert
        assert isinstance(rows[0][0], str)
        assert isinstance(rows[0][1], int)

    def test_run_duckdb_quotes_reserved_end_column(self):
        """Test run_duckdb quotes the reserved ``end`` column.

        Given:
            A table with the reserved-word ``end`` column.
        When:
            run_duckdb selects the quoted column.
        Then:
            It should return the end value without a parse error.
        """
        # Arrange / Act
        rows = run_duckdb(
            'SELECT "end" FROM t', {"t": [("chr1", 1, 99)]}, _INTERVAL_COLUMNS
        )

        # Assert
        assert rows == [(99,)]

    def test_run_duckdb_empty_table_returns_zero_rows(self):
        """Test run_duckdb tolerates an empty table via its ``if rows`` guard.

        Given:
            An empty interval table.
        When:
            run_duckdb selects from it.
        Then:
            It should return zero rows (the guard skips the INSERT).
        """
        # Arrange / Act
        rows = run_duckdb("SELECT chrom FROM t", {"t": []}, _INTERVAL_COLUMNS)

        # Assert
        assert rows == []

    def test_run_duckdb_custom_columns(self):
        """Test run_duckdb registers a custom column schema.

        Given:
            A table with a custom ``(name, score)`` schema.
        When:
            run_duckdb selects from it.
        Then:
            It should map the custom utf8/int64 columns and return the row.
        """
        # Arrange / Act
        rows = run_duckdb(
            "SELECT name, score FROM t",
            {"t": [("gene1", 42)]},
            (("name", "utf8"), ("score", "int64")),
        )

        # Assert
        assert rows == [("gene1", 42)]


class TestRunDataFusion:
    """`run_datafusion` schema, types, quoting, and empty handling."""

    def test_run_datafusion_returns_normalized_rows(self):
        """Test run_datafusion registers a table and returns normalized rows.

        Given:
            A two-row interval table and a plain SELECT.
        When:
            run_datafusion executes the SQL.
        Then:
            It should return both rows as sorted plain-tuple values.
        """
        # Arrange / Act
        rows = run_datafusion(
            'SELECT chrom, start, "end" FROM t',
            {"t": [("chr2", 5, 6), ("chr1", 1, 2)]},
            _INTERVAL_COLUMNS,
        )

        # Assert
        assert rows == [("chr1", 1, 2), ("chr2", 5, 6)]

    def test_run_datafusion_maps_utf8_and_int64_types(self):
        """Test run_datafusion maps utf8 to pyarrow utf8 and int64 to int64.

        Given:
            A row whose chrom is a string and coordinates are ints.
        When:
            run_datafusion executes a SELECT.
        Then:
            The returned values should be plain str and int after coercion.
        """
        # Arrange / Act
        rows = run_datafusion(
            "SELECT chrom, start FROM t", {"t": [("chr1", 1, 2)]}, _INTERVAL_COLUMNS
        )

        # Assert
        assert isinstance(rows[0][0], str)
        assert isinstance(rows[0][1], int)

    def test_run_datafusion_quotes_reserved_end_column(self):
        """Test run_datafusion quotes the reserved ``end`` column.

        Given:
            A table with the reserved-word ``end`` column.
        When:
            run_datafusion selects the quoted column.
        Then:
            It should return the end value without a parse error.
        """
        # Arrange / Act
        rows = run_datafusion(
            'SELECT "end" FROM t', {"t": [("chr1", 1, 99)]}, _INTERVAL_COLUMNS
        )

        # Assert
        assert rows == [(99,)]

    def test_run_datafusion_empty_table_returns_zero_rows(self):
        """Test run_datafusion handles an empty table via a synthesized batch.

        Given:
            An empty interval table (which would panic raw
            ``register_record_batches``).
        When:
            run_datafusion selects from it.
        Then:
            It should return zero rows, proving the loader synthesizes an empty
            record batch where DuckDB instead skips the INSERT.
        """
        # Arrange / Act
        rows = run_datafusion("SELECT chrom FROM t", {"t": []}, _INTERVAL_COLUMNS)

        # Assert
        assert rows == []

    def test_run_datafusion_custom_columns(self):
        """Test run_datafusion registers a custom column schema.

        Given:
            A table with a custom ``(name, score)`` schema.
        When:
            run_datafusion selects from it.
        Then:
            It should map the custom utf8/int64 columns and return the row.
        """
        # Arrange / Act
        rows = run_datafusion(
            "SELECT name, score FROM t",
            {"t": [("gene1", 42)]},
            (("name", "utf8"), ("score", "int64")),
        )

        # Assert
        assert rows == [("gene1", 42)]


class TestRawRegisterRecordBatchesEmpty:
    """Pin the raw DataFusion empty-table panic the runner guards against."""

    def test_raw_register_record_batches_panics_on_empty(self):
        """Test raw register_record_batches panics on a truly-empty partition.

        Given:
            An empty pyarrow table whose ``to_batches()`` yields no batches.
        When:
            ``SessionContext.register_record_batches`` is called with it
            directly (NOT via the oracle helper).
        Then:
            DataFusion should raise/panic on the empty partition, documenting
            why :func:`run_datafusion` synthesizes an empty batch instead.
        """
        # Arrange
        import pyarrow as pa
        from datafusion import SessionContext

        schema = arrow_schema(_INTERVAL_COLUMNS)
        ctx = SessionContext()
        empty_batches = pa.table(
            {"chrom": [], "start": [], "end": []}, schema=schema
        ).to_batches()
        assert empty_batches == []  # no batches is the panic trigger

        # Act / Assert
        with pytest.raises(BaseException):
            ctx.register_record_batches("t", [empty_batches])


class TestCrossRunnerParity:
    """Output-shape parity between the two runners on trivial SQL."""

    def test_runners_agree_on_shape_and_values(self):
        """Test both runners produce identical normalized output for one query.

        Given:
            The same table data and trivial SELECT for both engines.
        When:
            run_duckdb and run_datafusion each execute it.
        Then:
            Their normalized results should be identical -- the parity the
            oracle relies on for its differential comparison.
        """
        # Arrange
        data = {"t": [("chr1", 1, 2), ("chr2", 3, 4)]}
        sql = 'SELECT chrom, start, "end" FROM t'

        # Act
        ddb = run_duckdb(sql, data, _INTERVAL_COLUMNS)
        df = run_datafusion(sql, data, _INTERVAL_COLUMNS)

        # Assert
        assert ddb == df
