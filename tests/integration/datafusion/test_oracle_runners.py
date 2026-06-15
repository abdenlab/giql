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
    """`run_duckdb` empty-table handling (its ``if rows`` guard)."""

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


class TestRunDataFusion:
    """`run_datafusion` empty-table handling (the synthesized-batch guard)."""

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
        with pytest.raises(BaseException, match="index out of bounds"):
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

    def test_runners_agree_on_empty_table(self):
        """Test both runners produce identical output for an empty table.

        Given:
            An empty table and a trivial SELECT — DuckDB skips the INSERT while
            DataFusion synthesizes an empty record batch.
        When:
            run_duckdb and run_datafusion each execute it.
        Then:
            Both should return zero rows, proving the two divergent empty-table
            paths converge on identical normalized output.
        """
        # Arrange
        data = {"t": []}
        sql = "SELECT chrom FROM t"

        # Act
        ddb = run_duckdb(sql, data, _INTERVAL_COLUMNS)
        df = run_datafusion(sql, data, _INTERVAL_COLUMNS)

        # Assert
        assert ddb == df == []

    def test_runners_agree_on_null_round_trip(self):
        """Test both runners normalize a NULL-bearing column identically.

        Given:
            A table with a NULL ``start`` in one row — DataFusion promotes the
            int64 column to float64 around the NULL while DuckDB keeps ints.
        When:
            run_duckdb and run_datafusion each select it.
        Then:
            Their normalized output should be identical, proving the scalar
            coercion collapses the NULL surrogate and the int/float promotion.
        """
        # Arrange
        data = {"t": [("chr1", None, 2), ("chr2", 3, 4)]}
        sql = 'SELECT chrom, start, "end" FROM t'

        # Act
        ddb = run_duckdb(sql, data, _INTERVAL_COLUMNS)
        df = run_datafusion(sql, data, _INTERVAL_COLUMNS)

        # Assert
        assert ddb == df
