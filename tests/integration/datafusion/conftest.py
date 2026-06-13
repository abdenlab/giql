"""Fixtures for DataFusion target execution smoke tests (issue #132).

These tests prove that SQL produced via ``transpile(..., dialect="datafusion")``
parses and executes on a real DataFusion engine. The broad cross-target result
oracle and the full DataFusion integration lane are deferred to #139.
"""

import pytest

pytest.importorskip("datafusion")
pytest.importorskip("pyarrow")

pytestmark = pytest.mark.integration


@pytest.fixture
def datafusion_ctx():
    """Return a builder that registers interval tables into a SessionContext.

    The returned callable accepts ``name=[(chrom, start, end), ...]`` keyword
    arguments and yields a DataFusion ``SessionContext`` with each table
    registered under the default ``chrom``/``start``/``end`` schema (matching
    the default :class:`giql.Table` column mapping).
    """
    import pyarrow as pa
    from datafusion import SessionContext

    schema = pa.schema(
        [
            ("chrom", pa.utf8()),
            ("start", pa.int64()),
            ("end", pa.int64()),
        ]
    )

    def _build(**tables):
        ctx = SessionContext()
        for name, rows in tables.items():
            arrays = {
                "chrom": [r[0] for r in rows],
                "start": [r[1] for r in rows],
                "end": [r[2] for r in rows],
            }
            ctx.register_record_batches(
                name, [pa.table(arrays, schema=schema).to_batches()]
            )
        return ctx

    return _build
