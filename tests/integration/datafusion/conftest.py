"""Fixtures for DataFusion target execution smoke tests (issue #132).

These tests prove that SQL produced via ``transpile(..., dialect="datafusion")``
parses and executes on a real DataFusion engine. The broad cross-target result
oracle and the full DataFusion integration lane are deferred to #139.
"""

import pytest

pytest.importorskip("datafusion")
pytest.importorskip("pyarrow")

pytestmark = pytest.mark.integration

from .._oracle import arrow_schema  # noqa: E402
from .._oracle import register_record_batches  # noqa: E402

# The default GIQL interval schema, shared with the cross-target oracle.
_INTERVAL_COLUMNS = (("chrom", "utf8"), ("start", "int64"), ("end", "int64"))


@pytest.fixture
def datafusion_ctx():
    """Return a builder that registers interval tables into a SessionContext.

    The returned callable accepts ``name=[(chrom, start, end), ...]`` keyword
    arguments and yields a DataFusion ``SessionContext`` with each table
    registered under the default ``chrom``/``start``/``end`` schema (matching
    the default :class:`giql.Table` column mapping).

    The pyarrow loader dance is the one shared helper in
    :mod:`tests.integration._oracle` (Finding 8), so this fixture and the
    cross-target oracle register tables identically.
    """
    from datafusion import SessionContext

    schema = arrow_schema(_INTERVAL_COLUMNS)

    def _build(**tables):
        ctx = SessionContext()
        for name, rows in tables.items():
            register_record_batches(ctx, name, rows, schema, _INTERVAL_COLUMNS)
        return ctx

    return _build
