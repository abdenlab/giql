"""Fixtures for the ``dialect="datafusion-bio"`` execution lane (#77).

These tests prove that SQL produced via ``transpile(..., dialect="datafusion-bio")``
executes on polars-bio's DataFusion session and, where the shape allows, is
planned through its native ``IntervalJoinExec``. Row identity comes from the
cross-target oracle (``tests/integration/conftest.py``) with ``datafusion-bio``
opted in as a third engine; plan identity comes from :func:`polars_bio_explain`,
since row identity alone cannot tell the native join from the hash-join fallback.

polars-bio keeps one process-global session and prints a tqdm progress bar on
every ``collect()``; the bar is silenced here through tqdm's own environment
switch so the lane's output stays readable.
"""

import os

import pytest

pytest.importorskip("polars_bio")
pytest.importorskip("polars")

pytestmark = pytest.mark.integration

os.environ.setdefault("TQDM_DISABLE", "1")


@pytest.fixture
def polars_bio_explain():
    """Return a callable that renders polars-bio's ``EXPLAIN`` output as one string.

    ``pb.sql("EXPLAIN ...")`` yields a two-column frame (plan type, plan text);
    the callable flattens every cell so a test can substring-assert on operator
    names such as ``IntervalJoinExec`` or ``HashJoinExec``.
    """

    def _explain(sql: str) -> str:
        import polars_bio as pb

        frame = pb.sql("EXPLAIN " + sql).collect()
        return "\n".join(str(value) for row in frame.rows() for value in row)

    return _explain


@pytest.fixture
def polars_bio_register():
    """Return a callable that registers one polars frame under a table name.

    Accepts ``(name, columns, rows)`` where ``columns`` is a ``{name: polars
    dtype}`` mapping so a test can register tables with different schemas and
    integer widths side by side (the oracle fixture applies one column spec to
    every table). Registration overwrites any prior table of the same name.
    """

    def _register(name: str, columns: dict, rows: list[tuple]) -> None:
        import polars as pl
        import polars_bio as pb

        frame = pl.DataFrame([tuple(r) for r in rows], schema=columns, orient="row")
        pb.from_polars(name, frame)

    return _register
