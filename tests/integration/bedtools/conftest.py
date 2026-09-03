"""Pytest fixtures for bedtools integration tests."""

import shutil

import pytest

from giql import transpile

duckdb = pytest.importorskip("duckdb")
pytest.importorskip("pybedtools")

if not shutil.which("bedtools"):
    pytest.skip(
        "bedtools binary not found in PATH",
        allow_module_level=True,
    )

pytestmark = pytest.mark.integration

from .utils.data_models import GenomicInterval  # noqa: E402
from .utils.duckdb_loader import load_intervals  # noqa: E402


@pytest.fixture(scope="function")
def duckdb_connection():
    """Provide clean DuckDB connection for each test.

    Each test gets a fresh in-memory database with no shared state.
    """
    conn = duckdb.connect(":memory:")
    yield conn
    conn.close()


@pytest.fixture
def partial_overlap_intervals(duckdb_connection):
    """Load an A/B interval pair with partial overlaps on two chromosomes.

    Also carries one A interval with no overlap and one on a chromosome absent
    from B, so a dialect rewrite that drops unmatched contigs is visible. Yields
    the two tuple lists for handing to the bedtools wrapper.
    """
    intervals_a = [
        GenomicInterval("chr1", 100, 200, "a1", 100, "+"),
        GenomicInterval("chr1", 300, 400, "a2", 200, "-"),
        GenomicInterval("chr2", 50, 150, "a3", 300, "+"),
        GenomicInterval("chr1", 5000, 6000, "a_lonely", 400, "-"),
        GenomicInterval("chrX", 10, 20, "a_other_chrom", 500, "+"),
    ]
    intervals_b = [
        GenomicInterval("chr1", 150, 250, "b1", 100, "+"),
        GenomicInterval("chr1", 350, 450, "b2", 200, "+"),
        GenomicInterval("chr2", 100, 120, "b3", 300, "-"),
        GenomicInterval("chrY", 1, 5, "b_other_chrom", 400, "-"),
    ]
    rows_a = [i.to_tuple() for i in intervals_a]
    rows_b = [i.to_tuple() for i in intervals_b]
    load_intervals(duckdb_connection, "intervals_a", rows_a)
    load_intervals(duckdb_connection, "intervals_b", rows_b)
    yield rows_a, rows_b


@pytest.fixture(scope="function")
def giql_query(duckdb_connection):
    """Provide a helper that loads data, transpiles GIQL, and executes.

    Usage::

        result = giql_query(
            "SELECT * FROM t WHERE interval INTERSECTS 'chr1:1-100'",
            tables=["t"],
            t=[GenomicInterval("chr1", 50, 150, "x", 0, "+")],
        )
    """

    def _run(query: str, *, tables: list[str], **table_data):
        for name, intervals in table_data.items():
            load_intervals(
                duckdb_connection,
                name,
                [i.to_tuple() for i in intervals],
            )
        # Transpile for DuckDB: this helper always executes on the DuckDB
        # connection, and some operators (e.g. a star-projected CLUSTER, #184) emit
        # DuckDB-specific spellings such as ``* EXCLUDE`` that the portable generic
        # ``* EXCEPT`` form does not share.
        sql = transpile(query, tables=tables, dialect="duckdb")
        return duckdb_connection.execute(sql).fetchall()

    return _run
