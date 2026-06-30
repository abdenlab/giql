"""Cross-target coordinate-encoding sweep on the DataFusion engine (#145).

The base ``cross_target_oracle`` lane (``test_cross_target_oracle.py``) runs every
operator at the default 0-based half-open encoding. Issue #145 extends the
DataFusion coverage across all four ``(coordinate_system, interval_type)``
conventions and the custom-column / strand schema axes, so the capability-driven
canonicalization the canonicalizer and the NEAREST passthrough now perform
(``SELECT * REPLACE`` on DuckDB vs. the portable ``SELECT * EXCEPT`` form on the
generic / DataFusion family) is actually executed on the real DataFusion engine
for every encoding — not just proven byte-shaped in the unit suite.

For DISJOIN the full three-target oracle runs: the generic and datafusion targets
emit the portable ``* EXCEPT`` canonical wrapper and passthrough, both executed on
DataFusion, while duckdb runs the ``* REPLACE`` form on DuckDB.

For NEAREST the sweep restricts to ``("datafusion", "duckdb")``. A non-canonical
correlated NEAREST has no runnable generic target: the generic capability set is
lateral-capable (so it emits the ``LATERAL`` form, which DataFusion cannot plan)
yet not ``* REPLACE``-capable (so it emits ``* EXCEPT``, which DuckDB cannot run).
The datafusion target (decorrelated fallback + ``* EXCEPT``, on DataFusion) and the
duckdb target (``LATERAL`` + ``* REPLACE``, on DuckDB) still pin cross-engine
identity for every encoding.
"""

import pytest

pytest.importorskip("duckdb")
pytest.importorskip("datafusion")
pytest.importorskip("pyarrow")

from giql import Table  # noqa: E402

from ..coordinate_space.encodings import CONVENTIONS  # noqa: E402
from ..coordinate_space.encodings import encode  # noqa: E402
from ..coordinate_space.encodings import make_table  # noqa: E402

pytestmark = pytest.mark.integration


def _interval(canonical_start, canonical_end, coordinate_system, interval_type):
    """Encode a canonical interval into a ``(start, end)`` tuple for a convention."""
    return encode(canonical_start, canonical_end, coordinate_system, interval_type)


class TestDisjoinEncodingSweepOnDataFusion:
    """DISJOIN canonical-wrapper + passthrough across all encodings on DataFusion."""

    @pytest.mark.parametrize(("coordinate_system", "interval_type"), CONVENTIONS)
    def test_disjoin_sweep_agrees_across_targets(
        self, cross_target_oracle, coordinate_system, interval_type
    ):
        """Test DISJOIN agrees across targets for every coordinate encoding.

        Given:
            Two overlapping intervals stored under one of the four
            (coordinate_system, interval_type) conventions, declared on a Table
            with that encoding.
        When:
            The DISJOIN oracle runs all three targets — generic and datafusion
            emit the portable EXCEPT canonical wrapper and passthrough (executed
            on DataFusion), duckdb emits the REPLACE form (executed on DuckDB).
        Then:
            Every target should return the same sub-segments, with the
            passthrough and disjoin columns de-canonicalized back to the declared
            encoding, proving the EXCEPT canonicalization runs on DataFusion for
            every encoding.
        """
        # Arrange
        cs, it = coordinate_system, interval_type
        feats = [
            ("chr1",) + _interval(0, 100, cs, it),
            ("chr1",) + _interval(50, 150, cs, it),
        ]
        canonical_rows = [
            ((0, 100), (0, 50)),
            ((0, 100), (50, 100)),
            ((50, 150), (50, 100)),
            ((50, 150), (100, 150)),
        ]
        expected = [
            _interval(*parent, cs, it) + _interval(*sub, cs, it)
            for parent, sub in canonical_rows
        ]

        # Act & assert
        cross_target_oracle(
            'SELECT start, "end", disjoin_start, disjoin_end FROM DISJOIN(feats)',
            tables=[make_table("feats", cs, it)],
            feats=feats,
            expected=expected,
        )


class TestNearestEncodingSweepOnDataFusion:
    """NEAREST canonical-wrapper + passthrough across encodings on DataFusion."""

    @pytest.mark.parametrize(("coordinate_system", "interval_type"), CONVENTIONS)
    def test_correlated_nearest_non_canonical_target_agrees_datafusion_vs_duckdb(
        self, cross_target_oracle, coordinate_system, interval_type
    ):
        """Test correlated NEAREST over a non-canonical target agrees per encoding.

        Given:
            A canonical peaks table and a candidate genes table stored under one
            of the four conventions, declared with that encoding.
        When:
            A correlated NEAREST k=1 runs on the datafusion target (decorrelated
            fallback + portable EXCEPT passthrough, on DataFusion) and the duckdb
            target (LATERAL + REPLACE, on DuckDB).
        Then:
            Both should return the same nearest gene with its interval columns
            de-canonicalized back to the declared encoding and an
            encoding-invariant distance, proving the EXCEPT canonical wrapper and
            passthrough run on DataFusion for every encoding.
        """
        # Arrange
        cs, it = coordinate_system, interval_type
        genes = [
            ("chr1",) + _interval(1000, 1100, cs, it),
            ("chr1",) + _interval(50, 60, cs, it),
            ("chr1",) + _interval(280, 290, cs, it),
        ]
        nearest_start, nearest_end = _interval(280, 290, cs, it)

        # Act & assert
        cross_target_oracle(
            "SELECT b.start AS b_start, b.\"end\" AS b_end, distance "
            "FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b",
            tables=[Table("peaks"), make_table("genes", cs, it)],
            peaks=[("chr1", 200, 300)],
            genes=genes,
            expected=[(nearest_start, nearest_end, 0)],
            targets=("datafusion", "duckdb"),
        )


    def test_correlated_nearest_custom_columns_agree_datafusion_vs_duckdb(
        self, cross_target_oracle
    ):
        """Test correlated NEAREST over a non-canonical custom-column target agrees.

        Given:
            A peaks table and a non-canonical (1based/closed) genes table that both
            name their interval columns gc / gs / ge rather than chrom/start/end.
        When:
            A correlated NEAREST k=1 runs on the datafusion target (decorrelated
            fallback + EXCEPT passthrough over the custom names) and the duckdb
            target (LATERAL + REPLACE).
        Then:
            Both should return the same nearest gene with its custom-named interval
            columns de-canonicalized, proving the EXCEPT NEAREST passthrough
            resolves custom column names on DataFusion (the DISJOIN custom-column
            path's NEAREST counterpart).
        """
        # Arrange — both tables share one schema (the oracle registers all table
        # data under a single `columns` mapping); only genes is non-canonical.
        columns = (("gc", "utf8"), ("gs", "int64"), ("ge", "int64"))
        peaks = Table("peaks", chrom_col="gc", start_col="gs", end_col="ge")
        genes = Table(
            "genes",
            chrom_col="gc",
            start_col="gs",
            end_col="ge",
            coordinate_system="1based",
            interval_type="closed",
        )
        nearest_start, nearest_end = encode(280, 290, "1based", "closed")

        # Act & assert
        cross_target_oracle(
            "SELECT b.gs AS b_start, b.ge AS b_end, distance "
            "FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b",
            tables=[peaks, genes],
            columns=columns,
            peaks=[("chr1", 200, 300)],
            genes=[
                ("chr1",) + encode(1000, 1100, "1based", "closed"),
                ("chr1",) + encode(50, 60, "1based", "closed"),
                ("chr1",) + encode(280, 290, "1based", "closed"),
            ],
            expected=[(nearest_start, nearest_end, 0)],
            targets=("datafusion", "duckdb"),
        )


class TestDataFusionSchemaAxes:
    """Custom column names and strand schemas exercised on DataFusion (#145)."""

    def test_disjoin_custom_column_names_agree_across_targets(
        self, cross_target_oracle
    ):
        """Test DISJOIN over a custom-column table agrees across targets.

        Given:
            Two overlapping intervals in a table whose interval columns are named
            ch / s / e rather than the canonical chrom / start / end.
        When:
            The DISJOIN oracle runs all three targets.
        Then:
            Every target should resolve the custom column names and return the
            same sub-segments, proving custom-column resolution runs on DataFusion.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT s, e, disjoin_start, disjoin_end FROM DISJOIN(feats)",
            tables=[Table("feats", chrom_col="ch", start_col="s", end_col="e")],
            columns=(("ch", "utf8"), ("s", "int64"), ("e", "int64")),
            feats=[("chr1", 0, 100), ("chr1", 50, 150)],
            expected=[
                (0, 100, 0, 50),
                (0, 100, 50, 100),
                (50, 150, 50, 100),
                (50, 150, 100, 150),
            ],
        )

    def test_stranded_merge_agrees_across_targets(self, cross_target_oracle):
        """Test stranded MERGE aggregates within strand on every target.

        Given:
            Overlapping intervals on chr1 split across the + and - strands.
        When:
            A stranded MERGE runs on all three targets (generic and datafusion on
            DataFusion, duckdb on DuckDB).
        Then:
            Every target should merge only same-strand overlaps, yielding one
            merged interval per strand, and agree — exercising the strand schema
            axis on DataFusion.
        """
        # Arrange / Act / Assert — MERGE rewrites the whole SELECT, projecting
        # chrom, strand, MIN(start) AS start, MAX(end) AS end.
        cross_target_oracle(
            "SELECT MERGE(interval, stranded := true) FROM feats",
            tables=[Table("feats", strand_col="strand")],
            columns=(
                ("chrom", "utf8"),
                ("start", "int64"),
                ("end", "int64"),
                ("strand", "utf8"),
            ),
            feats=[
                ("chr1", 0, 100, "+"),
                ("chr1", 50, 150, "+"),
                ("chr1", 0, 100, "-"),
            ],
            expected=[
                ("chr1", "+", 0, 150),
                ("chr1", "-", 0, 100),
            ],
        )
