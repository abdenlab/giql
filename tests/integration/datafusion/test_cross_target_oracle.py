"""Cross-target result-identity tests via the ``cross_target_oracle`` fixture.

These exercise the reusable oracle (``tests/integration/conftest.py``) over the
operators that already emit identical generic SQL across Generic and DataFusion
and run correctly on DuckDB: INTERSECTS (literal + column-to-column join),
CONTAINS, WITHIN, and standalone NEAREST. No operator has been migrated to the
expander registry yet (epic #137), so this lane locks in the verification path
every later migration (#140-#144) will consume.

NEAREST's expansion uses a correlated ``LATERAL`` subquery, which DataFusion has
no physical plan for today; its case therefore restricts the oracle to the
DuckDB-executed targets and is the one documented cross-target gap.
"""

import pytest

pytest.importorskip("duckdb")
pytest.importorskip("datafusion")
pytest.importorskip("pyarrow")

pytestmark = pytest.mark.integration


class TestCrossTargetOracleIntersects:
    """INTERSECTS identity across Generic, DataFusion, and DuckDB."""

    def test_literal_intersects_returns_overlapping_rows(self, cross_target_oracle):
        """Test a literal INTERSECTS yields identical rows on every target.

        Given:
            A peaks table with one interval overlapping chr1:1000-2000, one
            non-overlapping interval on chr1, and one on chr2.
        When:
            A literal-range INTERSECTS query runs through every target on its
            engine.
        Then:
            Each target should return exactly the overlapping chr1 row and all
            targets should agree.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            'SELECT chrom, start, "end" FROM peaks '
            "WHERE interval INTERSECTS 'chr1:1000-2000'",
            peaks=[
                ("chr1", 1500, 1800),
                ("chr1", 5000, 6000),
                ("chr2", 1500, 1800),
            ],
            expected=[("chr1", 1500, 1800)],
        )

    def test_literal_intersects_no_overlap_returns_zero_rows(self, cross_target_oracle):
        """Test a non-overlapping literal INTERSECTS yields no rows on any target.

        Given:
            A peaks table whose intervals do not overlap the query range.
        When:
            The literal INTERSECTS query runs through every target.
        Then:
            Every target should return zero rows in agreement.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            'SELECT chrom, start, "end" FROM peaks '
            "WHERE interval INTERSECTS 'chr1:1000-2000'",
            peaks=[
                ("chr1", 5000, 6000),
                ("chr2", 1500, 1800),
            ],
            expected=[],
        )

    def test_column_to_column_intersects_join_agrees(self, cross_target_oracle):
        """Test a column-to-column INTERSECTS join agrees across all targets.

        Given:
            Peaks and genes tables where one peak overlaps one gene on chr1 and
            other intervals do not overlap or sit on another chromosome.
        When:
            A column-to-column INTERSECTS join runs — generic/datafusion via the
            binned equi-join on DataFusion and duckdb via the IEJoin on DuckDB.
        Then:
            Every target should return the single overlapping pair, proving the
            structurally different DuckDB IEJoin SQL yields identical rows.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.chrom, a.start AS a_start, b.start AS b_start "
            "FROM peaks a JOIN genes b ON a.interval INTERSECTS b.interval",
            peaks=[
                ("chr1", 100, 500),
                ("chr1", 1000, 2000),
                ("chr2", 100, 500),
            ],
            genes=[
                ("chr1", 300, 600),
                ("chr1", 5000, 6000),
                ("chr2", 9000, 9500),
            ],
            expected=[("chr1", 100, 300)],
        )


class TestCrossTargetOracleContainsWithin:
    """CONTAINS and WITHIN identity across all targets."""

    def test_contains_returns_enclosing_rows(self, cross_target_oracle):
        """Test CONTAINS yields identical enclosing rows across targets.

        Given:
            A peaks table with one interval enclosing chr1:1200-1300, one that
            only partially overlaps it, and one on chr2.
        When:
            A literal CONTAINS query runs through every target.
        Then:
            Every target should return only the fully-enclosing chr1 row.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            'SELECT chrom, start, "end" FROM peaks '
            "WHERE interval CONTAINS 'chr1:1200-1300'",
            peaks=[
                ("chr1", 1000, 2000),
                ("chr1", 1250, 1400),
                ("chr2", 1000, 2000),
            ],
            expected=[("chr1", 1000, 2000)],
        )

    def test_within_returns_enclosed_rows(self, cross_target_oracle):
        """Test WITHIN yields identical enclosed rows across targets.

        Given:
            A peaks table with one interval fully inside chr1:1000-2000, one that
            spills past it, and one on chr2.
        When:
            A literal WITHIN query runs through every target.
        Then:
            Every target should return only the fully-enclosed chr1 row.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            'SELECT chrom, start, "end" FROM peaks '
            "WHERE interval WITHIN 'chr1:1000-2000'",
            peaks=[
                ("chr1", 1200, 1800),
                ("chr1", 1500, 2500),
                ("chr2", 1200, 1800),
            ],
            expected=[("chr1", 1200, 1800)],
        )


class TestCrossTargetOracleNearest:
    """Standalone NEAREST identity on the targets whose engine supports LATERAL."""

    def test_standalone_nearest_k1_agrees_on_duckdb_targets(self, cross_target_oracle):
        """Test standalone NEAREST k=1 agrees across the LATERAL-capable targets.

        Given:
            A single-row peaks table and three candidate genes at varying
            distances on chr1.
        When:
            A correlated ``CROSS JOIN LATERAL NEAREST(..., k := 1)`` query runs
            on the generic and duckdb targets, both executed on DuckDB.
        Then:
            Both targets should return the single nearest gene and agree.

        The generic target is routed to DuckDB here (not its default DataFusion
        engine): the correlated LATERAL this expansion emits has no DataFusion
        physical plan — the one documented cross-target gap.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.chrom, a.start AS a_start, b.start AS b_start "
            "FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b",
            peaks=[("chr1", 200, 300)],
            genes=[
                ("chr1", 1000, 1100),
                ("chr1", 50, 60),
                ("chr1", 280, 290),
            ],
            expected=[("chr1", 200, 280)],
            targets=("generic", "duckdb"),
            engines={"generic": "duckdb"},
        )

    def test_nearest_on_datafusion_lateral_is_unsupported(self, cross_target_oracle):
        """Test the NEAREST LATERAL expansion has no DataFusion physical plan.

        Given:
            The same single-row NEAREST query.
        When:
            The oracle is asked to run only the datafusion target on DataFusion.
        Then:
            DataFusion should reject the correlated LATERAL, pinning the gap that
            justifies excluding it from the identity assertion above.
        """
        # Arrange / Act
        with pytest.raises(Exception) as excinfo:
            cross_target_oracle(
                "SELECT a.chrom, a.start AS a_start, b.start AS b_start "
                "FROM peaks a "
                "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b",
                peaks=[("chr1", 200, 300)],
                genes=[("chr1", 280, 290)],
                expected=[("chr1", 200, 280)],
                targets=("datafusion",),
            )

        # Assert
        assert "not implemented" in str(excinfo.value).lower()
