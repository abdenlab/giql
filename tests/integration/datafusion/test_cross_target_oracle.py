"""Cross-target result-identity tests via the ``cross_target_oracle`` fixture.

These exercise the reusable oracle (``tests/integration/conftest.py``) over the
operators that already emit identical generic SQL across Generic and DataFusion
and run correctly on DuckDB: INTERSECTS (literal + column-to-column join),
CONTAINS, WITHIN, and standalone NEAREST. No operator has been migrated to the
expander registry yet (epic #137), so this lane locks in the verification path
every later migration (#140-#144) will consume.

NEAREST's expansion uses a correlated ``LATERAL`` subquery, which DataFusion has
no physical plan for today; its generic-vs-duckdb equivalence case runs both on
DuckDB, and the full three-target oracle is pinned as a strict xfail (#142) that
auto-promotes when DataFusion gains correlated LATERAL — the one documented
cross-target gap.
"""

import pytest

pytest.importorskip("duckdb")
pytest.importorskip("datafusion")
pytest.importorskip("pyarrow")

from giql import Table  # noqa: E402

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

    def test_standalone_nearest_k1_agrees_generic_vs_duckdb_on_duckdb(
        self, cross_target_oracle
    ):
        """Test NEAREST k=1 generic and duckdb SQL agree when both run on DuckDB.

        Given:
            A single-row peaks table and three candidate genes at varying
            distances on chr1.
        When:
            A correlated ``CROSS JOIN LATERAL NEAREST(..., k := 1)`` query runs
            for the generic and duckdb targets, both executed on DuckDB (the
            generic target routed via ``engines={"generic": "duckdb"}``).
        Then:
            The generic and duckdb SQL should both return the single nearest
            gene and agree.

        This is a generic-vs-duckdb *equivalence* check, not a cross-target
        identity check: DataFusion is never run here. The generic target is
        routed to DuckDB because the correlated LATERAL this expansion emits has
        no DataFusion physical plan — the one documented cross-target gap.
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

    @pytest.mark.xfail(
        strict=True,
        raises=Exception,
        reason="DataFusion lacks correlated LATERAL (OuterReferenceColumn) — #142",
    )
    def test_nearest_full_oracle_xfails_on_datafusion_lateral(self, cross_target_oracle):
        """Test the full NEAREST oracle xfails on DataFusion's missing LATERAL.

        Given:
            The single-row NEAREST query and a candidate gene on chr1.
        When:
            The oracle runs all three targets — the datafusion target executes
            the correlated LATERAL on DataFusion, which has no physical plan.
        Then:
            DataFusion should raise its ``OuterReferenceColumn`` "not
            implemented" error, pinning the gap as a strict xfail that
            auto-promotes (XPASS -> fail) when #142 lands. An UNRELATED
            DataFusion error must still fail loudly, so the match is narrowed
            to the LATERAL signature before re-raising.
        """
        # Arrange / Act
        try:
            cross_target_oracle(
                "SELECT a.chrom, a.start AS a_start, b.start AS b_start "
                "FROM peaks a "
                "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b",
                peaks=[("chr1", 200, 300)],
                genes=[("chr1", 280, 290)],
                expected=[("chr1", 200, 280)],
            )
        except Exception as exc:  # noqa: BLE001
            message = str(exc).lower()
            # Assert: narrow the xfail to the documented LATERAL gap so any
            # unrelated DataFusion failure escapes and fails the test loudly.
            assert "outerreferencecolumn" in message or "not implemented" in message, (
                f"unexpected DataFusion error, not the LATERAL gap: {exc!r}"
            )
            raise


class TestCrossTargetOracleIntersectsAnyAll:
    """INTERSECTS ANY/ALL identity across all targets (T2)."""

    def test_intersects_any_returns_rows_matching_either_range(
        self, cross_target_oracle
    ):
        """Test INTERSECTS ANY returns rows overlapping at least one range.

        Given:
            Peaks overlapping the first range, the second range, neither, and a
            different chromosome.
        When:
            An INTERSECTS ANY query over two ranges runs on every target.
        Then:
            Every target should return exactly the two peaks overlapping a
            listed range, in agreement.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            'SELECT chrom, start, "end" FROM peaks '
            "WHERE interval INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')",
            peaks=[
                ("chr1", 1500, 1800),
                ("chr1", 5500, 5600),
                ("chr1", 3000, 3100),
                ("chr2", 1500, 1800),
            ],
            expected=[("chr1", 1500, 1800), ("chr1", 5500, 5600)],
        )

    def test_intersects_any_no_match_returns_zero_rows(self, cross_target_oracle):
        """Test INTERSECTS ANY with no overlapping peak returns zero rows.

        Given:
            A peak overlapping neither listed range.
        When:
            The INTERSECTS ANY query runs on every target.
        Then:
            Every target should return zero rows in agreement.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            'SELECT chrom, start, "end" FROM peaks '
            "WHERE interval INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')",
            peaks=[("chr1", 8000, 8100)],
            expected=[],
        )

    def test_intersects_all_returns_rows_matching_every_range(self, cross_target_oracle):
        """Test INTERSECTS ALL returns rows overlapping every listed range.

        Given:
            One peak overlapping both ranges and one overlapping neither.
        When:
            An INTERSECTS ALL query over two ranges runs on every target.
        Then:
            Every target should return only the peak overlapping both ranges.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            'SELECT chrom, start, "end" FROM peaks '
            "WHERE interval INTERSECTS ALL('chr1:1000-2000', 'chr1:1500-1700')",
            peaks=[("chr1", 1400, 1800), ("chr1", 100, 200)],
            expected=[("chr1", 1400, 1800)],
        )

    def test_intersects_all_no_full_match_returns_zero_rows(self, cross_target_oracle):
        """Test INTERSECTS ALL returns zero rows when no peak overlaps all ranges.

        Given:
            Peaks each overlapping only one of the two ranges.
        When:
            The INTERSECTS ALL query runs on every target.
        Then:
            Every target should return zero rows in agreement.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            'SELECT chrom, start, "end" FROM peaks '
            "WHERE interval INTERSECTS ALL('chr1:1000-1200', 'chr1:5000-6000')",
            peaks=[("chr1", 1050, 1100), ("chr1", 5050, 5100)],
            expected=[],
        )


class TestCrossTargetOracleDistance:
    """DISTANCE column-to-column identity across all targets (T2)."""

    def test_distance_column_to_column_filters_near_pairs(self, cross_target_oracle):
        """Test a column-to-column DISTANCE predicate agrees across targets.

        Given:
            A peak and two genes: one within the distance threshold and one far
            beyond it on the same chromosome.
        When:
            A DISTANCE(a, b) < threshold join runs on every target.
        Then:
            Every target should return only the near pair.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.start AS a_start, b.start AS b_start "
            "FROM peaks a JOIN genes b ON a.chrom = b.chrom "
            "WHERE DISTANCE(a.interval, b.interval) < 100",
            peaks=[("chr1", 100, 200)],
            genes=[("chr1", 250, 300), ("chr1", 5000, 5100)],
            expected=[(100, 250)],
        )

    def test_distance_no_pair_within_threshold_returns_zero_rows(
        self, cross_target_oracle
    ):
        """Test a column-to-column DISTANCE predicate with no near pair is empty.

        Given:
            A peak and a single gene beyond the distance threshold.
        When:
            The DISTANCE join runs on every target.
        Then:
            Every target should return zero rows in agreement.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.start AS a_start, b.start AS b_start "
            "FROM peaks a JOIN genes b ON a.chrom = b.chrom "
            "WHERE DISTANCE(a.interval, b.interval) < 5",
            peaks=[("chr1", 100, 200)],
            genes=[("chr1", 5000, 5100)],
            expected=[],
        )


class TestCrossTargetOracleCluster:
    """CLUSTER identity across all targets (T2)."""

    def test_cluster_assigns_shared_ids_to_overlapping_runs(self, cross_target_oracle):
        """Test CLUSTER assigns identical run ids across targets.

        Given:
            Two overlapping intervals and one isolated interval on chr1.
        When:
            A CLUSTER(interval) projection runs on every target.
        Then:
            Every target should assign cluster id 1 to the overlapping run and 2
            to the isolated interval, in agreement.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            'SELECT chrom, start, "end", CLUSTER(interval) AS cid FROM peaks',
            peaks=[
                ("chr1", 100, 200),
                ("chr1", 150, 300),
                ("chr1", 5000, 6000),
            ],
            expected=[
                ("chr1", 100, 200, 1),
                ("chr1", 150, 300, 1),
                ("chr1", 5000, 6000, 2),
            ],
        )

    def test_cluster_empty_input_returns_zero_rows(self, cross_target_oracle):
        """Test CLUSTER over an empty table returns zero rows on every target.

        Given:
            An empty peaks table.
        When:
            The CLUSTER projection runs on every target (DataFusion via the
            synthesized empty record batch).
        Then:
            Every target should return zero rows in agreement.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            'SELECT chrom, start, "end", CLUSTER(interval) AS cid FROM peaks',
            peaks=[],
            expected=[],
        )


class TestCrossTargetOracleMerge:
    """MERGE identity across all targets (T2)."""

    def test_merge_collapses_overlapping_intervals(self, cross_target_oracle):
        """Test MERGE collapses overlapping intervals identically across targets.

        Given:
            Two overlapping intervals and one isolated interval on chr1.
        When:
            A MERGE(interval) query runs on every target.
        Then:
            Every target should return the merged span and the isolated
            interval, in agreement.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT MERGE(interval) FROM peaks",
            tables=[Table("peaks")],
            peaks=[
                ("chr1", 100, 200),
                ("chr1", 150, 300),
                ("chr1", 5000, 6000),
            ],
            expected=[("chr1", 100, 300), ("chr1", 5000, 6000)],
        )

    def test_merge_empty_input_returns_zero_rows(self, cross_target_oracle):
        """Test MERGE over an empty table returns zero rows on every target.

        Given:
            An empty peaks table.
        When:
            The MERGE query runs on every target (DataFusion via the synthesized
            empty record batch).
        Then:
            Every target should return zero rows in agreement.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT MERGE(interval) FROM peaks",
            tables=[Table("peaks")],
            peaks=[],
            expected=[],
        )


class TestCrossTargetOracleDisjoin:
    """DISJOIN behaviour: a DuckDB-only case plus the DataFusion duplicate-end gap."""

    def test_disjoin_splits_overlaps_on_duckdb(self, cross_target_oracle):
        """Test DISJOIN splits overlapping intervals at breakpoints on DuckDB.

        Given:
            Two overlapping intervals on chr1.
        When:
            A DISJOIN query runs on the duckdb target only.
        Then:
            DuckDB should return each original interval paired with the
            sub-segments cut at every breakpoint.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            'SELECT chrom, start, "end", disjoin_start, disjoin_end FROM DISJOIN(peaks)',
            tables=[Table("peaks")],
            peaks=[("chr1", 0, 100), ("chr1", 50, 150)],
            expected=[
                ("chr1", 0, 100, 0, 50),
                ("chr1", 0, 100, 50, 100),
                ("chr1", 50, 150, 50, 100),
                ("chr1", 50, 150, 100, 150),
            ],
            targets=("duckdb",),
        )

    @pytest.mark.xfail(
        strict=True,
        raises=Exception,
        reason='DISJOIN CTE * passthrough emits duplicate t."end" columns — #145',
    )
    def test_disjoin_full_oracle_xfails_on_datafusion_duplicate_end(
        self, cross_target_oracle
    ):
        """Test the full DISJOIN oracle xfails on DataFusion's duplicate-end names.

        Given:
            The same two overlapping intervals.
        When:
            The oracle runs all three targets — the datafusion target executes
            DISJOIN's CTE ``*`` passthrough, which projects ``t."end"`` twice.
        Then:
            DataFusion should reject the projection for non-unique expression
            names (DuckDB tolerates it). This is a strict xfail that promotes
            when DISJOIN's passthrough is de-duplicated; an unrelated DataFusion
            error escapes and fails loudly. Note: the cause is the duplicate
            output-column name, NOT canonicalization or ``* REPLACE`` (no
            ``* REPLACE`` is emitted today).
        """
        # Arrange / Act
        try:
            cross_target_oracle(
                'SELECT chrom, start, "end", disjoin_start, disjoin_end '
                "FROM DISJOIN(peaks)",
                tables=[Table("peaks")],
                peaks=[("chr1", 0, 100), ("chr1", 50, 150)],
                expected=[
                    ("chr1", 0, 100, 0, 50),
                    ("chr1", 0, 100, 50, 100),
                    ("chr1", 50, 150, 50, 100),
                    ("chr1", 50, 150, 100, 150),
                ],
            )
        except Exception as exc:  # noqa: BLE001
            message = str(exc).lower()
            # Assert: narrow to the duplicate-column-name signature so unrelated
            # DataFusion failures escape and fail the test loudly.
            assert "unique" in message and "name" in message, (
                f"unexpected DataFusion error, not the duplicate-end gap: {exc!r}"
            )
            raise


class TestCrossTargetOracleDataShapes:
    """Data-shape coverage for the oracle (T3)."""

    def test_multi_chrom_join_does_not_cross_chromosomes(self, cross_target_oracle):
        """Test a multi-chromosome join keeps overlaps within each chromosome.

        Given:
            Peaks on chr1 and chr2 and genes on chr1, chr2, and chr3 where only
            same-chromosome pairs overlap by coordinate.
        When:
            A column-to-column INTERSECTS join runs on every target.
        Then:
            Every target should return only the same-chromosome overlaps and
            never a cross-chromosome pair.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.chrom AS chrom, a.start AS a_start, b.start AS b_start "
            "FROM peaks a JOIN genes b ON a.interval INTERSECTS b.interval",
            peaks=[("chr1", 100, 500), ("chr2", 100, 500)],
            genes=[("chr1", 300, 400), ("chr2", 300, 400), ("chr3", 300, 400)],
            expected=[("chr1", 100, 300), ("chr2", 100, 300)],
        )

    def test_join_path_with_no_overlap_returns_zero_rows(self, cross_target_oracle):
        """Test the JOIN path returns zero rows when nothing overlaps.

        Given:
            A peak and a gene that do not overlap.
        When:
            A column-to-column INTERSECTS join runs on every target.
        Then:
            Every target should return zero rows in agreement.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.start AS a_start, b.start AS b_start "
            "FROM peaks a JOIN genes b ON a.interval INTERSECTS b.interval",
            peaks=[("chr1", 100, 200)],
            genes=[("chr1", 5000, 6000)],
            expected=[],
        )

    def test_empty_input_table_returns_zero_rows(self, cross_target_oracle):
        """Test an empty input table yields zero rows on every target.

        Given:
            An empty peaks table.
        When:
            A literal INTERSECTS query runs on every target — DataFusion
            registers the empty table via the synthesized empty record batch.
        Then:
            Every target should return zero rows in agreement.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            'SELECT chrom, start, "end" FROM peaks '
            "WHERE interval INTERSECTS 'chr1:1000-2000'",
            peaks=[],
            expected=[],
        )

    def test_half_open_touching_intervals_do_not_overlap(self, cross_target_oracle):
        """Test half-open adjacent intervals (end == start) do not overlap.

        Given:
            A peak ending exactly where a gene begins (half-open adjacency).
        When:
            A column-to-column INTERSECTS join runs on every target.
        Then:
            Every target should return zero rows — touching is not overlap under
            half-open semantics.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.start AS a_start, b.start AS b_start "
            "FROM peaks a JOIN genes b ON a.interval INTERSECTS b.interval",
            peaks=[("chr1", 100, 200)],
            genes=[("chr1", 200, 300)],
            expected=[],
        )

    def test_one_base_pair_overlap_matches(self, cross_target_oracle):
        """Test a single-base-pair overlap is counted as an overlap.

        Given:
            A peak extending one base past a gene's start.
        When:
            A column-to-column INTERSECTS join runs on every target.
        Then:
            Every target should return the single overlapping pair.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.start AS a_start, b.start AS b_start "
            "FROM peaks a JOIN genes b ON a.interval INTERSECTS b.interval",
            peaks=[("chr1", 100, 201)],
            genes=[("chr1", 200, 300)],
            expected=[(100, 200)],
        )

    def test_large_interval_spanning_bins_is_not_duplicated(self, cross_target_oracle):
        """Test a multi-bin interval yields one pair, not duplicates.

        Given:
            A peak spanning many join bins and a single gene inside it.
        When:
            A column-to-column INTERSECTS join runs on every target.
        Then:
            Every target should return exactly one pair — the binned join must
            dedup the cross-bin candidates, the oracle's multiset compare being
            what would catch a duplicate-row regression.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.start AS a_start, b.start AS b_start "
            "FROM peaks a JOIN genes b ON a.interval INTERSECTS b.interval",
            peaks=[("chr1", 0, 500000)],
            genes=[("chr1", 250000, 250100)],
            expected=[(0, 250000)],
        )
