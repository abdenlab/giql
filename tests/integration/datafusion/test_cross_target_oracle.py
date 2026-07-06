"""Cross-target result-identity tests via the ``cross_target_oracle`` fixture.

These exercise the reusable oracle (``tests/integration/conftest.py``) over the
operators that already emit identical generic SQL across Generic and DataFusion
and run correctly on DuckDB: INTERSECTS (literal + column-to-column join),
CONTAINS, WITHIN, and standalone NEAREST. DISTANCE, the spatial predicates,
NEAREST, and DISJOIN have since been migrated to the expander registry (epic
#137); this lane locks in the verification path each migration consumes.

For the non-join operators (DISTANCE, CONTAINS, WITHIN, ANY/ALL, CLUSTER,
MERGE) the generic and datafusion targets emit byte-identical SQL and both run
on the DataFusion engine, so the load-bearing comparison there is
DataFusion-vs-DuckDB. Only the column-to-column INTERSECTS join produces
genuinely divergent SQL across targets (the DuckDB IEJoin vs. the binned
equi-join).

NEAREST's correlated expansion is capability-driven (#142): lateral-capable
targets (generic, duckdb) emit the portable ``LATERAL`` subquery, while
DataFusion — which has no correlated-LATERAL physical plan — gets a decorrelated
window-function fallback. For a single correlated NEAREST per query the two forms
return the same rows on every projection shape, so the full three-target identity
oracle runs on every target — including ``SELECT *`` / ``SELECT b.*``: the
fallback's reserved ``__giql_x_*`` rank/key columns, which its decorrelated join
must expose, are projected away by a statement finalizer that wraps the enclosing
SELECT in ``SELECT * EXCEPT (...)`` on the DataFusion path (#160), so no reserved
column leaks. The wrapper re-locates its target join by a reserved ``meta`` marker
and mints its reserved column names per fallback (#172), so the former residual
compositions no longer leak either: two correlated NEAREST fallbacks in one query
(distinct per-run names, each ``* EXCEPT`` independently) and a correlated NEAREST
re-surfaced by an enclosing ``SELECT *`` outside its own SELECT — e.g. a wrapping
CLUSTER, whose ``copy()`` + ``transplant`` the marker survives — are both covered
below. DISJOIN is migrated onto the expander registry (#143) and its previously
pending-#153 gap (duplicate ``end`` output names) is closed, so its full
three-target oracle now runs as a real cross-target identity test.
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

    def test_correlated_nearest_k1_agrees_across_all_targets(self, cross_target_oracle):
        """Test correlated NEAREST k=1 returns identical rows on every target.

        Given:
            A single-row peaks table and three candidate genes at varying
            distances on chr1.
        When:
            A correlated ``CROSS JOIN LATERAL NEAREST(..., k := 1)`` query runs
            for every target — the generic and duckdb targets emit the portable
            LATERAL form (executed on DuckDB, the lateral-capable engine), and
            the datafusion target emits the decorrelated window-function fallback
            the #142 expander produces (executed on DataFusion).
        Then:
            Every target should return the single nearest gene and agree.

        Promoted from the ``_unsupported_pending_142`` expected-failure pin:
        DataFusion now plans correlated NEAREST through the capability-driven
        window-function fallback, so the full three-target oracle is a real
        identity test rather than a ``pytest.raises`` placeholder. The generic
        target is routed to DuckDB because its portable SQL is the LATERAL form,
        which only the datafusion-specific fallback decorrelates for DataFusion.
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
            engines={"generic": "duckdb"},
        )

    def test_correlated_nearest_k2_returns_two_nearest_across_targets(
        self, cross_target_oracle
    ):
        """Test correlated NEAREST k=2 picks the two nearest on every target.

        Given:
            One peak and four candidate genes, more than k of them on the peak's
            chromosome at distinct distances.
        When:
            A correlated ``NEAREST(..., k := 2)`` runs on every target — DuckDB
            via the LATERAL form, DataFusion via the decorrelated window fallback.
        Then:
            Every target should return the two nearest genes and agree, pinning
            the top-k fan-out of the fallback against the LATERAL form.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.start AS a_start, b.start AS b_start "
            "FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 2) b",
            peaks=[("chr1", 200, 300)],
            genes=[
                ("chr1", 1000, 1100),
                ("chr1", 50, 60),
                ("chr1", 280, 290),
                ("chr1", 310, 320),
            ],
            expected=[(200, 280), (200, 310)],
            engines={"generic": "duckdb"},
        )

    def test_correlated_nearest_duplicate_reference_rows_fan_out(
        self, cross_target_oracle
    ):
        """Test correlated NEAREST fans the top-k out to duplicate reference rows.

        Given:
            Two identical peak rows and two candidate genes.
        When:
            A correlated ``NEAREST(..., k := 1)`` runs on every target.
        Then:
            Every target should return the nearest gene once per duplicate peak
            (two rows), pinning the fallback's DISTINCT-then-rejoin fan-out so a
            duplicate outer row is not collapsed.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.start AS a_start, b.start AS b_start "
            "FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b",
            peaks=[("chr1", 200, 300), ("chr1", 200, 300)],
            genes=[("chr1", 280, 290), ("chr1", 50, 60)],
            expected=[(200, 280), (200, 280)],
            engines={"generic": "duckdb"},
        )

    def test_correlated_nearest_partitions_by_chromosome(self, cross_target_oracle):
        """Test correlated NEAREST keys the nearest per outer chromosome.

        Given:
            Peaks on chr1 and chr2 and candidate genes on both chromosomes.
        When:
            A correlated ``NEAREST(..., k := 1)`` runs on every target.
        Then:
            Each peak should pair with the nearest gene on its own chromosome and
            all targets agree, pinning the fallback's PARTITION BY reference key
            across distinct outer keys.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.chrom AS chrom, a.start AS a_start, b.start AS b_start "
            "FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b",
            peaks=[("chr1", 200, 300), ("chr2", 200, 300)],
            genes=[("chr1", 280, 290), ("chr2", 500, 510), ("chr2", 205, 215)],
            expected=[("chr1", 200, 280), ("chr2", 200, 205)],
            engines={"generic": "duckdb"},
        )

    def test_correlated_nearest_max_distance_boundary(self, cross_target_oracle):
        """Test correlated NEAREST drops candidates beyond max_distance everywhere.

        Given:
            A peak and two genes, one just inside and one far beyond a
            ``max_distance`` threshold.
        When:
            A correlated ``NEAREST(..., k := 5, max_distance := 100)`` runs on
            every target.
        Then:
            Every target should return only the in-threshold gene, pinning the
            ``max_distance`` filter through both the LATERAL and fallback forms.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.start AS a_start, b.start AS b_start "
            "FROM peaks a "
            "CROSS JOIN LATERAL NEAREST("
            "genes, reference := a.interval, k := 5, max_distance := 100) b",
            peaks=[("chr1", 200, 300)],
            genes=[("chr1", 360, 400), ("chr1", 5000, 5100)],
            expected=[(200, 360)],
            engines={"generic": "duckdb"},
        )

    def test_correlated_nearest_stranded_matches_strand(self, cross_target_oracle):
        """Test stranded correlated NEAREST matches strand on every target.

        Given:
            A ``+`` peak and two genes — a slightly farther ``+`` gene and a
            nearer ``-`` gene.
        When:
            A correlated ``NEAREST(..., k := 1, stranded := true)`` runs on every
            target.
        Then:
            Every target should return the same-strand (``+``) gene even though
            the opposite-strand gene is nearer, in agreement.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.start AS a_start, b.start AS b_start "
            "FROM peaks a "
            "CROSS JOIN LATERAL NEAREST("
            "genes, reference := a.interval, k := 1, stranded := true) b",
            tables=[Table("peaks"), Table("genes")],
            columns=_STRANDED_COLUMNS,
            peaks=[("chr1", 200, 300, "+")],
            genes=[("chr1", 280, 290, "+"), ("chr1", 250, 260, "-")],
            expected=[(200, 280)],
            engines={"generic": "duckdb"},
        )

    def test_correlated_nearest_signed_distance_agrees(self, cross_target_oracle):
        """Test signed correlated NEAREST reports signed distances everywhere.

        Given:
            A peak with one upstream and one downstream candidate gene.
        When:
            A correlated ``NEAREST(..., k := 2, signed := true)`` projects the
            ``distance`` column on every target.
        Then:
            Every target should report a negative distance for the upstream gene
            and a positive one for the downstream gene, in agreement.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.start AS a_start, b.start AS b_start, b.distance AS d "
            "FROM peaks a "
            "CROSS JOIN LATERAL NEAREST("
            "genes, reference := a.interval, k := 2, signed := true) b",
            peaks=[("chr1", 200, 300)],
            genes=[("chr1", 50, 60), ("chr1", 360, 400)],
            expected=[(200, 50, -141), (200, 360, 61)],
            engines={"generic": "duckdb"},
        )

    def test_correlated_nearest_k_th_distance_tie_breaks_on_coordinates(
        self, cross_target_oracle
    ):
        """Test the (start, end) tiebreaker picks the same k-th candidate everywhere.

        Given:
            One peak and three genes where two candidates are tied at the k-th
            (k=1) distance — both 100 bp away, one upstream and one downstream of
            the peak — so only the ``(start, end)`` tiebreaker can choose between
            them (the lower ``(start, end)`` wins).
        When:
            A correlated ``NEAREST(..., k := 1)`` runs on every target — DuckDB via
            the LATERAL form's ``ORDER BY ABS(distance), start, end LIMIT 1`` and
            DataFusion via the fallback's matching ``ROW_NUMBER()`` ordering.
        Then:
            Every target should return the lower-coordinate tied candidate (the
            upstream gene), so the LATERAL and window forms agree on the tie rather
            than ordering it by engine-dependent chance.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.start AS a_start, b.start AS b_start "
            "FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b",
            peaks=[("chr1", 200, 300)],
            # Upstream gene ends at 100 (gap 100); downstream gene starts at 400
            # (gap 100). Both tie at distance 100; (start, end) breaks the tie in
            # favor of the upstream gene (start 50 < start 400).
            genes=[("chr1", 50, 100), ("chr1", 400, 450)],
            expected=[(200, 50)],
            engines={"generic": "duckdb"},
        )

    def test_correlated_nearest_stranded_opposite_strands_same_position(
        self, cross_target_oracle
    ):
        """Test stranded NEAREST keys per-strand for co-located opposite-strand rows.

        Given:
            Two peaks at the *same* position but on opposite strands (``+`` and
            ``-``), and one same-position ``+`` gene plus one ``-`` gene, so the
            strand-augmented reference key must keep each outer row's nearest
            strand-matched independently.
        When:
            A correlated ``NEAREST(..., k := 1, stranded := true)`` runs on every
            target — DuckDB via the LATERAL form, DataFusion via the fallback whose
            reference key includes strand.
        Then:
            The ``+`` peak should pair with the ``+`` gene and the ``-`` peak with
            the ``-`` gene, in agreement, proving the fan-out keys by strand so two
            co-located opposite-strand outer rows are not collapsed.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.strand AS a_strand, b.start AS b_start "
            "FROM peaks a "
            "CROSS JOIN LATERAL NEAREST("
            "genes, reference := a.interval, k := 1, stranded := true) b",
            tables=[Table("peaks"), Table("genes")],
            columns=_STRANDED_COLUMNS,
            peaks=[("chr1", 200, 300, "+"), ("chr1", 200, 300, "-")],
            genes=[("chr1", 280, 290, "+"), ("chr1", 250, 260, "-")],
            expected=[("+", 280), ("-", 250)],
            engines={"generic": "duckdb"},
        )

    def test_correlated_nearest_max_distance_keeps_k_survivors(
        self, cross_target_oracle
    ):
        """Test max_distance with k>1 keeps every in-threshold survivor everywhere.

        Given:
            One peak and four genes where three sit within a ``max_distance``
            threshold at distinct distances and one sits beyond it, with k larger
            than the survivor count.
        When:
            A correlated ``NEAREST(..., k := 3, max_distance := 200)`` runs on
            every target — DuckDB via the LATERAL form, DataFusion via the window
            fallback.
        Then:
            Every target should return exactly the three in-threshold genes (the
            beyond-threshold gene dropped), proving ``max_distance`` and the top-k
            survivor set interact identically across both forms.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.start AS a_start, b.start AS b_start "
            "FROM peaks a "
            "CROSS JOIN LATERAL NEAREST("
            "genes, reference := a.interval, k := 3, max_distance := 200) b",
            peaks=[("chr1", 200, 300)],
            genes=[
                ("chr1", 350, 360),
                ("chr1", 420, 430),
                ("chr1", 480, 490),
                ("chr1", 5000, 5100),
            ],
            expected=[(200, 350), (200, 420), (200, 480)],
            engines={"generic": "duckdb"},
        )

    def test_correlated_nearest_unaliased_lateral_agrees_across_targets(
        self, cross_target_oracle
    ):
        """Test an unaliased correlated NEAREST agrees across targets (B3 on-engine).

        Given:
            A correlated ``CROSS JOIN LATERAL NEAREST(...)`` written *without* a
            table alias — legitimate GIQL that, before B3, raised on DataFusion
            while running on DuckDB.
        When:
            The query runs on every target — DuckDB via the LATERAL form and
            DataFusion via the decorrelated fallback, which now synthesizes the
            missing alias instead of asserting one.
        Then:
            Every target should return the single nearest gene and agree, proving
            the synthesized-alias fallback both transpiles and executes on the real
            DataFusion engine.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.start AS a_start "
            "FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1)",
            peaks=[("chr1", 200, 300)],
            genes=[
                ("chr1", 1000, 1100),
                ("chr1", 50, 60),
                ("chr1", 280, 290),
            ],
            expected=[(200,)],
            engines={"generic": "duckdb"},
        )

    def test_correlated_nearest_star_projection_agrees_across_targets(
        self, cross_target_oracle
    ):
        """Test SELECT b.* over a correlated NEAREST agrees on every target (#160).

        Given:
            A single-row peak and three candidate genes on chr1.
        When:
            A correlated ``NEAREST(..., k := 1)`` query projects ``b.*`` on every
            target — DuckDB emits the LATERAL form (``genes.* + distance``) and
            DataFusion's fallback wraps its output in ``SELECT * EXCEPT (...)`` to
            hide the reserved ``__giql_x_<n>_rk_*`` / ``__giql_x_<n>_rn`` columns its
            decorrelated join must expose.
        Then:
            Every target should return the single nearest gene with the same
            columns (``genes.* + distance``) and agree — the reserved columns no
            longer leak on DataFusion.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT b.* "
            "FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b",
            peaks=[("chr1", 200, 300)],
            genes=[
                ("chr1", 1000, 1100),
                ("chr1", 50, 60),
                ("chr1", 280, 290),
            ],
            expected=[("chr1", 280, 290, 0)],
            engines={"generic": "duckdb"},
        )

    def test_correlated_nearest_unqualified_star_agrees_across_targets(
        self, cross_target_oracle
    ):
        """Test unqualified SELECT * over a correlated NEAREST agrees per target.

        Given:
            A single-row peak and three candidate genes on chr1.
        When:
            A correlated ``NEAREST(..., k := 1)`` query projects unqualified
            ``SELECT *`` on every target — pulling both the outer peak row and the
            nearest gene; DataFusion's fallback wraps the statement in
            ``SELECT * EXCEPT (...)`` to strip the reserved join columns (#160).
        Then:
            Every target should return the outer peak columns plus the nearest
            gene's ``genes.* + distance`` and agree, with no reserved columns
            leaking on DataFusion.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT * "
            "FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b",
            peaks=[("chr1", 200, 300)],
            genes=[
                ("chr1", 1000, 1100),
                ("chr1", 50, 60),
                ("chr1", 280, 290),
            ],
            expected=[("chr1", 200, 300, "chr1", 280, 290, 0)],
            engines={"generic": "duckdb"},
        )

    def test_correlated_nearest_stranded_star_projection_agrees_across_targets(
        self, cross_target_oracle
    ):
        """Test stranded SELECT b.* correlated NEAREST agrees per target (#160).

        Given:
            A ``+`` peak and two candidate genes — a same-strand ``+`` gene and a
            nearer opposite-strand ``-`` gene — on the stranded schema.
        When:
            A stranded correlated ``NEAREST(..., k := 1, stranded := true)`` projects
            ``b.*`` on every target — DataFusion's fallback additionally exposes the
            reserved ``__giql_x_<n>_rk_strand`` key column, which its ``* EXCEPT``
            wrapper must strip.
        Then:
            Every target should return the same-strand nearest gene with identical
            columns (``genes.* + distance``) and agree — no reserved column, strand
            key included, leaks on DataFusion.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT b.* "
            "FROM peaks a "
            "CROSS JOIN LATERAL NEAREST("
            "genes, reference := a.interval, k := 1, stranded := true) b",
            tables=[Table("peaks"), Table("genes")],
            columns=_STRANDED_COLUMNS,
            peaks=[("chr1", 200, 300, "+")],
            genes=[("chr1", 280, 290, "+"), ("chr1", 250, 260, "-")],
            expected=[("chr1", 280, 290, "+", 0)],
            engines={"generic": "duckdb"},
        )

    def test_correlated_nearest_star_projection_k_gt_1_agrees_across_targets(
        self, cross_target_oracle
    ):
        """Test SELECT b.* over a k>1 correlated NEAREST agrees per target (#160).

        Given:
            One peak and four candidate genes at distinct distances on chr1.
        When:
            A correlated ``NEAREST(..., k := 3)`` projects ``b.*`` on every target —
            DataFusion's fallback wraps its multi-row output in ``SELECT * EXCEPT
            (...)``.
        Then:
            Every target should return the three nearest genes with identical columns
            and agree — the reserved columns are stripped and the multi-row top-k
            fan-out survives the wrapper.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT b.* "
            "FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 3) b",
            peaks=[("chr1", 200, 300)],
            genes=[
                ("chr1", 280, 290),
                ("chr1", 350, 360),
                ("chr1", 420, 430),
                ("chr1", 5000, 5100),
            ],
            expected=[
                ("chr1", 280, 290, 0),
                ("chr1", 350, 360, 51),
                ("chr1", 420, 430, 121),
            ],
            engines={"generic": "duckdb"},
        )

    def test_cluster_over_correlated_nearest_star_agrees_across_targets(
        self, cross_target_oracle
    ):
        """Test SELECT * CLUSTER over a correlated NEAREST agrees per target (#172).

        Given:
            Two peaks whose nearest genes fall in two well-separated clusters, with
            the correlated ``NEAREST(..., k := 1)`` nested inside a subquery that a
            ``SELECT *``-projecting ``CLUSTER`` wraps.
        When:
            The query runs on every target — on DataFusion the CLUSTER's
            copy+transplant relocates the fallback join, and the finalizer must
            re-locate it by its reserved marker to wrap away the reserved columns.
        Then:
            Every target should return each nearest gene (``genes.* + distance``),
            CLUSTER's own ``__giql_is_new_cluster`` helper, and the cluster id, and
            agree. No reserved ``__giql_x_*`` NEAREST column leaks through the
            enclosing ``SELECT *`` (formerly a #172 residual) — a leak would surface
            four extra columns and break the cross-target column-count agreement.
        """
        # Arrange / Act / Assert
        # The 5th column is CLUSTER's own ``__giql_is_new_cluster`` flag, which a
        # ``SELECT *, CLUSTER(...)`` surfaces on every target — a separate known
        # leak family (#161-related), orthogonal to #172. The #172 point is that no
        # ``__giql_x_*`` reserved NEAREST column joins it and the targets agree.
        cross_target_oracle(
            "SELECT *, CLUSTER(interval) AS cid FROM ("
            "SELECT b.* FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b"
            ") sub",
            peaks=[("chr1", 200, 300), ("chr1", 1000, 1100)],
            genes=[("chr1", 280, 290), ("chr1", 1050, 1060)],
            expected=[
                ("chr1", 280, 290, 0, 1, 1),
                ("chr1", 1050, 1060, 0, 1, 2),
            ],
            engines={"generic": "duckdb"},
        )

    def test_two_correlated_nearest_star_agrees_across_targets(
        self, cross_target_oracle
    ):
        """Test two correlated NEAREST projected via star agree per target (#172).

        Given:
            One peak and two candidate genes, with two correlated ``NEAREST(...,
            k := 1)`` fallbacks (``b`` and ``c``) both projected via star
            (``SELECT b.*, c.*``).
        When:
            The query runs on every target — on DataFusion each fallback exposes its
            own reserved rank/key columns, minted under distinct per-fallback tags so
            two nested ``* EXCEPT`` wrappers can each strip their own set.
        Then:
            Every target should return the nearest gene twice over
            (``b.genes.* + distance`` and ``c.genes.* + distance``) and agree — no
            reserved column from either fallback leaks (formerly a #172 residual).
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT b.*, c.* FROM peaks a "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) b "
            "CROSS JOIN LATERAL NEAREST(genes, reference := a.interval, k := 1) c",
            peaks=[("chr1", 200, 300)],
            genes=[("chr1", 280, 290), ("chr1", 50, 60)],
            expected=[("chr1", 280, 290, 0, "chr1", 280, 290, 0)],
            engines={"generic": "duckdb"},
        )


#: A chrom/start/end/strand schema for the stranded NEAREST oracle cases (the
#: default oracle schema carries no strand column).
_STRANDED_COLUMNS = (
    ("chrom", "utf8"),
    ("start", "int64"),
    ("end", "int64"),
    ("strand", "utf8"),
)


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

    def test_distance_upstream_gap_b_precedes_a(self, cross_target_oracle):
        """Test DISTANCE measures the upstream gap when B precedes A.

        Given:
            A peak and a gene that ends before the peak begins on the same
            chromosome, exercising the ``ELSE`` (upstream-gap) CASE branch.
        When:
            A DISTANCE(a, b) < threshold join runs on every target.
        Then:
            Every target should return the pair, the gap being measured from B's
            end to A's start.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.start AS a_start, b.start AS b_start "
            "FROM peaks a JOIN genes b ON a.chrom = b.chrom "
            "WHERE DISTANCE(a.interval, b.interval) < 100",
            peaks=[("chr1", 250, 300)],
            genes=[("chr1", 100, 200)],
            expected=[(250, 100)],
        )

    def test_distance_overlapping_pair_is_zero(self, cross_target_oracle):
        """Test DISTANCE is zero for overlapping intervals.

        Given:
            A peak and a gene that overlap on the same chromosome, exercising the
            overlap (``THEN 0``) CASE branch.
        When:
            A DISTANCE(a, b) < threshold join runs on every target.
        Then:
            Every target should return the pair, the overlap yielding distance 0.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.start AS a_start, b.start AS b_start "
            "FROM peaks a JOIN genes b ON a.chrom = b.chrom "
            "WHERE DISTANCE(a.interval, b.interval) < 100",
            peaks=[("chr1", 100, 300)],
            genes=[("chr1", 200, 400)],
            expected=[(100, 200)],
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

    def test_cluster_over_cte_from_agrees_across_targets(self, cross_target_oracle):
        """Test CLUSTER over a CTE FROM agrees across targets (#174).

        Given:
            Two overlapping intervals and one isolated interval on chr1, wrapped in
            a pass-through CTE that the CLUSTER runs over.
        When:
            A CLUSTER projection over the CTE runs on every target.
        Then:
            Every target should preserve the enclosing WITH — keeping the emitted
            SQL executable — and assign identical cluster ids, in agreement, proving
            the CTE-scoping fix holds on DataFusion as well as DuckDB.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "WITH sub AS (SELECT * FROM peaks) "
            'SELECT chrom, start, "end", CLUSTER(interval) AS cid FROM sub',
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

    def test_cluster_with_limit_agrees_across_targets(self, cross_target_oracle):
        """Test CLUSTER with an outer LIMIT agrees across targets (#181).

        Given:
            Four intervals forming two clusters on chr1, ordered by start with an
            outer LIMIT of two rows.
        When:
            A CLUSTER projection with ORDER BY / LIMIT runs on every target.
        Then:
            Every target should preserve the LIMIT — returning the same two ordered
            rows rather than the whole table — and agree, proving the result-clause
            preservation holds on DataFusion as well as DuckDB.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            'SELECT chrom, start, "end", CLUSTER(interval) AS cid '
            "FROM peaks ORDER BY start LIMIT 2",
            peaks=[
                ("chr1", 10, 20),
                ("chr1", 15, 30),
                ("chr1", 100, 110),
                ("chr1", 105, 120),
            ],
            expected=[
                ("chr1", 10, 20, 1),
                ("chr1", 15, 30, 1),
            ],
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

    def test_merge_over_cte_from_agrees_across_targets(self, cross_target_oracle):
        """Test MERGE over a CTE FROM agrees across targets (#174).

        Given:
            Two overlapping intervals and one isolated interval on chr1, wrapped in
            a pass-through CTE that the MERGE runs over.
        When:
            A MERGE query over the CTE runs on every target.
        Then:
            Every target should preserve the enclosing WITH — keeping the emitted
            SQL executable — and return the merged span and the isolated interval,
            in agreement, proving the CTE-scoping fix holds on DataFusion as well as
            DuckDB.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "WITH sub AS (SELECT * FROM peaks) SELECT MERGE(interval) FROM sub",
            tables=[Table("peaks")],
            peaks=[
                ("chr1", 100, 200),
                ("chr1", 150, 300),
                ("chr1", 5000, 6000),
            ],
            expected=[("chr1", 100, 300), ("chr1", 5000, 6000)],
        )

    def test_merge_with_limit_agrees_across_targets(self, cross_target_oracle):
        """Test MERGE with an outer LIMIT agrees across targets (#181).

        Given:
            Two merged runs on chr1 with an outer LIMIT of one row (MERGE already
            emits ORDER BY chrom, start, so the window is deterministic).
        When:
            A MERGE query with a LIMIT runs on every target.
        Then:
            Every target should preserve the LIMIT — returning only the first merged
            span rather than both — and agree, proving the result-clause preservation
            holds on DataFusion as well as DuckDB.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT MERGE(interval) FROM peaks LIMIT 1",
            tables=[Table("peaks")],
            peaks=[
                ("chr1", 100, 200),
                ("chr1", 150, 300),
                ("chr1", 5000, 6000),
            ],
            expected=[("chr1", 100, 300)],
        )

    def test_merge_with_user_order_by_limit_agrees_across_targets(
        self, cross_target_oracle
    ):
        """Test MERGE honors a user ORDER BY under LIMIT across targets (#181).

        Given:
            Two merged runs on chr1 with a user ORDER BY on the merged end descending
            and an outer LIMIT of one row.
        When:
            A MERGE query with ORDER BY / LIMIT runs on every target.
        Then:
            Every target should order by the user's ORDER BY (not MERGE's fixed
            chrom, start default) and return the largest-end merged span, in
            agreement — proving MERGE honors the user ordering so the preserved LIMIT
            slices the rows the user asked for.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            'SELECT MERGE(interval) FROM peaks ORDER BY "end" DESC LIMIT 1',
            tables=[Table("peaks")],
            peaks=[
                ("chr1", 100, 200),
                ("chr1", 150, 300),
                ("chr1", 5000, 6000),
            ],
            expected=[("chr1", 5000, 6000)],
        )


class TestCrossTargetOracleDisjoin:
    """DISJOIN cross-target identity across Generic, DataFusion, and DuckDB (#143).

    DISJOIN is migrated onto the expander registry (#143), and the duplicate
    unaliased ``t."end"`` columns in the ``__giql_dj_cuts`` CTE UNION branches are
    aliased in every branch (#153), so DataFusion no longer rejects the projection
    for non-unique names. The full three-target oracle therefore runs and agrees.
    """

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

    def test_disjoin_agrees_across_all_targets(self, cross_target_oracle):
        """Test the full DISJOIN oracle returns identical rows on every target.

        Given:
            Two overlapping intervals on chr1.
        When:
            The oracle runs all three targets — generic and datafusion execute the
            registry-expanded DISJOIN (with every ``__giql_dj_cuts`` UNION branch
            aliased, #153) on DataFusion, and duckdb runs on DuckDB.
        Then:
            Every target should return the same sub-segments, proving DISJOIN now
            executes on DataFusion (the #153 duplicate-name gap is closed) and
            agrees with DuckDB.
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
        )

    def test_disjoin_non_canonical_passthrough_agrees_across_targets(
        self, cross_target_oracle
    ):
        """Test the non-canonical EXCEPT passthrough agrees with DuckDB REPLACE.

        Given:
            Two overlapping intervals on chr1 in a table declared **1-based
            closed** (non-canonical), each carrying an extra ``name`` passthrough
            column.
        When:
            The oracle runs all three targets — generic and datafusion exercise
            the portable ``SELECT * EXCEPT (start, end), ...`` passthrough on
            DataFusion (the #143 headline path, which runs on no engine without
            this case), and duckdb runs the in-place ``* REPLACE`` form on DuckDB.
        Then:
            Every target should return the same rows: the de-canonicalized
            (1-based closed) interval columns, the passthrough ``name``, and the
            disjoin sub-segments. The projection lists **explicit** columns so the
            row comparison is order-agnostic across the EXCEPT-re-append vs.
            REPLACE-in-place column orders the two forms produce.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            'SELECT name, start, "end", disjoin_start, disjoin_end FROM DISJOIN(feats)',
            tables=[
                Table(
                    "feats",
                    coordinate_system="1based",
                    interval_type="closed",
                )
            ],
            columns=(
                ("chrom", "utf8"),
                ("start", "int64"),
                ("end", "int64"),
                ("name", "utf8"),
            ),
            feats=[("chr1", 1, 100, "a"), ("chr1", 50, 150, "b")],
            expected=[
                ("a", 1, 100, 1, 49),
                ("a", 1, 100, 50, 100),
                ("b", 50, 150, 50, 100),
                ("b", 50, 150, 101, 150),
            ],
        )


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
        """Test two intervals sharing many bins yield one pair, not duplicates.

        Given:
            A peak (0-500000) and a gene (5000-495000) that co-occupy dozens of
            shared join bins (the bin width is 10000), so the binned candidate
            join produces the same pair once per shared bin.
        When:
            A column-to-column INTERSECTS join runs on every target.
        Then:
            Every target should return exactly one pair — the binned join's
            ``SELECT DISTINCT`` must collapse the duplicate cross-bin candidates,
            and the oracle's multiset compare is what would catch a stripped
            DISTINCT as a duplicate-row regression.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            "SELECT a.start AS a_start, b.start AS b_start "
            "FROM peaks a JOIN genes b ON a.interval INTERSECTS b.interval",
            peaks=[("chr1", 0, 500000)],
            genes=[("chr1", 5000, 495000)],
            expected=[(0, 5000)],
        )
