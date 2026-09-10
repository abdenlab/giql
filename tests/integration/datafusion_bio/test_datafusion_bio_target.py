"""Execution and plan tests for the ``dialect="datafusion-bio"`` target (#77).

Every oracle case here runs the vanilla ``datafusion`` target, the ``duckdb``
target and the ``datafusion-bio`` target and asserts row identity across all
three; the cases that matter for performance additionally assert on polars-bio's
``EXPLAIN`` output, because the whole point of the target is which physical
operator the join lands on. The downstream cCRE insertion-matrix script's exact
shape (LEFT JOIN + ``COUNT(right.start)`` + GROUP BY over a custom-column table)
is covered directly against polars-bio on both integer widths.
"""

import pytest

pytest.importorskip("polars_bio")
pytest.importorskip("polars")

import polars as pl  # noqa: E402
import polars_bio as pb  # noqa: E402

from giql import Table  # noqa: E402
from giql import transpile  # noqa: E402

pytestmark = pytest.mark.integration

_TARGETS = ("datafusion", "duckdb", "datafusion-bio")
_INT32_COLUMNS = (("chrom", "utf8"), ("start", "int32"), ("end", "int32"))
_STRANDED_COLUMNS = (
    ("chrom", "utf8"),
    ("start", "int64"),
    ("end", "int64"),
    ("strand", "utf8"),
)

# Three cCRE-like regions and four 1-bp insertion sites: two inside the first
# region (one at its last base), one book-ended at its exclusive end (must not
# count), and one on chr2 outside every region.
_CCRES = [("chr1", 100, 200), ("chr1", 300, 400), ("chr2", 100, 200)]
_INSERTIONS = [
    ("chr1", 150, 151),
    ("chr1", 199, 200),
    ("chr1", 200, 201),
    ("chr2", 50, 51),
]
_SCRIPT_QUERY = (
    "SELECT c.chrom, c.start, COUNT(r.start) AS n FROM ccres c "
    "LEFT JOIN insertions r ON c.interval INTERSECTS r.interval "
    "GROUP BY c.chrom, c.start"
)
_SCRIPT_EXPECTED = [("chr1", 100, 2), ("chr1", 300, 0), ("chr2", 100, 0)]


class TestScriptShape:
    """The cCRE insertion-matrix query shape on every target and both widths."""

    def test_left_join_count_should_agree_and_plan_interval_join_when_int64(
        self, cross_target_oracle, polars_bio_explain
    ):
        """Test the script shape on Int64 coordinates.

        Given:
            cCRE regions and 1-bp insertion sites on Int64 coordinates, with one
            region matching two sites, one site book-ended at a region's end,
            and two regions matching nothing.
        When:
            The LEFT JOIN + COUNT(r.start) + GROUP BY query runs on every target.
        Then:
            It should count 2, 0 and 0 on every engine, and polars-bio should plan
            the datafusion-bio SQL's inner count through IntervalJoinExec (the
            case the strict predicate fails at execution on) with the zero-fill
            riding a hash join on the keys.
        """
        # Arrange / Act
        cross_target_oracle(
            _SCRIPT_QUERY,
            ccres=_CCRES,
            insertions=_INSERTIONS,
            expected=_SCRIPT_EXPECTED,
            targets=_TARGETS,
        )
        plan = polars_bio_explain(cross_target_oracle.sql_by_target["datafusion-bio"])

        # Assert
        assert "IntervalJoinExec" in plan

    def test_left_join_count_should_agree_and_plan_interval_join_when_int32(
        self, cross_target_oracle, polars_bio_explain
    ):
        """Test the script shape on Int32 coordinates.

        Given:
            The same regions and sites declared with 32-bit coordinate columns.
        When:
            The LEFT JOIN + COUNT(r.start) + GROUP BY query runs on every target.
        Then:
            It should return the same counts and still plan through
            IntervalJoinExec, so the closed form costs nothing on the width the
            strict form already handled.
        """
        # Arrange / Act
        cross_target_oracle(
            _SCRIPT_QUERY,
            ccres=_CCRES,
            insertions=_INSERTIONS,
            expected=_SCRIPT_EXPECTED,
            targets=_TARGETS,
            columns=_INT32_COLUMNS,
            tables=["ccres", "insertions"],
        )
        plan = polars_bio_explain(cross_target_oracle.sql_by_target["datafusion-bio"])

        # Assert
        assert "IntervalJoinExec" in plan

    def test_left_join_count_should_return_zeros_when_right_table_empty(
        self, cross_target_oracle
    ):
        """Test the script shape against an empty insertions table.

        Given:
            cCRE regions and no insertion sites at all.
        When:
            The LEFT JOIN + COUNT(r.start) + GROUP BY query runs on every target.
        Then:
            It should keep every region with a count of 0 on every engine.
        """
        # Arrange / Act / Assert
        cross_target_oracle(
            _SCRIPT_QUERY,
            ccres=_CCRES,
            insertions=[],
            expected=[("chr1", 100, 0), ("chr1", 300, 0), ("chr2", 100, 0)],
            targets=_TARGETS,
        )

    def test_script_shape_should_execute_when_custom_columns_and_mixed_widths(
        self, polars_bio_register, polars_bio_explain
    ):
        """Test the script's exact table mapping against polars-bio directly.

        Given:
            A BED-like ccres table (Int64, with a ccre_id) and a BAM-derived
            insertions table whose chromosome column is ``rname`` and whose
            coordinates are Int32, as oxbow produces them.
        When:
            The script's query is transpiled with dialect="datafusion-bio" and
            executed on polars-bio.
        Then:
            It should return the per-cCRE counts, and the plan should use
            IntervalJoinExec despite the mixed integer widths across sides.
        """
        # Arrange
        polars_bio_register(
            "ccres",
            {"chrom": pl.Utf8, "start": pl.Int64, "end": pl.Int64, "ccre_id": pl.Utf8},
            [
                ("chr1", 100, 200, "E1"),
                ("chr1", 300, 400, "E2"),
                ("chr2", 100, 200, "E3"),
            ],
        )
        polars_bio_register(
            "insertions",
            {"rname": pl.Utf8, "start": pl.Int32, "end": pl.Int32},
            _INSERTIONS,
        )
        sql = transpile(
            "SELECT c.ccre_id, COUNT(r.start) AS n FROM ccres c "
            "LEFT JOIN insertions r ON c.interval INTERSECTS r.interval "
            "GROUP BY c.ccre_id",
            tables=[
                "ccres",
                Table(
                    "insertions",
                    chrom_col="rname",
                    start_col="start",
                    end_col="end",
                    strand_col=None,
                ),
            ],
            dialect="datafusion-bio",
        )

        # Act
        rows = sorted(pb.sql(sql).collect().rows())
        plan = polars_bio_explain(sql)

        # Assert
        assert rows == [("E1", 2), ("E2", 0), ("E3", 0)]
        assert "IntervalJoinExec" in plan


class TestJoinShapes:
    """Join shapes the interval-join rule accepts, plus the documented fallback."""

    def test_inner_join_should_agree_and_plan_interval_join(
        self, cross_target_oracle, polars_bio_explain
    ):
        """Test an inner column-to-column INTERSECTS join.

        Given:
            Two interval tables with one overlapping pair and one book-ended pair.
        When:
            An inner INTERSECTS join projecting both sides runs on every target.
        Then:
            It should return only the overlapping pair everywhere and plan through
            IntervalJoinExec on polars-bio.
        """
        # Arrange / Act
        cross_target_oracle(
            "SELECT a.start, b.start AS bs FROM a JOIN b "
            "ON a.interval INTERSECTS b.interval",
            a=[("chr1", 10, 20), ("chr1", 100, 110)],
            b=[("chr1", 15, 25), ("chr1", 20, 30), ("chr2", 10, 20)],
            expected=[(10, 15)],
            targets=_TARGETS,
        )
        plan = polars_bio_explain(cross_target_oracle.sql_by_target["datafusion-bio"])

        # Assert
        assert "IntervalJoinExec" in plan

    @pytest.mark.parametrize(
        ("query", "expected"),
        [
            (
                "SELECT a.start, b.start AS bs FROM a LEFT JOIN b "
                "ON a.interval INTERSECTS b.interval",
                [(10, 15), (10, 15), (100, None)],
            ),
            (
                "SELECT a.start, b.start AS bs FROM a RIGHT JOIN b "
                "ON a.interval INTERSECTS b.interval",
                [(10, 15), (10, 15), (None, 500)],
            ),
            (
                "SELECT a.start, b.start AS bs FROM a FULL OUTER JOIN b "
                "ON a.interval INTERSECTS b.interval",
                [(10, 15), (10, 15), (100, None), (None, 500)],
            ),
            (
                "SELECT a.start FROM a WHERE EXISTS "
                "(SELECT 1 FROM b WHERE a.interval INTERSECTS b.interval)",
                [(10,)],
            ),
            (
                "SELECT a.start FROM a WHERE NOT EXISTS "
                "(SELECT 1 FROM b WHERE a.interval INTERSECTS b.interval)",
                [(100,)],
            ),
        ],
        ids=["left", "right", "full", "exists", "not-exists"],
    )
    def test_non_inner_join_should_agree_and_avoid_interval_join(
        self, cross_target_oracle, polars_bio_explain, query, expected
    ):
        """Test outer, semi and anti INTERSECTS shapes stay correct on polars-bio.

        Given:
            A left table with one interval matching two right intervals and one
            matching nothing, and a right table with one unmatched interval.
        When:
            An outer, EXISTS or NOT EXISTS INTERSECTS query runs on every target.
        Then:
            It should return the same rows everywhere, and the polars-bio plan
            should NOT use IntervalJoinExec: that operator returns the inner pairs
            for every join type upstream, so the target deliberately emits a
            rule-defeating predicate and takes the hash join instead.
        """
        # Arrange / Act
        cross_target_oracle(
            query,
            a=[("chr1", 10, 20), ("chr1", 100, 110)],
            b=[("chr1", 15, 25), ("chr1", 15, 18), ("chr1", 500, 510)],
            expected=expected,
            targets=_TARGETS,
        )
        plan = polars_bio_explain(cross_target_oracle.sql_by_target["datafusion-bio"])

        # Assert
        assert "IntervalJoinExec" not in plan
        assert "HashJoinExec" in plan

    @pytest.mark.parametrize("operator", ["CONTAINS", "WITHIN"])
    def test_left_containment_join_should_agree_and_avoid_interval_join(
        self, cross_target_oracle, polars_bio_explain, operator
    ):
        """Test LEFT CONTAINS / WITHIN joins stay correct on polars-bio.

        Given:
            A left table whose first interval nests with a right interval and
            whose second matches nothing.
        When:
            A LEFT CONTAINS or WITHIN join runs on every target.
        Then:
            It should preserve the unmatched left row everywhere and avoid
            IntervalJoinExec on polars-bio, since the generic non-strict form
            would otherwise reach the join-type-blind operator.
        """
        # Arrange / Act
        cross_target_oracle(
            "SELECT a.start, b.start AS bs FROM a LEFT JOIN b "
            f"ON a.interval {operator} b.interval",
            a=[("chr1", 10, 20), ("chr1", 100, 110)],
            b=[("chr1", 10, 20)],
            expected=[(10, 10), (100, None)],
            targets=_TARGETS,
        )
        plan = polars_bio_explain(cross_target_oracle.sql_by_target["datafusion-bio"])

        # Assert
        assert "IntervalJoinExec" not in plan

    @pytest.mark.parametrize(
        ("operator", "expected"),
        [("CONTAINS", [(10, 12)]), ("WITHIN", [(12, 10)])],
    )
    def test_containment_join_should_agree_and_plan_interval_join(
        self, cross_target_oracle, polars_bio_explain, operator, expected
    ):
        """Test CONTAINS / WITHIN column-to-column joins.

        Given:
            A wide interval and a narrow one nested inside it, on both tables.
        When:
            A CONTAINS or WITHIN join runs on every target.
        Then:
            It should return the nested pair in the operator's orientation
            everywhere and plan through IntervalJoinExec on polars-bio without
            any datafusion-bio-specific rewrite.
        """
        # Arrange / Act
        cross_target_oracle(
            "SELECT a.start, b.start AS bs FROM a JOIN b "
            f"ON a.interval {operator} b.interval",
            a=[("chr1", 10, 20), ("chr1", 12, 15)],
            b=[("chr1", 10, 20), ("chr1", 12, 15)],
            expected=expected + [(10, 10), (12, 12)],
            targets=_TARGETS,
        )
        plan = polars_bio_explain(cross_target_oracle.sql_by_target["datafusion-bio"])

        # Assert
        assert "IntervalJoinExec" in plan

    def test_grouped_count_star_should_agree_and_plan_interval_join(
        self, cross_target_oracle, polars_bio_explain
    ):
        """Test a grouped COUNT(*) over an inner INTERSECTS join.

        Given:
            Two left intervals, one overlapping two right intervals and one
            overlapping one.
        When:
            A COUNT(*) ... GROUP BY left-keys query runs on every target.
        Then:
            It should return the per-key counts everywhere and plan through
            IntervalJoinExec on polars-bio, with COUNT(*) left untouched because
            the grouped keys already carry a column through the join.
        """
        # Arrange / Act
        cross_target_oracle(
            "SELECT a.chrom, a.start, COUNT(*) AS n FROM a JOIN b "
            "ON a.interval INTERSECTS b.interval GROUP BY a.chrom, a.start",
            a=[("chr1", 10, 20), ("chr1", 100, 110)],
            b=[("chr1", 15, 25), ("chr1", 12, 18), ("chr1", 105, 106)],
            expected=[("chr1", 10, 2), ("chr1", 100, 1)],
            targets=_TARGETS,
        )
        sql = cross_target_oracle.sql_by_target["datafusion-bio"]
        plan = polars_bio_explain(sql)

        # Assert
        assert "COUNT(*)" in sql
        assert "IntervalJoinExec" in plan

    def test_bare_count_star_should_execute_when_projection_carried(
        self, cross_target_oracle, polars_bio_explain
    ):
        """Test the bare COUNT(*) carried-projection workaround end to end.

        Given:
            Two interval tables with three overlapping pairs.
        When:
            SELECT COUNT(*) with no other projected column runs on every target.
        Then:
            It should return 3 everywhere; the datafusion-bio SQL should count the
            FROM side's start column instead, and still plan through
            IntervalJoinExec, since an empty-projection interval join fails
            upstream.
        """
        # Arrange / Act
        cross_target_oracle(
            "SELECT COUNT(*) AS n FROM a JOIN b ON a.interval INTERSECTS b.interval",
            a=[("chr1", 10, 20), ("chr1", 100, 110)],
            b=[("chr1", 15, 25), ("chr1", 12, 18), ("chr1", 105, 106)],
            expected=[(3,)],
            targets=_TARGETS,
        )
        sql = cross_target_oracle.sql_by_target["datafusion-bio"]
        plan = polars_bio_explain(sql)

        # Assert
        assert 'COUNT(a."start") AS n' in sql
        assert "IntervalJoinExec" in plan

    def test_stranded_equality_should_agree_and_plan_interval_join(
        self, cross_target_oracle, polars_bio_explain
    ):
        """Test INTERSECTS plus a same-strand equality conjunct.

        Given:
            Stranded interval tables where one overlapping pair shares a strand
            and another does not.
        When:
            An INTERSECTS AND a.strand = b.strand join runs on every target.
        Then:
            It should return only the same-strand pair everywhere and keep the
            IntervalJoinExec plan, the equality riding as a hash key.
        """
        # Arrange / Act
        cross_target_oracle(
            "SELECT a.start, b.start AS bs FROM a JOIN b "
            "ON a.interval INTERSECTS b.interval AND a.strand = b.strand",
            a=[("chr1", 10, 20, "+"), ("chr1", 100, 110, "-")],
            b=[("chr1", 15, 25, "+"), ("chr1", 105, 115, "+")],
            expected=[(10, 15)],
            targets=_TARGETS,
            columns=_STRANDED_COLUMNS,
            tables=[Table("a"), Table("b")],
        )
        plan = polars_bio_explain(cross_target_oracle.sql_by_target["datafusion-bio"])

        # Assert
        assert "IntervalJoinExec" in plan

    def test_two_sided_residual_should_agree_but_fall_back_to_hash_join(
        self, cross_target_oracle, polars_bio_explain
    ):
        """Test the documented fallback for a two-sided non-equi residual.

        Given:
            Stranded interval tables where one overlapping pair shares a strand
            and another does not.
        When:
            An INTERSECTS AND a.strand <> b.strand join runs on every target.
        Then:
            It should return only the opposite-strand pair everywhere, but the
            polars-bio plan should NOT use IntervalJoinExec: the third non-equi
            conjunct defeats the rule and the join falls back to a hash join with
            a residual filter, as the performance guide documents.
        """
        # Arrange / Act
        cross_target_oracle(
            "SELECT a.start, b.start AS bs FROM a JOIN b "
            "ON a.interval INTERSECTS b.interval AND a.strand <> b.strand",
            a=[("chr1", 10, 20, "+"), ("chr1", 100, 110, "-")],
            b=[("chr1", 15, 25, "+"), ("chr1", 105, 115, "+")],
            expected=[(100, 105)],
            targets=_TARGETS,
            columns=_STRANDED_COLUMNS,
            tables=[Table("a"), Table("b")],
        )
        plan = polars_bio_explain(cross_target_oracle.sql_by_target["datafusion-bio"])

        # Assert
        assert "IntervalJoinExec" not in plan
        assert "HashJoinExec" in plan
