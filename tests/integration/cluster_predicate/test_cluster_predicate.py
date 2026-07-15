"""Integration tests for the CLUSTER/MERGE predicate argument.

These execute the generated SQL against DuckDB and assert the functional
behavior of ``predicate :=``: run-length encoding of equal-valued runs,
reconstruction of Bioconductor ``disjoin()``-style depth-annotated output via a
DISJOIN -> depth-aggregation -> predicate-MERGE pipeline, and the documented
single-linkage drift for non-equivalence predicates.
"""

import pytest

pytestmark = pytest.mark.integration


def _ival(chrom, start, end, name, score, strand="+"):
    """Build a 6-tuple for ``load_intervals``."""
    return (chrom, start, end, name, score, strand)


class TestMergePredicate:
    """Functional behavior of MERGE(interval, predicate := ...)."""

    def test_merge_without_predicate_coalesces_adjacent_run(self, giql_query):
        """Test that MERGE without a predicate collapses an adjacent run.

        Given:
            Five abutting intervals spanning [0, 50) with mixed scores.
        When:
            MERGE(interval) runs with no predicate.
        Then:
            It should collapse the whole adjacent run into a single interval.
        """
        # Arrange & act
        rows = giql_query(
            "SELECT MERGE(interval) FROM intervals",
            tables=["intervals"],
            intervals=[
                _ival("chr1", 0, 10, "a", 5),
                _ival("chr1", 10, 20, "b", 5),
                _ival("chr1", 20, 30, "c", 3),
                _ival("chr1", 30, 40, "d", 3),
                _ival("chr1", 40, 50, "e", 5),
            ],
        )

        # Assert
        assert rows == [("chr1", 0, 50)]

    def test_merge_with_equality_predicate_run_length_encodes(self, giql_query):
        """Test that an equality predicate run-length-encodes adjacent runs.

        Given:
            Five abutting intervals whose scores form the runs 5,5 | 3,3 | 5.
        When:
            MERGE(interval, predicate := score = PREV(score)) runs.
        Then:
            It should emit one merged interval per maximal equal-score run.
        """
        # Arrange & act
        rows = giql_query(
            "SELECT MERGE(interval, predicate := score = PREV(score)) FROM intervals",
            tables=["intervals"],
            intervals=[
                _ival("chr1", 0, 10, "a", 5),
                _ival("chr1", 10, 20, "b", 5),
                _ival("chr1", 20, 30, "c", 3),
                _ival("chr1", 30, 40, "d", 3),
                _ival("chr1", 40, 50, "e", 5),
            ],
        )

        # Assert
        assert rows == [
            ("chr1", 0, 20),
            ("chr1", 20, 40),
            ("chr1", 40, 50),
        ]

    def test_merge_predicate_drifts_under_single_linkage(self, giql_query):
        """Test that a non-equivalence predicate exhibits single-linkage drift.

        Given:
            Three abutting intervals with scores 10, 13, 16, where each
            consecutive pair differs by 3 but the extremes differ by 6.
        When:
            MERGE(interval, predicate := ABS(score - PREV(score)) < 5) runs.
        Then:
            It should merge the entire run into one interval even though the
            cluster's extremes violate the predicate (documented drift).
        """
        # Arrange & act
        rows = giql_query(
            "SELECT MERGE(interval, predicate := ABS(score - PREV(score)) < 5) "
            "FROM intervals",
            tables=["intervals"],
            intervals=[
                _ival("chr1", 0, 10, "a", 10),
                _ival("chr1", 10, 20, "b", 13),
                _ival("chr1", 20, 30, "c", 16),
            ],
        )

        # Assert
        assert rows == [("chr1", 0, 30)]

    def test_merge_compound_predicate_breaks_on_either_column(self, giql_query):
        """Test that a multi-column AND predicate breaks on any column change.

        Given:
            Four abutting intervals where the (name, strand) pair changes first
            on strand and then on name partway through the run.
        When:
            MERGE(interval, predicate := strand = PREV(strand) AND
            name = PREV(name)) runs.
        Then:
            It should break the merge wherever either column changes, yielding
            one merged interval per homogeneous (name, strand) run.
        """
        # Arrange & act
        rows = giql_query(
            "SELECT MERGE(interval, predicate := strand = PREV(strand) "
            "AND name = PREV(name)) FROM intervals",
            tables=["intervals"],
            intervals=[
                _ival("chr1", 0, 10, "g", 5, "+"),
                _ival("chr1", 10, 20, "g", 5, "+"),
                _ival("chr1", 20, 30, "g", 5, "-"),
                _ival("chr1", 30, 40, "h", 5, "-"),
            ],
        )

        # Assert
        assert rows == [
            ("chr1", 0, 20),
            ("chr1", 20, 30),
            ("chr1", 30, 40),
        ]

    def test_merge_stranded_predicate_evaluates_within_strand(self, giql_query):
        """Test that the predicate is evaluated within each strand partition.

        Given:
            Equal-name abutting intervals on both the + and - strands.
        When:
            MERGE(interval, stranded := true, predicate := name = PREV(name))
            runs.
        Then:
            It should merge each strand's equal-name run independently, emitting
            one interval per strand rather than collapsing across strands.
        """
        # Arrange & act
        rows = giql_query(
            "SELECT MERGE(interval, stranded := true, predicate := name = PREV(name)) "
            "FROM intervals",
            tables=["intervals"],
            intervals=[
                _ival("chr1", 0, 10, "a", 5, "+"),
                _ival("chr1", 10, 20, "a", 5, "+"),
                _ival("chr1", 0, 10, "a", 5, "-"),
                _ival("chr1", 10, 20, "a", 5, "-"),
            ],
        )

        # Assert
        # The emitted MERGE SQL orders only by (chrom, start); these two rows
        # tie on (chr1, 0) and differ only by strand, so the row order is not
        # deterministic. Compare order-independently.
        assert sorted(rows) == sorted(
            [
                ("chr1", "+", 0, 20),
                ("chr1", "-", 0, 20),
            ]
        )


class TestClusterPredicate:
    """Functional behavior of CLUSTER(interval, predicate := ...)."""

    def test_cluster_with_equality_predicate_assigns_run_ids(self, giql_query):
        """Test that an equality predicate assigns one cluster id per run.

        Given:
            Five abutting intervals whose scores form the runs 5,5 | 3,3 | 5.
        When:
            CLUSTER(interval, predicate := score = PREV(score)) runs.
        Then:
            It should assign a distinct cluster id to each maximal equal-score
            run, comparing each row only to its immediate predecessor.
        """
        # Arrange & act
        rows = giql_query(
            "SELECT name, CLUSTER(interval, predicate := score = PREV(score)) AS cid "
            "FROM intervals",
            tables=["intervals"],
            intervals=[
                _ival("chr1", 0, 10, "a", 5),
                _ival("chr1", 10, 20, "b", 5),
                _ival("chr1", 20, 30, "c", 3),
                _ival("chr1", 30, 40, "d", 3),
                _ival("chr1", 40, 50, "e", 5),
            ],
        )

        # Assert
        cid = dict(rows)
        assert cid["a"] == cid["b"]
        assert cid["c"] == cid["d"]
        assert len({cid["a"], cid["c"], cid["e"]}) == 3

    def test_cluster_distance_and_predicate_compose(self, giql_query):
        """Test that the distance gate and predicate gate both apply.

        Given:
            Three equal-score intervals separated by a 50bp gap then a 150bp
            gap, clustered with CLUSTER(interval, 100, predicate := score =
            PREV(score)).
        When:
            The generated SQL runs in DuckDB.
        Then:
            It should keep the 50bp-gap pair together but start a new cluster at
            the 150bp gap even though the scores match, proving distance and
            predicate compose.
        """
        # Arrange & act
        rows = giql_query(
            "SELECT name, CLUSTER(interval, 100, predicate := score = PREV(score)) "
            "AS cid FROM intervals",
            tables=["intervals"],
            intervals=[
                _ival("chr1", 100, 200, "i1", 5),
                _ival("chr1", 250, 350, "i2", 5),  # 50bp gap, score matches
                _ival("chr1", 500, 600, "i3", 5),  # 150bp gap, score matches
            ],
        )

        # Assert
        cid = dict(rows)
        assert cid["i1"] == cid["i2"]
        assert cid["i3"] != cid["i1"]

    def test_cluster_predicate_resets_at_chromosome_boundary(self, giql_query):
        """Test that the predicate does not carry across chromosomes.

        Given:
            Equal-score abutting intervals on chr1 and chr2, clustered with
            CLUSTER(interval, predicate := score = PREV(score)).
        When:
            The generated SQL runs in DuckDB.
        Then:
            It should begin a fresh cluster for the first row of each
            chromosome, since the predecessor LAG resets at the per-chromosome
            partition boundary.
        """
        # Arrange & act
        rows = giql_query(
            "SELECT chrom, name, CLUSTER(interval, predicate := score = PREV(score)) "
            "AS cid FROM intervals",
            tables=["intervals"],
            intervals=[
                _ival("chr1", 0, 10, "a", 5),
                _ival("chr1", 10, 20, "b", 5),
                _ival("chr2", 0, 10, "c", 5),
                _ival("chr2", 10, 20, "d", 5),
            ],
        )

        # Assert
        by_name = {name: (chrom, cid) for chrom, name, cid in rows}
        # Same chromosome, equal score, adjacent -> same cluster.
        assert by_name["a"][1] == by_name["b"][1]
        assert by_name["c"][1] == by_name["d"][1]
        # chr2's first row starts its own cluster (ids are per-partition).
        assert {by_name["a"], by_name["c"]} == {("chr1", 1), ("chr2", 1)}

    def test_cluster_reserved_word_predicate_executes(self, giql_query):
        """Test that a predicate over reserved-word columns executes.

        Given:
            Three abutting intervals and a predicate over the reserved-word
            genomic columns start and end (start = PREV(end)).
        When:
            The generated SQL runs in DuckDB.
        Then:
            It should emit valid quoted SQL and merge the abutting run, proving
            reserved-word predicate columns are handled.
        """
        # Arrange & act
        rows = giql_query(
            "SELECT MERGE(interval, predicate := start = PREV(end)) FROM intervals",
            tables=["intervals"],
            intervals=[
                _ival("chr1", 0, 10, "a", 5),
                _ival("chr1", 10, 20, "b", 5),
                _ival("chr1", 20, 30, "c", 5),
            ],
        )

        # Assert
        assert rows == [("chr1", 0, 30)]

    def test_cluster_predicate_column_absent_from_projection_executes(self, giql_query):
        """Test that a predicate column need not appear in the projection.

        Given:
            A CLUSTER query whose explicit projection omits the predicate
            column score (selecting only chrom, start, end, and the cluster id).
        When:
            The generated SQL runs in DuckDB.
        Then:
            It should still cluster by the score runs, confirming the predicate
            column is projected into the intermediate CTE.
        """
        # Arrange & act
        rows = giql_query(
            'SELECT chrom, start, "end", '
            "CLUSTER(interval, predicate := score = PREV(score)) AS cid "
            "FROM intervals",
            tables=["intervals"],
            intervals=[
                _ival("chr1", 0, 10, "a", 5),
                _ival("chr1", 10, 20, "b", 5),
                _ival("chr1", 20, 30, "c", 3),
            ],
        )

        # Assert
        cids = [row[-1] for row in rows]
        assert cids[0] == cids[1]
        assert cids[2] != cids[0]


class TestDisjoinDepthMergePipeline:
    """DISJOIN -> depth-aggregation -> predicate-MERGE reconstructs disjoin()."""

    _FEATURES = [
        _ival("chr1", 0, 20, "a", 0),
        _ival("chr1", 0, 20, "b", 0),
        _ival("chr1", 20, 40, "c", 0),
        _ival("chr1", 20, 40, "d", 0),
        _ival("chr1", 10, 30, "e", 0),
    ]
    # Per-breakpoint coverage depth: [0,10)=2, [10,20)=3, [20,30)=3, [30,40)=2.

    _DEPTH_SEGMENTS_QUERY = """
        SELECT disjoin_chrom AS chrom,
               disjoin_start AS start,
               disjoin_end AS "end",
               COUNT(*) AS depth
        FROM DISJOIN(features)
        GROUP BY disjoin_chrom, disjoin_start, disjoin_end
    """

    def test_disjoin_depth_segments_form_expected_partition(self, giql_query):
        """Test that DISJOIN + depth aggregation yields per-breakpoint depths.

        Given:
            Two doubled abutting regions overlaid by a spanning interval,
            producing per-breakpoint coverage depths 2, 3, 3, 2.
        When:
            DISJOIN(features) is aggregated to per-segment coverage depth.
        Then:
            It should yield four disjoint segments carrying those depths.
        """
        # Arrange & act
        rows = giql_query(
            self._DEPTH_SEGMENTS_QUERY + ' ORDER BY "start"',
            tables=["features"],
            features=self._FEATURES,
        )

        # Assert
        assert rows == [
            ("chr1", 0, 10, 2),
            ("chr1", 10, 20, 3),
            ("chr1", 20, 30, 3),
            ("chr1", 30, 40, 2),
        ]

    def test_predicate_merge_remerges_equal_depth_segments(self, giql_query):
        """Test that predicate-MERGE re-clusters runs of equal coverage depth.

        Given:
            The depth-annotated disjoint segments (depths 2, 3, 3, 2).
        When:
            MERGE(interval, predicate := depth = PREV(depth)) runs over them.
        Then:
            It should coalesce the adjacent depth-3 run into [10, 30) while
            keeping the depth-2 flanks distinct, reproducing the re-clustered
            partition Bioconductor disjoin() emits.
        """
        # Arrange & act
        rows = giql_query(
            "SELECT MERGE(interval, predicate := depth = PREV(depth)) "
            f"FROM ({self._DEPTH_SEGMENTS_QUERY}) AS segments",
            tables=["features"],
            features=self._FEATURES,
        )

        # Assert
        assert rows == [
            ("chr1", 0, 10),
            ("chr1", 10, 30),
            ("chr1", 30, 40),
        ]

    def test_predicate_merge_differs_from_unconditioned_merge(self, giql_query):
        """Test that the predicate is what preserves the depth boundaries.

        Given:
            The same depth-annotated disjoint segments.
        When:
            MERGE(interval) runs over them with no predicate.
        Then:
            It should collapse every adjacent segment into one interval,
            confirming the predicate alone preserves the coverage structure.
        """
        # Arrange & act
        rows = giql_query(
            f"SELECT MERGE(interval) FROM ({self._DEPTH_SEGMENTS_QUERY}) AS segments",
            tables=["features"],
            features=self._FEATURES,
        )

        # Assert
        assert rows == [("chr1", 0, 40)]
