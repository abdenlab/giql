"""Integration tests for the DataFusion execution engine."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

try:
    import datafusion  # noqa: F401

    HAS_DATAFUSION = True
except ImportError:
    HAS_DATAFUSION = False

pytestmark = pytest.mark.skipif(
    not HAS_DATAFUSION, reason="datafusion not installed"
)

from giql.engines.datafusion import DataFusionEngine
from giql.table import Table


@pytest.fixture
def engine():
    return DataFusionEngine()


@pytest.fixture
def peaks_table():
    return pa.table(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr2", "chr2"],
            "start": [100, 200, 500, 100, 300],
            "end": [200, 400, 600, 250, 500],
            "name": ["p1", "p2", "p3", "p4", "p5"],
            "strand": ["+", "-", "+", "+", "-"],
        }
    )


@pytest.fixture
def engine_with_peaks(engine, peaks_table):
    engine.register_arrow("peaks", peaks_table)
    return engine


class TestDataFusionRegistration:
    """Tests for table registration."""

    def test_register_arrow_table(self, engine, peaks_table):
        """
        GIVEN a DataFusionEngine and an Arrow table
        WHEN registering the table with register_arrow
        THEN it should be queryable via raw SQL
        """
        engine.register_arrow("peaks", peaks_table)
        result = engine.sql("SELECT COUNT(*) AS cnt FROM peaks")
        assert result.column("cnt").to_pylist() == [5]

    def test_register_arrow_record_batch(self, engine, peaks_table):
        """
        GIVEN a DataFusionEngine and an Arrow RecordBatch
        WHEN registering the batch with register_arrow
        THEN it should be queryable via raw SQL
        """
        batch = peaks_table.to_batches()[0]
        engine.register_arrow("peaks", batch)
        result = engine.sql("SELECT COUNT(*) AS cnt FROM peaks")
        assert result.column("cnt").to_pylist() == [5]

    def test_register_parquet(self, engine, peaks_table, tmp_path):
        """
        GIVEN a DataFusionEngine and a Parquet file
        WHEN registering the file with register_parquet
        THEN it should be queryable via raw SQL
        """
        path = tmp_path / "peaks.parquet"
        pq.write_table(peaks_table, path)
        engine.register_parquet("peaks", path)
        result = engine.sql("SELECT COUNT(*) AS cnt FROM peaks")
        assert result.column("cnt").to_pylist() == [5]

    def test_register_with_custom_table_config(self, engine):
        """
        GIVEN a DataFusionEngine and Arrow data with non-default columns
        WHEN registering with a custom Table config
        THEN GIQL queries should use the custom column mappings
        """
        data = pa.table(
            {
                "chr": ["chr1", "chr1"],
                "pos_start": [100, 500],
                "pos_end": [200, 600],
            }
        )
        table_config = Table(
            "variants",
            chrom_col="chr",
            start_col="pos_start",
            end_col="pos_end",
            strand_col=None,
        )
        engine.register_arrow("variants", data, table=table_config)
        result = engine.query(
            "SELECT * FROM variants"
            " WHERE interval INTERSECTS 'chr1:150-550'"
        )
        assert result.num_rows == 2

    def test_context_property(self, engine):
        """
        GIVEN a DataFusionEngine
        WHEN accessing the context property
        THEN it should return a DataFusion SessionContext
        """
        import datafusion

        assert isinstance(engine.context, datafusion.SessionContext)


class TestDataFusionSpatialOps:
    """Tests for spatial operations via DataFusion."""

    def test_intersects_literal(self, engine_with_peaks):
        """
        GIVEN peaks registered in DataFusion
        WHEN querying with INTERSECTS against a literal range
        THEN only overlapping intervals should be returned
        """
        result = engine_with_peaks.query(
            "SELECT name FROM peaks"
            " WHERE interval INTERSECTS 'chr1:150-250'"
        )
        names = sorted(result.column("name").to_pylist())
        assert names == ["p1", "p2"]

    def test_contains_literal(self, engine_with_peaks):
        """
        GIVEN peaks registered in DataFusion
        WHEN querying with CONTAINS against a literal point
        THEN only intervals containing the point should be returned
        """
        result = engine_with_peaks.query(
            "SELECT name FROM peaks"
            " WHERE interval CONTAINS 'chr1:300'"
        )
        names = result.column("name").to_pylist()
        assert names == ["p2"]

    def test_within_literal(self, engine_with_peaks):
        """
        GIVEN peaks registered in DataFusion
        WHEN querying with WITHIN against a literal range
        THEN only intervals fully within the range should be returned
        """
        result = engine_with_peaks.query(
            "SELECT name FROM peaks"
            " WHERE interval WITHIN 'chr1:50-650'"
        )
        names = sorted(result.column("name").to_pylist())
        assert "p1" in names
        assert "p3" in names

    def test_intersects_cross_chrom(self, engine_with_peaks):
        """
        GIVEN peaks on multiple chromosomes registered in DataFusion
        WHEN querying with INTERSECTS on chr2
        THEN only chr2 intervals should be returned
        """
        result = engine_with_peaks.query(
            "SELECT name FROM peaks"
            " WHERE interval INTERSECTS 'chr2:200-350'"
        )
        names = sorted(result.column("name").to_pylist())
        assert names == ["p4", "p5"]


class TestDataFusionAggregation:
    """Tests for CLUSTER and MERGE operations."""

    def test_cluster(self, engine):
        """
        GIVEN contiguous/overlapping intervals registered in DataFusion
        WHEN querying with CLUSTER
        THEN each row should get a cluster_id grouping overlapping intervals
        """
        data = pa.table(
            {
                "chrom": ["chr1", "chr1", "chr1", "chr1"],
                "start": [100, 150, 500, 550],
                "end": [200, 300, 600, 700],
                "strand": ["+", "+", "+", "+"],
            }
        )
        engine.register_arrow("intervals", data)
        result = engine.query(
            "SELECT *, CLUSTER(interval) AS cluster_id FROM intervals"
        )
        cluster_ids = result.column("cluster_id").to_pylist()
        # First two overlap -> same cluster; last two overlap -> same cluster
        assert cluster_ids[0] == cluster_ids[1]
        assert cluster_ids[2] == cluster_ids[3]
        assert cluster_ids[0] != cluster_ids[2]

    def test_merge(self, engine):
        """
        GIVEN overlapping intervals registered in DataFusion
        WHEN querying with MERGE
        THEN overlapping intervals should be merged into single intervals
        """
        data = pa.table(
            {
                "chrom": ["chr1", "chr1", "chr1", "chr1"],
                "start": [100, 150, 500, 550],
                "end": [200, 300, 600, 700],
                "strand": ["+", "+", "+", "+"],
            }
        )
        engine.register_arrow("intervals", data)
        result = engine.query(
            "SELECT MERGE(interval) FROM intervals"
        )
        assert result.num_rows == 2
        starts = sorted(result.column("start").to_pylist())
        ends = sorted(result.column("end").to_pylist())
        assert starts == [100, 500]
        assert ends == [300, 700]


class TestDataFusionDistance:
    """Tests for DISTANCE function."""

    def test_distance_gap(self, engine):
        """
        GIVEN two tables with non-overlapping intervals on the same chromosome
        WHEN computing DISTANCE between them
        THEN the distance should be the gap between the intervals
        """
        left = pa.table(
            {
                "chrom": ["chr1"],
                "start": [100],
                "end": [200],
                "strand": ["+"],
            }
        )
        right = pa.table(
            {
                "chrom": ["chr1"],
                "start": [300],
                "end": [400],
                "strand": ["+"],
            }
        )
        engine.register_arrow("a", left)
        engine.register_arrow("b", right)
        result = engine.query(
            "SELECT DISTANCE(a.interval, b.interval) AS dist"
            " FROM a CROSS JOIN b"
        )
        dist = result.column("dist").to_pylist()[0]
        assert dist == 100

    def test_distance_overlapping(self, engine):
        """
        GIVEN two tables with overlapping intervals on the same chromosome
        WHEN computing DISTANCE between them
        THEN the distance should be zero (overlapping intervals have no gap)
        """
        left = pa.table(
            {
                "chrom": ["chr1"],
                "start": [100],
                "end": [300],
                "strand": ["+"],
            }
        )
        right = pa.table(
            {
                "chrom": ["chr1"],
                "start": [200],
                "end": [400],
                "strand": ["+"],
            }
        )
        engine.register_arrow("a", left)
        engine.register_arrow("b", right)
        result = engine.query(
            "SELECT DISTANCE(a.interval, b.interval) AS dist"
            " FROM a CROSS JOIN b"
        )
        dist = result.column("dist").to_pylist()[0]
        assert dist == 0

    def test_distance_cross_chrom(self, engine):
        """
        GIVEN two tables with intervals on different chromosomes
        WHEN computing DISTANCE between them
        THEN the distance should be NULL
        """
        left = pa.table(
            {
                "chrom": ["chr1"],
                "start": [100],
                "end": [200],
                "strand": ["+"],
            }
        )
        right = pa.table(
            {
                "chrom": ["chr2"],
                "start": [100],
                "end": [200],
                "strand": ["+"],
            }
        )
        engine.register_arrow("a", left)
        engine.register_arrow("b", right)
        result = engine.query(
            "SELECT DISTANCE(a.interval, b.interval) AS dist"
            " FROM a CROSS JOIN b"
        )
        dist = result.column("dist").to_pylist()[0]
        assert dist is None


class TestDataFusionNearest:
    """Tests for standalone NEAREST (ORDER BY + LIMIT mode)."""

    @pytest.mark.xfail(
        reason="DataFusion planner bug with complex ORDER BY expressions in NEAREST",
        raises=Exception,
    )
    def test_nearest_basic(self, engine):
        """
        GIVEN a targets table and a query interval
        WHEN using NEAREST to find the closest interval
        THEN the nearest interval by distance should be returned
        """
        targets = pa.table(
            {
                "chrom": ["chr1", "chr1", "chr1"],
                "start": [100, 500, 900],
                "end": [200, 600, 1000],
                "name": ["t1", "t2", "t3"],
                "strand": ["+", "+", "+"],
            }
        )
        engine.register_arrow("targets", targets)
        result = engine.query(
            "SELECT * FROM NEAREST(targets, reference='chr1:450-470', k=1)"
        )
        names = result.column("name").to_pylist()
        assert names == ["t2"]
