"""Tests for CLUSTER and MERGE operations."""

import pytest

from giql import GIQLEngine


@pytest.fixture
def cluster_test_data_csv(tmp_path):
    """Create sample data for cluster testing."""
    csv_content = """
    id,chromosome,start_pos,end_pos,name
    1,chr1,100,200,f1
    2,chr1,180,250,f2
    3,chr1,250,500,f3
    4,chr1,501,1000,f4
    5,chr2,100,200,f5
    6,chr2,300,400,f6
    """
    csv_path = tmp_path / "features.csv"
    csv_path.write_text(csv_content.strip())
    return str(csv_path)


@pytest.fixture
def stranded_test_data_csv(tmp_path):
    """Create stranded data for cluster testing."""
    csv_content = """
    id,chromosome,start_pos,end_pos,strand,name
    1,chr1,100,200,+,f1
    2,chr1,180,250,+,f2
    3,chr1,200,300,-,f3
    4,chr1,250,350,-,f4
    5,chr1,400,500,+,f5
    """
    csv_path = tmp_path / "stranded_features.csv"
    csv_path.write_text(csv_content.strip())
    return str(csv_path)


@pytest.fixture
def duckdb_cluster_engine(cluster_test_data_csv):
    """DuckDB engine with cluster test data loaded."""
    engine = GIQLEngine(target_dialect="duckdb", verbose=True)
    engine.load_csv("features", cluster_test_data_csv)
    engine.register_table_schema(
        "features",
        {
            "id": "INTEGER",
            "chromosome": "VARCHAR",
            "start_pos": "BIGINT",
            "end_pos": "BIGINT",
            "name": "VARCHAR",
        },
        genomic_column="position",
    )
    yield engine
    engine.close()


@pytest.fixture
def duckdb_stranded_engine(stranded_test_data_csv):
    """DuckDB engine with stranded test data loaded."""
    engine = GIQLEngine(target_dialect="duckdb", verbose=True)
    engine.load_csv("stranded_features", stranded_test_data_csv)
    engine.register_table_schema(
        "stranded_features",
        {
            "id": "INTEGER",
            "chromosome": "VARCHAR",
            "start_pos": "BIGINT",
            "end_pos": "BIGINT",
            "strand": "VARCHAR",
            "name": "VARCHAR",
        },
        genomic_column="position",
        strand_col="strand",
    )
    yield engine
    engine.close()


class TestCluster:
    """Tests for CLUSTER window function."""

    def test_basic_cluster(self, duckdb_cluster_engine):
        """Test basic CLUSTER operation."""
        result = duckdb_cluster_engine.query("""
            SELECT
                id,
                chromosome,
                start_pos,
                end_pos,
                name,
                CLUSTER(position) AS cluster_id
            FROM features
            ORDER BY chromosome, start_pos
        """)

        # Expected clusters:
        # chr1: features 1,2,3 are cluster 1 (overlapping/bookended)
        # chr1: feature 4 is cluster 2 (gap at 501)
        # chr2: feature 5 is cluster 1
        # chr2: feature 6 is cluster 2 (gap at 300)

        assert len(result) == 6

        # Check chr1 clusters
        chr1_results = result[result["chromosome"] == "chr1"]
        assert chr1_results.iloc[0]["cluster_id"] == chr1_results.iloc[1]["cluster_id"]
        assert chr1_results.iloc[1]["cluster_id"] == chr1_results.iloc[2]["cluster_id"]
        assert chr1_results.iloc[2]["cluster_id"] != chr1_results.iloc[3]["cluster_id"]

        # Check chr2 clusters
        chr2_results = result[result["chromosome"] == "chr2"]
        assert chr2_results.iloc[0]["cluster_id"] != chr2_results.iloc[1]["cluster_id"]

    def test_cluster_with_distance(self, duckdb_cluster_engine):
        """Test CLUSTER with distance parameter."""
        result = duckdb_cluster_engine.query("""
            SELECT
                id,
                chromosome,
                start_pos,
                end_pos,
                name,
                CLUSTER(position, 100) AS cluster_id
            FROM features
            ORDER BY chromosome, start_pos
        """)

        # With distance=100, chr1 features 1,2,3,4 should all be in same cluster
        # (gap of 1bp at position 501 is within 100bp tolerance)
        chr1_results = result[result["chromosome"] == "chr1"]
        cluster_ids = chr1_results["cluster_id"].tolist()
        assert len(set(cluster_ids)) == 1  # All in same cluster

    def test_stranded_cluster(self, duckdb_stranded_engine):
        """Test CLUSTER with stranded=true."""
        result = duckdb_stranded_engine.query("""
            SELECT
                id,
                chromosome,
                start_pos,
                end_pos,
                strand,
                name,
                CLUSTER(position, stranded=true) AS cluster_id
            FROM stranded_features
            ORDER BY chromosome, start_pos
        """)

        # Features should cluster only within the same strand:
        # + strand: f1,f2 overlap -> cluster 1, f5 is separate -> cluster 2
        # - strand: f3,f4 overlap -> cluster 1
        # Note: cluster_id numbering restarts for each partition (strand)

        assert len(result) == 5

        # Extract features
        f1 = result[result["id"] == 1].iloc[0]
        f2 = result[result["id"] == 2].iloc[0]
        f3 = result[result["id"] == 3].iloc[0]
        f4 = result[result["id"] == 4].iloc[0]
        f5 = result[result["id"] == 5].iloc[0]

        # Check that f1 and f2 (both +, overlapping) have same cluster_id
        assert f1["cluster_id"] == f2["cluster_id"]
        assert f1["strand"] == "+"
        assert f2["strand"] == "+"

        # Check that f3 and f4 (both -, overlapping) have same cluster_id
        assert f3["cluster_id"] == f4["cluster_id"]
        assert f3["strand"] == "-"
        assert f4["strand"] == "-"

        # Check that f5 (+ strand, separated) has different cluster from f1/f2
        assert f5["cluster_id"] != f1["cluster_id"]
        assert f5["strand"] == "+"

        # Verify stranded clustering works: compare with non-stranded
        result_nonstranded = duckdb_stranded_engine.query("""
            SELECT
                id,
                CLUSTER(position) AS cluster_id
            FROM stranded_features
            ORDER BY id
        """)

        # Without stranded, f1-f4 should all be in same cluster (overlapping)
        ns_f1 = result_nonstranded[result_nonstranded["id"] == 1].iloc[0]
        ns_f2 = result_nonstranded[result_nonstranded["id"] == 2].iloc[0]
        ns_f3 = result_nonstranded[result_nonstranded["id"] == 3].iloc[0]
        ns_f4 = result_nonstranded[result_nonstranded["id"] == 4].iloc[0]

        assert ns_f1["cluster_id"] == ns_f2["cluster_id"]
        assert ns_f2["cluster_id"] == ns_f3["cluster_id"]
        assert ns_f3["cluster_id"] == ns_f4["cluster_id"]


class TestMerge:
    """Tests for MERGE aggregate function."""

    def test_basic_merge(self, duckdb_cluster_engine):
        """Test basic MERGE operation."""
        result = duckdb_cluster_engine.query("""
            SELECT MERGE(position)
            FROM features
        """)

        # Expected merged intervals:
        # chr1: features 1,2,3 merge into [100, 500]
        # chr1: feature 4 stays as [501, 1000]
        # chr2: feature 5 stays as [100, 200]
        # chr2: feature 6 stays as [300, 400]

        assert len(result) == 4

        # Check chr1 merged intervals
        chr1_results = result[result["chromosome"] == "chr1"].sort_values("start_pos")
        assert len(chr1_results) == 2
        assert chr1_results.iloc[0]["start_pos"] == 100
        assert chr1_results.iloc[0]["end_pos"] == 500
        assert chr1_results.iloc[1]["start_pos"] == 501
        assert chr1_results.iloc[1]["end_pos"] == 1000

        # Check chr2 stays separate
        chr2_results = result[result["chromosome"] == "chr2"].sort_values("start_pos")
        assert len(chr2_results) == 2
        assert chr2_results.iloc[0]["start_pos"] == 100
        assert chr2_results.iloc[0]["end_pos"] == 200
        assert chr2_results.iloc[1]["start_pos"] == 300
        assert chr2_results.iloc[1]["end_pos"] == 400

    def test_merge_with_distance(self, duckdb_cluster_engine):
        """Test MERGE with distance parameter."""
        result = duckdb_cluster_engine.query("""
            SELECT MERGE(position, 100)
            FROM features
        """)

        # With distance=100, chr1 features 1-4 should merge into one interval
        chr1_results = result[result["chromosome"] == "chr1"]
        assert len(chr1_results) == 1
        assert chr1_results.iloc[0]["start_pos"] == 100
        assert chr1_results.iloc[0]["end_pos"] == 1000

    def test_merge_with_aggregation(self, duckdb_cluster_engine):
        """Test MERGE with additional aggregation columns."""
        result = duckdb_cluster_engine.query("""
            SELECT MERGE(position), COUNT(*) as feature_count
            FROM features
        """)

        # chr1 should have 2 merged intervals with counts
        chr1_results = result[result["chromosome"] == "chr1"].sort_values("start_pos")
        assert len(chr1_results) == 2
        assert chr1_results.iloc[0]["feature_count"] == 3  # f1, f2, f3 merged
        assert chr1_results.iloc[1]["feature_count"] == 1  # f4 alone

    def test_stranded_merge(self, duckdb_stranded_engine):
        """Test MERGE with stranded=true."""
        result = duckdb_stranded_engine.query("""
            SELECT MERGE(position, stranded=true)
            FROM stranded_features
        """)

        # + strand: f1,f2 merge -> [100,250], f5 stays -> [400,500]
        # - strand: f3,f4 merge -> [200,350]
        assert len(result) == 3

        plus_strand = result[result["strand"] == "+"].sort_values("start_pos")
        assert len(plus_strand) == 2
        assert plus_strand.iloc[0]["start_pos"] == 100
        assert plus_strand.iloc[0]["end_pos"] == 250
        assert plus_strand.iloc[1]["start_pos"] == 400
        assert plus_strand.iloc[1]["end_pos"] == 500

        minus_strand = result[result["strand"] == "-"]
        assert len(minus_strand) == 1
        assert minus_strand.iloc[0]["start_pos"] == 200
        assert minus_strand.iloc[0]["end_pos"] == 350
