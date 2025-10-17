"""
Pytest fixtures for integration tests.
"""

import pandas as pd
import pytest

from giql import GIQLEngine


@pytest.fixture(scope="session")
def to_df():
    """Fixture providing a helper to convert cursors to DataFrames.

    Returns a function that materializes cursor results for testing.
    Session-scoped since it's a pure function with no state.

    Usage:
        result = to_df(engine.execute("SELECT ..."))
    """

    def _to_df(cursor):
        if cursor.description:
            columns = [desc[0] for desc in cursor.description]
            return pd.DataFrame(cursor.fetchall(), columns=columns)
        return pd.DataFrame()

    return _to_df


@pytest.fixture
def sample_variants_csv(tmp_path):
    """Create sample variants CSV."""
    csv_content = """
    id,chromosome,start_pos,end_pos,ref,alt,quality
    1,chr1,1500,1600,A,T,30.0
    2,chr1,10500,10600,G,C,40.0
    3,chr1,15000,15100,T,A,25.0
    4,chr2,500,600,C,G,35.0
    5,chr2,5500,5600,A,T,20.0
    6,chr1,25000,25100,G,A,35.0
    7,chr2,15000,15100,T,C,28.0
    8,chr3,1000,1100,A,G,32.0
    """
    csv_path = tmp_path / "variants.csv"
    csv_path.write_text(csv_content.strip())
    return str(csv_path)


@pytest.fixture
def sample_genes_csv(tmp_path):
    """Create sample genes CSV."""
    csv_content = """
    gene_id,name,chromosome,start_pos,end_pos,strand
    1,GENE1,chr1,1000,2000,+
    2,GENE2,chr1,10000,11000,-
    3,GENE3,chr1,14000,16000,+
    4,GENE4,chr2,400,700,+
    5,GENE5,chr2,5000,6000,-
    """
    csv_path = tmp_path / "genes.csv"
    csv_path.write_text(csv_content.strip())
    return str(csv_path)


@pytest.fixture(params=["duckdb", "sqlite"])
def engine_with_variants(request, sample_variants_csv):
    """Create engine with loaded variants data for different dialects."""
    dialect = request.param

    engine = GIQLEngine(target_dialect=dialect, verbose=False)
    engine.load_csv("variants", sample_variants_csv)
    engine.register_table_schema(
        "variants",
        {
            "id": "INTEGER",
            "chromosome": "VARCHAR",
            "start_pos": "BIGINT",
            "end_pos": "BIGINT",
            "ref": "VARCHAR",
            "alt": "VARCHAR",
            "quality": "FLOAT",
        },
        genomic_column="position",
    )

    yield engine
    engine.close()


@pytest.fixture
def duckdb_engine_with_data(sample_variants_csv, sample_genes_csv):
    """DuckDB engine with both variants and genes loaded."""
    engine = GIQLEngine(target_dialect="duckdb", verbose=False)
    engine.load_csv("variants", sample_variants_csv)
    engine.load_csv("genes", sample_genes_csv)

    engine.register_table_schema(
        "variants",
        {
            "id": "INTEGER",
            "chromosome": "VARCHAR",
            "start_pos": "BIGINT",
            "end_pos": "BIGINT",
            "ref": "VARCHAR",
            "alt": "VARCHAR",
            "quality": "FLOAT",
        },
        genomic_column="position",
    )

    engine.register_table_schema(
        "genes",
        {
            "gene_id": "INTEGER",
            "name": "VARCHAR",
            "chromosome": "VARCHAR",
            "start_pos": "BIGINT",
            "end_pos": "BIGINT",
            "strand": "VARCHAR",
        },
        genomic_column="position",
    )

    yield engine
    engine.close()
