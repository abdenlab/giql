"""Per-operation timing functions for DataFusion, polars-bio, and bedtools.

All backend imports are lazy (inside function bodies).
"""

from __future__ import annotations

from pathlib import Path

from benchmarks.data import GenomicDataset
from benchmarks.operations import QUERY_REGION, QUERY_REGION_BED


def make_datafusion_engine(
    primary: GenomicDataset,
    secondary: GenomicDataset | None = None,
):
    """Create and configure a DataFusionEngine. Not timed.

    Registers primary as "peaks" and secondary (if given) as "b".
    """
    from giql.engines.datafusion import DataFusionEngine
    from giql.table import Table

    engine = DataFusionEngine()
    engine.register_parquet("peaks", primary.parquet_path, Table("peaks", strand_col=None))
    if secondary is not None:
        engine.register_parquet("b", secondary.parquet_path, Table("b", strand_col=None))
    return engine


def run_datafusion_intersect_filter(engine) -> None:
    engine.query(f"SELECT * FROM peaks WHERE interval INTERSECTS '{QUERY_REGION}'")


def run_datafusion_intersect_join(engine) -> None:
    engine.query(
        "SELECT DISTINCT peaks.* FROM peaks JOIN b ON peaks.interval INTERSECTS b.interval"
    )


def run_datafusion_merge(engine) -> None:
    engine.query("SELECT MERGE(interval) FROM peaks")


def run_datafusion_merge_unsorted(engine) -> None:
    engine.query(
        "SELECT MERGE(interval) FROM (SELECT * FROM peaks ORDER BY chrom, start) AS t"
    )


def run_datafusion_cluster(engine) -> None:
    engine.query("SELECT *, CLUSTER(interval) AS cid FROM peaks")


def run_datafusion_cluster_unsorted(engine) -> None:
    engine.query(
        "SELECT *, CLUSTER(interval) AS cid "
        "FROM (SELECT * FROM peaks ORDER BY chrom, start) AS t"
    )


def setup_polarsbio() -> None:
    """Configure polars-bio for 0-based half-open coordinates. Call once."""
    import polars_bio as pb

    pb.set_option("datafusion.bio.coordinate_system_zero_based", True)
    pb.set_option("datafusion.bio.show_progress", False)


def run_polarsbio_intersect_filter(parquet_path: Path, query_parquet: Path) -> None:
    import polars_bio as pb

    pb.overlap(
        str(parquet_path), str(query_parquet),
    ).select(["chrom_1", "start_1", "end_1"]).unique().collect()


def run_polarsbio_intersect_join(parquet_a: Path, parquet_b: Path) -> None:
    import polars_bio as pb

    pb.overlap(
        str(parquet_a), str(parquet_b),
    ).select(["chrom_1", "start_1", "end_1"]).unique().collect()


def run_polarsbio_merge(parquet_path: Path) -> None:
    import polars_bio as pb

    pb.merge(str(parquet_path)).collect()


def run_polarsbio_cluster(parquet_path: Path) -> None:
    import polars_bio as pb

    pb.cluster(str(parquet_path)).collect()


def run_bedtools_intersect_filter(bed_path: Path, query_bed_path: Path) -> None:
    import pybedtools

    a = pybedtools.BedTool(str(bed_path))
    b = pybedtools.BedTool(str(query_bed_path))
    a.intersect(b, u=True).saveas()


def run_bedtools_intersect_join(bed_a: Path, bed_b: Path) -> None:
    import pybedtools

    a = pybedtools.BedTool(str(bed_a))
    b = pybedtools.BedTool(str(bed_b))
    a.intersect(b, u=True).saveas()


def run_bedtools_merge(bed_path: Path) -> None:
    import pybedtools

    pybedtools.BedTool(str(bed_path)).merge().saveas()


def run_bedtools_merge_unsorted(bed_path: Path) -> None:
    import pybedtools

    pybedtools.BedTool(str(bed_path)).sort().merge().saveas()


def run_bedtools_cluster(bed_path: Path) -> None:
    import pybedtools

    pybedtools.BedTool(str(bed_path)).cluster().saveas()


def run_bedtools_cluster_unsorted(bed_path: Path) -> None:
    import pybedtools

    pybedtools.BedTool(str(bed_path)).sort().cluster().saveas()
