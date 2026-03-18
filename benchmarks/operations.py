"""Operation registry for GIQL vs bedtools benchmarks."""

from __future__ import annotations

from dataclasses import dataclass

QUERY_REGION = "chr1:1000000-2000000"
QUERY_REGION_BED = ("chr1", 1_000_000, 2_000_000)


@dataclass(frozen=True)
class Operation:
    name: str
    label: str
    needs_secondary: bool
    datafusion_supported: bool


ALL_OPS: dict[str, Operation] = {
    "intersect_filter": Operation(
        name="intersect_filter",
        label="INTERSECTS filter",
        needs_secondary=False,
        datafusion_supported=True,
    ),
    "intersect_join": Operation(
        name="intersect_join",
        label="INTERSECTS join",
        needs_secondary=True,
        datafusion_supported=True,
    ),
    "intersect_pairs": Operation(
        name="intersect_pairs",
        label="INTERSECTS pairs",
        needs_secondary=True,
        datafusion_supported=True,
    ),
    "merge": Operation(
        name="merge",
        label="MERGE",
        needs_secondary=False,
        datafusion_supported=True,
    ),
    "cluster": Operation(
        name="cluster",
        label="CLUSTER",
        needs_secondary=False,
        datafusion_supported=True,
    ),
}
