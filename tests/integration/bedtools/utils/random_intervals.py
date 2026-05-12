"""Deterministic random-interval generator for bedtools integration tests."""

import random

from .data_models import GenomicInterval


def generate_random_intervals(
    *,
    seed: int,
    prefix: str,
    count_per_chrom: int = 30,
    n_chroms: int = 3,
    start_max: int = 100_000,
    min_size: int = 100,
    max_size: int = 1000,
    strand: str = "+",
) -> list[GenomicInterval]:
    """Generate a deterministic list of GenomicInterval samples.

    Used by scale tests to produce realistic multi-chromosome input
    sets without duplicating the same random-loop boilerplate. The
    seed determines the exact sample — callers expecting identical
    bedtools and GIQL outputs must pass the same seed to both sides.
    """
    rng = random.Random(seed)
    intervals: list[GenomicInterval] = []
    for chrom_num in range(1, n_chroms + 1):
        for i in range(count_per_chrom):
            start = rng.randint(0, start_max)
            size = rng.randint(min_size, max_size)
            intervals.append(
                GenomicInterval(
                    f"chr{chrom_num}",
                    start,
                    start + size,
                    f"{prefix}_{chrom_num}_{i}",
                    0,
                    strand,
                )
            )
    return intervals
