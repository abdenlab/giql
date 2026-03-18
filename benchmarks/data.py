"""Dataset generation and BED/Parquet conversion for benchmarks."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

HG38_CHROM_SIZES: dict[str, int] = {
    "chr1": 248_956_422,
    "chr2": 242_193_529,
    "chr3": 198_295_559,
    "chr4": 190_214_555,
    "chr5": 181_538_259,
    "chr6": 170_805_979,
    "chr7": 159_345_973,
    "chr8": 145_138_636,
    "chr9": 138_394_717,
    "chr10": 133_797_422,
    "chr11": 135_086_622,
    "chr12": 133_275_309,
    "chr13": 114_364_328,
    "chr14": 107_043_718,
    "chr15": 101_991_189,
    "chr16": 90_338_345,
    "chr17": 83_257_441,
    "chr18": 80_373_285,
    "chr19": 58_617_616,
    "chr20": 64_444_167,
    "chr21": 46_709_983,
    "chr22": 50_818_468,
    "chrX": 156_040_895,
}


@dataclass(frozen=True)
class GenomicDataset:
    name: str
    bed_path: Path
    parquet_path: Path
    n_intervals: int


def _chrom_sort_key(chrom: str) -> tuple[int, str]:
    """Natural chromosome ordering: chr1→1, chrX→23."""
    body = chrom[3:]  # strip "chr"
    try:
        return (int(body), "")
    except ValueError:
        return (23, body)


def _sort_records(
    chroms: list[str],
    starts: list[int],
    ends: list[int],
) -> tuple[list[str], list[int], list[int]]:
    """Sort genomic records by natural chrom order then start."""
    order = sorted(
        range(len(chroms)),
        key=lambda i: (_chrom_sort_key(chroms[i]), starts[i]),
    )
    return (
        [chroms[i] for i in order],
        [starts[i] for i in order],
        [ends[i] for i in order],
    )


def _shuffle_within_chroms(
    chroms: list[str],
    starts: list[int],
    ends: list[int],
    rng: np.random.Generator,
) -> tuple[list[str], list[int], list[int]]:
    """Group by chromosome in natural order but shuffle intervals within each.

    Chromosomes remain contiguous; only the start-position ordering within
    each chromosome is destroyed.
    """
    from collections import defaultdict

    groups: dict[str, list[int]] = defaultdict(list)
    for i, c in enumerate(chroms):
        groups[c].append(i)

    result_chroms: list[str] = []
    result_starts: list[int] = []
    result_ends: list[int] = []
    for c in sorted(groups, key=_chrom_sort_key):
        indices = groups[c]
        rng.shuffle(indices)
        for i in indices:
            result_chroms.append(chroms[i])
            result_starts.append(starts[i])
            result_ends.append(ends[i])

    return result_chroms, result_starts, result_ends


def _write_bed(bed_path: Path, chroms: list[str], starts: list[int], ends: list[int]) -> None:
    with open(bed_path, "w") as f:
        for c, s, e in zip(chroms, starts, ends):
            f.write(f"{c}\t{s}\t{e}\n")


def _write_parquet(
    parquet_path: Path,
    chroms: list[str],
    starts: list[int],
    ends: list[int],
    *,
    row_group_per_chrom: bool = False,
) -> None:
    """Write a Parquet file.

    Data must already be sorted by chromosome (natural order) then start.

    When row_group_per_chrom is True, each chromosome becomes its own row
    group so that DataFusion can use min/max statistics for row-group
    pruning. When False (default), all data goes into a single row group
    — worst case for pruning but representative of naive file layouts.
    """
    schema = pa.schema([
        ("chrom", pa.string()),
        ("start", pa.int32()),
        ("end", pa.int32()),
    ])

    if not row_group_per_chrom:
        table = pa.table(
            {
                "chrom": pa.array(chroms, type=pa.string()),
                "start": pa.array(starts, type=pa.int32()),
                "end": pa.array(ends, type=pa.int32()),
            }
        )
        pq.write_table(table, parquet_path)
        return

    writer = pq.ParquetWriter(parquet_path, schema)
    n = len(chroms)
    rg_start = 0
    while rg_start < n:
        current_chrom = chroms[rg_start]
        rg_end = rg_start + 1
        while rg_end < n and chroms[rg_end] == current_chrom:
            rg_end += 1
        batch = pa.table(
            {
                "chrom": pa.array(chroms[rg_start:rg_end], type=pa.string()),
                "start": pa.array(starts[rg_start:rg_end], type=pa.int32()),
                "end": pa.array(ends[rg_start:rg_end], type=pa.int32()),
            }
        )
        writer.write_table(batch)
        rg_start = rg_end
    writer.close()


def generate_dataset(
    n: int,
    workdir: Path,
    name: str | None = None,
    seed: int = 42,
    *,
    row_group_per_chrom: bool = False,
    unsorted: bool = False,
) -> GenomicDataset:
    """Generate a synthetic genomic dataset with n intervals.

    Intervals are distributed proportionally to chromosome size, with
    log-normal lengths (median ~500bp, clipped to [50, 50_000]).

    When unsorted is False (default), intervals are fully sorted by
    chromosome then start. When True, chromosomes remain contiguous
    but intervals within each chromosome are shuffled.
    """
    rng = np.random.default_rng(seed)
    if name is None:
        name = f"synth_{n}"

    chrom_names = list(HG38_CHROM_SIZES.keys())
    chrom_sizes = np.array([HG38_CHROM_SIZES[c] for c in chrom_names], dtype=np.float64)
    probs = chrom_sizes / chrom_sizes.sum()

    chrom_indices = rng.choice(len(chrom_names), size=n, p=probs)
    lengths = np.clip(
        rng.lognormal(mean=math.log(500), sigma=1.0, size=n).astype(int),
        50,
        50_000,
    )

    chroms: list[str] = []
    starts: list[int] = []
    ends: list[int] = []
    for idx, length in zip(chrom_indices, lengths):
        c = chrom_names[idx]
        size = int(chrom_sizes[idx])
        max_start = max(0, size - length)
        start = int(rng.integers(0, max_start + 1))
        chroms.append(c)
        starts.append(start)
        ends.append(start + length)

    # Always sort first so _shuffle_within_chroms has grouped chromosomes
    chroms, starts, ends = _sort_records(chroms, starts, ends)
    if unsorted:
        chroms, starts, ends = _shuffle_within_chroms(chroms, starts, ends, rng)

    bed_path = workdir / f"{name}.bed"
    parquet_path = workdir / f"{name}.parquet"
    _write_bed(bed_path, chroms, starts, ends)
    _write_parquet(parquet_path, chroms, starts, ends, row_group_per_chrom=row_group_per_chrom)

    return GenomicDataset(name=name, bed_path=bed_path, parquet_path=parquet_path, n_intervals=n)


def load_from_bed(
    bed_path: Path,
    workdir: Path,
    name: str | None = None,
    *,
    row_group_per_chrom: bool = False,
) -> GenomicDataset:
    """Load a BED file, sort it, and write sorted BED + Parquet to workdir."""
    if name is None:
        name = bed_path.stem

    chroms: list[str] = []
    starts: list[int] = []
    ends: list[int] = []
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            chroms.append(parts[0])
            starts.append(int(parts[1]))
            ends.append(int(parts[2]))

    chroms, starts, ends = _sort_records(chroms, starts, ends)

    sorted_bed = workdir / f"{name}.bed"
    parquet_path = workdir / f"{name}.parquet"
    _write_bed(sorted_bed, chroms, starts, ends)
    _write_parquet(parquet_path, chroms, starts, ends, row_group_per_chrom=row_group_per_chrom)

    return GenomicDataset(
        name=name,
        bed_path=sorted_bed,
        parquet_path=parquet_path,
        n_intervals=len(chroms),
    )


def load_from_parquet(
    parquet_path: Path,
    workdir: Path,
    name: str | None = None,
    *,
    row_group_per_chrom: bool = False,
) -> GenomicDataset:
    """Load a Parquet file, sort by chrom/start, write sorted BED + Parquet to workdir."""
    if name is None:
        name = parquet_path.stem

    table = pq.read_table(parquet_path, columns=["chrom", "start", "end"])
    chroms = table["chrom"].to_pylist()
    starts = table["start"].to_pylist()
    ends = table["end"].to_pylist()

    chroms, starts, ends = _sort_records(chroms, starts, ends)

    sorted_bed = workdir / f"{name}.bed"
    sorted_parquet = workdir / f"{name}.parquet"
    _write_bed(sorted_bed, chroms, starts, ends)
    _write_parquet(sorted_parquet, chroms, starts, ends, row_group_per_chrom=row_group_per_chrom)

    return GenomicDataset(
        name=name,
        bed_path=sorted_bed,
        parquet_path=sorted_parquet,
        n_intervals=len(chroms),
    )
