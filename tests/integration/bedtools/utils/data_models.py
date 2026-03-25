"""Data models for bedtools integration testing.

This module defines the core data structures used throughout the test suite:
- GenomicInterval: Represents a single genomic interval
- SimulatedDataset: Collection of intervals for testing
- ComparisonResult: Result of comparing GIQL vs bedtools outputs
- IntervalGeneratorConfig: Configuration for dataset generation
"""

from dataclasses import dataclass
from dataclasses import field
from typing import List


@dataclass
class GenomicInterval:
    """Represents a single genomic interval with all BED file fields.

    Attributes:
        chrom: Chromosome name (e.g., "chr1", "chr2", "chrX")
        start: Start position (0-based, inclusive)
        end: End position (0-based, exclusive)
        name: Optional interval name/identifier
        score: Optional score value (0-1000)
        strand: Optional strand ("+", "-", or ".")
    """

    chrom: str
    start: int
    end: int
    name: str | None = None
    score: int | None = None
    strand: str | None = None

    def __post_init__(self):
        if self.start >= self.end:
            raise ValueError(
                f"Invalid interval: start ({self.start}) >= end ({self.end})"
            )
        if self.start < 0:
            raise ValueError(f"Invalid interval: start ({self.start}) < 0")
        if self.strand and self.strand not in ["+", "-", "."]:
            raise ValueError(f"Invalid strand: {self.strand}")
        if self.score is not None and not (0 <= self.score <= 1000):
            raise ValueError(f"Invalid score: {self.score}")

    def to_tuple(self) -> tuple:
        return (
            self.chrom,
            self.start,
            self.end,
            self.name,
            self.score,
            self.strand,
        )


@dataclass
class SimulatedDataset:
    """Collection of genomic intervals with controlled properties for testing."""

    name: str
    intervals: List[GenomicInterval]
    scenario_type: str
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        if len(self.intervals) == 0:
            raise ValueError("Dataset must contain at least one interval")


@dataclass
class ComparisonResult:
    """Result of comparing GIQL and bedtools outputs."""

    match: bool
    giql_row_count: int
    bedtools_row_count: int
    differences: List[str] = field(default_factory=list)
    comparison_metadata: dict = field(default_factory=dict)

    def __bool__(self) -> bool:
        return self.match

    def failure_message(self) -> str:
        if self.match:
            return "Results match"

        msg = [
            "Results do not match",
            f"  GIQL rows: {self.giql_row_count}",
            f"  Bedtools rows: {self.bedtools_row_count}",
        ]

        if self.differences:
            msg.append("  Differences:")
            for diff in self.differences[:10]:
                msg.append(f"    - {diff}")
            if len(self.differences) > 10:
                msg.append(f"    ... and {len(self.differences) - 10} more")

        return "\n".join(msg)


@dataclass
class IntervalGeneratorConfig:
    """Configuration for simulated dataset generation."""

    chromosome_count: int = 3
    intervals_per_chromosome: int = 100
    min_interval_size: int = 100
    max_interval_size: int = 1000
    overlap_probability: float = 0.3
    strand_distribution: dict = field(
        default_factory=lambda: {"+": 0.45, "-": 0.45, ".": 0.1}
    )
    seed: int = 42

    def __post_init__(self):
        if self.chromosome_count <= 0:
            raise ValueError("chromosome_count must be > 0")
        if self.intervals_per_chromosome <= 0:
            raise ValueError("intervals_per_chromosome must be > 0")
        if self.min_interval_size < 1:
            raise ValueError("min_interval_size must be >= 1")
        if self.max_interval_size < self.min_interval_size:
            raise ValueError("max_interval_size must be >= min_interval_size")
        if not (0.0 <= self.overlap_probability <= 1.0):
            raise ValueError("overlap_probability must be in [0.0, 1.0]")
        if abs(sum(self.strand_distribution.values()) - 1.0) > 1e-6:
            raise ValueError("strand_distribution must sum to 1.0")
