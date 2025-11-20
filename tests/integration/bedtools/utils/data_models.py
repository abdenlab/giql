"""Data models for bedtools integration testing.

This module defines the core data structures used throughout the test suite:
- GenomicInterval: Represents a single genomic interval
- SimulatedDataset: Collection of intervals for testing
- ComparisonResult: Result of comparing GIQL vs bedtools outputs
- IntervalGeneratorConfig: Configuration for dataset generation
- BedtoolsVersion: Bedtools version information
"""

import re
from dataclasses import dataclass
from dataclasses import field
from pathlib import Path
from typing import Dict
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
        """Validate interval fields."""
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

    def to_bed_line(self, format="bed6") -> str:
        """Convert to BED format line.

        Args:
            format: Output format ('bed3' or 'bed6')

        Returns:
            Tab-separated BED format string
        """
        if format == "bed3":
            return f"{self.chrom}\t{self.start}\t{self.end}"
        elif format == "bed6":
            name = self.name or "."
            score = self.score if self.score is not None else 0
            strand = self.strand or "."
            return f"{self.chrom}\t{self.start}\t{self.end}\t{name}\t{score}\t{strand}"
        else:
            raise ValueError(f"Unsupported BED format: {format}")


@dataclass
class SimulatedDataset:
    """Collection of genomic intervals with controlled properties for testing.

    Attributes:
        name: Dataset identifier (e.g., "intervals_a", "intervals_b")
        intervals: List of genomic intervals
        scenario_type: Scenario descriptor (e.g., "overlapping", "adjacent")
        metadata: Generation parameters (seed, chromosome_count, etc.)
    """

    name: str
    intervals: List[GenomicInterval]
    scenario_type: str
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        """Validate dataset has at least one interval."""
        if len(self.intervals) == 0:
            raise ValueError("Dataset must contain at least one interval")

    def to_bed_file(self, path: Path, format="bed6"):
        """Export to BED file.

        Args:
            path: Output file path
            format: BED format ('bed3' or 'bed6')
        """
        with open(path, "w") as f:
            for interval in self.intervals:
                f.write(interval.to_bed_line(format) + "\n")

    def to_duckdb_table(self, conn, table_name: str):
        """Load into DuckDB table.

        Args:
            conn: DuckDB connection
            table_name: Name of table to create
        """
        rows = [
            (i.chrom, i.start, i.end, i.name, i.score, i.strand) for i in self.intervals
        ]
        conn.execute(f"""
            CREATE TABLE {table_name} (
                chrom VARCHAR,
                start INTEGER,
                end INTEGER,
                name VARCHAR,
                score INTEGER,
                strand VARCHAR
            )
        """)
        conn.executemany(f"INSERT INTO {table_name} VALUES (?,?,?,?,?,?)", rows)


@dataclass
class ComparisonResult:
    """Result of comparing GIQL and bedtools outputs.

    Attributes:
        match: Whether results match
        giql_row_count: Number of rows from GIQL query
        bedtools_row_count: Number of rows from bedtools output
        differences: Specific differences found (if match=False)
        comparison_metadata: Epsilon used, sort order, etc.
    """

    match: bool
    giql_row_count: int
    bedtools_row_count: int
    differences: List[str] = field(default_factory=list)
    comparison_metadata: dict = field(default_factory=dict)

    def __bool__(self) -> bool:
        """Allow direct boolean evaluation in assertions."""
        return self.match

    def failure_message(self) -> str:
        """Generate detailed failure message for test output.

        Returns:
            Formatted failure message with differences
        """
        if self.match:
            return "✓ Results match"

        msg = [
            f"✗ Results do not match",
            f"  GIQL rows: {self.giql_row_count}",
            f"  Bedtools rows: {self.bedtools_row_count}",
        ]

        if self.differences:
            msg.append("  Differences:")
            for diff in self.differences[:10]:  # Limit to first 10
                msg.append(f"    - {diff}")
            if len(self.differences) > 10:
                msg.append(f"    ... and {len(self.differences) - 10} more")

        return "\n".join(msg)


@dataclass
class IntervalGeneratorConfig:
    """Configuration for simulated dataset generation.

    Attributes:
        chromosome_count: Number of chromosomes to generate
        intervals_per_chromosome: Intervals per chromosome
        min_interval_size: Minimum interval length
        max_interval_size: Maximum interval length
        overlap_probability: Probability of overlap (0.0-1.0)
        strand_distribution: Proportions of +/-/. strands
        seed: Random seed for reproducibility
    """

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
        """Validate configuration parameters."""
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


@dataclass
class BedtoolsVersion:
    """Represents bedtools version information.

    Attributes:
        major: Major version number
        minor: Minor version number
        patch: Patch version number
        raw_version_string: Original version string from bedtools
    """

    major: int
    minor: int
    patch: int
    raw_version_string: str

    def is_compatible(self) -> bool:
        """Check if version meets minimum requirement (2.30.0).

        Returns:
            True if version >= 2.30.0
        """
        return (self.major, self.minor, self.patch) >= (2, 30, 0)

    def __str__(self) -> str:
        """Return version as string."""
        return f"{self.major}.{self.minor}.{self.patch}"

    @classmethod
    def from_string(cls, version_str: str) -> "BedtoolsVersion":
        """Parse version from bedtools --version output.

        Args:
            version_str: Version string from bedtools (e.g., "bedtools v2.30.0")

        Returns:
            BedtoolsVersion instance

        Raises:
            ValueError: If version string cannot be parsed
        """
        match = re.search(r"v?(\d+)\.(\d+)\.(\d+)", version_str)
        if not match:
            raise ValueError(f"Could not parse version from: {version_str}")
        major, minor, patch = map(int, match.groups())
        return cls(major, minor, patch, version_str)
