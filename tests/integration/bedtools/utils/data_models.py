"""Data models for bedtools integration testing.

- GenomicInterval: Represents a single genomic interval
- ComparisonResult: Result of comparing GIQL vs bedtools outputs
"""

from dataclasses import dataclass
from dataclasses import field


@dataclass
class GenomicInterval:
    """Represents a single genomic interval with all BED file fields."""

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
class ComparisonResult:
    """Result of comparing GIQL and bedtools outputs."""

    match: bool
    giql_row_count: int
    bedtools_row_count: int
    differences: list[str] = field(default_factory=list)
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
