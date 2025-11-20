"""
Test Interface Contracts for Bedtools Integration Test Suite

This file defines the public interfaces (protocols/abstract base classes) that
components in the test suite must implement. These contracts ensure consistency
and enable testing of test infrastructure itself.
"""

from abc import ABC
from abc import abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Protocol
from typing import Tuple

# ============================================================================
# Data Generation Contracts
# ============================================================================


class IntervalGeneratorProtocol(Protocol):
    """Protocol for generating simulated genomic datasets."""

    def generate_overlapping(
        self, chromosome: str, count: int
    ) -> List[Tuple[str, int, int, str, int, str]]:
        """Generate intervals with controlled overlap.

        Args:
            chromosome: Chromosome name
            count: Number of intervals to generate

        Returns:
            List of (chrom, start, end, name, score, strand) tuples
        """
        ...

    def generate_adjacent(
        self, chromosome: str, count: int
    ) -> List[Tuple[str, int, int, str, int, str]]:
        """Generate adjacent non-overlapping intervals.

        Args:
            chromosome: Chromosome name
            count: Number of intervals to generate

        Returns:
            List of (chrom, start, end, name, score, strand) tuples
        """
        ...

    def generate_separated(
        self, chromosome: str, count: int, min_gap: int
    ) -> List[Tuple[str, int, int, str, int, str]]:
        """Generate intervals with minimum gap between them.

        Args:
            chromosome: Chromosome name
            count: Number of intervals to generate
            min_gap: Minimum gap size between intervals

        Returns:
            List of (chrom, start, end, name, score, strand) tuples
        """
        ...

    def generate_multi_chromosome(
        self, chromosome_count: int, intervals_per_chrom: int
    ) -> List[Tuple[str, int, int, str, int, str]]:
        """Generate intervals across multiple chromosomes.

        Args:
            chromosome_count: Number of chromosomes
            intervals_per_chrom: Intervals per chromosome

        Returns:
            List of (chrom, start, end, name, score, strand) tuples
        """
        ...


# ============================================================================
# Bedtools Integration Contracts
# ============================================================================


class BedtoolsWrapperProtocol(Protocol):
    """Protocol for executing bedtools commands."""

    def check_version(self) -> Tuple[int, int, int]:
        """Check installed bedtools version.

        Returns:
            (major, minor, patch) version tuple

        Raises:
            BedtoolsNotFoundError: If bedtools not in PATH
        """
        ...

    def intersect(
        self,
        file_a: Path,
        file_b: Path,
        same_strand: bool = False,
        opposite_strand: bool = False,
    ) -> str:
        """Execute bedtools intersect command.

        Args:
            file_a: First BED file path
            file_b: Second BED file path
            same_strand: Only intersect same strand (-s flag)
            opposite_strand: Only intersect opposite strand (-S flag)

        Returns:
            Bedtools stdout output

        Raises:
            BedtoolsExecutionError: If command fails
        """
        ...

    def merge(self, file: Path, strand_specific: bool = False) -> str:
        """Execute bedtools merge command.

        Args:
            file: BED file path
            strand_specific: Merge per-strand (-s flag)

        Returns:
            Bedtools stdout output

        Raises:
            BedtoolsExecutionError: If command fails
        """
        ...

    def closest(
        self,
        file_a: Path,
        file_b: Path,
        report_distance: bool = False,
        same_strand: bool = False,
    ) -> str:
        """Execute bedtools closest command.

        Args:
            file_a: First BED file path
            file_b: Second BED file path
            report_distance: Report distance (-d flag)
            same_strand: Only report same strand (-s flag)

        Returns:
            Bedtools stdout output

        Raises:
            BedtoolsExecutionError: If command fails
        """
        ...


# ============================================================================
# Result Comparison Contracts
# ============================================================================


@dataclass
class ComparisonResult:
    """Result of comparing GIQL and bedtools outputs."""

    match: bool
    giql_row_count: int
    bedtools_row_count: int
    differences: List[str]
    epsilon: float = 1e-9

    def __bool__(self) -> bool:
        """Allow direct boolean evaluation."""
        return self.match

    def failure_message(self) -> str:
        """Generate detailed failure message."""
        if self.match:
            return "✓ Results match"

        msg = [
            f"✗ Results do not match",
            f"  GIQL rows: {self.giql_row_count}",
            f"  Bedtools rows: {self.bedtools_row_count}",
            f"  Differences:",
        ]
        for diff in self.differences[:10]:  # Limit to first 10 differences
            msg.append(f"    - {diff}")
        if len(self.differences) > 10:
            msg.append(f"    ... and {len(self.differences) - 10} more")
        return "\n".join(msg)


class ResultComparatorProtocol(Protocol):
    """Protocol for comparing GIQL and bedtools results."""

    def compare(
        self, giql_rows: List[Tuple], bedtools_rows: List[Tuple], epsilon: float = 1e-9
    ) -> ComparisonResult:
        """Compare result sets with appropriate tolerance.

        Args:
            giql_rows: Rows from GIQL query execution
            bedtools_rows: Rows from bedtools output
            epsilon: Tolerance for floating-point comparisons

        Returns:
            ComparisonResult with match status and differences
        """
        ...

    def parse_bedtools_output(self, output: str) -> List[Tuple]:
        """Parse bedtools command output into structured rows.

        Args:
            output: Raw bedtools stdout

        Returns:
            List of tuples matching BED format
        """
        ...


# ============================================================================
# Database Fixture Contracts
# ============================================================================


class DatabaseFixtureProtocol(Protocol):
    """Protocol for database test fixtures."""

    def get_connection(self) -> Any:
        """Get clean DuckDB connection.

        Returns:
            DuckDB connection to in-memory database
        """
        ...

    def load_intervals(
        self,
        connection: Any,
        table_name: str,
        intervals: List[Tuple[str, int, int, str, int, str]],
    ) -> None:
        """Load intervals into database table.

        Args:
            connection: DuckDB connection
            table_name: Table name to create
            intervals: List of (chrom, start, end, name, score, strand) tuples
        """
        ...

    def export_to_bed(
        self, connection: Any, table_name: str, output_path: Path, format: str = "bed6"
    ) -> None:
        """Export table to BED file.

        Args:
            connection: DuckDB connection
            table_name: Table to export
            output_path: BED file output path
            format: BED format ('bed3' or 'bed6')
        """
        ...


# ============================================================================
# Test Case Contracts
# ============================================================================


class TestScenarioProtocol(Protocol):
    """Protocol for test scenario definitions."""

    @property
    def name(self) -> str:
        """Scenario name for test identification."""
        ...

    @property
    def description(self) -> str:
        """Human-readable scenario description."""
        ...

    def setup_data(self, generator: IntervalGeneratorProtocol) -> Dict[str, List[Tuple]]:
        """Generate test data for this scenario.

        Args:
            generator: Interval generator instance

        Returns:
            Dict mapping table names to interval lists
        """
        ...

    def get_giql_query(self) -> str:
        """Get GIQL query for this scenario.

        Returns:
            GIQL query string
        """
        ...

    def get_bedtools_command(self, bed_files: Dict[str, Path]) -> List[str]:
        """Get bedtools command for this scenario.

        Args:
            bed_files: Dict mapping table names to BED file paths

        Returns:
            Bedtools command as list of arguments
        """
        ...


# ============================================================================
# Error Handling Contracts
# ============================================================================


class BedtoolsNotFoundError(Exception):
    """Raised when bedtools binary not found in PATH."""

    pass


class BedtoolsVersionError(Exception):
    """Raised when bedtools version is incompatible."""

    pass


class BedtoolsExecutionError(Exception):
    """Raised when bedtools command execution fails."""

    pass


class InvalidIntervalError(Exception):
    """Raised when interval validation fails."""

    pass


class ComparisonError(Exception):
    """Raised when result comparison encounters unexpected state."""

    pass
