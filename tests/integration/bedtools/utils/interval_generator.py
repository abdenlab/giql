"""Interval generator for creating simulated genomic datasets.

This module provides the IntervalGenerator class for creating test datasets
with controlled properties (overlap density, strand distribution, etc.).
"""

import random
from typing import List
from typing import Tuple

from .data_models import GenomicInterval
from .data_models import IntervalGeneratorConfig
from .data_models import SimulatedDataset


class IntervalGenerator:
    """Generate simulated genomic intervals for testing.

    Provides methods for generating intervals with various patterns:
    - Overlapping intervals
    - Adjacent intervals
    - Separated intervals
    - Multi-chromosome datasets
    - Strand-specific datasets
    """

    def __init__(self, config: IntervalGeneratorConfig | None = None):
        """Initialize interval generator.

        Args:
            config: Generator configuration (uses defaults if None)
        """
        self.config = config or IntervalGeneratorConfig()
        self.rng = random.Random(self.config.seed)

    def _choose_strand(self) -> str:
        """Choose strand based on configured distribution.

        Returns:
            Strand ('+', '-', or '.')
        """
        r = self.rng.random()
        cumulative = 0.0
        for strand, prob in self.config.strand_distribution.items():
            cumulative += prob
            if r <= cumulative:
                return strand
        return "."  # Fallback

    def _generate_interval_size(self) -> int:
        """Generate random interval size within configured range.

        Returns:
            Interval size in base pairs
        """
        return self.rng.randint(
            self.config.min_interval_size, self.config.max_interval_size
        )

    def generate_basic(
        self, chromosome: str, count: int, max_position: int = 1000000
    ) -> List[GenomicInterval]:
        """Generate basic random intervals.

        Args:
            chromosome: Chromosome name
            count: Number of intervals to generate
            max_position: Maximum chromosome position

        Returns:
            List of genomic intervals
        """
        intervals = []
        for i in range(count):
            size = self._generate_interval_size()
            start = self.rng.randint(0, max_position - size)
            end = start + size
            strand = self._choose_strand()

            intervals.append(
                GenomicInterval(
                    chrom=chromosome,
                    start=start,
                    end=end,
                    name=f"interval_{i}",
                    score=self.rng.randint(0, 1000),
                    strand=strand,
                )
            )

        return intervals

    def generate_dataset(
        self,
        name: str,
        scenario_type: str = "basic",
        chromosome_count: int | None = None,
        intervals_per_chrom: int | None = None,
    ) -> SimulatedDataset:
        """Generate a complete simulated dataset.

        Args:
            name: Dataset identifier
            scenario_type: Type of scenario ("basic", "overlapping", etc.)
            chromosome_count: Number of chromosomes (uses config default if None)
            intervals_per_chrom: Intervals per chromosome (uses config default if None)

        Returns:
            SimulatedDataset with generated intervals
        """
        chrom_count = chromosome_count or self.config.chromosome_count
        interval_count = intervals_per_chrom or self.config.intervals_per_chromosome

        all_intervals = []
        for i in range(chrom_count):
            chrom_name = f"chr{i + 1}"
            intervals = self.generate_basic(chrom_name, interval_count)
            all_intervals.extend(intervals)

        return SimulatedDataset(
            name=name,
            intervals=all_intervals,
            scenario_type=scenario_type,
            metadata={
                "chromosome_count": chrom_count,
                "intervals_per_chromosome": interval_count,
                "seed": self.config.seed,
                "total_intervals": len(all_intervals),
            },
        )

    def generate_overlapping_scenarios(
        self, chromosome: str, count: int, overlap_size: int = 50
    ) -> List[GenomicInterval]:
        """Generate overlapping intervals with controlled overlap.

        Args:
            chromosome: Chromosome name
            count: Number of intervals to generate
            overlap_size: Size of overlap between adjacent intervals

        Returns:
            List of overlapping genomic intervals
        """
        intervals = []
        base_size = self.config.min_interval_size
        current_start = 100

        for i in range(count):
            start = current_start
            end = start + base_size
            strand = self._choose_strand()

            intervals.append(
                GenomicInterval(
                    chrom=chromosome,
                    start=start,
                    end=end,
                    name=f"overlap_{i}",
                    score=self.rng.randint(0, 1000),
                    strand=strand,
                )
            )

            # Next interval starts before current ends (creating overlap)
            current_start = end - overlap_size

        return intervals

    def generate_adjacent_scenarios(
        self, chromosome: str, count: int
    ) -> List[GenomicInterval]:
        """Generate adjacent intervals (touching but not overlapping).

        Args:
            chromosome: Chromosome name
            count: Number of intervals to generate

        Returns:
            List of adjacent genomic intervals
        """
        intervals = []
        base_size = self.config.min_interval_size
        current_start = 100

        for i in range(count):
            start = current_start
            end = start + base_size
            strand = self._choose_strand()

            intervals.append(
                GenomicInterval(
                    chrom=chromosome,
                    start=start,
                    end=end,
                    name=f"adjacent_{i}",
                    score=self.rng.randint(0, 1000),
                    strand=strand,
                )
            )

            # Next interval starts exactly where current ends
            current_start = end

        return intervals

    def generate_separated_scenarios(
        self, chromosome: str, count: int, gap_size: int = 100
    ) -> List[GenomicInterval]:
        """Generate separated intervals with gaps between them.

        Args:
            chromosome: Chromosome name
            count: Number of intervals to generate
            gap_size: Size of gap between intervals

        Returns:
            List of separated genomic intervals
        """
        intervals = []
        base_size = self.config.min_interval_size
        current_start = 100

        for i in range(count):
            start = current_start
            end = start + base_size
            strand = self._choose_strand()

            intervals.append(
                GenomicInterval(
                    chrom=chromosome,
                    start=start,
                    end=end,
                    name=f"separated_{i}",
                    score=self.rng.randint(0, 1000),
                    strand=strand,
                )
            )

            # Next interval starts after a gap
            current_start = end + gap_size

        return intervals

    def generate_multi_chromosome_scenarios(
        self,
        chromosome_count: int,
        intervals_per_chrom: int,
        scenario_func: str = "basic",
    ) -> List[GenomicInterval]:
        """Generate intervals across multiple chromosomes.

        Args:
            chromosome_count: Number of chromosomes
            intervals_per_chrom: Number of intervals per chromosome
            scenario_func: Scenario type ("basic", "overlapping", "adjacent",
                          "separated")

        Returns:
            List of genomic intervals across multiple chromosomes
        """
        all_intervals = []

        for i in range(chromosome_count):
            chrom_name = f"chr{i + 1}"

            if scenario_func == "overlapping":
                intervals = self.generate_overlapping_scenarios(
                    chrom_name, intervals_per_chrom
                )
            elif scenario_func == "adjacent":
                intervals = self.generate_adjacent_scenarios(
                    chrom_name, intervals_per_chrom
                )
            elif scenario_func == "separated":
                intervals = self.generate_separated_scenarios(
                    chrom_name, intervals_per_chrom
                )
            else:  # basic
                intervals = self.generate_basic(chrom_name, intervals_per_chrom)

            all_intervals.extend(intervals)

        return all_intervals

    def generate_same_strand_pairs(
        self, chromosome: str, pair_count: int, strand: str = "+"
    ) -> Tuple[List[GenomicInterval], List[GenomicInterval]]:
        """Generate two sets of intervals on the same strand.

        Args:
            chromosome: Chromosome name
            pair_count: Number of interval pairs to generate
            strand: Strand to use for all intervals ('+' or '-')

        Returns:
            Tuple of (intervals_a, intervals_b) on same strand
        """
        intervals_a = []
        intervals_b = []
        base_size = self.config.min_interval_size
        current_start = 100

        for i in range(pair_count):
            # Interval A
            start_a = current_start
            end_a = start_a + base_size
            intervals_a.append(
                GenomicInterval(
                    chrom=chromosome,
                    start=start_a,
                    end=end_a,
                    name=f"a{i}",
                    score=self.rng.randint(0, 1000),
                    strand=strand,
                )
            )

            # Interval B - overlaps A, same strand
            start_b = start_a + (base_size // 2)
            end_b = start_b + base_size
            intervals_b.append(
                GenomicInterval(
                    chrom=chromosome,
                    start=start_b,
                    end=end_b,
                    name=f"b{i}",
                    score=self.rng.randint(0, 1000),
                    strand=strand,
                )
            )

            # Move to next region
            current_start = end_b + 100

        return intervals_a, intervals_b

    def generate_opposite_strand_pairs(
        self, chromosome: str, pair_count: int
    ) -> Tuple[List[GenomicInterval], List[GenomicInterval]]:
        """Generate two sets of intervals on opposite strands.

        Args:
            chromosome: Chromosome name
            pair_count: Number of interval pairs to generate

        Returns:
            Tuple of (intervals_a, intervals_b) on opposite strands
        """
        intervals_a = []
        intervals_b = []
        base_size = self.config.min_interval_size
        current_start = 100

        for i in range(pair_count):
            # Interval A on + strand
            start_a = current_start
            end_a = start_a + base_size
            intervals_a.append(
                GenomicInterval(
                    chrom=chromosome,
                    start=start_a,
                    end=end_a,
                    name=f"a{i}",
                    score=self.rng.randint(0, 1000),
                    strand="+",
                )
            )

            # Interval B - overlaps A, opposite strand (-)
            start_b = start_a + (base_size // 2)
            end_b = start_b + base_size
            intervals_b.append(
                GenomicInterval(
                    chrom=chromosome,
                    start=start_b,
                    end=end_b,
                    name=f"b{i}",
                    score=self.rng.randint(0, 1000),
                    strand="-",
                )
            )

            # Move to next region
            current_start = end_b + 100

        return intervals_a, intervals_b

    def generate_mixed_strand_intervals(
        self, chromosome: str, count: int
    ) -> List[GenomicInterval]:
        """Generate intervals with mixed strand assignments.

        Args:
            chromosome: Chromosome name
            count: Number of intervals to generate

        Returns:
            List of intervals with randomly assigned strands (+, -, .)
        """
        intervals = []
        base_size = self.config.min_interval_size
        strands = ["+", "-", "."]
        current_start = 100

        for i in range(count):
            start = current_start
            end = start + base_size
            # Randomly choose strand from +, -, .
            strand = self.rng.choice(strands)

            intervals.append(
                GenomicInterval(
                    chrom=chromosome,
                    start=start,
                    end=end,
                    name=f"mixed_{i}",
                    score=self.rng.randint(0, 1000),
                    strand=strand,
                )
            )

            current_start = end + 50  # Small gap

        return intervals
