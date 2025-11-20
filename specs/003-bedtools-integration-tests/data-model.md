# Data Model: Bedtools Integration Test Suite

**Date**: 2025-11-17
**Feature**: 003-bedtools-integration-tests

## Overview

This document defines the data structures and entities used in the bedtools integration test suite. Since this is a testing infrastructure feature, the "data model" primarily describes test data structures, configuration objects, and result representations rather than persistent database schemas.

## Core Entities

### 1. GenomicInterval

Represents a single genomic interval with all BED file fields.

**Fields**:
- `chrom` (str): Chromosome name (e.g., "chr1", "chr2", "chrX")
- `start` (int): Start position (0-based, inclusive)
- `end` (int): End position (0-based, exclusive)
- `name` (str, optional): Interval name/identifier
- `score` (int, optional): Score value (0-1000)
- `strand` (str, optional): Strand ("+", "-", or ".")

**Validation Rules**:
- `start` must be < `end` (per clarifications: exclude malformed data)
- `start` must be >= 0 (non-negative positions)
- `end` must be > `start`
- `strand` must be one of: "+", "-", ".", or null
- `score` must be in range [0, 1000] if present

**Relationships**:
- Member of `SimulatedDataset`
- Compared in `ComparisonResult`

**Python Representation**:
```python
@dataclass
class GenomicInterval:
    chrom: str
    start: int
    end: int
    name: str | None = None
    score: int | None = None
    strand: str | None = None

    def __post_init__(self):
        if self.start >= self.end:
            raise ValueError(f"Invalid interval: start ({self.start}) >= end ({self.end})")
        if self.start < 0:
            raise ValueError(f"Invalid interval: start ({self.start}) < 0")
        if self.strand and self.strand not in ['+', '-', '.']:
            raise ValueError(f"Invalid strand: {self.strand}")
        if self.score is not None and not (0 <= self.score <= 1000):
            raise ValueError(f"Invalid score: {self.score}")

    def to_bed_line(self, format='bed6') -> str:
        """Convert to BED format line."""
        if format == 'bed3':
            return f"{self.chrom}\t{self.start}\t{self.end}"
        elif format == 'bed6':
            name = self.name or '.'
            score = self.score if self.score is not None else 0
            strand = self.strand or '.'
            return f"{self.chrom}\t{self.start}\t{self.end}\t{name}\t{score}\t{strand}"
```

### 2. SimulatedDataset

Collection of genomic intervals with controlled properties for testing.

**Fields**:
- `name` (str): Dataset identifier (e.g., "intervals_a", "intervals_b")
- `intervals` (List[GenomicInterval]): List of intervals
- `scenario_type` (str): Scenario descriptor (e.g., "overlapping", "adjacent", "separated")
- `metadata` (dict): Generation parameters (seed, chromosome_count, etc.)

**Validation Rules**:
- Must contain at least 1 interval (per clarifications: no empty datasets)
- All intervals must be well-formed (validated by GenomicInterval)
- Intervals within dataset need not be sorted (sorting is operation-specific)

**Relationships**:
- Contains multiple `GenomicInterval` instances
- Used in `TestCase` as input data

**Python Representation**:
```python
@dataclass
class SimulatedDataset:
    name: str
    intervals: List[GenomicInterval]
    scenario_type: str
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        if len(self.intervals) == 0:
            raise ValueError("Dataset must contain at least one interval")

    def to_bed_file(self, path: Path, format='bed6'):
        """Export to BED file."""
        with open(path, 'w') as f:
            for interval in self.intervals:
                f.write(interval.to_bed_line(format) + '\n')

    def to_duckdb_table(self, conn: duckdb.Connection, table_name: str):
        """Load into DuckDB table."""
        rows = [(i.chrom, i.start, i.end, i.name, i.score, i.strand)
                for i in self.intervals]
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
```

### 3. TestCase

Links a GIQL query, bedtools command, input datasets, and validation logic.

**Fields**:
- `test_name` (str): Descriptive test name
- `giql_query` (str): GIQL query to execute
- `bedtools_command` (List[str]): Bedtools command arguments
- `input_datasets` (Dict[str, SimulatedDataset]): Named input datasets
- `expected_behavior` (str): Description of expected outcome
- `operation` (str): Operation being tested (e.g., "intersect", "merge")

**Validation Rules**:
- GIQL query must be non-empty
- Bedtools command must start with valid operation
- Input datasets must match GIQL query and bedtools inputs

**Relationships**:
- References multiple `SimulatedDataset` instances
- Produces `ComparisonResult`

**Python Representation**:
```python
@dataclass
class TestCase:
    test_name: str
    giql_query: str
    bedtools_command: List[str]
    input_datasets: Dict[str, SimulatedDataset]
    expected_behavior: str
    operation: str

    def execute(self, conn: duckdb.Connection, temp_dir: Path) -> 'ComparisonResult':
        """Execute test case and return comparison result."""
        # Load data into DuckDB
        for name, dataset in self.input_datasets.items():
            dataset.to_duckdb_table(conn, name)

        # Execute GIQL query
        giql_result = conn.execute(self.giql_query).fetchall()

        # Export datasets to BED files
        bed_files = {}
        for name, dataset in self.input_datasets.items():
            bed_path = temp_dir / f"{name}.bed"
            dataset.to_bed_file(bed_path)
            bed_files[name] = bed_path

        # Execute bedtools command (replace placeholders with actual paths)
        cmd_with_paths = self._substitute_paths(self.bedtools_command, bed_files)
        bedtools_output = run_bedtools_command(cmd_with_paths)
        bedtools_result = parse_bedtools_output(bedtools_output)

        # Compare results
        return compare_results(giql_result, bedtools_result)
```

### 4. ComparisonResult

Result of comparing GIQL and bedtools outputs.

**Fields**:
- `match` (bool): Whether results match
- `giql_row_count` (int): Number of rows from GIQL query
- `bedtools_row_count` (int): Number of rows from bedtools output
- `differences` (List[str]): Specific differences found (if match=False)
- `comparison_metadata` (dict): Epsilon used, sort order, etc.

**Validation Rules**:
- If `match` is True, `differences` should be empty
- If `match` is False, `differences` should contain at least one item

**Relationships**:
- Produced by `TestCase.execute()`
- Used in test assertions

**Python Representation**:
```python
@dataclass
class ComparisonResult:
    match: bool
    giql_row_count: int
    bedtools_row_count: int
    differences: List[str] = field(default_factory=list)
    comparison_metadata: dict = field(default_factory=dict)

    def __bool__(self) -> bool:
        """Allow direct boolean evaluation in assertions."""
        return self.match

    def failure_message(self) -> str:
        """Generate detailed failure message for test output."""
        if self.match:
            return "Results match"
        msg = [f"Results do not match (GIQL: {self.giql_row_count} rows, bedtools: {self.bedtools_row_count} rows)"]
        msg.extend([f"  - {diff}" for diff in self.differences])
        return '\n'.join(msg)
```

### 5. IntervalGeneratorConfig

Configuration for simulated dataset generation.

**Fields**:
- `chromosome_count` (int): Number of chromosomes to generate
- `intervals_per_chromosome` (int): Intervals per chromosome
- `min_interval_size` (int): Minimum interval length
- `max_interval_size` (int): Maximum interval length
- `overlap_probability` (float): Probability of overlap (0.0-1.0)
- `strand_distribution` (dict): Proportions of +/-/. strands
- `seed` (int): Random seed for reproducibility

**Validation Rules**:
- `chromosome_count` > 0
- `intervals_per_chromosome` > 0
- `min_interval_size` >= 1
- `max_interval_size` >= `min_interval_size`
- `overlap_probability` in [0.0, 1.0]
- `strand_distribution` values sum to 1.0

**Python Representation**:
```python
@dataclass
class IntervalGeneratorConfig:
    chromosome_count: int = 3
    intervals_per_chromosome: int = 100
    min_interval_size: int = 100
    max_interval_size: int = 1000
    overlap_probability: float = 0.3
    strand_distribution: dict = field(default_factory=lambda: {'+': 0.45, '-': 0.45, '.': 0.1})
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
```

### 6. BedtoolsVersion

Represents bedtools version information.

**Fields**:
- `major` (int): Major version number
- `minor` (int): Minor version number
- `patch` (int): Patch version number
- `raw_version_string` (str): Original version string from bedtools

**Validation Rules**:
- Version must be >= 2.30.0 (per FR-010)

**Python Representation**:
```python
@dataclass
class BedtoolsVersion:
    major: int
    minor: int
    patch: int
    raw_version_string: str

    def is_compatible(self) -> bool:
        """Check if version meets minimum requirement (2.30.0)."""
        return (self.major, self.minor, self.patch) >= (2, 30, 0)

    def __str__(self) -> str:
        return f"{self.major}.{self.minor}.{self.patch}"

    @classmethod
    def from_string(cls, version_str: str) -> 'BedtoolsVersion':
        """Parse version from bedtools --version output."""
        # Example: "bedtools v2.30.0"
        match = re.search(r'v?(\d+)\.(\d+)\.(\d+)', version_str)
        if not match:
            raise ValueError(f"Could not parse version from: {version_str}")
        major, minor, patch = map(int, match.groups())
        return cls(major, minor, patch, version_str)
```

## Entity Relationships

```text
IntervalGeneratorConfig → IntervalGenerator
                              ↓
                        SimulatedDataset
                              ↓ (contains)
                        GenomicInterval

TestCase → SimulatedDataset (input_datasets)
         → GIQL Query (string)
         → Bedtools Command (list)
         ↓
    ComparisonResult

BedtoolsVersion → (validated in pytest session fixture)
```

## State Transitions

### GenomicInterval
- Immutable once created (dataclass frozen=True recommended)
- No state transitions

### SimulatedDataset
- Created by IntervalGenerator
- Exported to BED file (one-way, no state change)
- Loaded into DuckDB table (one-way, no state change)

### TestCase
- Created with configuration
- Executed once → produces ComparisonResult
- No state transitions (stateless execution)

### ComparisonResult
- Created from comparison
- Immutable
- No state transitions

## Data Flow

1. **Test Setup**:
   - `IntervalGeneratorConfig` → `IntervalGenerator.generate()` → `SimulatedDataset`
   - `SimulatedDataset` → `to_duckdb_table()` → DuckDB tables
   - `SimulatedDataset` → `to_bed_file()` → BED files

2. **Test Execution**:
   - GIQL query on DuckDB → result rows
   - Bedtools command on BED files → output rows
   - Both results → `compare_results()` → `ComparisonResult`

3. **Test Assertion**:
   - `ComparisonResult.match` → pytest assertion
   - If failed: `ComparisonResult.failure_message()` → test output

## Storage Considerations

- **DuckDB**: In-memory databases (`:memory:` connection string), no persistent storage
- **BED files**: Temporary files in temp directory, cleaned up after each test
- **Test results**: pytest captures, no persistent storage needed
- **Configuration**: Defined in Python code, no external config files

## Performance Considerations

- **Dataset size**: Target 100-1000 intervals per dataset for fast tests
- **Large interval tests**: Separate tests with larger datasets (up to 10,000 intervals per SC-005)
- **Memory**: DuckDB in-memory keeps all data in RAM, limit dataset sizes accordingly
- **File I/O**: BED file export/import is I/O bound, use tempfs if available in CI

## Validation Summary

All entities include validation at construction time using `__post_init__` in dataclasses. This ensures:
- No malformed intervals reach the test execution stage
- Configuration errors caught early with clear messages
- Type safety via Python type hints and runtime checks
