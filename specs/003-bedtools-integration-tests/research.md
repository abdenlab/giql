# Research: Bedtools Integration Test Suite

**Date**: 2025-11-17
**Feature**: 003-bedtools-integration-tests

## Overview

This document consolidates research findings for implementing an integration test suite that validates GIQL query results against bedtools command outputs. Research focused on test data generation strategies, result comparison approaches, pytest fixture patterns, and bedtools integration techniques.

## Research Areas

### 1. Simulated Genomic Dataset Generation

**Decision**: Use programmatic generation with controlled randomness (seeded random for reproducibility)

**Rationale**:
- Provides fine-grained control over test scenarios (overlap density, chromosome distribution, strand patterns)
- Enables comprehensive coverage of edge cases (boundary conditions, large intervals, multiple chromosomes)
- Faster than loading static files and more flexible for parameterized tests
- Deterministic (seeded) generation ensures reproducible test failures

**Alternatives Considered**:
- **Static BED files**: Rejected because inflexible, doesn't scale to comprehensive scenario coverage, harder to maintain
- **Real genomic data**: Rejected because too large, unpredictable, adds unnecessary complexity, harder to isolate failure causes
- **Hypothesis property-based testing**: Considered as complement, not replacement, for targeted scenario testing

**Implementation Approach**:
- Create `IntervalGenerator` class with configurable properties:
  - `chromosome_count`: Number of chromosomes (default: 3-5 for tests)
  - `interval_count_per_chrom`: Intervals per chromosome
  - `overlap_probability`: Control overlap density (0.0-1.0)
  - `strand_distribution`: Proportions of +/- /unstranded intervals
  - `min_interval_size`, `max_interval_size`: Interval length bounds
  - `seed`: Random seed for reproducibility
- Generate scenarios: overlapping, adjacent, separated, boundary-spanning, cross-chromosome
- Output formats: In-memory list of dicts, DuckDB table, BED file

### 2. DuckDB to BED File Export

**Decision**: Use Python string formatting with TSV writer for BED export

**Rationale**:
- BED format is simple TSV: chrom, start, end, [name, score, strand]
- No complex parsing required, just column mapping and formatting
- DuckDB query results → list of tuples → write to temp file
- Maintain column order: chrom, chromStart, chromEnd, name, score, strand (BED6 format)

**Alternatives Considered**:
- **pandas DataFrame**: Rejected as unnecessary dependency, adds overhead
- **pybedtools**: Rejected to avoid dependency on bedtools Python wrapper (we're testing against bedtools binary)

**Implementation Approach**:
```python
def export_to_bed(rows, bed_path, format='bed6'):
    """Export query results to BED file.

    Args:
        rows: List of (chrom, start, end, name, score, strand) tuples
        bed_path: Output BED file path
        format: 'bed3', 'bed6' (default bed6 for full info)
    """
    with open(bed_path, 'w') as f:
        for row in rows:
            if format == 'bed3':
                line = f"{row[0]}\t{row[1]}\t{row[2]}\n"
            elif format == 'bed6':
                line = f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\t{row[5]}\n"
            f.write(line)
```

### 3. Bedtools Command Execution

**Decision**: Use subprocess with explicit command construction and output parsing

**Rationale**:
- Direct control over bedtools CLI arguments
- Capture stdout/stderr for debugging test failures
- Version checking via `bedtools --version`
- Timeout protection for hanging commands
- Clear error messages when bedtools unavailable

**Alternatives Considered**:
- **pybedtools wrapper**: Rejected to maintain independence from Python wrappers (test against canonical bedtools binary)
- **Docker container**: Rejected as adds complexity, not needed for dev/CI environments with bedtools installed

**Implementation Approach**:
```python
def run_bedtools_command(command_args, input_files=None, timeout=30):
    """Execute bedtools command and return output.

    Args:
        command_args: List of command arguments ['intersect', '-a', 'file1.bed', '-b', 'file2.bed']
        input_files: Optional dict of stdin inputs
        timeout: Command timeout in seconds

    Returns:
        stdout as string

    Raises:
        BedtoolsNotFoundError: bedtools not in PATH
        BedtoolsVersionError: bedtools version < 2.30.0
        BedtoolsExecutionError: Command failed
    """
    cmd = ['bedtools'] + command_args
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=True
    )
    return result.stdout
```

### 4. Result Comparison Strategy

**Decision**: Multi-level comparison with exact integer matching, epsilon float comparison, and order-independent set comparison

**Rationale**:
- Integer positions/counts must match exactly (per clarifications)
- Floating-point values need epsilon tolerance (standard practice: 1e-9 or 1e-10)
- SQL query results may have different row ordering than bedtools
- Clear diff output showing mismatches for debugging

**Alternatives Considered**:
- **Byte-for-byte comparison**: Rejected because row ordering differs
- **Sorted comparison only**: Insufficient, doesn't handle floating-point precision
- **pytest-approx for all**: Rejected because positions must be exact

**Implementation Approach**:
```python
def compare_results(giql_rows, bedtools_rows, epsilon=1e-9):
    """Compare GIQL and bedtools results with appropriate tolerance.

    Args:
        giql_rows: List of result tuples from GIQL query
        bedtools_rows: List of result tuples from bedtools output
        epsilon: Tolerance for floating-point comparisons

    Returns:
        ComparisonResult(match=bool, differences=list)
    """
    # Sort both result sets for order-independent comparison
    giql_sorted = sorted(giql_rows, key=sort_key)
    bt_sorted = sorted(bedtools_rows, key=sort_key)

    # Compare counts
    if len(giql_sorted) != len(bt_sorted):
        return ComparisonResult(False, f"Row count mismatch: {len(giql_sorted)} vs {len(bt_sorted)}")

    # Compare each row
    for i, (g_row, bt_row) in enumerate(zip(giql_sorted, bt_sorted)):
        for j, (g_val, bt_val) in enumerate(zip(g_row, bt_row)):
            if isinstance(g_val, float):
                if abs(g_val - bt_val) > epsilon:
                    return ComparisonResult(False, f"Row {i}, col {j}: {g_val} != {bt_val} (epsilon={epsilon})")
            else:
                if g_val != bt_val:
                    return ComparisonResult(False, f"Row {i}, col {j}: {g_val} != {bt_val}")

    return ComparisonResult(True, None)
```

### 5. Pytest Fixture Design

**Decision**: Hierarchical fixtures with function scope for isolation

**Rationale**:
- **Function-scoped database fixture**: Each test gets clean DuckDB instance (maximum isolation per constitution)
- **Session-scoped bedtools version check**: Verify once per test run, fail fast
- **Function-scoped temp directory**: Clean workspace for each test's BED files
- **Parameterized scenario fixtures**: Reuse dataset generators across tests

**Alternatives Considered**:
- **Module-scoped database**: Rejected due to shared state risk
- **Class-scoped fixtures**: Not needed, tests are independent functions
- **Global temp directory**: Rejected due to cleanup and isolation concerns

**Implementation Approach**:
```python
# conftest.py

import pytest
import duckdb
import tempfile
import shutil
from pathlib import Path

@pytest.fixture(scope="session")
def bedtools_version():
    """Verify bedtools version 2.30.0+ is installed."""
    result = subprocess.run(['bedtools', '--version'], capture_output=True, text=True)
    version = parse_version(result.stdout)
    if version < (2, 30, 0):
        pytest.fail(f"Bedtools version 2.30.0+ required, found {version}")
    return version

@pytest.fixture(scope="function")
def duckdb_connection():
    """Provide clean DuckDB connection for each test."""
    conn = duckdb.connect(':memory:')
    yield conn
    conn.close()

@pytest.fixture(scope="function")
def temp_bed_dir():
    """Provide temporary directory for BED files."""
    temp_dir = tempfile.mkdtemp(prefix='giql_bedtools_test_')
    yield Path(temp_dir)
    shutil.rmtree(temp_dir)

@pytest.fixture
def interval_generator():
    """Provide configured interval generator."""
    return IntervalGenerator(seed=42)  # Deterministic for reproducibility
```

### 6. Test Organization Strategy

**Decision**: One test file per bedtools operation, parameterized tests for scenarios

**Rationale**:
- Clear mapping: `test_intersect.py` → bedtools intersect command
- Parameterized tests reduce code duplication for similar scenarios
- Easy to run subsets: `pytest tests/integration/bedtools/test_intersect.py`
- Follows constitution test naming: `test_<operation>_<scenario>`

**Alternatives Considered**:
- **Single monolithic test file**: Rejected as unmaintainable for comprehensive coverage
- **Test class hierarchy**: Rejected as unnecessary complexity, functions sufficient
- **Scenario-based files**: Rejected as harder to map to bedtools operations

**Implementation Approach**:
```python
# test_intersect.py

@pytest.mark.parametrize("overlap_type", ["full", "partial", "adjacent", "none"])
def test_intersect_overlap_types(duckdb_connection, temp_bed_dir, interval_generator, overlap_type):
    """Test INTERSECTS operator with various overlap patterns.

    Given:
        Two tables with intervals of different overlap patterns
    When:
        GIQL query uses INTERSECTS operator
    Then:
        Results match bedtools intersect output exactly
    """
    # Arrange
    intervals_a, intervals_b = interval_generator.generate_overlap_scenario(overlap_type)
    # ... test implementation
```

### 7. Bedtools Operations Coverage

**Decision**: Prioritize P1 operations (intersect, merge, nearest), then P2 (distance, strand-aware)

**Rationale**:
- Aligns with user story priorities in spec
- P1 operations are foundational, must be correct first
- Incremental test development matches implementation priorities
- Easy to track progress against 80% coverage goal (SC-001)

**Operations to Test** (in priority order):
1. **P1**: `bedtools intersect` → GIQL INTERSECTS
2. **P1**: `bedtools merge` → GIQL MERGE
3. **P1**: `bedtools closest` → GIQL NEAREST
4. **P2**: `bedtools closest -d` → GIQL DISTANCE
5. **P2**: `bedtools intersect -s` → GIQL strand-aware INTERSECTS
6. **P2**: `bedtools closest -s` → GIQL strand-aware NEAREST
7. **P3**: Pipeline workflows (intersect | merge, etc.)

**Command-Line Flags to Test**:
- `-a`, `-b`: Input files (all operations)
- `-s`: Same strand (intersect, closest)
- `-S`: Opposite strand (intersect, closest)
- `-d`: Report distance (closest)
- No flags: Default behavior (all)

## Summary

Key technical decisions finalized:
- **Data generation**: Programmatic with seeded random for reproducibility
- **DuckDB export**: Simple TSV formatting to BED files
- **Bedtools execution**: Direct subprocess with version checking
- **Result comparison**: Multi-level with exact int, epsilon float, order-independent
- **Pytest fixtures**: Function-scoped for isolation, session-scoped for version check
- **Test organization**: One file per operation, parameterized scenarios
- **Coverage priority**: P1 operations first (intersect, merge, nearest)

All technical unknowns from plan.md resolved. Ready for Phase 1 design.
