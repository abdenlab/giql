# Quickstart: Bedtools Integration Test Suite

**Feature**: 003-bedtools-integration-tests
**Date**: 2025-11-17

## Overview

This guide helps you quickly get started with the bedtools integration test suite, whether you're adding new tests, debugging failures, or understanding the test architecture.

## Prerequisites

Before running the tests, ensure you have:

1. **Python 3.11+** installed
2. **Bedtools 2.30.0+** installed and in PATH
   ```bash
   bedtools --version
   # Should output: bedtools v2.30.0 or higher
   ```
3. **pytest** installed
   ```bash
   pip install pytest pytest-xdist
   ```
4. **DuckDB Python driver** installed
   ```bash
   pip install duckdb
   ```

## Quick Start

### Run All Bedtools Integration Tests

```bash
# From repository root
pytest tests/integration/bedtools/ -v
```

### Run Specific Operation Tests

```bash
# Test only intersect operations
pytest tests/integration/bedtools/test_intersect.py -v

# Test only merge operations
pytest tests/integration/bedtools/test_merge.py -v

# Test only nearest operations
pytest tests/integration/bedtools/test_nearest.py -v
```

### Run Tests in Parallel

```bash
# Use pytest-xdist for parallel execution
pytest tests/integration/bedtools/ -n auto
```

### Run with Detailed Output

```bash
# Show stdout and detailed failure messages
pytest tests/integration/bedtools/ -v -s

# Show only failed test details
pytest tests/integration/bedtools/ -v --tb=short
```

## Test Structure Overview

```text
tests/integration/bedtools/
├── conftest.py              # Shared pytest fixtures
├── test_intersect.py        # Intersect operation tests
├── test_merge.py            # Merge operation tests
├── test_nearest.py          # Nearest/closest operation tests
├── test_distance.py         # Distance calculation tests
├── test_strand_aware.py     # Strand-specific tests
├── test_workflows.py        # Multi-operation workflow tests
└── utils/
    ├── bedtools_wrapper.py  # Bedtools command execution
    ├── comparison.py        # Result comparison logic
    └── bed_export.py        # BED file export utilities
```

## Writing Your First Test

### Example: Test Intersect Operation

```python
# tests/integration/bedtools/test_intersect.py

import pytest

def test_intersect_basic_overlap(duckdb_connection, temp_bed_dir, interval_generator):
    """Test INTERSECTS predicate finds overlapping intervals.

    Given:
        Two tables with genomic intervals where some intervals overlap
    When:
        A GIQL query uses INTERSECTS predicate in WHERE clause
    Then:
        Results match bedtools intersect output exactly
    """
    # Arrange: Generate overlapping intervals
    intervals_a = interval_generator.generate_overlapping('chr1', count=10)
    intervals_b = interval_generator.generate_overlapping('chr1', count=10)

    # Load into DuckDB
    load_intervals(duckdb_connection, 'intervals_a', intervals_a)
    load_intervals(duckdb_connection, 'intervals_b', intervals_b)

    # Export to BED files
    bed_a = temp_bed_dir / 'a.bed'
    bed_b = temp_bed_dir / 'b.bed'
    export_to_bed(duckdb_connection, 'intervals_a', bed_a)
    export_to_bed(duckdb_connection, 'intervals_b', bed_b)

    # Act: Execute GIQL query
    giql_query = """
        SELECT a.chrom, a.start, a.end, a.name, a.score, a.strand
        FROM intervals_a a, intervals_b b
        WHERE a.position INTERSECTS b.position
    """
    giql_result = duckdb_connection.execute(giql_query).fetchall()

    # Execute bedtools command
    bedtools_output = run_bedtools_command(['intersect', '-a', str(bed_a), '-b', str(bed_b)])
    bedtools_result = parse_bedtools_output(bedtools_output)

    # Assert: Compare results
    comparison = compare_results(giql_result, bedtools_result)
    assert comparison, comparison.failure_message()
```

## Common Fixtures

### `duckdb_connection`

Provides a clean DuckDB in-memory connection for each test.

```python
def test_example(duckdb_connection):
    # Fresh database, no shared state
    duckdb_connection.execute("CREATE TABLE test_table (...)")
```

### `temp_bed_dir`

Provides a temporary directory for BED files, automatically cleaned up.

```python
def test_example(temp_bed_dir):
    bed_file = temp_bed_dir / 'test.bed'
    # Write BED data
    # Directory cleaned up after test
```

### `interval_generator`

Provides configured interval generator with deterministic seed.

```python
def test_example(interval_generator):
    intervals = interval_generator.generate_overlapping('chr1', count=100)
    # Same intervals every test run (deterministic)
```

### `bedtools_version`

Session-scoped fixture that verifies bedtools is installed and compatible.

```python
def test_example(bedtools_version):
    # Test automatically skipped if bedtools not available or version < 2.30.0
```

## Debugging Test Failures

### View Detailed Comparison Output

When tests fail, the `ComparisonResult.failure_message()` provides detailed output:

```text
✗ Results do not match
  GIQL rows: 15
  Bedtools rows: 17
  Differences:
    - Row count mismatch
    - Missing in GIQL: chr1:1000-2000 (score=100, strand=+)
    - Extra in GIQL: chr1:1500-2500 (score=200, strand=-)
```

### Inspect Generated Data

Add `-s` flag to see print statements:

```python
def test_debug_example(interval_generator):
    intervals = interval_generator.generate_overlapping('chr1', 10)
    print(f"Generated {len(intervals)} intervals:")
    for interval in intervals:
        print(f"  {interval}")
```

```bash
pytest tests/integration/bedtools/test_intersect.py::test_debug_example -v -s
```

### Compare BED Files Manually

Preserve temp directory to inspect BED files:

```python
def test_manual_inspection(temp_bed_dir):
    print(f"BED files in: {temp_bed_dir}")
    # ... test code ...
    import time
    time.sleep(60)  # Keep files for 60 seconds
```

Or use pytest's `--basetemp`:

```bash
pytest tests/integration/bedtools/ --basetemp=/tmp/giql_test_debug
```

### Run Single Test

```bash
pytest tests/integration/bedtools/test_intersect.py::test_intersect_basic_overlap -v -s
```

## Test Scenarios

The test suite validates GIQL against bedtools across these scenarios:

### Overlap Patterns
- **Full overlap**: Intervals completely contained
- **Partial overlap**: Intervals partially intersect
- **Adjacent**: Intervals touch but don't overlap
- **Separated**: Intervals have gaps between them

### Chromosome Coverage
- **Single chromosome**: All intervals on chr1
- **Multiple chromosomes**: Intervals across chr1, chr2, chr3
- **Boundary cases**: Intervals at chromosome start/end

### Strand Configurations
- **Same strand**: Both intervals on + or both on -
- **Opposite strands**: One on +, one on -
- **Unstranded**: Strand = "."
- **Mixed**: Combination of stranded and unstranded

### Interval Sizes
- **Small intervals**: 100-500bp
- **Medium intervals**: 500-5000bp
- **Large intervals**: 5000-50000bp
- **Very large**: Entire chromosome arms

## Performance Benchmarks

Expected test execution times:

- **Single test**: < 100ms
- **Test file** (e.g., test_intersect.py): < 5 seconds
- **Full suite**: < 5 minutes (per SC-002)

If tests run slower, check:

1. Dataset sizes (should be 100-1000 intervals for most tests)
2. Bedtools execution (ensure bedtools installed locally, not via network mount)
3. Temp directory location (use local filesystem, not network drive)

## Continuous Integration

### GitHub Actions Example

```yaml
name: Bedtools Integration Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3

      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.11'

      - name: Install bedtools
        run: |
          sudo apt-get update
          sudo apt-get install -y bedtools
          bedtools --version

      - name: Install dependencies
        run: |
          pip install pytest pytest-xdist duckdb

      - name: Run tests
        run: |
          pytest tests/integration/bedtools/ -v -n auto
```

## Extending the Test Suite

### Adding a New Operation Test

1. **Create test file**: `tests/integration/bedtools/test_newoperation.py`
2. **Follow naming convention**: `test_<operation>_<scenario>`
3. **Use Given-When-Then docstrings** (per constitution)
4. **Parameterize scenarios** where applicable:

```python
@pytest.mark.parametrize("overlap_type", ["full", "partial", "none"])
def test_newop_overlap_types(overlap_type, duckdb_connection, ...):
    """Test new operation with various overlap patterns.

    Given:
        Intervals with {overlap_type} overlap
    When:
        GIQL NEWOP operator is used
    Then:
        Results match bedtools newop output
    """
```

### Adding a New Scenario

1. **Update `interval_generator`** with new generation method
2. **Add scenario to existing test files** via parameterization
3. **Document expected behavior** in test docstring

## Troubleshooting

### Bedtools Not Found

```
BedtoolsNotFoundError: bedtools binary not found in PATH
```

**Solution**: Install bedtools or add to PATH:
```bash
# macOS
brew install bedtools

# Ubuntu/Debian
sudo apt-get install bedtools

# Verify
which bedtools
bedtools --version
```

### Version Incompatible

```
BedtoolsVersionError: bedtools version 2.28.0 found, 2.30.0+ required
```

**Solution**: Update bedtools:
```bash
# macOS
brew upgrade bedtools

# Ubuntu/Debian - may need to install from source
```

### DuckDB Import Error

```
ModuleNotFoundError: No module named 'duckdb'
```

**Solution**:
```bash
pip install duckdb
```

### Test Failures After GIQL Changes

If tests fail after modifying GIQL query operators:

1. **Check if GIQL semantics changed** - may be expected
2. **Verify bedtools behavior** - run bedtools command manually
3. **Update expected results** if GIQL intentionally changed
4. **Add new tests** for new edge cases discovered

## Next Steps

- **Read**: [data-model.md](./data-model.md) for entity definitions
- **Read**: [contracts/README.md](./contracts/README.md) for interface contracts
- **Read**: [research.md](./research.md) for design decisions
- **Implement**: Use `/speckit.tasks` to generate implementation tasks
