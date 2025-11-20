# Implementation Plan: Bedtools Integration Test Suite

**Branch**: `003-bedtools-integration-tests` | **Date**: 2025-11-17 | **Spec**: [spec.md](./spec.md)
**Input**: Feature specification from `/specs/003-bedtools-integration-tests/spec.md`

**Note**: This template is filled in by the `/speckit.plan` command. See `.specify/templates/commands/plan.md` for the execution workflow.

## Summary

This feature implements a comprehensive integration test suite that validates GIQL query results against bedtools command outputs using simulated genomic datasets. The test suite executes GIQL queries against a DuckDB database and compares results to bedtools operations on equivalent BED files, ensuring GIQL correctly replicates bedtools functionality for core genomic interval operations (intersect, merge, nearest, distance calculations). Tests use pytest fixtures for database isolation and include comprehensive scenario coverage with varied overlap patterns, multiple chromosomes, strand combinations, and boundary cases.

## Technical Context

**Language/Version**: Python 3.11+
**Primary Dependencies**: pytest, DuckDB Python driver, bedtools (external binary 2.30.0+), GIQL library (project under test)
**Storage**: DuckDB in-memory databases (via pytest fixtures)
**Testing**: pytest with fixtures, parameterized tests, subprocess for bedtools execution
**Target Platform**: Development/CI environments (Linux, macOS)
**Project Type**: Testing infrastructure (single project structure under tests/)
**Performance Goals**: Complete test suite execution in <5 minutes, dataset generation <10 seconds for 10k intervals
**Constraints**: Requires bedtools 2.30.0+ installed in PATH, temp directory write access
**Scale/Scope**: Validate 80%+ of bedtools operations available in GIQL, comprehensive scenario coverage

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

### Core Principles Alignment

**I. Expressive** - ✅ PASS (N/A for testing infrastructure)
- This feature tests GIQL expressiveness but doesn't directly implement query operators
- Test suite validates that GIQL operators correctly express genomic operations

**II. Readable** - ✅ PASS
- Test structure follows AAA pattern with clear Given-When-Then docstrings per constitution
- Test names descriptive: `test_<operation>_<scenario>` format
- Comparison logic clearly documents what constitutes a "match" between GIQL and bedtools

**III. Portable** - ✅ PASS
- Tests validate GIQL portability by ensuring consistent semantics
- DuckDB chosen as test database (already supported per constitution)
- Test suite itself is portable across Linux/macOS development environments

**IV. Canonical** - ✅ PASS
- Directly validates canonical behavior by comparing GIQL to bedtools (the reference standard)
- Constitution states "bedtools behavior is the reference standard" - this test suite enforces that

### Testing Requirements Compliance

**Coverage** - ✅ PLANNED
- Integration tests for GIQL operators (intersect, merge, nearest, distance)
- Comprehensive edge case coverage per FR-005
- Each test verifies distinct behavior (no overlap)

**Test Structure** - ✅ COMPLIANT
- AAA pattern required by constitution
- Given-When-Then docstrings for all test functions
- Mirror source structure: tests/integration/bedtools/ (to be created)
- Test naming: `test_<operation>_<scenario>` format

**Test Isolation** - ✅ COMPLIANT
- Each test independent via pytest fixtures (FR-004, FR-013)
- DuckDB in-memory databases for fast, isolated tests
- No shared state between tests
- Deterministic test data generation (FR-005)

**Database Testing** - ✅ COMPLIANT
- In-memory DuckDB per constitution guidance
- Fresh database instances via fixtures (no shared state)
- No external database servers required

### Gate Status (Pre-Research)

✅ **PASS** - All constitution requirements met. Proceeded to Phase 0 research.

### Gate Re-Evaluation (Post-Design)

**I. Expressive** - ✅ PASS (N/A for testing infrastructure)
- Design maintains focus on validating GIQL expressiveness
- No changes to assessment

**II. Readable** - ✅ PASS
- Contracts use clear protocols with descriptive names
- Data model entities have self-documenting fields
- Quickstart provides clear examples with Given-When-Then structure
- All design artifacts maintain readability focus

**III. Portable** - ✅ PASS
- Design uses portable technologies (Python stdlib, pytest, DuckDB)
- No platform-specific dependencies beyond bedtools (which is cross-platform)
- Test suite runs on Linux/macOS as documented

**IV. Canonical** - ✅ PASS
- Design enforces canonical behavior through bedtools comparison
- ComparisonResult ensures exact matching per bedtools reference
- No deviations from bedtools semantics

**Testing Requirements** - ✅ PASS
- Data model validates all test structure requirements
- Contracts enforce AAA pattern and Given-When-Then docstrings
- Fixtures provide function-scoped isolation per constitution
- quickstart.md demonstrates compliant test structure

✅ **FINAL GATE PASS** - Design maintains all constitution compliance. Ready for implementation (`/speckit.tasks`).

## Project Structure

### Documentation (this feature)

```text
specs/[###-feature]/
├── plan.md              # This file (/speckit.plan command output)
├── research.md          # Phase 0 output (/speckit.plan command)
├── data-model.md        # Phase 1 output (/speckit.plan command)
├── quickstart.md        # Phase 1 output (/speckit.plan command)
├── contracts/           # Phase 1 output (/speckit.plan command)
└── tasks.md             # Phase 2 output (/speckit.tasks command - NOT created by /speckit.plan)
```

### Source Code (repository root)

```text
tests/
├── integration/
│   └── bedtools/                     # New: Bedtools integration tests
│       ├── conftest.py               # Pytest fixtures for DuckDB, bedtools
│       ├── test_data_generator.py    # Simulated dataset generation
│       ├── test_intersect.py         # Intersect operation tests
│       ├── test_merge.py             # Merge operation tests
│       ├── test_nearest.py           # Nearest/closest operation tests
│       ├── test_distance.py          # Distance calculation tests
│       ├── test_strand_aware.py      # Strand-specific operation tests
│       ├── test_workflows.py         # Multi-operation workflow tests
│       └── utils/
│           ├── bedtools_wrapper.py   # Bedtools command execution
│           ├── comparison.py         # Result comparison logic
│           └── bed_export.py         # Export DuckDB data to BED files
└── unit/                             # Existing unit tests (unchanged)
```

**Structure Decision**: Single project structure (Option 1) with new integration test directory. Tests are added under `tests/integration/bedtools/` to maintain clear separation from existing unit tests. This follows the constitution requirement to mirror source structure and keep tests organized by test type (contract/integration/unit).

## Complexity Tracking

No constitution violations. This section is not applicable.
