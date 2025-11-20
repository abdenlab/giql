# Feature Specification: Bedtools Integration Test Suite

**Feature Branch**: `003-bedtools-integration-tests`
**Created**: 2025-11-17
**Status**: Draft
**Input**: User description: "Implement an integration test suite that evaluates GIQL queries against bedtools for accuracy using a simulated dataset. Replicate as much of bedtools functionality as is feasible using GIQL. Make the DB engine to use a pytest fixture."

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Basic Interval Operations Validation (Priority: P1)

As a developer, I want to validate that GIQL's core interval operations produce identical results to bedtools, so that I can trust GIQL as a replacement for bedtools workflows in genomic data analysis.

**Why this priority**: Core interval operations (intersect, merge, closest/nearest) are the foundation of most genomic analyses. These must be correct before any other functionality can be trusted. This represents the minimum viable test suite.

**Independent Test**: Can be fully tested by running a subset of bedtools commands (intersect, merge, nearest) against the same test dataset and comparing outputs, delivering confidence in basic operations.

**Acceptance Scenarios**:

1. **Given** two BED files with overlapping genomic intervals in the test database, **When** a GIQL query performs an intersection operation, **Then** the results match bedtools intersect output exactly (same intervals, same counts)
2. **Given** a BED file with adjacent intervals, **When** a GIQL query performs a merge operation, **Then** the merged intervals match bedtools merge output
3. **Given** two BED files with non-overlapping intervals, **When** a GIQL query finds nearest neighbors, **Then** the results match bedtools closest/nearest output

---

### User Story 2 - Strand-Aware Operations (Priority: P2)

As a developer, I want to verify that GIQL correctly handles strand-specific interval operations, so that strand-sensitive genomic analyses produce accurate results.

**Why this priority**: Many genomic features are strand-specific (genes, transcripts), and incorrect strand handling leads to biologically meaningless results. This is critical but builds on basic operations.

**Independent Test**: Can be tested independently by running bedtools commands with -s and -S flags against strand-annotated test data and comparing outputs.

**Acceptance Scenarios**:

1. **Given** BED files with strand information (+/-), **When** a GIQL query performs strand-specific intersection, **Then** results match bedtools intersect -s output
2. **Given** intervals on opposite strands, **When** a GIQL query finds nearest neighbors with strand consideration, **Then** results match bedtools closest with appropriate strand flags
3. **Given** intervals with mixed strands, **When** a GIQL query ignores strand information, **Then** results match bedtools default (strand-ignorant) behavior

---

### User Story 3 - Distance Calculations (Priority: P2)

As a developer, I want to validate that GIQL's distance calculations match bedtools, so that proximity-based analyses are accurate.

**Why this priority**: Distance calculations are fundamental to many genomic analyses (regulatory element discovery, TAD analysis). This builds on basic operations but is critical for scientific accuracy.

**Independent Test**: Can be tested independently by running bedtools with distance calculation flags and comparing numeric outputs.

**Acceptance Scenarios**:

1. **Given** two non-overlapping intervals, **When** a GIQL query calculates the distance between them, **Then** the distance value matches bedtools output
2. **Given** overlapping intervals, **When** a GIQL query calculates distance (overlap), **Then** the distance value matches bedtools (typically 0 or negative for overlap)
3. **Given** intervals on different chromosomes, **When** a GIQL query attempts distance calculation, **Then** the behavior matches bedtools (typically undefined/null)

---

### User Story 4 - Complex Multi-Operation Workflows (Priority: P3)

As a developer, I want to validate that chained GIQL operations produce results matching bedtools pipelines, so that complex genomic workflows can be trusted.

**Why this priority**: While important for real-world usage, this depends on all basic operations working correctly first. This is the integration validation of multiple features.

**Independent Test**: Can be tested by creating bedtools pipeline commands (using pipes or intermediate files) and comparing final outputs to equivalent GIQL query chains.

**Acceptance Scenarios**:

1. **Given** a workflow that intersects two files then merges the results, **When** executed as a GIQL query, **Then** the final output matches a bedtools pipeline (intersect | merge)
2. **Given** a workflow that finds nearest neighbors then filters by distance, **When** executed as a GIQL query, **Then** results match equivalent bedtools commands
3. **Given** a complex workflow with multiple steps, **When** intermediate results are compared, **Then** each step matches the corresponding bedtools operation

---

### Edge Cases

- What happens when test data contains non-standard chromosome names (chrM, chr1_random)?

**Coverage**: Test datasets include varied overlap patterns (overlapping, adjacent, separated), multiple chromosomes, all strand combinations (+/-, mixed, unstranded), chromosome boundary cases, and very large genomic intervals.

**Note**: Test suite only validates well-formed BED files with non-empty datasets. Malformed data (intervals with start > end, zero-length intervals) and empty input files/tables are excluded from test scope.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: Test suite MUST execute GIQL queries against simulated genomic datasets stored in a database
- **FR-002**: Test suite MUST execute equivalent bedtools commands against the same source data exported to BED format files
- **FR-003**: Test suite MUST compare GIQL query results to bedtools command outputs for equivalence using exact match for positions/counts, epsilon tolerance for floating-point values, and order-independent comparison
- **FR-004**: Test suite MUST use a test fixture to provide isolated database instances for each test
- **FR-005**: Test suite MUST generate simulated genomic datasets with comprehensive coverage including varied overlap patterns (overlapping, adjacent, separated), multiple chromosomes, all strand combinations (+/-, mixed, unstranded), and boundary cases (chromosome edges, very large intervals)
- **FR-006**: Test suite MUST test the following bedtools operations: intersect, merge, closest/nearest, distance calculations
- **FR-007**: Test suite MUST test strand-aware operations (same-strand, opposite-strand, ignore-strand)
- **FR-008**: Test suite MUST validate both the presence/absence of intervals and their exact positions in results
- **FR-009**: Test suite MUST provide clear failure messages indicating which bedtools operation failed and how GIQL output differed
- **FR-010**: Test suite MUST verify bedtools version 2.30.0+ is installed and fail with clear error message if bedtools is unavailable or version is incompatible
- **FR-011**: Test suite MUST use temporary directories for bedtools input/output files that are cleaned up after tests
- **FR-012**: Test suite MUST validate that GIQL operations preserve all required BED file columns (chrom, start, end, name, score, strand)
- **FR-013**: Database fixture MUST provide a clean database instance for each test with isolated state
- **FR-014**: Test suite MUST document which bedtools features are NOT replicated in GIQL and why

### Key Entities

- **Genomic Interval**: Represents a region on a genome with chromosome, start position, end position, optional name, optional score, and optional strand
- **Simulated Dataset**: A collection of genomic intervals with comprehensive scenario coverage including varied overlap patterns, multiple chromosomes, strand combinations, and boundary cases
- **Test Case**: Links a GIQL query, equivalent bedtools command, input dataset, and expected output validation logic
- **Database Fixture**: Provides a database instance configured for testing with clean state for each test

## Clarifications

### Session 2025-11-17

- Q: For malformed input data (e.g., intervals with start > end, zero-length intervals), how should the test suite handle these? → A: Exclude malformed data entirely - test suite only validates well-formed BED files
- Q: When comparing GIQL results to bedtools output, what comparison tolerance should be used for determining a match? → A: Exact match for positions/counts, epsilon for floating-point, order-independent comparison
- Q: How should the test suite validate empty dataset scenarios? → A: Skip empty datasets - focus only on non-empty test cases
- Q: For handling bedtools version differences, what compatibility strategy should the test suite use? → A: Pin to specific version - require bedtools 2.30.0+ and fail if version mismatches
- Q: For simulated test dataset diversity, what range of genomic scenarios should be covered? → A: Comprehensive - varied overlap patterns, multiple chromosomes, strand combinations, boundary cases

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Test suite validates at least 80% of common bedtools operations available in GIQL
- **SC-002**: All integration tests complete execution in under 5 minutes on standard development hardware
- **SC-003**: Test suite achieves 100% pass rate when GIQL operations are correct (no false failures)
- **SC-004**: Test failures provide actionable diagnostics showing exact differences between GIQL and bedtools outputs
- **SC-005**: Simulated datasets are generated in under 10 seconds for datasets with up to 10,000 intervals
- **SC-006**: Test suite can be executed with a single command without additional setup steps (assuming bedtools and database are installed)
- **SC-007**: Documentation clearly lists which bedtools features are tested and which are not yet supported

## Assumptions

- Bedtools is installed and available in the system PATH when running integration tests
- Test environment has write access to temporary directories for bedtools file I/O
- Database system is properly configured and accessible
- GIQL query syntax is stable enough to write deterministic test queries
- Bedtools version is 2.30.0 or higher (for consistent output formats)
- Simulated datasets do not need to represent real biological data, only valid BED format
- Integer positions and counts require exact match; floating-point values use epsilon tolerance; result row ordering may differ between GIQL and bedtools
- Tests run in a controlled environment where bedtools behavior is deterministic
- DuckDB will be used as the database engine for implementation
