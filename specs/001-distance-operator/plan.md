# Implementation Plan: DISTANCE UDF Operator

**Branch**: `001-distance-operator` | **Date**: 2025-11-06 | **Spec**: [spec.md](./spec.md)
**Input**: Feature specification from `/specs/001-distance-operator/spec.md`

**Note**: This template is filled in by the `/speckit.plan` command. See `.specify/templates/commands/plan.md` for the execution workflow.

## Summary

Implement a DISTANCE() UDF (User-Defined Function) operator that calculates genomic distances between intervals, replicating bedtools' `closest -d` functionality. The operator accepts flexible inputs (position columns or literal genomic ranges) and supports optional `stranded` and `signed` parameters for strand-specific and directional distance calculations. The implementation will transpile GIQL queries to SQL with UDF calls across DuckDB, SQLite, and PostgreSQL backends, maintaining portability while providing a more expressive and readable syntax than raw SQL.

## Technical Context

**Language/Version**: Python 3.11+
**Primary Dependencies**: sqlglot (>=20.0.0), duckdb (>=1.4.0), psycopg2-binary (>=2.9.10)
**Storage**: N/A (operates on user-provided database connections)
**Testing**: pytest (>=7.0.0), pytest-cov (>=4.0.0), hypothesis (property-based testing)
**Target Platform**: Cross-platform (Linux, macOS, Windows) - Python library + CLI
**Project Type**: Single project (Python library with CLI)
**Performance Goals**: Performance is database-dependent; GIQL adds minimal overhead for query transpilation and UDF dispatch
**Constraints**:
- Must work on DuckDB, SQLite, and PostgreSQL without requiring database-specific features
- Generated SQL must be human-readable
- No pre-sorting requirement (unlike bedtools)
**Scale/Scope**:
- Single new UDF operator with 2 optional parameters
- ~500-800 LOC for implementation (UDF + SQL generator + tests)
- 3 database dialect implementations

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

### Expressive ✅
- **Requirement**: Common genomic operations should have dedicated operators
- **Compliance**: DISTANCE() provides a concise, dedicated operator for a common genomic analysis pattern (finding nearest features)
- **Status**: PASS

### Readable ✅
- **Requirement**: Operator names must be self-describing and queries should express intent
- **Compliance**: `DISTANCE(a.position, b.position)` clearly communicates distance calculation; optional named parameters (`stranded=true`, `signed=true`) are self-documenting
- **Status**: PASS

### Portable ✅
- **Requirement**: Must work across DuckDB, SQLite, PostgreSQL without database-specific features
- **Compliance**: UDF pattern works across all three backends; implementation uses standard SQL types and logic
- **Status**: PASS

### Canonical ✅
- **Requirement**: Must match established tools (bedtools) semantics
- **Compliance**: Explicitly designed to replicate bedtools `closest -d/-D/-s` behavior with explicit mappings documented
- **Status**: PASS

### Code Quality & Testing ✅
- **Requirement**: 100% public API coverage, behavior-focused tests, Given-When-Then docstrings
- **Compliance**: Implementation plan includes:
  - Unit tests for UDF logic (all edge cases)
  - Cross-dialect transpilation tests
  - Integration tests for "find closest" workflows
  - Property-based tests with Hypothesis
- **Status**: PASS (pending implementation)

**Gate Status**: ✅ PASS - No violations, proceed to Phase 0

## Project Structure

### Documentation (this feature)

```text
specs/001-distance-operator/
├── plan.md              # This file (/speckit.plan command output)
├── spec.md              # Feature specification (copied from .specify/specs/)
├── research.md          # Phase 0 output - UDF implementation patterns
├── data-model.md        # Phase 1 output - DISTANCE function interface design
├── quickstart.md        # Phase 1 output - Usage examples
└── tasks.md             # Phase 2 output (/speckit.tasks command - NOT created by /speckit.plan)
```

### Source Code (repository root)

```text
src/giql/
├── expressions.py           # Add GIQLDistance expression node
├── dialect.py               # Add DISTANCE parsing logic
├── transformer.py           # (May need updates for position resolution)
├── generators/
│   ├── base.py             # Add distance_sql() generation method
│   ├── duckdb.py           # DuckDB-specific UDF registration
│   ├── sqlite.py           # SQLite-specific UDF registration
│   └── postgres.py         # PostgreSQL-specific UDF registration
├── udfs/
│   └── distance.py         # NEW: Core distance calculation logic
└── engine.py               # Register DISTANCE UDF on connection

tests/
├── unit/
│   ├── test_distance_udf.py           # NEW: Unit tests for distance calculation logic
│   ├── test_distance_parsing.py       # NEW: Parser tests for DISTANCE syntax
│   └── test_distance_transpilation.py # NEW: SQL generation tests
├── integration/
│   └── test_distance_queries.py       # NEW: End-to-end query tests
└── property/
    └── test_distance_properties.py    # NEW: Hypothesis-based property tests
```

**Structure Decision**: Single project structure (existing GIQL library). New code integrates into existing modules following established patterns:
- `expressions.py` for AST nodes (like `GIQLCluster`, `GIQLMerge`)
- `dialect.py` for parsing
- `generators/` for SQL generation per dialect
- New `udfs/` module for UDF implementations
- Tests mirror source structure per constitution

## Complexity Tracking

> **Fill ONLY if Constitution Check has violations that must be justified**

N/A - No constitution violations

## Phase 0: Research & Unknowns

### Research Topics

1. **UDF Registration Patterns Across Dialects**
   - **Question**: What are the exact API differences for registering Python UDFs in DuckDB vs SQLite vs PostgreSQL?
   - **Needed For**: Implementing database-specific UDF registration in generators

2. **Handling Named Parameters in UDFs**
   - **Question**: How do we pass named parameters (`stranded=true`, `signed=true`) through SQL to Python UDFs in each dialect?
   - **Needed For**: Ensuring `stranded` and `signed` parameters work correctly across backends

3. **NULL Handling in Window Function ORDER BY**
   - **Question**: How do NULL distance values (from different chromosomes/strands) behave in `ORDER BY DISTANCE(...)` clauses across different SQL dialects?
   - **Needed For**: Ensuring consistent behavior when finding closest features

4. **Strand Information Resolution**
   - **Question**: How does the existing codebase resolve position columns to chromosome/start/end? Does it already handle strand columns?
   - **Needed For**: Extending position resolution to include strand information

### Research Output Location

See `research.md` (generated in Phase 0)

## Phase 1: Design Artifacts

### Data Model

**Entity**: DISTANCE Function

**Inputs**:
- `interval_a`: GenomicInterval (position column or literal range string)
- `interval_b`: GenomicInterval (position column or literal range string)
- `stranded`: Boolean (default=false)
- `signed`: Boolean (default=false)

**Output**:
- Integer | NULL

**Resolution Pipeline**:
1. GIQL Parser → `GIQLDistance` AST node
2. SQL Generator → Resolve position columns to (chrom, start, end, strand)
3. SQL Generator → Generate UDF call with resolved columns + boolean flags
4. Database → Execute UDF, return distance value

See `data-model.md` for full details (generated in Phase 1)

### API Contracts

**GIQL Syntax Contract**:
```sql
DISTANCE(interval_a, interval_b [, stranded=false] [, signed=false]) -> INTEGER | NULL
```

**SQL UDF Signature Contract** (internal, post-transpilation):
```sql
DISTANCE(chrom_a TEXT, start_a INT, end_a INT, strand_a TEXT,
         chrom_b TEXT, start_b INT, end_b INT, strand_b TEXT,
         stranded BOOLEAN, signed BOOLEAN) -> INTEGER | NULL
```

See `contracts/` directory for OpenAPI-style documentation (generated in Phase 1)

### Quickstart

Basic usage examples for users:
1. Calculate distance between two tables
2. Find closest feature using window functions
3. Strand-specific distance queries
4. Signed distance (directional) queries

See `quickstart.md` (generated in Phase 1)

## Phase 2: Task Breakdown

*This phase is handled by `/speckit.tasks` command - not part of `/speckit.plan` output*

Tasks will be generated based on:
- Implementation of `GIQLDistance` expression node
- Parser updates for DISTANCE syntax
- SQL generator updates for each dialect
- UDF implementation and registration
- Comprehensive test suite (unit, integration, property-based)
- Documentation and examples

## Implementation Notes

### Key Design Decisions

1. **UDF vs Special Syntax**: Chose UDF approach (like MERGE/CLUSTER) rather than special predicate syntax (like INTERSECTS) because DISTANCE is a scalar function returning a value, not a boolean predicate

2. **Position Resolution**: Extends existing position column resolution pattern from INTERSECTS/WITHIN to include strand information

3. **Signed Distance Orientation**: When `signed=true`, directionality is always relative to reference genome coordinates (FR-013), not strand-specific orientation

4. **NULL Semantics**: Returns NULL for:
   - Different chromosomes (always)
   - Different strands (when `stranded=true`)
   - Any NULL input (chromosome, start, end, or strand)

### Integration Points

- **Parser (`dialect.py`)**: Add DISTANCE function parsing similar to CLUSTER/MERGE
- **Expressions (`expressions.py`)**: Add `GIQLDistance` class with `from_arg_list()` for named parameters
- **SQL Generation (`generators/base.py`)**: Add `distance_sql()` method to resolve positions and generate UDF call
- **UDF Registration (`engine.py`)**: Register DISTANCE UDF when creating database connections

### Testing Strategy

1. **Unit Tests** (test behavior, not implementation):
   - Distance calculation algorithm (overlaps → 0, gaps → positive, NULL cases)
   - Signed distance logic (upstream → negative, downstream → positive)
   - Stranded mode logic (different strands → NULL when enabled)

2. **Parser Tests**:
   - Valid DISTANCE syntax variations
   - Named parameter parsing
   - Error messages for invalid syntax

3. **Transpilation Tests**:
   - GIQL query → SQL output for each dialect
   - Position column resolution (including strand)
   - Literal range handling

4. **Integration Tests**:
   - End-to-end "find closest feature" workflow
   - Cross-dialect consistency
   - Performance with realistic datasets

5. **Property-Based Tests** (Hypothesis):
   - Distance calculation invariants (symmetric except for signed mode)
   - Overlap detection correctness
   - NULL handling across all input combinations

### Performance Considerations

- **No Pre-Sorting**: Unlike bedtools, GIQL DISTANCE doesn't require sorted input (documented as a feature, but users should be aware CROSS JOINs can be expensive)
- **Chromosome Pre-Filter**: Examples should demonstrate `WHERE a.chromosome = b.chromosome` for performance
- **Database-Dependent**: Actual query performance depends on the underlying database engine and dataset characteristics

### Documentation Requirements

1. Update README.md with DISTANCE examples
2. Add bedtools equivalence table
3. Document performance characteristics and best practices
4. Add "Finding Closest Features" tutorial section
5. Update CLI help text to mention DISTANCE

## Next Steps

After plan approval:
1. Run `/speckit.tasks` to generate detailed task breakdown
2. Begin Phase 0 research (UDF patterns, named parameters, etc.)
3. Create Phase 1 design artifacts (data-model.md, quickstart.md)
4. Implement in task order with test-first approach per constitution
