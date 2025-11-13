# Implementation Plan: NEAREST Operator

**Branch**: `002-nearest-operator` | **Date**: 2025-11-07 | **Spec**: [spec.md](./spec.md)
**Input**: Feature specification from `/specs/002-nearest-operator/spec.md`

**Note**: This template is filled in by the `/speckit.plan` command. See `.specify/templates/commands/plan.md` for the execution workflow.

## Summary

The NEAREST operator provides a table-valued function for finding k-nearest genomic features, eliminating the need for users to write complex window functions. It supports two modes: correlated mode for per-row k-NN queries (`CROSS JOIN LATERAL NEAREST(genes, k=3)`) and standalone mode for literal reference points (`FROM NEAREST(genes, reference='chr1:1000-2000', k=3)`). Implementation will transpile NEAREST to SQL using RANK window functions over the existing DISTANCE() UDF, with support for distance constraints, strand-specific, and directional queries across all supported dialects (DuckDB, SQLite, PostgreSQL).

## Technical Context

**Language/Version**: Python 3.11+
**Primary Dependencies**: sqlglot (>=20.0.0), duckdb (>=1.4.0), psycopg2-binary (>=2.9.10)
**Storage**: N/A (operates on existing database connections)
**Testing**: pytest with property-based testing (Hypothesis)
**Target Platform**: Cross-platform (Linux, macOS, Windows) - Python library + CLI
**Project Type**: Single project (SQL transpiler library)
**Performance Goals**: Transpilation <100ms for typical queries; execution performance depends on database backend query optimizer
**Constraints**: Must transpile to valid SQL for DuckDB, SQLite, PostgreSQL; no database-specific features in core operator; LATERAL JOIN support required
**Scale/Scope**: Small feature addition (~5-8 new files, ~1000-1500 LOC total including tests)

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

### Core Principles Alignment

✅ **I. Expressive**
- NEAREST provides dedicated operator for common k-NN pattern
- Reduces SQL boilerplate (eliminates manual window functions, RANK, PARTITION BY)
- Composes naturally with standard SQL (WHERE, ORDER BY on result set)
- **Gate status**: PASS - Follows GIQL's mandate to provide dedicated operators for common genomic operations

✅ **II. Readable**
- Self-describing name: `NEAREST(genes, k=3)` clearly expresses intent
- Mirrors logical flow: "find nearest genes"
- Generated SQL will be human-readable (CTEs with descriptive names)
- Documentation will include biological context examples
- **Gate status**: PASS - Query communicates "what" (find nearest) not "how" (window functions)

✅ **III. Portable**
- Transpiles to standard SQL (window functions: RANK, PARTITION BY, ORDER BY)
- No database-specific features required
- Explicit testing for all three supported dialects (DuckDB, SQLite, PostgreSQL)
- LATERAL JOIN is supported by PostgreSQL natively; other dialects will need transpilation strategy
- **Gate status**: PASS with RESEARCH NEEDED - Must verify LATERAL JOIN compatibility/transpilation for SQLite and DuckDB

✅ **IV. Canonical**
- Extends bedtools `closest` semantics (k=1 matches bedtools behavior)
- Follows spatial database conventions (PostGIS LATERAL pattern, SQL Server APPLY)
- Tie handling matches bedtools `-t all` (RANK vs ROW_NUMBER)
- Edge cases align with DISTANCE operator (NULL for different chromosomes, strand handling)
- **Gate status**: PASS - Strong alignment with both genomics tools (bedtools) and spatial databases (PostGIS)

### Testing Requirements Alignment

✅ **Test Coverage**
- Will include comprehensive unit tests for transpilation logic
- Integration tests for all supported dialects
- Property-based tests (Hypothesis) for edge cases
- Tests for both standalone and correlated modes
- **Gate status**: PASS - Follows constitution testing requirements

✅ **Test Structure**
- AAA pattern (Arrange-Act-Assert)
- Given-When-Then docstrings for all test functions
- Tests will verify behavior (transpiled SQL + execution results), not implementation
- **Gate status**: PASS - Aligned with "Test Behavior, Not Implementation" philosophy

### Research Requirements

The following areas need research before design phase:

1. **LATERAL JOIN Portability** (NEEDS CLARIFICATION)
   - PostgreSQL has native LATERAL support
   - DuckDB LATERAL support status?
   - SQLite LATERAL support status (likely NO - need transpilation strategy)
   - **Alternative**: Transpile LATERAL to correlated subqueries for SQLite

2. **Window Function Optimization** (NEEDS CLARIFICATION)
   - Best practices for RANK vs ROW_NUMBER performance across dialects
   - Impact of PARTITION BY on large datasets
   - Query optimizer behavior with DISTANCE() in ORDER BY

3. **AST Integration** (NEEDS CLARIFICATION)
   - How to represent NEAREST in sqlglot AST
   - Custom expression class vs function call node
   - Parameter handling for named arguments (k, reference, max_distance, stranded, signed)

4. **Schema Resolution** (NEEDS CLARIFICATION)
   - How to resolve target table name to position columns
   - Implicit `.position` column convention enforcement
   - Error handling when target table/column not found

## Project Structure

### Documentation (this feature)

```text
specs/002-nearest-operator/
├── plan.md              # This file
├── research.md          # Phase 0: LATERAL portability, AST design, schema resolution
├── data-model.md        # Phase 1: NEAREST expression AST nodes, transpilation data flow
├── quickstart.md        # Phase 1: User guide with examples
├── contracts/           # Phase 1: N/A (library, not API)
├── checklists/
│   └── requirements.md  # Spec validation checklist (complete)
└── spec.md              # Feature specification (complete)
```

### Source Code (repository root)

```text
src/
├── giql/
│   ├── expressions.py              # NEW: GIQLNearest expression class
│   ├── transformer.py              # MODIFY: Add NEAREST parsing logic
│   ├── range_parser.py             # REUSE: Parse literal genomic ranges
│   ├── schema.py                   # MODIFY: Table/column resolution for NEAREST
│   ├── generators/
│   │   ├── base.py                 # MODIFY: Add nearest_sql() method
│   │   ├── duckdb.py               # MODIFY: DuckDB-specific NEAREST (if needed)
│   │   └── __init__.py             # (existing)
│   ├── engine.py                   # REUSE: Existing transpilation/execution
│   └── __init__.py                 # (existing)

tests/
├── test_nearest_parsing.py         # NEW: Parse NEAREST syntax
├── test_nearest_transpilation.py   # NEW: GIQL→SQL for all dialects
├── test_nearest_integration.py     # NEW: End-to-end execution tests
├── test_nearest_edge_cases.py      # NEW: Ties, k=0, empty results
├── conftest.py                      # MODIFY: Add NEAREST test fixtures
└── (existing test files)           # (unchanged)
```

**Structure Decision**: Single project structure. NEAREST is a new operator in the existing GIQL transpiler. Core implementation involves: (1) new GIQLNearest expression class, (2) parsing logic in transformer, (3) SQL generation in generators, (4) comprehensive test suite. Estimated ~5-8 new files, ~1000-1500 LOC including tests.

## Complexity Tracking

> **Fill ONLY if Constitution Check has violations that must be justified**

No violations identified. All constitution principles pass initial gate check. Research phase will resolve NEEDS CLARIFICATION items (LATERAL portability, AST design patterns).

## Phase 0: Research

**Goal**: Resolve all NEEDS CLARIFICATION items from Technical Context

### Research Tasks

1. **LATERAL JOIN Cross-Dialect Support**
   - **Question**: How to implement LATERAL semantics in SQLite (lacks native support)?
   - **Approach**:
     - Test PostgreSQL LATERAL with DISTANCE()
     - Test DuckDB LATERAL support and syntax
     - Research SQLite alternatives: correlated subqueries, CTEs
   - **Deliverable**: Transpilation strategy for each dialect

2. **Window Function Performance Patterns**
   - **Question**: What's the optimal window function strategy for k-NN queries?
   - **Approach**:
     - Compare RANK vs ROW_NUMBER performance
     - Analyze PARTITION BY cardinality impact
     - Review DuckDB/PostgreSQL query optimizer docs
   - **Deliverable**: Best practices for window function generation

3. **sqlglot AST Extension Patterns**
   - **Question**: How to add custom NEAREST expression to sqlglot AST?
   - **Approach**:
     - Review existing custom expressions (Intersects, Contains, GIQLDistance)
     - Identify pattern for function-like table expressions
     - Design parameter handling for named arguments
   - **Deliverable**: AST design for GIQLNearest expression

4. **Schema Resolution Strategy**
   - **Question**: How to resolve target table → position columns at transpilation time?
   - **Approach**:
     - Review existing SchemaInfo usage in base generator
     - Design table name resolution logic
     - Define error handling for missing tables/columns
   - **Deliverable**: Schema resolution algorithm

**Output**: `research.md` with findings and decisions

## Phase 1: Design & Contracts

**Prerequisites**: research.md complete with all decisions finalized

### 1. Data Model (data-model.md)

**Entities**:

- **GIQLNearest Expression** (AST Node)
  - `target_table`: Table name to search (string)
  - `reference`: Optional position reference (column ref or literal range)
  - `k`: Maximum neighbors to return (integer, default=1)
  - `max_distance`: Optional distance threshold (integer or None)
  - `stranded`: Boolean flag for strand-specific search (default=False)
  - `signed`: Boolean flag for directional distance (default=False)
  - Parent class: `exp.Expression` (sqlglot base)

- **TranspilationContext** (Internal State)
  - `mode`: "standalone" | "correlated" (detected from query structure)
  - `outer_table`: Table name from LATERAL context (for implicit reference resolution)
  - `target_schema`: Resolved position columns for target table
  - `reference_schema`: Resolved position columns for reference (if column ref)

- **Generated SQL Structure** (Output)
  - **Correlated mode**: CTE with RANK window function partitioned by query interval
  - **Standalone mode**: Simple SELECT with ORDER BY + LIMIT
  - Both include: DISTANCE() UDF calls, chromosome pre-filter, distance constraints

### 2. API Contracts

N/A - GIQL is a transpiler library, not a REST/GraphQL API. User-facing interface is:
- **Python API**: `engine.execute("SELECT * FROM NEAREST(...)")`
- **CLI**: `giql "SELECT * FROM NEAREST(...)"`

Contract is the GIQL syntax itself, documented in quickstart.md

### 3. Quickstart Guide (quickstart.md)

**Sections**:
1. Introduction: What is NEAREST and why use it?
2. Basic Usage: Correlated mode (per-row k-NN)
3. Standalone Mode: Literal reference points
4. Advanced Features: max_distance, stranded, signed
5. Filtering Results: Combining with WHERE clauses
6. Performance Tips: Chromosome pre-filtering, distance constraints
7. Comparison with bedtools: Equivalence table

### 4. Agent Context Update

Run `.specify/scripts/bash/update-agent-context.sh claude` to add:
- New technology: (none - reusing existing stack)
- New commands: (none - using existing `pytest`, `ruff check .`)

## Phase 2: Task Generation

**NOT executed by /speckit.plan** - This is handled by `/speckit.tasks` command

## Notes

- **Dependency**: Requires DISTANCE() UDF from feature 001-distance-operator (already implemented)
- **Backward Compatibility**: New operator, no breaking changes
- **Documentation Impact**: Update README.md, add examples to docs/examples.rst
- **Version Impact**: MINOR version bump (new operator)
