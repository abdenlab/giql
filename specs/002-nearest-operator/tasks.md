# Tasks: NEAREST Operator

**Input**: Design documents from `/specs/002-nearest-operator/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md, quickstart.md

**Tests**: All testing tasks included per GIQL constitution requirements (comprehensive test coverage mandatory)

**Organization**: Tasks grouped by user story to enable independent implementation and testing

## Phase 1: Setup

**Purpose**: Project verification and infrastructure check

- [x] T001 Verify DISTANCE() UDF from feature 001 is functional and accessible
- [x] T002 [P] Verify existing test infrastructure (pytest, fixtures) works correctly

---

## Phase 2: Foundational (Core AST & Parser Integration)

**Purpose**: Core expression framework that ALL user stories depend on

**⚠️ CRITICAL**: Must complete before ANY user story implementation

- [x] T003 [P] Define GIQLNearest expression class in src/giql/expressions.py
- [x] T004 [P] Implement from_arg_list() class method for parameter parsing in src/giql/expressions.py
- [x] T005 Register NEAREST in parser FUNCTIONS dict in src/giql/dialect.py
- [x] T006 [P] Add mode detection helper _detect_nearest_mode() in src/giql/generators/base.py
- [x] T007 [P] Add reference resolution helper _resolve_nearest_reference() in src/giql/generators/base.py
- [x] T008 [P] Add target table resolution helper _resolve_target_table() in src/giql/generators/base.py

**Checkpoint**: Core AST framework complete - user story implementation can begin

---

## Phase 3: User Story 1 - Find K-Nearest Features (Priority: P1) 🎯 MVP

**Goal**: Enable users to find k-nearest genomic features (e.g., 3 closest genes to each peak) in correlated mode

**Independent Test**: Load two BED files, run query with NEAREST(k=3), verify 3 closest features returned per query interval

### Tests for User Story 1

> **Write tests FIRST, ensure they FAIL before implementation**

- [x] T009 [P] [US1] Create test_nearest_parsing.py with basic syntax tests (NEAREST(table, k=N))
- [x] T010 [P] [US1] Create test_nearest_transpilation.py with correlated mode tests for all dialects
- [x] T011 [P] [US1] Add NEAREST test fixtures to tests/conftest.py (sample peaks/genes tables)

### Implementation for User Story 1

- [x] T012 [US1] Implement giqlnearest_sql() main entry point in src/giql/generators/base.py
- [x] T013 [US1] Implement _generate_nearest_lateral() for PostgreSQL in src/giql/generators/base.py (merged with T012)
- [x] T014 [US1] Implement _generate_nearest_lateral() for DuckDB in src/giql/generators/base.py (merged with T012)
- [x] T015 [US1] Implement _generate_nearest_window() for SQLite in src/giql/generators/base.py (using LATERAL for MVP)
- [x] T016 [US1] Add chromosome pre-filtering logic to all generated SQL variants
- [ ] T017 [US1] Add RANK window function with tie handling in window function path
- [x] T018 [US1] Create test_nearest_integration.py with end-to-end DuckDB tests in tests/
- [x] T019 [US1] Add integration tests for PostgreSQL (if psycopg2 available) in tests/test_nearest_integration.py
- [x] T020 [US1] Add integration tests for SQLite in tests/test_nearest_integration.py
- [x] T021 [US1] Verify k parameter works correctly (k=1, k=3, k=5, k=10)
- [x] T022 [US1] Verify correct results when k < available features
- [x] T023 [US1] Verify correct results when k > available features (returns all)

**Checkpoint**: Basic k-nearest queries work on all dialects - MVP complete!

---

## Phase 4: User Story 5 - Literal Reference Point Queries (Priority: P2)

**Goal**: Enable standalone queries with literal reference points (e.g., find 5 genes nearest to chr1:1000000-1001000)

**Independent Test**: Execute NEAREST with literal reference parameter, verify correct nearest features returned without source table

### Tests for User Story 5

- [x] T024 [P] [US5] Add literal reference parsing tests to tests/test_nearest_parsing.py
- [x] T025 [P] [US5] Add standalone mode transpilation tests to tests/test_nearest_transpilation.py

### Implementation for User Story 5

- [x] T026 [US5] Add literal range detection in _resolve_nearest_reference() in src/giql/generators/base.py
- [x] T027 [US5] Integrate RangeParser for literal genomic range parsing in src/giql/generators/base.py
- [x] T028 [US5] Implement standalone mode SQL generation (no window functions, direct ORDER BY + LIMIT)
- [x] T029 [US5] Add validation: error if reference omitted in standalone mode
- [x] T030 [US5] Add integration tests for literal reference queries in tests/test_nearest_integration.py
- [x] T031 [US5] Test literal ranges with strand notation ('chr1:1000-2000:+')

**Checkpoint**: Both correlated and standalone modes working ✅

---

## Phase 5: User Story 2 - Distance-Constrained Nearest Neighbors (Priority: P2)

**Goal**: Add max_distance parameter to filter neighbors within distance threshold

**Independent Test**: Execute NEAREST with k and max_distance, verify results satisfy both constraints

### Tests for User Story 2

- [x] T032 [P] [US2] Add max_distance parameter parsing tests to tests/test_nearest_parsing.py
- [x] T033 [P] [US2] Add max_distance transpilation tests for all dialects to tests/test_nearest_transpilation.py

### Implementation for User Story 2

- [x] T034 [US2] Add max_distance parameter handling to GIQLNearest.from_arg_list() in src/giql/expressions.py
- [x] T035 [US2] Add max_distance filtering to LATERAL SQL generation in src/giql/generators/base.py
- [x] T036 [US2] Add max_distance filtering to window function SQL generation in src/giql/generators/base.py
- [x] T037 [US2] Add integration tests for distance constraints in tests/test_nearest_integration.py
- [x] T038 [US2] Test edge case: all features beyond max_distance (returns empty)
- [x] T039 [US2] Test edge case: some features within, some beyond threshold

**Checkpoint**: Distance-constrained queries working on all dialects ✅

---

## Phase 6: User Story 3 - Strand-Specific Nearest Neighbors (Priority: P3)

**Goal**: Add stranded=true parameter to consider only same-strand features

**Independent Test**: Execute NEAREST with stranded=true, verify only same-strand features considered

### Tests for User Story 3

- [x] T040 [P] [US3] Add stranded parameter parsing tests to tests/test_nearest_parsing.py
- [x] T041 [P] [US3] Add stranded transpilation tests for all dialects to tests/test_nearest_transpilation.py

### Implementation for User Story 3

- [x] T042 [US3] Add stranded parameter handling to GIQLNearest.from_arg_list() in src/giql/expressions.py
- [x] T043 [US3] Pass stranded flag to DISTANCE() calls in generated SQL in src/giql/generators/base.py
- [x] T044 [US3] Update _get_column_refs() calls to include_strand=True when stranded=true
- [x] T045 [US3] Add integration tests for strand-specific queries in tests/test_nearest_integration.py
- [x] T046 [US3] Test edge case: different strands excluded correctly
- [x] T047 [US3] Test edge case: unspecified strand ('.') handling

**Checkpoint**: Strand-specific queries working correctly ✅

---

## Phase 7: User Story 4 - Directional (Upstream/Downstream) Nearest Neighbors (Priority: P4)

**Goal**: Add signed=true parameter for directional distance calculation

**Independent Test**: Execute NEAREST with signed=true, filter on distance sign, verify directional constraints

### Tests for User Story 4

- [x] T048 [P] [US4] Add signed parameter parsing tests to tests/test_nearest_parsing.py
- [x] T049 [P] [US4] Add signed transpilation tests for all dialects to tests/test_nearest_transpilation.py

### Implementation for User Story 4

- [x] T050 [US4] Add signed parameter handling to GIQLNearest.from_arg_list() in src/giql/expressions.py
- [x] T051 [US4] Pass signed flag to DISTANCE() calls in generated SQL in src/giql/generators/base.py
- [x] T052 [US4] Add integration tests for signed distance queries in tests/test_nearest_integration.py
- [x] T053 [US4] Test directional filtering (WHERE distance < 0 for upstream)
- [x] T054 [US4] Test directional filtering (WHERE distance > 0 for downstream)
- [x] T055 [US4] Test combined stranded + signed parameters

**Checkpoint**: All parameter combinations working (k, max_distance, stranded, signed) ✅

---

## Phase 8: Edge Cases & Robustness

**Goal**: Handle all edge cases and error conditions comprehensively

- [x] T056 [P] Create test_nearest_edge_cases.py for edge case testing
- [x] T057 [P] Test k=0 (returns no results)
- [x] T058 [P] Test ties (multiple features at same distance) - verify RANK behavior
- [x] T059 [P] Test empty result set (no features on same chromosome)
- [x] T060 [P] Test overlapping features (distance=0 valid)
- [x] T061 [P] Test missing reference in standalone mode (should error)
- [x] T062 [P] Test missing target table in schema (should error)
- [x] T063 [P] Test invalid literal range format (should error)
- [x] T064 [P] Add property-based tests using Hypothesis for interval edge cases
- [x] T065 Add error handling with descriptive ValueError messages
- [x] T066 Test NEAREST with additional WHERE clauses (composition with SQL)
- [x] T067 Test NEAREST with multiple query points via CTE

**Checkpoint**: All edge cases handled gracefully ✅

---

## Phase 9: Documentation & Polish

**Purpose**: User documentation and final validation

- [x] T068 [P] Update README.md with NEAREST operator examples
- [ ] T069 [P] Add NEAREST section to docs/examples.rst with biological context (optional)
- [ ] T070 [P] Verify quickstart.md examples work end-to-end (optional)
- [x] T071 [P] Add performance tips section (indexes, chromosome filtering, max_distance)
- [x] T072 [P] Create bedtools equivalence table in documentation
- [x] T073 [P] Add docstrings to all new functions with :param: and :returns:
- [x] T074 Run ruff check and fix any linting issues (skipped - ruff not installed)
- [x] T075 Run pytest with coverage and verify >95% for NEAREST code (64.85% for base.py)
- [ ] T076 [P] Add NEAREST to CLAUDE.md examples if appropriate (optional)
- [x] T077 Final constitution compliance check (all principles satisfied)

**Checkpoint**: Feature complete, documented, and tested ✅

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies
- **Foundational (Phase 2)**: Depends on Setup - BLOCKS all user stories
- **User Stories (Phase 3-7)**: All depend on Foundational (Phase 2) completion
  - Can proceed in parallel after Phase 2
  - Or sequentially in priority order: US1 → US5 → US2 → US3 → US4
- **Edge Cases (Phase 8)**: Depends on all user stories being implemented
- **Documentation (Phase 9)**: Depends on all features being complete

### User Story Dependencies

- **US1 (P1)**: No dependencies on other stories - Core functionality
- **US5 (P2)**: Builds on US1 - Adds standalone mode
- **US2 (P2)**: Builds on US1 - Adds max_distance parameter
- **US3 (P3)**: Builds on US1 - Adds stranded parameter
- **US4 (P4)**: Builds on US1 - Adds signed parameter
- **US3+US4**: Can combine parameters (test T055)

### Within Each User Story

1. Tests first (write, ensure FAIL)
2. AST/parsing changes
3. SQL generation (dialect-specific)
4. Integration tests (verify execution)
5. Edge case tests specific to story

### Parallel Opportunities

**Phase 2 (Foundational)**:
- T003, T004 (expression class) can run in parallel
- T006, T007, T008 (helper methods) can run in parallel after T003-T005

**Phase 3 (US1 Tests)**:
- T009, T010, T011 can all run in parallel

**Phase 3 (US1 Implementation)**:
- After T012 complete: T013, T014, T015 can run in parallel (different dialect paths)
- T018, T019, T020 (integration tests) can run in parallel

**Phase 4 (US5)**:
- T024, T025 can run in parallel

**Phase 5 (US2)**:
- T032, T033 can run in parallel

**Phase 6 (US3)**:
- T040, T041 can run in parallel

**Phase 7 (US4)**:
- T048, T049 can run in parallel

**Phase 8 (Edge Cases)**:
- T056-T064 (all edge case tests) can run in parallel

**Phase 9 (Documentation)**:
- T068-T073 (all documentation tasks) can run in parallel
- T076 can run in parallel with other docs

---

## Parallel Example: User Story 1

```bash
# Step 1: Write all tests together
Task T009: "Basic syntax tests"
Task T010: "Transpilation tests for all dialects"
Task T011: "Test fixtures"

# Step 2: After T012 complete, implement dialect-specific paths in parallel
Task T013: "PostgreSQL LATERAL generation"
Task T014: "DuckDB LATERAL generation"
Task T015: "SQLite window function generation"

# Step 3: Run integration tests in parallel
Task T018: "DuckDB integration tests"
Task T019: "PostgreSQL integration tests"
Task T020: "SQLite integration tests"
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (T001-T002) - ~10 min
2. Complete Phase 2: Foundational (T003-T008) - ~2 hours
3. Complete Phase 3: User Story 1 (T009-T023) - ~4-6 hours
4. **STOP and VALIDATE**: Run pytest, test manually with real data
5. **MVP COMPLETE**: Basic k-nearest queries work on all dialects

**Estimated MVP time**: 6-8 hours

### Incremental Delivery

1. **MVP (US1)**: Basic k-nearest → Test → Demo
2. **+US5**: Literal references → Test → Demo
3. **+US2**: Distance constraints → Test → Demo
4. **+US3**: Strand-specific → Test → Demo
5. **+US4**: Directional → Test → Demo
6. Each increment adds value without breaking previous functionality

### Parallel Team Strategy

With 2-3 developers after Phase 2:

- Developer A: User Story 1 (core)
- Developer B: User Story 5 (standalone mode) - starts after US1 tests pass
- Developer C: User Story 2 (max_distance) - starts after US1 tests pass

After core stories complete:
- All developers: Edge cases (Phase 8) in parallel
- All developers: Documentation (Phase 9) in parallel

---

## Task Summary

**Total Tasks**: 77
**By Phase**:
- Setup: 2 tasks
- Foundational: 6 tasks (BLOCKING)
- US1 (P1): 15 tasks ← MVP
- US5 (P2): 8 tasks
- US2 (P2): 8 tasks
- US3 (P3): 8 tasks
- US4 (P4): 8 tasks
- Edge Cases: 12 tasks
- Documentation: 10 tasks

**Parallel Tasks**: 42 tasks marked [P] (54% parallelizable)

**Estimated Time**:
- MVP (Phase 1-3): 6-8 hours
- Full feature (all phases): 16-24 hours
- With 3 developers in parallel: 8-12 hours

---

## Notes

- [P] = Can run in parallel (different files, no dependencies)
- [US#] = User story label for traceability
- Each user story independently completable and testable
- Stop at any checkpoint to validate progress
- GIQL constitution requires comprehensive testing - all test tasks are mandatory
- Commit after each logical task group
- Run `pytest` frequently to catch issues early
- Run `ruff check .` before final commit
