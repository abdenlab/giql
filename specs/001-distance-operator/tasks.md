# Tasks: DISTANCE UDF Operator

**Input**: Design documents from `/specs/001-distance-operator/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md, quickstart.md

**Tests**: Tests are included per constitution requirements (100% public API coverage, behavior-focused)

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3, US4)
- Include exact file paths in descriptions

## Path Conventions

- **Single project**: `src/giql/`, `tests/` at repository root
- Paths follow existing GIQL structure

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Project initialization and file structure

- [ ] T001 Create `src/giql/udfs/` directory for distance calculation logic
- [ ] T002 Create test directories: `tests/unit/`, `tests/integration/`, `tests/property/`
- [ ] T003 [P] Install hypothesis for property-based testing if not already present

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Core infrastructure that MUST be complete before ANY user story can be implemented

**⚠️ CRITICAL**: No user story work can begin until this phase is complete

- [ ] T004 Add `GIQLDistance` expression class to `src/giql/expressions.py`
- [ ] T005 [P] Extend `_get_column_refs()` in `src/giql/generators/base.py` to optionally return strand column
- [ ] T006 [P] Add DISTANCE function parsing to `src/giql/dialect.py`
- [ ] T007 Register `GIQLDistance` in `src/giql/dialect.py` parser token mappings

**Checkpoint**: Foundation ready - user story implementation can now begin in parallel

---

## Phase 3: User Story 1 - Basic Distance Calculation (Priority: P1) 🎯 MVP

**Goal**: Calculate genomic distance between two intervals, returning 0 for overlaps, positive integers for gaps, NULL for different chromosomes

**Independent Test**: Load two BED files, run query with DISTANCE(), verify correct distances (0 for overlaps, positive for gaps, NULL for different chromosomes)

### Tests for User Story 1

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [ ] T008 [P] [US1] Unit test for overlapping intervals (distance=0) in `tests/unit/test_distance_udf.py`
- [ ] T009 [P] [US1] Unit test for non-overlapping intervals (positive distance) in `tests/unit/test_distance_udf.py`
- [ ] T010 [P] [US1] Unit test for different chromosomes (NULL) in `tests/unit/test_distance_udf.py`
- [ ] T011 [P] [US1] Unit test for adjacent/bookended intervals (distance=0) in `tests/unit/test_distance_udf.py`
- [ ] T012 [P] [US1] Unit test for zero-width intervals (point features) in `tests/unit/test_distance_udf.py`
- [ ] T013 [P] [US1] Parser test for `DISTANCE(a.position, b.position)` syntax in `tests/unit/test_distance_parsing.py`
- [ ] T014 [P] [US1] Parser test for literal range `DISTANCE(a.pos, 'chr1:100-200')` in `tests/unit/test_distance_parsing.py`
- [ ] T015 [P] [US1] Transpilation test for DuckDB SQL generation in `tests/unit/test_distance_transpilation.py`
- [ ] T016 [P] [US1] Transpilation test for SQLite SQL generation in `tests/unit/test_distance_transpilation.py`
- [ ] T017 [P] [US1] Transpilation test for PostgreSQL SQL generation in `tests/unit/test_distance_transpilation.py`
- [ ] T018 [US1] Integration test for basic DISTANCE query execution in `tests/integration/test_distance_queries.py`

### Implementation for User Story 1

- [ ] T019 [US1] Implement `from_arg_list()` classmethod in `GIQLDistance` to parse positional arguments
- [ ] T020 [US1] Implement `giqldistance_sql()` method in `src/giql/generators/base.py` for SQL CASE expression generation
- [ ] T021 [US1] Implement `_generate_distance_case()` helper method in `src/giql/generators/base.py` for SQL CASE logic (different chromosomes, overlaps, gaps)
- [ ] T022 [US1] Handle column-to-column distance in `giqldistance_sql()` (resolve both position columns)
- [ ] T023 [US1] Handle literal range distance in `giqldistance_sql()` (parse range, resolve one position column)
- [ ] T024 [US1] Add NULL handling for different chromosomes in SQL CASE expression
- [ ] T025 [US1] Add overlap detection logic (return 0) in SQL CASE expression
- [ ] T026 [US1] Add gap distance calculation logic in SQL CASE expression
- [ ] T027 [US1] Add error handling for invalid range literals with clear error messages
- [ ] T028 [US1] Verify all T008-T018 tests now pass

**Checkpoint**: At this point, User Story 1 should be fully functional - basic DISTANCE operator works for overlapping/non-overlapping intervals

---

## Phase 4: User Story 2 - Finding Nearest Features (Priority: P2)

**Goal**: Enable queries that combine DISTANCE() with window functions to find closest features, replicating `bedtools closest -d` functionality

**Independent Test**: Create query combining DISTANCE() with ROW_NUMBER() window function to select minimum distance feature per query interval

### Tests for User Story 2

- [ ] T029 [P] [US2] Integration test for "find closest feature" query pattern in `tests/integration/test_distance_queries.py`
- [ ] T030 [P] [US2] Integration test for handling ties (multiple features at same distance) in `tests/integration/test_distance_queries.py`
- [ ] T031 [P] [US2] Integration test for NULL distance handling in ORDER BY with `NULLS LAST` in `tests/integration/test_distance_queries.py`
- [ ] T032 [P] [US2] Property-based test for distance calculation invariants with Hypothesis in `tests/property/test_distance_properties.py`

### Implementation for User Story 2

- [ ] T033 [US2] Add documentation comment to `giqldistance_sql()` explaining window function usage pattern
- [ ] T034 [US2] Add example queries to quickstart.md showing "find closest feature" pattern (already in quickstart.md, verify accuracy)
- [ ] T035 [US2] Add example showing k-nearest neighbors pattern to quickstart.md (already in quickstart.md, verify accuracy)
- [ ] T036 [US2] Verify all T029-T032 tests pass

**Checkpoint**: At this point, User Stories 1 AND 2 should both work - users can find closest features using window functions

---

## Phase 5: User Story 3 - Strand-Specific Distance (Priority: P3)

**Goal**: Calculate distances only between features on the same strand using `stranded=true` parameter, mirroring bedtools `-s` flag

**Independent Test**: Pass `stranded=true` to DISTANCE() and verify it returns NULL for intervals on different strands

### Tests for User Story 3

- [ ] T037 [P] [US3] Unit test for same-strand distance (returns value) in `tests/unit/test_distance_udf.py`
- [ ] T038 [P] [US3] Unit test for different-strand distance (returns NULL) in `tests/unit/test_distance_udf.py`
- [ ] T039 [P] [US3] Unit test for unspecified strand handling (strand='.') in `tests/unit/test_distance_udf.py`
- [ ] T040 [P] [US3] Parser test for named parameter `stranded=true` in `tests/unit/test_distance_parsing.py`
- [ ] T041 [P] [US3] Transpilation test including strand comparison in SQL CASE in `tests/unit/test_distance_transpilation.py`
- [ ] T042 [US3] Integration test for stranded DISTANCE queries in `tests/integration/test_distance_queries.py`

### Implementation for User Story 3

- [ ] T043 [US3] Update `from_arg_list()` in `GIQLDistance` to parse `stranded` named parameter
- [ ] T044 [US3] Update `giqldistance_sql()` to evaluate `stranded` parameter (boolean AST node to bool value)
- [ ] T045 [US3] Call `_get_column_refs()` with `include_strand=True` when `stranded=true`
- [ ] T046 [US3] Add strand comparison clause to SQL CASE expression when `stranded=true`
- [ ] T047 [US3] Handle missing strand column error with clear error message
- [ ] T048 [US3] Add documentation for stranded mode to quickstart.md (already present, verify accuracy)
- [ ] T049 [US3] Verify all T037-T042 tests pass

**Checkpoint**: All three user stories work independently - basic distance, finding nearest, and strand-specific distance

---

## Phase 6: User Story 4 - Signed Distance (Priority: P4)

**Goal**: Return directional distance (negative=upstream, positive=downstream) using `signed=true` parameter, mirroring bedtools `-D ref` flag

**Independent Test**: Pass `signed=true` to DISTANCE() and verify negative values for upstream features, positive for downstream

### Tests for User Story 4

- [ ] T050 [P] [US4] Unit test for upstream distance (negative value) in `tests/unit/test_distance_udf.py`
- [ ] T051 [P] [US4] Unit test for downstream distance (positive value) in `tests/unit/test_distance_udf.py`
- [ ] T052 [P] [US4] Unit test for overlapping intervals with signed=true (still returns 0) in `tests/unit/test_distance_udf.py`
- [ ] T053 [P] [US4] Unit test for combined `stranded=true, signed=true` in `tests/unit/test_distance_udf.py`
- [ ] T054 [P] [US4] Parser test for named parameter `signed=true` in `tests/unit/test_distance_parsing.py`
- [ ] T055 [P] [US4] Parser test for both parameters `stranded=true, signed=true` in `tests/unit/test_distance_parsing.py`
- [ ] T056 [P] [US4] Transpilation test for signed distance SQL in `tests/unit/test_distance_transpilation.py`
- [ ] T057 [US4] Integration test for signed DISTANCE queries in `tests/integration/test_distance_queries.py`
- [ ] T058 [US4] Integration test for combined stranded+signed queries in `tests/integration/test_distance_queries.py`

### Implementation for User Story 4

- [ ] T059 [US4] Update `from_arg_list()` in `GIQLDistance` to parse `signed` named parameter
- [ ] T060 [US4] Update `giqldistance_sql()` to evaluate `signed` parameter
- [ ] T061 [US4] Modify SQL CASE expression to negate distance when A is after B and signed=true
- [ ] T062 [US4] Ensure overlapping intervals still return 0 regardless of signed parameter
- [ ] T063 [US4] Implement FR-013: signed distance always relative to reference genome coordinates (not strand-specific)
- [ ] T064 [US4] Add documentation for signed mode to quickstart.md (already present, verify accuracy)
- [ ] T065 [US4] Add bedtools equivalence mapping to README showing `closest -D ref` mapping
- [ ] T066 [US4] Verify all T050-T058 tests pass

**Checkpoint**: All user stories complete and independently functional

---

## Phase 7: Polish & Cross-Cutting Concerns

**Purpose**: Improvements that affect multiple user stories, documentation, and finalization

- [ ] T067 [P] Add DISTANCE operator to README.md with basic examples
- [ ] T068 [P] Add bedtools `closest` equivalence table to README.md
- [ ] T069 [P] Update CLI help text if DISTANCE needs to be mentioned
- [ ] T070 [P] Add DISTANCE examples to demo.ipynb (create cells showing US1-US4 patterns)
- [ ] T071 [P] Add performance notes to README about pre-filtering by chromosome
- [ ] T072 [P] Add note about no pre-sorting requirement (unlike bedtools) to README
- [ ] T073 Code review: Check all SQL generation follows readable formatting
- [ ] T074 Code review: Verify error messages are clear and actionable
- [ ] T075 Code review: Ensure constitution compliance (expressive, readable, portable, canonical)
- [ ] T076 Run all tests across all three dialects (DuckDB, SQLite, PostgreSQL)
- [ ] T077 Run ruff linting and formatting on all modified files
- [ ] T078 Verify test coverage reaches 100% for public DISTANCE API
- [ ] T079 [P] Property-based test for symmetric distance (unsigned) in `tests/property/test_distance_properties.py`
- [ ] T080 [P] Property-based test for overlap detection correctness in `tests/property/test_distance_properties.py`
- [ ] T081 [P] Property-based test for NULL handling across all input combinations in `tests/property/test_distance_properties.py`
- [ ] T082 Run quickstart.md validation (execute all example queries manually or via script)
- [ ] T083 Update CHANGELOG or release notes with DISTANCE operator addition

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies - can start immediately
- **Foundational (Phase 2)**: Depends on Setup completion - BLOCKS all user stories
- **User Stories (Phase 3-6)**: All depend on Foundational phase completion
  - User Story 1 (P1): Basic distance - No dependencies on other stories
  - User Story 2 (P2): Finding nearest - Depends on US1 (builds on basic DISTANCE)
  - User Story 3 (P3): Stranded distance - Independent of US2, depends on US1
  - User Story 4 (P4): Signed distance - Independent of US2/US3, depends on US1
- **Polish (Phase 7)**: Depends on all desired user stories being complete

### User Story Dependencies

```
Foundation (Phase 2)
    ↓
User Story 1 (P1) - Basic DISTANCE
    ↓
    ├─→ User Story 2 (P2) - Finding Nearest (builds on US1)
    ├─→ User Story 3 (P3) - Stranded (extends US1)
    └─→ User Story 4 (P4) - Signed (extends US1)
```

### Within Each User Story

1. Tests MUST be written FIRST and FAIL before implementation
2. Parser tests before SQL generation tests
3. Unit tests before integration tests
4. Core logic before edge cases
5. Story complete before moving to next priority

### Parallel Opportunities

**Setup Phase**:
- All 3 tasks can run in parallel (directory creation is independent)

**Foundational Phase**:
- T005, T006 can run in parallel (different files)
- T004 must complete before T007

**User Story 1**:
- Tests T008-T017 can all run in parallel (different test files/functions)
- Implementation T022-T023 can run in parallel (handle both input types)

**User Story 2**:
- Tests T029-T032 can run in parallel

**User Story 3**:
- Tests T037-T041 can run in parallel

**User Story 4**:
- Tests T050-T058 can run in parallel

**Polish Phase**:
- Documentation tasks T067-T072 can run in parallel
- Property-based tests T079-T081 can run in parallel

**Cross-Story Parallelism**:
- Once US1 completes, US2, US3, and US4 can proceed in parallel (if team capacity allows)
- US3 and US4 are fully independent - can be worked on simultaneously

---

## Parallel Example: User Story 1

```bash
# Launch all unit tests for User Story 1 together:
Task T008: "Unit test for overlapping intervals in tests/unit/test_distance_udf.py"
Task T009: "Unit test for non-overlapping intervals in tests/unit/test_distance_udf.py"
Task T010: "Unit test for different chromosomes in tests/unit/test_distance_udf.py"
Task T011: "Unit test for adjacent intervals in tests/unit/test_distance_udf.py"
Task T012: "Unit test for zero-width intervals in tests/unit/test_distance_udf.py"

# Launch parser tests in parallel:
Task T013: "Parser test for DISTANCE(a.position, b.position) in tests/unit/test_distance_parsing.py"
Task T014: "Parser test for literal range in tests/unit/test_distance_parsing.py"

# Launch transpilation tests for all dialects in parallel:
Task T015: "Transpilation test for DuckDB in tests/unit/test_distance_transpilation.py"
Task T016: "Transpilation test for SQLite in tests/unit/test_distance_transpilation.py"
Task T017: "Transpilation test for PostgreSQL in tests/unit/test_distance_transpilation.py"
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (T001-T003)
2. Complete Phase 2: Foundational (T004-T007) - CRITICAL BLOCKER
3. Complete Phase 3: User Story 1 (T008-T028)
4. **STOP and VALIDATE**: Run all US1 tests, execute example queries, verify basic DISTANCE works
5. Deploy/demo if ready - **This is a functional MVP!**

### Incremental Delivery

1. Foundation (Phase 1-2) → Core DISTANCE infrastructure ready
2. Add US1 → Test independently → **Deploy/Demo (MVP!)**
3. Add US2 → Test independently → Deploy/Demo (now supports "find closest" workflow)
4. Add US3 → Test independently → Deploy/Demo (now supports strand-specific queries)
5. Add US4 → Test independently → Deploy/Demo (now supports directional distance)
6. Each story adds value without breaking previous stories

### Parallel Team Strategy

With multiple developers:

1. Team completes Setup + Foundational together (Phase 1-2)
2. Once Foundational is done:
   - **Developer A**: User Story 1 (T008-T028)
   - Wait for US1 to complete (it's a dependency)
3. After US1 completes:
   - **Developer A**: User Story 2 (T029-T036)
   - **Developer B**: User Story 3 (T037-T049)
   - **Developer C**: User Story 4 (T050-T066)
4. Stories complete and integrate independently

---

## Testing Strategy Summary

### Test Coverage Requirements (per Constitution)

- **100% coverage** of public API (DISTANCE function)
- **Behavior-focused** tests (not implementation details)
- **Given-When-Then** docstrings for all test functions
- **AAA pattern** (Arrange-Act-Assert) with clear visual separation

### Test Types by Phase

| Phase | Unit Tests | Integration Tests | Property Tests |
|-------|------------|-------------------|----------------|
| US1 - Basic | 5 (overlap, gaps, NULL, edge cases) | 1 (basic query) | 0 |
| US2 - Nearest | 0 | 4 (window functions, ties, NULL ordering) | 1 (invariants) |
| US3 - Stranded | 3 (same/diff strand, unspecified) | 1 (stranded queries) | 0 |
| US4 - Signed | 5 (upstream/downstream/combined) | 2 (signed, stranded+signed) | 0 |
| Polish | 0 | 0 | 3 (symmetric, overlap, NULL) |

### Total Test Count

- **Unit Tests**: 18 (distance logic, parsing, transpilation)
- **Integration Tests**: 8 (end-to-end query execution)
- **Property-Based Tests**: 4 (Hypothesis-generated edge cases)
- **Total**: 30 tests

---

## Notes

- [P] tasks = different files, no dependencies - can run in parallel
- [Story] label maps task to specific user story for traceability
- Each user story should be independently completable and testable
- **TDD approach**: Write tests first, verify they fail, then implement
- Commit after each task or logical group
- Stop at any checkpoint to validate story independently
- Avoid: vague tasks, same file conflicts, cross-story dependencies that break independence
- **Constitution alignment**: All tasks designed to maintain expressive, readable, portable, canonical principles

---

## Quick Reference

**Total Tasks**: 83
**MVP Tasks** (US1 only): 28 (Setup + Foundation + US1)
**Estimated LOC**: ~500-800 (implementation + tests)

**Task Breakdown by Phase**:
- Phase 1 (Setup): 3 tasks
- Phase 2 (Foundational): 4 tasks
- Phase 3 (US1 - Basic): 21 tasks (11 tests + 10 implementation)
- Phase 4 (US2 - Nearest): 8 tasks (4 tests + 4 implementation)
- Phase 5 (US3 - Stranded): 13 tasks (6 tests + 7 implementation)
- Phase 6 (US4 - Signed): 17 tasks (9 tests + 8 implementation)
- Phase 7 (Polish): 17 tasks (documentation + property tests + validation)

**Key Files to Modify**:
- `src/giql/expressions.py` - Add GIQLDistance class
- `src/giql/dialect.py` - Add DISTANCE parsing
- `src/giql/generators/base.py` - Add SQL generation logic
- `tests/*` - Comprehensive test suite (30+ tests)
- `README.md`, `demo.ipynb` - Documentation and examples
