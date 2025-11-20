# Tasks: Bedtools Integration Test Suite

**Input**: Design documents from `/specs/003-bedtools-integration-tests/`
**Prerequisites**: plan.md (required), spec.md (required for user stories), research.md, data-model.md, contracts/

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3, US4)
- Include exact file paths in descriptions

## Path Conventions

All paths are relative to repository root. Test files under `tests/integration/bedtools/`.

---

## Phase 1: Setup (Test Infrastructure)

**Purpose**: Initialize test directory structure and install dependencies

- [x] T001 Create test directory structure: tests/integration/bedtools/ and tests/integration/bedtools/utils/
- [x] T002 Create empty __init__.py files in tests/integration/bedtools/ and tests/integration/bedtools/utils/
- [x] T003 Verify pytest is installed and accessible (pytest --version) - NOTE: Not installed, tests require `pytest` package
- [x] T004 Verify DuckDB Python driver is installed (pip list | grep duckdb) - NOTE: Not installed, tests require `duckdb` package
- [x] T005 Verify bedtools 2.30.0+ is installed and in PATH (bedtools --version) - NOTE: Not installed, tests require bedtools 2.30.0+

---

## Phase 2: Foundational (Core Test Infrastructure)

**Purpose**: Core test utilities and fixtures that ALL user stories depend on

**⚠️ CRITICAL**: No user story test implementation can begin until this phase is complete

- [x] T006 [P] Implement GenomicInterval dataclass in tests/integration/bedtools/utils/data_models.py
- [x] T007 [P] Implement SimulatedDataset dataclass in tests/integration/bedtools/utils/data_models.py
- [x] T008 [P] Implement ComparisonResult dataclass in tests/integration/bedtools/utils/data_models.py
- [x] T009 [P] Implement IntervalGeneratorConfig dataclass in tests/integration/bedtools/utils/data_models.py
- [x] T010 [P] Implement BedtoolsVersion dataclass in tests/integration/bedtools/utils/data_models.py
- [x] T011 Implement IntervalGenerator class with seeded random generation in tests/integration/bedtools/utils/interval_generator.py
- [x] T012 [P] Implement BED file export function (export_to_bed) in tests/integration/bedtools/utils/bed_export.py
- [x] T013 [P] Implement DuckDB table loading function (load_intervals) in tests/integration/bedtools/utils/bed_export.py
- [x] T014 Implement bedtools version checking in tests/integration/bedtools/utils/bedtools_wrapper.py
- [x] T015 Implement bedtools command execution (run_bedtools_command) in tests/integration/bedtools/utils/bedtools_wrapper.py
- [x] T016 Implement bedtools output parsing (parse_bedtools_output) in tests/integration/bedtools/utils/bedtools_wrapper.py
- [x] T017 Implement result comparison logic (compare_results) in tests/integration/bedtools/utils/comparison.py
- [x] T018 Implement order-independent row sorting in tests/integration/bedtools/utils/comparison.py
- [x] T019 Implement epsilon float comparison in tests/integration/bedtools/utils/comparison.py
- [x] T020 Create bedtools_version session-scoped fixture in tests/integration/bedtools/conftest.py
- [x] T021 Create duckdb_connection function-scoped fixture in tests/integration/bedtools/conftest.py
- [x] T022 Create temp_bed_dir function-scoped fixture in tests/integration/bedtools/conftest.py
- [x] T023 Create interval_generator function-scoped fixture in tests/integration/bedtools/conftest.py

**Checkpoint**: Foundation ready - user story test implementation can now begin in parallel

---

## Phase 3: User Story 1 - Basic Interval Operations Validation (Priority: P1) 🎯 MVP

**Goal**: Validate that GIQL's core interval operations (intersect, merge, nearest) produce identical results to bedtools

**Independent Test**: Run pytest tests/integration/bedtools/test_intersect.py tests/integration/bedtools/test_merge.py tests/integration/bedtools/test_nearest.py and verify all pass

### Implementation for User Story 1

- [x] T024 [P] [US1] Implement test_intersect_basic_overlap in tests/integration/bedtools/test_intersect.py
- [x] T025 [P] [US1] Implement test_intersect_partial_overlap in tests/integration/bedtools/test_intersect.py
- [x] T026 [P] [US1] Implement test_intersect_no_overlap in tests/integration/bedtools/test_intersect.py
- [x] T027 [P] [US1] Implement test_intersect_adjacent_intervals in tests/integration/bedtools/test_intersect.py
- [x] T028 [P] [US1] Implement test_intersect_multiple_chromosomes in tests/integration/bedtools/test_intersect.py
- [x] T029 [P] [US1] Implement test_merge_adjacent_intervals in tests/integration/bedtools/test_merge.py
- [x] T030 [P] [US1] Implement test_merge_overlapping_intervals in tests/integration/bedtools/test_merge.py
- [x] T031 [P] [US1] Implement test_merge_separated_intervals in tests/integration/bedtools/test_merge.py
- [x] T032 [P] [US1] Implement test_merge_multiple_chromosomes in tests/integration/bedtools/test_merge.py
- [x] T033 [P] [US1] Implement test_nearest_non_overlapping in tests/integration/bedtools/test_nearest.py
- [x] T034 [P] [US1] Implement test_nearest_multiple_candidates in tests/integration/bedtools/test_nearest.py
- [x] T035 [P] [US1] Implement test_nearest_cross_chromosome in tests/integration/bedtools/test_nearest.py
- [x] T036 [P] [US1] Implement test_nearest_boundary_cases in tests/integration/bedtools/test_nearest.py
- [x] T037 [US1] Add interval generator methods for overlapping scenarios in tests/integration/bedtools/utils/interval_generator.py
- [x] T038 [US1] Add interval generator methods for adjacent scenarios in tests/integration/bedtools/utils/interval_generator.py
- [x] T039 [US1] Add interval generator methods for separated scenarios in tests/integration/bedtools/utils/interval_generator.py
- [x] T040 [US1] Add interval generator methods for multi-chromosome scenarios in tests/integration/bedtools/utils/interval_generator.py

**Checkpoint**: At this point, basic interval operations (intersect, merge, nearest) are validated against bedtools. This represents the minimum viable test suite (MVP).

---

## Phase 4: User Story 2 - Strand-Aware Operations (Priority: P2)

**Goal**: Verify that GIQL correctly handles strand-specific interval operations

**Independent Test**: Run pytest tests/integration/bedtools/test_strand_aware.py and verify all pass

### Implementation for User Story 2

- [x] T041 [P] [US2] Implement test_intersect_same_strand in tests/integration/bedtools/test_strand_aware.py
- [x] T042 [P] [US2] Implement test_intersect_opposite_strand in tests/integration/bedtools/test_strand_aware.py
- [x] T043 [P] [US2] Implement test_intersect_ignore_strand in tests/integration/bedtools/test_strand_aware.py
- [x] T044 [P] [US2] Implement test_intersect_mixed_strands in tests/integration/bedtools/test_strand_aware.py
- [x] T045 [P] [US2] Implement test_nearest_same_strand in tests/integration/bedtools/test_strand_aware.py
- [x] T046 [P] [US2] Implement test_nearest_opposite_strand in tests/integration/bedtools/test_strand_aware.py
- [x] T047 [P] [US2] Implement test_nearest_ignore_strand in tests/integration/bedtools/test_strand_aware.py
- [x] T048 [P] [US2] Implement test_merge_strand_specific in tests/integration/bedtools/test_strand_aware.py
- [x] T049 [US2] Add interval generator methods for strand combinations in tests/integration/bedtools/utils/interval_generator.py
- [x] T050 [US2] Add bedtools wrapper methods for strand flags (-s, -S) in tests/integration/bedtools/utils/bedtools_wrapper.py

**Checkpoint**: At this point, strand-aware operations are validated. Test suite now covers both basic operations (US1) and strand-specific behaviors (US2).

---

## Phase 5: User Story 3 - Distance Calculations (Priority: P2)

**Goal**: Validate that GIQL's distance calculations match bedtools

**Independent Test**: Run pytest tests/integration/bedtools/test_distance.py and verify all pass

### Implementation for User Story 3

- [ ] T051 [P] [US3] Implement test_distance_non_overlapping_intervals in tests/integration/bedtools/test_distance.py
- [ ] T052 [P] [US3] Implement test_distance_overlapping_intervals in tests/integration/bedtools/test_distance.py
- [ ] T053 [P] [US3] Implement test_distance_adjacent_intervals in tests/integration/bedtools/test_distance.py
- [ ] T054 [P] [US3] Implement test_distance_cross_chromosome in tests/integration/bedtools/test_distance.py
- [ ] T055 [P] [US3] Implement test_distance_large_gaps in tests/integration/bedtools/test_distance.py
- [ ] T056 [P] [US3] Implement test_distance_with_strand_consideration in tests/integration/bedtools/test_distance.py
- [ ] T057 [US3] Add bedtools wrapper method for distance reporting (-d flag) in tests/integration/bedtools/utils/bedtools_wrapper.py
- [ ] T058 [US3] Add floating-point distance comparison handling in tests/integration/bedtools/utils/comparison.py

**Checkpoint**: At this point, distance calculations are validated. Test suite now covers basic operations (US1), strand-aware (US2), and distance calculations (US3).

---

## Phase 6: User Story 4 - Complex Multi-Operation Workflows (Priority: P3)

**Goal**: Validate that chained GIQL operations produce results matching bedtools pipelines

**Independent Test**: Run pytest tests/integration/bedtools/test_workflows.py and verify all pass

### Implementation for User Story 4

- [ ] T059 [P] [US4] Implement test_workflow_intersect_then_merge in tests/integration/bedtools/test_workflows.py
- [ ] T060 [P] [US4] Implement test_workflow_nearest_then_filter_distance in tests/integration/bedtools/test_workflows.py
- [ ] T061 [P] [US4] Implement test_workflow_intersect_merge_nearest in tests/integration/bedtools/test_workflows.py
- [ ] T062 [P] [US4] Implement test_workflow_strand_specific_pipeline in tests/integration/bedtools/test_workflows.py
- [ ] T063 [US4] Add bedtools pipeline execution support (piped commands) in tests/integration/bedtools/utils/bedtools_wrapper.py
- [ ] T064 [US4] Add intermediate result validation in tests/integration/bedtools/utils/comparison.py

**Checkpoint**: All user stories are now complete. Test suite validates basic operations, strand-aware operations, distance calculations, and complex workflows.

---

## Phase 7: Edge Cases & Comprehensive Coverage

**Purpose**: Add edge case coverage and validate comprehensive scenarios

- [ ] T065 [P] Implement test for non-standard chromosome names (chrM, chr1_random) in tests/integration/bedtools/test_edge_cases.py
- [ ] T066 [P] Implement test for very large genomic intervals in tests/integration/bedtools/test_edge_cases.py
- [ ] T067 [P] Implement test for chromosome boundary cases in tests/integration/bedtools/test_edge_cases.py
- [ ] T068 [P] Implement test for high-density overlapping intervals in tests/integration/bedtools/test_edge_cases.py
- [ ] T069 Add interval generator method for boundary cases in tests/integration/bedtools/utils/interval_generator.py
- [ ] T070 Add interval generator method for very large intervals in tests/integration/bedtools/utils/interval_generator.py

---

## Phase 8: Polish & Documentation

**Purpose**: Improvements that affect multiple user stories and documentation

- [ ] T071 [P] Add comprehensive docstrings to all test functions with Given-When-Then format
- [ ] T072 [P] Add type hints to all utility functions and classes
- [ ] T073 [P] Create test_data_generator.py with reusable dataset generation patterns
- [ ] T074 Document bedtools features NOT replicated in GIQL in tests/integration/bedtools/README.md
- [ ] T075 Add test execution instructions to tests/integration/bedtools/README.md
- [ ] T076 Add pytest markers for slow tests (large datasets) in conftest.py
- [ ] T077 Add pytest parameterization for different dataset sizes in conftest.py
- [ ] T078 Validate all tests complete in <5 minutes (run pytest with timing)
- [ ] T079 Validate dataset generation completes in <10 seconds for 10k intervals
- [ ] T080 Run full test suite and verify 100% pass rate

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies - can start immediately
- **Foundational (Phase 2)**: Depends on Setup completion - BLOCKS all user stories
- **User Stories (Phase 3-6)**: All depend on Foundational phase completion
  - User stories can then proceed in parallel (if staffed)
  - Or sequentially in priority order (P1 → P2 → P3)
- **Edge Cases (Phase 7)**: Can start after US1 completion, benefits from all user stories
- **Polish (Phase 8)**: Depends on all user stories being complete

### User Story Dependencies

- **User Story 1 (P1)**: Can start after Foundational (Phase 2) - No dependencies on other stories
- **User Story 2 (P2)**: Can start after Foundational (Phase 2) - Independent, but builds on US1 test patterns
- **User Story 3 (P2)**: Can start after Foundational (Phase 2) - Independent, but builds on US1 test patterns
- **User Story 4 (P3)**: Can start after Foundational (Phase 2) - May reference US1/US2/US3 tests but is independently testable

### Within Each User Story

- All test files marked [P] can be implemented in parallel
- Interval generator enhancements should follow test implementation to ensure tests drive requirements
- Bedtools wrapper enhancements should follow test implementation

### Parallel Opportunities

- **Phase 1**: All setup tasks can run in parallel
- **Phase 2**: All data model classes (T006-T010) can be implemented in parallel
- **Phase 2**: Export/wrapper/comparison utilities (T012-T019) can be implemented in parallel after data models
- **Phase 2**: All fixtures (T020-T023) can be implemented in parallel after utilities
- **Phase 3**: All US1 test files (T024-T036) can be implemented in parallel
- **Phase 4**: All US2 tests (T041-T048) can be implemented in parallel
- **Phase 5**: All US3 tests (T051-T056) can be implemented in parallel
- **Phase 6**: All US4 tests (T059-T062) can be implemented in parallel
- **Once Foundational phase completes, all user stories (Phase 3-6) can start in parallel if team capacity allows**

---

## Parallel Example: User Story 1

```bash
# Launch all test files for User Story 1 together:
Task T024: "Implement test_intersect_basic_overlap in tests/integration/bedtools/test_intersect.py"
Task T025: "Implement test_intersect_partial_overlap in tests/integration/bedtools/test_intersect.py"
Task T026: "Implement test_intersect_no_overlap in tests/integration/bedtools/test_intersect.py"
Task T027: "Implement test_intersect_adjacent_intervals in tests/integration/bedtools/test_intersect.py"
Task T028: "Implement test_intersect_multiple_chromosomes in tests/integration/bedtools/test_intersect.py"
Task T029: "Implement test_merge_adjacent_intervals in tests/integration/bedtools/test_merge.py"
Task T030: "Implement test_merge_overlapping_intervals in tests/integration/bedtools/test_merge.py"
Task T031: "Implement test_merge_separated_intervals in tests/integration/bedtools/test_merge.py"
Task T032: "Implement test_merge_multiple_chromosomes in tests/integration/bedtools/test_merge.py"
Task T033: "Implement test_nearest_non_overlapping in tests/integration/bedtools/test_nearest.py"
Task T034: "Implement test_nearest_multiple_candidates in tests/integration/bedtools/test_nearest.py"
Task T035: "Implement test_nearest_cross_chromosome in tests/integration/bedtools/test_nearest.py"
Task T036: "Implement test_nearest_boundary_cases in tests/integration/bedtools/test_nearest.py"

# Then implement supporting generator methods (depends on tests defining requirements):
Task T037: "Add interval generator methods for overlapping scenarios"
Task T038: "Add interval generator methods for adjacent scenarios"
Task T039: "Add interval generator methods for separated scenarios"
Task T040: "Add interval generator methods for multi-chromosome scenarios"
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (T001-T005)
2. Complete Phase 2: Foundational (T006-T023) - CRITICAL, blocks all stories
3. Complete Phase 3: User Story 1 (T024-T040)
4. **STOP and VALIDATE**: Run pytest tests/integration/bedtools/test_intersect.py tests/integration/bedtools/test_merge.py tests/integration/bedtools/test_nearest.py
5. Verify all US1 tests pass independently

**Result**: Minimum viable test suite that validates core interval operations (intersect, merge, nearest) against bedtools

### Incremental Delivery

1. Complete Setup + Foundational → Foundation ready
2. Add User Story 1 → Test independently → MVP ready (covers 60% of bedtools operations)
3. Add User Story 2 → Test independently → Strand-aware coverage added
4. Add User Story 3 → Test independently → Distance calculations added
5. Add User Story 4 → Test independently → Workflow validation complete
6. Add Edge Cases → Comprehensive coverage (80%+ goal met)
7. Polish & Document → Production ready

Each increment adds value and maintains all previous functionality.

### Parallel Team Strategy

With multiple developers:

1. Team completes Setup + Foundational together (T001-T023)
2. Once Foundational is done:
   - Developer A: User Story 1 (T024-T040) - MVP priority
   - Developer B: User Story 2 (T041-T050)
   - Developer C: User Story 3 (T051-T058)
   - Developer D: User Story 4 (T059-T064)
3. Stories complete and integrate independently
4. Team collaborates on Edge Cases (T065-T070)
5. Team collaborates on Polish (T071-T080)

---

## Task Summary

**Total Tasks**: 80

**By Phase**:
- Phase 1 (Setup): 5 tasks
- Phase 2 (Foundational): 18 tasks (CRITICAL - blocks all user stories)
- Phase 3 (US1 - Basic Operations): 17 tasks
- Phase 4 (US2 - Strand-Aware): 10 tasks
- Phase 5 (US3 - Distance): 8 tasks
- Phase 6 (US4 - Workflows): 6 tasks
- Phase 7 (Edge Cases): 6 tasks
- Phase 8 (Polish): 10 tasks

**By User Story**:
- US1 (Basic Operations): 17 tasks - MVP
- US2 (Strand-Aware): 10 tasks
- US3 (Distance): 8 tasks
- US4 (Workflows): 6 tasks
- Infrastructure: 39 tasks (Setup + Foundational + Edge Cases + Polish)

**Parallel Opportunities**: 58 tasks marked [P] can run in parallel (within their phase constraints)

**Independent Test Criteria**:
- US1: pytest tests/integration/bedtools/test_{intersect,merge,nearest}.py
- US2: pytest tests/integration/bedtools/test_strand_aware.py
- US3: pytest tests/integration/bedtools/test_distance.py
- US4: pytest tests/integration/bedtools/test_workflows.py

**MVP Scope**: Phase 1 + Phase 2 + Phase 3 (User Story 1) = 40 tasks
- Delivers core interval operation validation (intersect, merge, nearest)
- Provides foundation for all subsequent enhancements
- Independently testable and valuable

---

## Format Validation

✅ All tasks follow checklist format: `- [ ] [ID] [P?] [Story?] Description with file path`
✅ All tasks have sequential IDs (T001-T080)
✅ All user story tasks have [US#] labels
✅ All parallel tasks have [P] markers
✅ All tasks include specific file paths
✅ Tasks organized by user story for independent implementation
