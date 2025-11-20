<!--
SYNC IMPACT REPORT
==================
Version Change: N/A → 1.0.0
Change Type: MAJOR - Initial constitution creation for GIQL project
Date: 2025-11-03

Modified Principles:
- Initial creation with 4 core principles: Expressive, Readable, Portable, Canonical

Added Sections:
- Core Principles (4 principles)
- Code Quality & Testing (philosophy, requirements, test structure, property-based testing)
- Release & Versioning (software versioning, backwards compatibility)
- Governance (amendment process, constitution versioning, compliance review)

Removed Sections: N/A (initial creation)

Templates Requiring Updates:
- ✅ .specify/templates/* - Will reference these principles

Follow-up TODOs: None
-->

# GIQL Project Constitution

## Core Principles

### I. Expressive
GIQL should make common genomic operations clear, concise and natural to express. It SHOULD reduce SQL boilerplate for the operations and predicates that genomics practitioners use most frequently.

**Rules:**
- Common genomic operations (intersect, merge, cluster, within, contains, etc.) MUST have dedicated operators
- Genomic interval syntax MUST support standard formats (e.g., 'chr1:1000-2000')
- Complex genomic workflows SHOULD be expressible in a single query
- Operators MUST compose naturally with standard SQL (WHERE, JOIN, GROUP BY, etc.)

**Rationale:** Genomic analysis involves repetitive patterns that are cumbersome in standard SQL. GIQL exists to reduce boilerplate and make genomic queries more natural. If users still need to write complex SQL for common operations, GIQL fails its purpose.

### II. Readable
Like Python, GIQL queries should clearly express their intent even to someone who isn't proficient in SQL. Queries should communicate "what" they do rather than "how" they do it. A bioinformatician should be able to understand a query's purpose without being a database expert.

**Rules:**
- Operator names MUST be self-describing (INTERSECTS, WITHIN, CONTAINS, MERGE, CLUSTER)
- Query structure SHOULD mirror the logical flow of genomic analysis
- Generated SQL SHOULD be human-readable when transpiled (for debugging and learning)
- Documentation MUST include real-world examples with biological context
- Error messages MUST be clear and actionable, not database internals

**Rationale:** A domain-specific query language that requires deep SQL expertise creates barriers to adoption. Readable queries reduce cognitive load, improve collaboration, and make code review more effective.

### III. Portable
GIQL should strive to be usable with any database query engine that supports basic SQL. GIQL does not prescribe which underlying database engine users should adopt—users should be free to choose the best database for their needs without rewriting analysis code.

**Rules:**
- GIQL MUST transpile to standard SQL compatible with multiple backends
- Database-specific features MUST NOT be required for core operators
- Supported dialects MUST include at least: DuckDB, SQLite, PostgreSQL
- Transpilation MUST preserve query semantics across all supported backends
- New features MUST work on all supported backends or be clearly marked as dialect-specific
- Dependencies on specific database features SHOULD be minimized

**Rationale:** The genomics community uses diverse database systems depending on scale, infrastructure, and existing tooling. A portable query language protects against vendor lock-in and ensures longevity. Users' investment in learning GIQL should not be tied to a specific database vendor.

### IV. Canonical
GIQL should implement canonical operators and adhere to established conventions from widely-used genomics tools. Where multiple conventions exist, GIQL should follow the most prevalent or well-established approach.

**Rules:**
- Operator semantics MUST match established tools (bedtools, bioframe, etc.) when possible
- Naming conventions SHOULD align with community standards
- Edge cases (e.g., zero-width intervals, chromosome boundaries) MUST be handled consistently with established tools
- Deviations from conventions MUST be explicitly documented with rationale
- When in doubt, bedtools behavior is the reference standard

**Rationale:** Genomics practitioners already have mental models from existing tools. Following established conventions reduces the learning curve and makes GIQL predictable. Canonical implementations also ensure correctness through proven semantics that the community has validated over years of use.

## Code Quality & Testing

### Testing Philosophy: Test Behavior, Not Implementation
Tests MUST verify what code does, not how it does it. Developers should be able to completely refactor internal implementation (e.g., AST transformation logic, SQL generation) while tests continue to pass, as long as external behavior remains unchanged.

**Key Principles:**
- Test GIQL query → SQL transpilation output, not transformation internals
- Test query execution results, not intermediate data structures
- Test operator semantics (does INTERSECTS correctly find overlaps?), not implementation details
- Refactoring should not require rewriting tests

### Testing Requirements

**Coverage:**
- MUST achieve 100% test coverage of public APIs
- SHOULD minimize test overlap—each test should verify distinct behavior
- Tests SHOULD NOT touch internal implementation details (e.g., AST node types)
- Every pull request MUST maintain or improve test coverage
- Each GIQL operator MUST have comprehensive tests covering edge cases

**Test structure:**
- MUST follow AAA (Arrange-Act-Assert) pattern with clear visual separation
- MUST use Given-When-Then docstring format for all test functions/methods
- MUST mirror source directory structure in tests/ directory
- Test function naming: `test_<function_name>_<scenario>`
- Test file naming: `test_<module_name>.py`

**Example test structure:**
```python
def test_intersects_basic_overlap():
    """Test INTERSECTS predicate finds overlapping intervals.

    Given:
        Two tables with genomic intervals where some intervals overlap
    When:
        A GIQL query uses INTERSECTS predicate in WHERE clause
    Then:
        Only overlapping intervals should be returned
    """
    # Arrange
    engine = GIQLEngine(target_dialect="duckdb")
    engine.register_table_schema(...)

    # Act
    result = engine.execute("SELECT * FROM a WHERE a.position INTERSECTS b.position")

    # Assert
    assert len(result) == expected_count
    assert all(intervals_actually_overlap(row) for row in result)
```

### Test Isolation and Independence

**Isolation principles:**
- Each test MUST be independent (no test should depend on state from another)
- Tests MUST be deterministic (same input → same output, always)
- SHOULD mock only at system boundaries (database connections, file I/O)
- SHOULD NOT mock internal methods—test the complete unit of behavior

**Database testing:**
- Use in-memory databases (DuckDB, SQLite) for fast, isolated tests
- Tests SHOULD create fresh database instances (avoid shared state)
- Tests SHOULD NOT require external database servers
- Integration tests with real databases SHOULD be clearly marked and optional

**When to use dependencies:**
- SHOULD use `@pytest.mark.dependency` to declare test dependencies and minimize failure overlaps (dependencies prevent cascade failures from masking root causes)

**Mocking guidelines:**
- SHOULD mock where used, not where defined (`mocker.patch('module_using_it.object')`)
- SHOULD use `autospec=True` to prevent mock drift and catch interface mismatches
- SHOULD avoid over-mocking (only mock essential boundaries)
- SHOULD use descriptive names: `mock_<name>` for mocks, `spy_<name>` for spies

### Property-Based Testing with Hypothesis

Use Hypothesis for:
- Testing query transpilation invariants (valid GIQL → valid SQL)
- Roundtrip operations (parse → generate → parse should be idempotent)
- Testing operator semantics across diverse inputs
- Generating edge cases for genomic intervals (empty, zero-width, chromosome boundaries)

**Example:**
```python
from hypothesis import given
from hypothesis import strategies as st

@given(st.lists(genomic_interval_strategy(), min_size=2))
def test_merge_produces_non_overlapping_intervals(intervals):
    """Test that MERGE operator produces non-overlapping output.

    Given:
        Any list of genomic intervals
    When:
        MERGE operator is applied
    Then:
        Resulting intervals should not overlap each other
    """
    result = execute_merge_query(intervals)
    assert all(not overlaps(a, b) for a, b in consecutive_pairs(result))
```

### Transpilation Testing

**Critical tests for transpilation:**
- Input GIQL query → Output SQL for each supported dialect
- Verify SQL is valid (can be parsed by target database)
- Verify SQL produces semantically correct results
- Test operator composition (INTERSECTS + GROUP BY + ORDER BY, etc.)
- Test edge cases: NULL handling, empty results, chromosome name variations

**SQL validation:**
- Generated SQL MUST be syntactically valid for target dialect
- Generated SQL SHOULD be formatted for readability
- Generated SQL SHOULD NOT include unnecessary complexity (no redundant subqueries)

### Pytest Fixtures and Test Data

**Fixture usage:**
- Use fixtures for reusable setup/teardown logic
- Place shared fixtures in `conftest.py` files
- Choose appropriate scope (function, class, module, session)
- Chain fixtures to build complex test environments incrementally

**Fixture scopes:**
- **Function** (default): Maximum isolation, new instance per test
- **Class**: Shared across test class methods
- **Module**: Shared across all tests in a file
- **Session**: Shared across entire test run (use sparingly)

**Test data management:**
- Define simple data inline within tests
- Use parameterized tests (`@pytest.mark.parametrize`) for multiple scenarios
- Create custom Hypothesis strategies for domain objects

### Test Performance

**Speed principles:**
- Unit tests should execute rapidly (milliseconds, not seconds)
- Avoid `await asyncio.sleep()` or `time.sleep()` when possible
- Fast tests encourage frequent execution during development

### Pull Request Requirements
All contributions MUST meet these standards before merge:
- **Test coverage**: New features MUST include comprehensive unit tests
- **Existing tests**: All existing tests MUST pass (CI/CD enforced)
- **Test quality**: Proper isolation, appropriate fixtures, naming conventions
- **Documentation**: Public APIs MUST have complete docstrings with examples
- **Transpilation tests**: New operators MUST include tests for all supported dialects
- **Code review**: At least one approving review required
- **Constitution compliance**: Changes must align with core principles

### Code Style
- Use Ruff for automatic formatting and linting
- Follow PEP 8 conventions
- Prefer simple, readable code over clever optimizations
- Avoid deep nesting (max 3-4 levels)
- Type hints MUST be complete for public APIs
- Docstrings MUST use reST format with `:param:`, `:returns:`, `:raises:`

## Release & Versioning

### Software Versioning
GIQL follows semantic versioning (MAJOR.MINOR.PATCH) for software releases:
- **MAJOR**: Breaking changes to GIQL syntax or query semantics
- **MINOR**: New operators, new supported dialects, backwards-compatible features
- **PATCH**: Bug fixes, transpilation improvements, documentation updates

**Examples:**
- Adding new operator (OVERLAPS): MINOR version bump
- Changing INTERSECTS semantics: MAJOR version bump
- Fixing bug in MERGE transpilation: PATCH version bump
- Adding PostgreSQL support: MINOR version bump

### Backwards Compatibility
- **Minor and patch versions**: MUST maintain backwards compatibility within major version
- **Query compatibility**: Existing GIQL queries MUST continue to work across MINOR/PATCH updates
- **Breaking changes**: MUST be communicated through deprecation warnings at least one minor version before removal
- **Deprecation cycle**: Deprecated syntax MUST remain functional (with warnings) until next major revision
- **Transpilation stability**: SQL output may change across versions as long as semantics are preserved

**Compatibility guarantees:**
- GIQL query syntax (what users write)
- Query semantics (what results are returned)
- Public Python API (`GIQLEngine`, `execute()`, `transpile()`, etc.)

**NOT guaranteed across versions:**
- Exact SQL output (transpilation may be optimized)
- Internal AST structure
- Performance characteristics (may improve or regress)

### Database Dialect Support
- New dialect support is a MINOR version change
- Dropping dialect support is a MAJOR version change
- Dialect-specific features MUST be clearly documented
- Core operators MUST work on all supported dialects

## Governance

### Amendment Process
This constitution supersedes all other practices. Amendments to this constitution require:
1. Documentation of proposed change with rationale
2. Approval through pull request review process
3. Version increment per constitution versioning rules (below)
4. Update to SYNC IMPACT REPORT at top of this file

### Constitution Versioning
The constitution itself follows semantic versioning (MAJOR.MINOR.PATCH):
- **MAJOR**: Backward incompatible governance/principle removals or redefinitions
- **MINOR**: New principle/section added or materially expanded guidance
- **PATCH**: Clarifications, wording, typo fixes, non-semantic refinements

### Compliance Review
All pull requests and code reviews MUST verify compliance with constitution principles:
- Constitution principles take precedence over convenience or convention
- Violations MUST be explicitly justified in pull request description
- Reviewers MUST check alignment with Expressive, Readable, Portable, and Canonical principles
- New operators MUST demonstrate alignment with all four core principles
- Use CLAUDE.md for agent-specific runtime development guidance

**Version**: 1.0.0 | **Ratified**: 2025-11-03 | **Last Amended**: 2025-11-03
