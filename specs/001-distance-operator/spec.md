# Feature Specification: DISTANCE UDF Operator

**Feature Branch**: `001-distance-operator`
**Created**: 2025-11-03
**Status**: Draft
**Input**: User description: "I'd like to add a DISTANCE UDF-based operator that mirrors bedtools' 'closest -d' command"

## Constitution Alignment

This feature aligns with GIQL's core principles:

- **Expressive**: Provides a concise way to find nearest genomic features and calculate distances, avoiding complex SQL window functions and self-joins
- **Readable**: `DISTANCE(a.position, b.position)` clearly expresses intent to calculate genomic distance between intervals
- **Portable**: Can be implemented as a UDF (User-Defined Function) across all supported SQL dialects (DuckDB, SQLite, PostgreSQL)
- **Canonical**: Mirrors bedtools' `closest -d` behavior, following established genomics tool conventions

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Basic Distance Calculation (Priority: P1)

A genomics researcher wants to find the closest feature in dataset B for each feature in dataset A and calculate the genomic distance between them.

**Why this priority**: Core use case that delivers immediate value. This is the most common "closest feature" query pattern.

**Independent Test**: Can be fully tested by loading two BED files, running a query with DISTANCE(), and verifying correct distances are calculated (0 for overlaps, positive integers for gaps).

**Acceptance Scenarios**:

1. **Given** two genomic features that overlap, **When** DISTANCE() is calculated, **Then** it returns 0
2. **Given** two genomic features on the same chromosome separated by a gap, **When** DISTANCE() is calculated, **Then** it returns the number of base pairs between them
3. **Given** two genomic features on different chromosomes, **When** DISTANCE() is calculated, **Then** it returns NULL or a special value indicating no distance is meaningful

---

### User Story 2 - Finding Nearest Features (Priority: P2)

A researcher wants to identify the single closest feature in dataset B for each feature in dataset A, replicating `bedtools closest -d` functionality.

**Why this priority**: Builds on P1 to provide the complete "closest feature" workflow. Common pattern in regulatory element analysis (e.g., finding nearest gene to each ChIP-seq peak).

**Independent Test**: Can be tested by creating a query that combines DISTANCE() with window functions (ROW_NUMBER) to select the minimum distance feature per query interval.

**Acceptance Scenarios**:

1. **Given** a feature in A with multiple features in B, **When** querying for the closest feature, **Then** the feature with minimum distance is returned
2. **Given** a feature in A with tied closest features (same distance), **When** querying for the closest feature, **Then** all tied features are returned (matching bedtools default behavior)
3. **Given** a feature in A with no features in B on the same chromosome, **When** querying for the closest feature, **Then** NULL or appropriate indicator is returned

---

### User Story 3 - Strand-Specific Distance (Priority: P3)

A researcher wants to calculate distances only between features on the same strand, mirroring bedtools' `-s` flag.

**Why this priority**: Important for strand-specific analyses but less common than basic distance calculations. Follows the same pattern as MERGE and CLUSTER operators.

**Independent Test**: Can be tested by passing `stranded=true` parameter to DISTANCE() and verifying it returns NULL for intervals on different strands (similar to different chromosomes).

**Acceptance Scenarios**:

1. **Given** two features on the same strand, **When** DISTANCE(..., stranded=true) is calculated, **Then** it returns the distance value
2. **Given** two features on different strands, **When** DISTANCE(..., stranded=true) is calculated, **Then** it returns NULL
3. **Given** two features with unspecified strand (e.g., '.'), **When** DISTANCE(..., stranded=true) is calculated, **Then** behavior matches bedtools (treat as separate strands)

---

### User Story 4 - Signed Distance (Priority: P4)

A researcher wants to know not just the distance but also the directionality (upstream vs downstream) between features, mirroring bedtools' `-D` flag.

**Why this priority**: Advanced use case for directional analyses. Less commonly needed and adds complexity.

**Independent Test**: Can be tested by passing `signed=true` parameter to DISTANCE() function and verifying negative values for upstream features and positive for downstream.

**Acceptance Scenarios**:

1. **Given** feature B is upstream of feature A, **When** DISTANCE(..., signed=true) is calculated, **Then** it returns a negative distance
2. **Given** feature B is downstream of feature A, **When** DISTANCE(..., signed=true) is calculated, **Then** it returns a positive distance
3. **Given** features overlap, **When** DISTANCE(..., signed=true) is calculated, **Then** it returns 0

---

### Edge Cases

- What happens when features are on different chromosomes? (Return NULL)
- What happens when features are on different strands with `stranded=true`? (Return NULL, analogous to different chromosomes)
- What happens when features are on different strands with `stranded=false`? (Calculate distance normally, ignoring strand)
- What happens with unspecified strand values (e.g., '.')? (Treat as distinct from '+' and '-'; return NULL if comparing '.' to '+' with `stranded=true`)
- What happens when one feature is zero-width (start == end)? (Calculate distance from point)
- What happens when intervals overlap partially? (Return 0, matching bedtools behavior)
- What happens when intervals are adjacent (bookended, e.g., end_a == start_b)? (Return 0, matching bedtools half-open interval convention)
- How are chromosome name comparisons handled (chr1 vs 1)? (Rely on user to ensure consistent naming)
- What happens with NULL chromosome, start, end, or strand values? (Return NULL)
- What happens with `signed=true` for overlapping intervals? (Return 0, no directionality for overlaps)

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST provide a DISTANCE() UDF that calculates genomic distance between two intervals
- **FR-002**: DISTANCE() MUST return 0 for overlapping intervals (matching bedtools semantics)
- **FR-003**: DISTANCE() MUST return the minimum distance (in base pairs) between two non-overlapping intervals on the same chromosome
- **FR-004**: DISTANCE() MUST return NULL when intervals are on different chromosomes
- **FR-005**: DISTANCE() MUST be implemented as a UDF for DuckDB, SQLite, and PostgreSQL
- **FR-006**: System MUST handle zero-width intervals (point features) correctly
- **FR-007**: DISTANCE() MUST accept flexible genomic interval inputs (matching INTERSECTS pattern):
  - **Position column references**: `DISTANCE(a.position, b.position)`
  - **Literal genomic ranges**: `DISTANCE(a.position, 'chr2:3000-4000')` OR `DISTANCE('chr1:1000-2000', b.position)`
  - **Two literal ranges**: `DISTANCE('chr1:1000-2000', 'chr2:3000-4000')`
- **FR-008**: DISTANCE() MUST accept optional named parameter `stranded` (boolean, default=false) for strand-specific distance calculation
- **FR-009**: When `stranded=true`, DISTANCE() MUST return NULL if intervals are on different strands (analogous to different chromosomes)
- **FR-010**: DISTANCE() MUST accept optional named parameter `signed` (boolean, default=false) for directional distance calculation
- **FR-011**: When `signed=true`, DISTANCE() MUST return negative values for upstream intervals and positive for downstream intervals
- **FR-012**: Both `stranded` and `signed` parameters MAY be used together
- **FR-013**: When both `signed=true` AND `stranded=true`, directionality MUST be relative to reference genome coordinates (lower coordinates = upstream) regardless of feature strand orientation
- **FR-014**: SQL generator MUST resolve position column references to physical chromosome/start/end/strand columns at transpilation time
- **FR-015**: Documentation MUST include examples showing how to combine DISTANCE() with window functions to replicate `bedtools closest -d`
- **FR-016**: Documentation MUST include examples showing how to use `stranded=true` parameter to replicate `bedtools closest -d -s`
- **FR-017**: Documentation MUST include examples showing how to use `signed=true` parameter to replicate `bedtools closest -D ref`
- **FR-018**: Documentation MUST include a note that data does not need to be presorted for correctness (unlike bedtools), though sorting may impact performance

### Implementation Considerations

- **IC-001**: DISTANCE() should NOT require pre-sorting of data (unlike bedtools closest)
- **IC-002**: For "find closest feature" queries, users will combine DISTANCE() with SQL window functions (ROW_NUMBER, RANK) or aggregations (MIN)
- **IC-003**: DISTANCE() follows the same input pattern as INTERSECTS predicate - accepting position columns and/or literal genomic ranges
- **IC-004**: The SQL generator resolves position column references at transpilation time, passing resolved chromosome/start/end/strand values to the underlying UDF
- **IC-005**: The actual UDF implementation receives explicit column values (chrom_a, start_a, end_a, strand_a, chrom_b, start_b, end_b, strand_b, stranded, signed) after SQL generator resolution
- **IC-006**: When `stranded=true`, the UDF needs strand information; literal ranges without strand specifications default to unstranded (strand='.')

### Key Entities

- **GenomicInterval**: Represents a genomic interval with chromosome, start, end, and strand
- **Distance**: Scalar integer representing base pairs between two intervals
  - 0 for overlapping intervals
  - Positive integer for gap between non-overlapping intervals (unsigned mode)
  - Positive/negative integer for directional distance (signed mode: negative=upstream, positive=downstream)
  - NULL for incompatible intervals (different chromosomes, or different strands when `stranded=true`)

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Users can calculate distance between any two genomic intervals using DISTANCE() UDF
- **SC-002**: DISTANCE() returns identical results to bedtools closest -d for overlapping intervals (0) and non-overlapping intervals (positive integer)
- **SC-003**: Users can replicate bedtools closest -d functionality by combining DISTANCE() with window functions in a single SQL query
- **SC-004**: DISTANCE() UDF works correctly on DuckDB, SQLite, and PostgreSQL backends
- **SC-005**: Documentation includes complete working examples showing DISTANCE() usage patterns

### Test Coverage

- Unit tests for DISTANCE() UDF with various interval configurations (overlapping, adjacent, separated, different chromosomes)
- Integration tests showing complete "find closest feature" query pattern
- Property-based tests using Hypothesis for edge cases (zero-width, boundary conditions)
- Cross-dialect tests ensuring consistent behavior across all supported databases

## API Design

### GIQL Function Signature

```sql
-- User-facing syntax (what users write in GIQL)
DISTANCE(interval_a, interval_b, stranded=false, signed=false) -> INTEGER | NULL

-- Interval arguments can be:
-- - Position column references: a.position, b.position
-- - Literal genomic ranges: 'chr1:1000-2000' or 'chr1:1000-2000:+'

-- Usage examples:
-- Basic distance (ignores strand)
SELECT
    a.*,
    b.*,
    DISTANCE(a.position, b.position) AS distance
FROM features_a a
CROSS JOIN features_b b
WHERE a.chromosome = b.chromosome  -- Pre-filter for performance

-- Strand-specific distance (NULL if different strands)
SELECT
    a.*,
    b.*,
    DISTANCE(a.position, b.position, stranded=true) AS distance
FROM features_a a
CROSS JOIN features_b b
WHERE a.chromosome = b.chromosome

-- Signed distance (directional: negative=upstream, positive=downstream)
SELECT
    a.*,
    b.*,
    DISTANCE(a.position, b.position, signed=true) AS signed_distance
FROM features_a a
CROSS JOIN features_b b
WHERE a.chromosome = b.chromosome

-- Both strand-specific AND signed
SELECT
    a.*,
    b.*,
    DISTANCE(a.position, b.position, stranded=true, signed=true) AS signed_distance
FROM features_a a
CROSS JOIN features_b b
WHERE a.chromosome = b.chromosome

-- Mixed: position column + literal range
SELECT
    a.*,
    DISTANCE(a.position, 'chr1:5000-10000') AS distance_to_region
FROM features_a a
```

### Finding Closest Features (bedtools closest -d)

```sql
-- Example: Find closest feature in B for each feature in A
WITH distances AS (
    SELECT
        a.*,
        b.chromosome AS b_chromosome,
        b.start_pos AS b_start_pos,
        b.end_pos AS b_end_pos,
        DISTANCE(a.position, b.position) AS distance,
        ROW_NUMBER() OVER (
            PARTITION BY a.chromosome, a.start_pos, a.end_pos
            ORDER BY DISTANCE(a.position, b.position)
        ) AS rank
    FROM features_a a
    CROSS JOIN features_b b
    WHERE a.chromosome = b.chromosome
)
SELECT *
FROM distances
WHERE rank = 1
ORDER BY chromosome, start_pos
```

## Distance Calculation Algorithm

The DISTANCE function follows bedtools semantics:

1. **Different chromosomes**: Return NULL
2. **Overlapping intervals** (start_a < end_b AND end_a > start_b): Return 0
3. **Non-overlapping intervals**:
   - If a is completely before b: distance = start_b - end_a
   - If b is completely before a: distance = start_a - end_b
   - Return the minimum distance

### Pseudocode

```python
def distance(chrom_a, start_a, end_a, strand_a,
             chrom_b, start_b, end_b, strand_b,
             stranded=False, signed=False):
    # Different chromosomes
    if chrom_a != chrom_b:
        return NULL

    # Different strands (when stranded mode is enabled)
    if stranded and strand_a != strand_b:
        return NULL

    # Check for overlap
    if start_a < end_b and end_a > start_b:
        return 0

    # Calculate gap distance
    if end_a <= start_b:
        # a is before b
        distance = start_b - end_a
        if signed:
            return distance  # positive = downstream
        return distance
    else:
        # b is before a
        distance = start_a - end_b
        if signed:
            return -distance  # negative = upstream
        return distance
```

## Bedtools Compatibility Notes

This implementation replicates the core behavior of `bedtools closest` with these mappings:

**GIQL to Bedtools Flag Mapping**:
- `DISTANCE(a.position, b.position)` → `bedtools closest -d`
- `DISTANCE(a.position, b.position, stranded=true)` → `bedtools closest -d -s`
- `DISTANCE(a.position, b.position, signed=true)` → `bedtools closest -D ref`
- `DISTANCE(a.position, b.position, stranded=true, signed=true)` → `bedtools closest -D ref -s`

**Implemented**:
- Basic distance calculation (equivalent to `-d` flag)
- Strand-specific distance (equivalent to `-s` flag via `stranded=true`)
- Signed distance (equivalent to `-D ref` flag via `signed=true`)
- Zero distance for overlaps
- NULL for different chromosomes (and different strands when `stranded=true`)

**User must implement via SQL**:
- Finding closest feature (use window functions: `ROW_NUMBER() OVER ... ORDER BY DISTANCE(...)`)
- Finding k-closest features (use window functions with `WHERE rank <= k`)
- Tie-breaking strategies (use RANK instead of ROW_NUMBER for all ties, equivalent to `-t all`)
- Ignoring overlaps (filter `WHERE distance > 0`, equivalent to `-io` flag)
- Multiple B file handling (use UNION ALL in SQL)

**Differences from bedtools**:
- **No pre-sorting required**: GIQL DISTANCE() does not require presorted input data (bedtools closest does). This simplifies usage but may have performance implications for very large datasets.
- **SQL-based workflows**: GIQL leverages SQL for advanced features (k-closest, tie-breaking) rather than command-line flags

## Implementation Notes

### Database-Specific UDF Registration

**DuckDB**:
```python
def register_distance_udf(conn):
    # UDF receives: chrom_a, start_a, end_a, strand_a,
    #               chrom_b, start_b, end_b, strand_b,
    #               stranded (bool), signed (bool)
    conn.create_function(
        "DISTANCE",
        distance_func,
        [str, int, int, str, str, int, int, str, bool, bool],
        int
    )
```

**SQLite**:
```python
def register_distance_udf(conn):
    # SQLite doesn't have native bool type, use integers (0/1)
    conn.create_function("DISTANCE", 10, distance_func)
```

**PostgreSQL**:
```sql
CREATE FUNCTION DISTANCE(
    chrom_a TEXT, start_a INTEGER, end_a INTEGER, strand_a TEXT,
    chrom_b TEXT, start_b INTEGER, end_b INTEGER, strand_b TEXT,
    stranded BOOLEAN DEFAULT FALSE,
    signed BOOLEAN DEFAULT FALSE
) RETURNS INTEGER AS $$
BEGIN
    -- Different chromosomes
    IF chrom_a != chrom_b THEN
        RETURN NULL;
    END IF;

    -- Different strands (when stranded mode enabled)
    IF stranded AND strand_a != strand_b THEN
        RETURN NULL;
    END IF;

    -- Check for overlap
    IF start_a < end_b AND end_a > start_b THEN
        RETURN 0;
    END IF;

    -- Calculate distance
    IF end_a <= start_b THEN
        IF signed THEN
            RETURN start_b - end_a;  -- positive = downstream
        ELSE
            RETURN start_b - end_a;
        END IF;
    ELSE
        IF signed THEN
            RETURN -(start_a - end_b);  -- negative = upstream
        ELSE
            RETURN start_a - end_b;
        END IF;
    END IF;
END;
$$ LANGUAGE plpgsql IMMUTABLE;
```

## Documentation Requirements

- Add DISTANCE() to README.md with examples showing:
  - Basic distance calculation: `DISTANCE(a.position, b.position)`
  - Strand-specific distance: `DISTANCE(a.position, b.position, stranded=true)`
  - Signed distance: `DISTANCE(a.position, b.position, signed=true)`
  - Combined stranded + signed
- Update docs/examples.rst with "Finding Closest Features" section
- Add bedtools closest equivalence table showing GIQL → bedtools flag mappings
- Include performance notes:
  - Data does NOT need to be presorted (unlike bedtools)
  - CROSS JOIN can be expensive on large datasets
  - Pre-filter by chromosome for better performance: `WHERE a.chromosome = b.chromosome`
- Add examples replicating common bedtools closest patterns:
  - `bedtools closest -d` → basic DISTANCE query
  - `bedtools closest -d -s` → DISTANCE with `stranded=true`
  - `bedtools closest -D ref` → DISTANCE with `signed=true`
  - `bedtools closest -d -io` → DISTANCE query with `WHERE distance > 0`

## Migration Path

1. Implement DISTANCE() UDF for all supported dialects
2. Add comprehensive unit tests
3. Document usage patterns in README
4. Add integration test showing full "closest feature" workflow
5. Update demo.ipynb with DISTANCE() examples
6. (Future) Consider adding syntactic sugar: `SELECT CLOSEST(a.position FROM b.position) AS nearest_feature`
