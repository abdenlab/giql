# Feature Specification: NEAREST Operator

**Feature Branch**: `002-nearest-operator`
**Created**: 2025-11-07
**Status**: Draft
**Input**: User description: "Implement a NEAREST operator that can be used to find the nearest k neighbors. Help me evaluate the best way to implement features like "select k nearest neighbors where distance < d"."

## Constitution Alignment

This feature aligns with GIQL's core principles:

- **Expressive**: Provides a concise way to find k-nearest genomic features, eliminating complex window functions and self-joins for this common pattern
- **Readable**: `NEAREST(genes, k=3)` clearly expresses intent to find the 3 closest features, making queries more intuitive
- **Portable**: Can be implemented across all supported SQL dialects using standard SQL constructs (window functions, subqueries)
- **Canonical**: Extends bedtools' `closest` behavior with SQL-native k-nearest neighbor semantics, a natural evolution of the DISTANCE operator

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Find K-Nearest Features (Priority: P1)

A genomics researcher wants to identify the k-nearest features in dataset B for each feature in dataset A (e.g., find the 3 closest genes to each ChIP-seq peak).

**Why this priority**: Core use case that delivers immediate value. Finding multiple nearest features is more common in practice than finding just the single closest feature (e.g., regulatory analysis often requires examining multiple nearby genes).

**Independent Test**: Can be fully tested by loading two BED files, running a query with NEAREST(k=3), and verifying the correct 3 closest features are returned for each query interval, ordered by distance.

**Acceptance Scenarios**:

1. **Given** a feature in A with 5+ features in B on the same chromosome, **When** NEAREST(k=3) is executed, **Then** exactly 3 features with the smallest distances are returned
2. **Given** a feature in A with only 2 features in B on the same chromosome, **When** NEAREST(k=3) is executed, **Then** only 2 features are returned (all available features)
3. **Given** multiple features have identical distances (ties), **When** NEAREST(k=3) is executed, **Then** all tied features are included (may return more than k features)

---

### User Story 2 - Distance-Constrained Nearest Neighbors (Priority: P2)

A researcher wants to find the k-nearest features within a maximum distance threshold (e.g., find up to 5 closest genes, but only those within 100kb).

**Why this priority**: Common constraint in regulatory genomics where spatial proximity has biological meaning only within certain distances. Builds on P1 by adding filtering capability.

**Independent Test**: Can be tested by executing NEAREST with both k and max_distance parameters and verifying that results satisfy both constraints (≤k features, all within max_distance).

**Acceptance Scenarios**:

1. **Given** a feature with 10 features in B, where 3 are within 50kb and 7 are beyond 100kb, **When** NEAREST(k=5, max_distance=100000) is executed, **Then** only the 3 features within 100kb are returned
2. **Given** a feature with 8 features in B all within 50kb, **When** NEAREST(k=5, max_distance=100000) is executed, **Then** the 5 closest features are returned
3. **Given** a feature with no features in B within 100kb, **When** NEAREST(k=5, max_distance=100000) is executed, **Then** no results are returned for that query feature

---

### User Story 3 - Strand-Specific Nearest Neighbors (Priority: P3)

A researcher wants to find k-nearest features considering only features on the same strand, similar to DISTANCE(stranded=true).

**Why this priority**: Important for strand-specific analyses but less common than basic k-nearest queries. Follows the same pattern as DISTANCE and other genomic operators.

**Independent Test**: Can be tested by executing NEAREST with stranded=true and verifying only same-strand features are considered.

**Acceptance Scenarios**:

1. **Given** a feature on '+' strand with 5 features in B (3 on '+', 2 on '-'), **When** NEAREST(k=3, stranded=true) is executed, **Then** only the 3 features on '+' strand are considered
2. **Given** a feature on '+' strand with only 1 feature in B on '+' strand, **When** NEAREST(k=3, stranded=true) is executed, **Then** only that 1 feature is returned
3. **Given** a feature with unspecified strand ('.'), **When** NEAREST(k=3, stranded=true) is executed, **Then** only features with strand='.' are considered

---

### User Story 4 - Directional (Upstream/Downstream) Nearest Neighbors (Priority: P4)

A researcher wants to find k-nearest features in a specific direction (upstream or downstream), using signed distances to filter results.

**Why this priority**: Advanced use case for directional analyses (e.g., finding nearest upstream promoters). Adds complexity and is less frequently needed.

**Independent Test**: Can be tested by executing NEAREST with signed=true and filtering on distance sign to verify directional constraints are applied.

**Acceptance Scenarios**:

1. **Given** a feature with 3 upstream and 4 downstream features in B, **When** NEAREST(k=3, signed=true) with WHERE distance < 0 is executed, **Then** the 3 nearest upstream features are returned
2. **Given** a feature with 2 upstream and 5 downstream features in B, **When** NEAREST(k=3, signed=true) with WHERE distance > 0 is executed, **Then** the 3 nearest downstream features are returned
3. **Given** a feature with overlapping features, **When** NEAREST(k=3, signed=true) is executed, **Then** overlaps (distance=0) are included in the k-nearest results

---

### User Story 5 - Literal Reference Point Queries (Priority: P2)

A researcher wants to find k-nearest features to a specific genomic location without needing a source table (e.g., find the 5 closest genes to chr1:1000000-1001000).

**Why this priority**: Common for ad-hoc queries and exploratory analysis. Complements the per-row pattern with standalone queries.

**Independent Test**: Can be tested by executing NEAREST with a literal reference parameter and verifying correct nearest features are returned.

**Acceptance Scenarios**:

1. **Given** a literal genomic position 'chr1:1000-2000', **When** NEAREST(genes, reference='chr1:1000-2000', k=5) is executed, **Then** the 5 closest genes are returned
2. **Given** a literal position with strand 'chr1:1000-2000:+', **When** NEAREST(genes, reference='chr1:1000-2000:+', k=3, stranded=true) is executed, **Then** only genes on '+' strand are considered
3. **Given** multiple literal positions via CTE, **When** NEAREST is used with CROSS JOIN LATERAL, **Then** k-nearest features are returned for each query position

---

### Edge Cases

- What happens when k is larger than the number of available features? (Return all available features, no padding)
- What happens when there are no features in B on the same chromosome? (Return no results for that query feature)
- What happens when multiple features have identical distances (ties)? (Include all tied features, may return >k features, matching bedtools `-t all` behavior)
- What happens when k=0? (Return no results, valid edge case)
- What happens with overlapping features (distance=0)? (Include in nearest results, distance=0 is valid)
- How are strand-specific queries handled when features have unspecified strand ('.')? (Treat '.' as distinct from '+' and '-', same as DISTANCE behavior)
- What happens when max_distance constraint eliminates all candidates? (Return empty result for that query feature)
- Can NEAREST be combined with other predicates (e.g., WHERE gene_type='protein_coding' on result)? (Yes, standard SQL composition applies)
- What happens when reference parameter is omitted in standalone mode? (Error - reference is required in standalone mode)
- What happens when reference parameter is omitted in LATERAL mode? (Use outer table's .position column by convention)
- What column is used for position when reference is omitted? (By convention, use `.position` column from outer table)

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST provide a NEAREST() table-valued function that returns k-nearest features from a target table
- **FR-002**: NEAREST() MUST accept a target table name as the first parameter
- **FR-003**: NEAREST() MUST support a `k` parameter specifying the maximum number of nearest features to return per query interval (default=1)
- **FR-004**: NEAREST() MUST support a `reference` parameter accepting:
  - Literal genomic ranges: `'chr1:1000-2000'` or `'chr1:1000-2000:+'`
  - Column references: `table.position` or `table.column_name`
  - When omitted in LATERAL context, defaults to outer table's `.position` column
- **FR-005**: NEAREST() MUST return results ordered by distance (ascending) for each query interval
- **FR-006**: When multiple features have identical distances, NEAREST() MUST return all tied features (may exceed k, matching bedtools `-t all` behavior)
- **FR-007**: NEAREST() MUST support an optional `max_distance` parameter to filter results by maximum distance threshold
- **FR-008**: When `max_distance` is specified, NEAREST() MUST exclude features beyond the distance threshold even if fewer than k features remain
- **FR-009**: NEAREST() MUST support an optional `stranded` parameter (boolean, default=false) for strand-specific neighbor finding
- **FR-010**: When `stranded=true`, NEAREST() MUST consider only features on the same strand as the query interval
- **FR-011**: NEAREST() MUST support an optional `signed` parameter (boolean, default=false) for directional distance calculation
- **FR-012**: When `signed=true`, NEAREST() MUST calculate signed distances (negative=upstream, positive=downstream)
- **FR-013**: NEAREST() MUST return all columns from both query and target tables in the result set
- **FR-014**: NEAREST() MUST include a computed `distance` column in the result set showing the distance to each neighbor
- **FR-015**: NEAREST() MUST handle cases where k exceeds available features by returning all available features
- **FR-016**: NEAREST() MUST handle cases where no features exist on the same chromosome by returning empty results for that query interval
- **FR-017**: NEAREST() MUST be implemented using SQL transpilation (not as a database UDF) for maximum portability
- **FR-018**: NEAREST() MUST transpile to valid SQL for DuckDB, SQLite, and PostgreSQL
- **FR-019**: Users MUST be able to combine NEAREST() with standard SQL clauses (WHERE, ORDER BY, GROUP BY) on the result set
- **FR-020**: NEAREST() MUST work in standalone mode (FROM NEAREST(...)) and correlated mode (CROSS JOIN LATERAL NEAREST(...))
- **FR-021**: In standalone mode with literal reference, NEAREST() MUST NOT require LATERAL JOIN
- **FR-022**: In correlated mode without explicit reference, NEAREST() MUST use the outer table's `.position` column by convention
- **FR-023**: NEAREST() MUST error if reference parameter is omitted in standalone mode (no outer table context)
- **FR-024**: Documentation MUST include examples showing common usage patterns (basic k-nearest, with distance constraints, stranded, directional, literal references)

### Implementation Considerations

- **IC-001**: NEAREST() should be implemented as SQL transpilation generating window functions (RANK) over DISTANCE() calculations
- **IC-002**: The transpilation should generate a CTE pattern: calculate distances → rank by distance → filter by k/max_distance
- **IC-003**: Consider performance implications of CROSS JOIN for large datasets (same as DISTANCE queries)
- **IC-004**: The generated SQL should leverage existing DISTANCE() UDF to avoid code duplication
- **IC-005**: Tie handling via RANK (not ROW_NUMBER) ensures all tied features are included (matching bedtools behavior)
- **IC-006**: The SQL generator must detect standalone vs correlated context to generate appropriate SQL
- **IC-007**: In standalone mode, generate simple query with WHERE/ORDER/LIMIT
- **IC-008**: In correlated mode, generate window function partitioned by query interval
- **IC-009**: When reference is omitted in LATERAL context, resolve to `outer_table.position` at transpilation time
- **IC-010**: The target table parameter must be resolved to actual table schema to access position columns

### Key Entities

- **GenomicInterval**: Represents a genomic interval with chromosome, start, end, and strand (same as DISTANCE)
- **NeighborSet**: Collection of k-nearest features for a single query interval, each with:
  - All query interval columns (in correlated mode)
  - All target interval columns
  - Computed distance value
  - Implicit rank/ordering by distance

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Users can find k-nearest genomic features using a single NEAREST() call without writing window functions
- **SC-002**: NEAREST(k=1) produces identical results to bedtools closest for single-nearest queries
- **SC-003**: NEAREST() queries execute successfully on all supported database backends
- **SC-004**: Users can combine NEAREST() with distance constraints in a single query
- **SC-005**: Users can query k-nearest features for a literal genomic position without creating temporary tables
- **SC-006**: Documentation includes complete working examples for all major use cases
- **SC-007**: 90% of users can express k-nearest neighbor queries without consulting SQL window function documentation

## API Design

### GIQL Syntax

#### Correlated Mode (Per-Row K-Nearest)

```sql
-- Basic k-nearest: find 3 closest genes for each peak
FROM peaks
CROSS JOIN LATERAL NEAREST(genes, k=3)
-- Implicitly uses peaks.position as reference

-- Explicit reference column
FROM peaks
CROSS JOIN LATERAL NEAREST(genes, reference=peaks.position, k=3)

-- With distance constraint: find up to 5 closest genes within 100kb
FROM peaks
CROSS JOIN LATERAL NEAREST(genes, k=5, max_distance=100000)

-- Strand-specific: find 3 closest same-strand features
FROM peaks
CROSS JOIN LATERAL NEAREST(genes, k=3, stranded=true)

-- Directional: find 3 closest upstream genes
FROM peaks
CROSS JOIN LATERAL NEAREST(genes, k=3, signed=true)
WHERE distance < 0  -- Negative = upstream

-- Combined constraints
FROM peaks
CROSS JOIN LATERAL NEAREST(
    genes,
    k=5,
    max_distance=50000,
    stranded=true,
    signed=true
)

-- Additional filtering on results
FROM peaks
CROSS JOIN LATERAL NEAREST(genes, k=5)
WHERE distance < 10000 AND gene_type = 'protein_coding'
```

#### Standalone Mode (Literal Reference Point)

```sql
-- Find 3 genes nearest to chr1:1000-2000
FROM NEAREST(genes, reference='chr1:1000-2000', k=3)

-- With distance constraint
FROM NEAREST(
    genes,
    reference='chr1:1000-2000',
    k=5,
    max_distance=100000
)

-- Strand-specific with strand in literal
FROM NEAREST(
    genes,
    reference='chr1:1000-2000:+',
    k=3,
    stranded=true
)

-- Multiple query points via CTE
WITH query_points AS (
    SELECT 'chr1:1000-2000' AS position UNION ALL
    SELECT 'chr2:5000-6000' UNION ALL
    SELECT 'chr3:10000-11000'
)
FROM query_points
CROSS JOIN LATERAL NEAREST(genes, reference=query_points.position, k=3)
```

### Expected Result Schema

```sql
-- Correlated mode: returns columns from both tables plus distance
| peak_chr | peak_start | peak_end | peak_name | gene_chr | gene_start | gene_end | gene_name | distance |
|----------|------------|----------|-----------|----------|------------|----------|-----------|----------|
| chr1     | 1000       | 2000     | peak1     | chr1     | 1500       | 2500     | gene1     | 0        |
| chr1     | 1000       | 2000     | peak1     | chr1     | 3000       | 4000     | gene2     | 1000     |
| chr1     | 1000       | 2000     | peak1     | chr1     | 5000       | 6000     | gene3     | 3000     |

-- Standalone mode: returns only target table columns plus distance
| gene_chr | gene_start | gene_end | gene_name | distance |
|----------|------------|----------|-----------|----------|
| chr1     | 1500       | 2500     | gene1     | 500      |
| chr1     | 3000       | 4000     | gene2     | 1000     |
| chr1     | 5000       | 6000     | gene3     | 3000     |
```

## Implementation Strategy

### Transpilation Pattern for Correlated Mode

Generate SQL using window functions and DISTANCE():

```sql
-- GIQL input:
-- FROM peaks CROSS JOIN LATERAL NEAREST(genes, k=3)

-- Transpiled SQL:
WITH ranked_distances AS (
    SELECT
        peaks.*,
        genes.*,
        DISTANCE(
            peaks.chromosome, peaks.start_pos, peaks.end_pos, peaks.strand,
            genes.chromosome, genes.start_pos, genes.end_pos, genes.strand,
            false,  -- stranded
            false   -- signed
        ) AS distance,
        RANK() OVER (
            PARTITION BY peaks.chromosome, peaks.start_pos, peaks.end_pos, peaks.strand
            ORDER BY DISTANCE(
                peaks.chromosome, peaks.start_pos, peaks.end_pos, peaks.strand,
                genes.chromosome, genes.start_pos, genes.end_pos, genes.strand,
                false, false
            ) ASC
        ) AS rank
    FROM peaks
    CROSS JOIN genes
    WHERE peaks.chromosome = genes.chromosome  -- Pre-filter for performance
)
SELECT
    peak_chr, peak_start, peak_end, peak_name,
    gene_chr, gene_start, gene_end, gene_name,
    distance
FROM ranked_distances
WHERE rank <= 3
ORDER BY peak_chr, peak_start, rank
```

### Transpilation Pattern for Standalone Mode

```sql
-- GIQL input:
-- FROM NEAREST(genes, reference='chr1:1000-2000', k=3)

-- Transpiled SQL:
SELECT
    genes.*,
    DISTANCE(
        'chr1', 1000, 2000, '.',
        genes.chromosome, genes.start_pos, genes.end_pos, genes.strand,
        false, false
    ) AS distance
FROM genes
WHERE genes.chromosome = 'chr1'
ORDER BY distance ASC
LIMIT 3
```

### With max_distance Constraint

```sql
-- GIQL input:
-- FROM peaks CROSS JOIN LATERAL NEAREST(genes, k=5, max_distance=100000)

-- Transpiled SQL (add WHERE clause on distance):
WITH ranked_distances AS (
    SELECT
        peaks.*,
        genes.*,
        DISTANCE(...) AS distance,
        RANK() OVER (PARTITION BY ... ORDER BY DISTANCE(...)) AS rank
    FROM peaks
    CROSS JOIN genes
    WHERE peaks.chromosome = genes.chromosome
      AND DISTANCE(...) <= 100000  -- Distance constraint
)
SELECT *
FROM ranked_distances
WHERE rank <= 5
ORDER BY peak_chr, peak_start, rank
```

## Comparison with Spatial Database Patterns

This design aligns with established spatial database conventions:

| System | K-Nearest Pattern | GIQL Equivalent |
|--------|------------------|-----------------|
| PostGIS | `CROSS JOIN LATERAL (... ORDER BY <-> LIMIT k)` | `CROSS JOIN LATERAL NEAREST(table, k=...)` |
| SQL Server | `CROSS APPLY (... TOP k ORDER BY STDistance)` | Same as above |
| Oracle | Nested SDO_NN queries | Less common pattern |

**Key alignment**: All major spatial databases use LATERAL/APPLY for per-row k-NN. GIQL's NEAREST operator provides syntactic sugar over this established pattern.

## Documentation Requirements

- Add NEAREST() to README.md with examples showing:
  - Correlated mode: `CROSS JOIN LATERAL NEAREST(genes, k=3)`
  - Standalone mode: `FROM NEAREST(genes, reference='chr1:1000-2000', k=3)`
  - With distance constraint: `NEAREST(genes, k=5, max_distance=100000)`
  - Strand-specific: `NEAREST(genes, k=3, stranded=true)`
  - Directional: `NEAREST(genes, k=3, signed=true)` with WHERE distance < 0
- Update docs/examples.rst with "Finding K-Nearest Features" section
- Include bedtools equivalence examples (bedtools closest → NEAREST with k=1)
- Add performance notes:
  - Pre-filter by chromosome happens automatically in transpiled SQL
  - Consider using max_distance to reduce candidate set
  - CROSS JOIN can be expensive on very large datasets
- Document tie-handling behavior (RANK vs ROW_NUMBER)
- Show how to combine NEAREST with additional filters (gene type, score thresholds)
- Explain standalone vs correlated modes with clear examples
- Document the `.position` column convention for implicit reference

## Assumptions

- Target tables have a `.position` column or explicit column references are provided
- Users understand LATERAL JOIN semantics (or will learn from documentation)
- Database query optimizers can efficiently handle window functions over large result sets
- DISTANCE() UDF is already implemented (dependency on feature 001-distance-operator)
- Standard SQL window functions (RANK, PARTITION BY) are available in all target dialects
- LATERAL JOIN (PostgreSQL) or equivalent CROSS APPLY (SQL Server) is available in target dialects
- The most common use case is k=1-10, not k=1000+ (affects performance considerations)
- Tie handling should match bedtools `-t all` behavior (include all ties) rather than arbitrary tie-breaking
