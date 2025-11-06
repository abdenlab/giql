# Quickstart: DISTANCE Operator

**Date**: 2025-11-06
**Audience**: GIQL Users

## Overview

The `DISTANCE()` operator calculates genomic distances between intervals, making it easy to find the nearest features, compute proximity metrics, and perform distance-based analyses.

## Basic Usage

### Calculate Distance Between Two Tables

```sql
SELECT
    a.name AS feature_a,
    b.name AS feature_b,
    DISTANCE(a.position, b.position) AS distance
FROM peaks a
CROSS JOIN genes b
WHERE a.chromosome = b.chromosome
LIMIT 10;
```

**Results**:
```
feature_a | feature_b | distance
----------+-----------+---------
peak_1    | gene_A    | 0        (overlapping)
peak_2    | gene_B    | 1500     (1.5kb apart)
peak_3    | gene_C    | NULL     (different chromosomes)
```

### Find Closest Feature

Find the single closest gene for each peak:

```sql
WITH ranked_distances AS (
    SELECT
        a.name AS peak,
        b.name AS nearest_gene,
        DISTANCE(a.position, b.position) AS distance,
        ROW_NUMBER() OVER (
            PARTITION BY a.chromosome, a.start_pos, a.end_pos
            ORDER BY DISTANCE(a.position, b.position) NULLS LAST
        ) AS rank
    FROM peaks a
    CROSS JOIN genes b
    WHERE a.chromosome = b.chromosome
)
SELECT peak, nearest_gene, distance
FROM ranked_distances
WHERE rank = 1
ORDER BY peak;
```

**bedtools equivalent**: `bedtools closest -a peaks.bed -b genes.bed -d`

## Advanced Usage

### Strand-Specific Distance

Calculate distances only between features on the same strand:

```sql
SELECT
    a.name,
    b.name,
    DISTANCE(a.position, b.position, stranded=true) AS distance
FROM forward_peaks a
CROSS JOIN forward_genes b
WHERE a.chromosome = b.chromosome;
```

**Notes**:
- Returns `NULL` if features are on different strands
- Useful for strand-specific regulatory analyses

**bedtools equivalent**: `bedtools closest -a peaks.bed -b genes.bed -d -s`

### Signed Distance (Directional)

Get directional distance showing upstream vs downstream:

```sql
SELECT
    a.name AS peak,
    b.name AS gene,
    DISTANCE(a.position, b.position, signed=true) AS signed_distance
FROM peaks a
CROSS JOIN genes b
WHERE a.chromosome = b.chromosome;
```

**Interpretation**:
- **Negative**: Gene is upstream of peak (lower coordinates)
- **Positive**: Gene is downstream of peak (higher coordinates)
- **0**: Overlapping

**bedtools equivalent**: `bedtools closest -a peaks.bed -b genes.bed -D ref`

### Combined: Strand-Specific AND Signed

```sql
SELECT
    a.name,
    b.name,
    DISTANCE(a.position, b.position, stranded=true, signed=true) AS signed_distance
FROM peaks a
CROSS JOIN genes b
WHERE a.chromosome = b.chromosome
  AND a.strand = b.strand;  -- Pre-filter for performance
```

**bedtools equivalent**: `bedtools closest -a peaks.bed -b genes.bed -D ref -s`

## Common Patterns

### Find K Nearest Neighbors

Find the 3 closest genes for each peak:

```sql
WITH ranked_distances AS (
    SELECT
        a.name AS peak,
        b.name AS gene,
        DISTANCE(a.position, b.position) AS distance,
        ROW_NUMBER() OVER (
            PARTITION BY a.chromosome, a.start_pos, a.end_pos
            ORDER BY DISTANCE(a.position, b.position) NULLS LAST
        ) AS rank
    FROM peaks a
    CROSS JOIN genes b
    WHERE a.chromosome = b.chromosome
)
SELECT peak, gene, distance, rank
FROM ranked_distances
WHERE rank <= 3
ORDER BY peak, rank;
```

### Exclude Overlaps (Closest Non-Overlapping Feature)

Find nearest gene that doesn't overlap the peak:

```sql
WITH ranked_distances AS (
    SELECT
        a.name AS peak,
        b.name AS gene,
        DISTANCE(a.position, b.position) AS distance,
        ROW_NUMBER() OVER (
            PARTITION BY a.chromosome, a.start_pos, a.end_pos
            ORDER BY DISTANCE(a.position, b.position) NULLS LAST
        ) AS rank
    FROM peaks a
    CROSS JOIN genes b
    WHERE a.chromosome = b.chromosome
      AND DISTANCE(a.position, b.position) > 0  -- Exclude overlaps
)
SELECT peak, gene, distance
FROM ranked_distances
WHERE rank = 1;
```

**bedtools equivalent**: `bedtools closest -a peaks.bed -b genes.bed -d -io`

### Distance to Specific Region

Calculate distance from each peak to a specific genomic region:

```sql
SELECT
    name,
    DISTANCE(position, 'chr1:1000000-1001000') AS dist_to_region
FROM peaks
WHERE chromosome = 'chr1';
```

### Aggregate Distance Statistics

Find average distance to nearest gene per chromosome:

```sql
WITH nearest_distances AS (
    SELECT
        a.chromosome,
        a.name AS peak,
        MIN(DISTANCE(a.position, b.position)) AS nearest_distance
    FROM peaks a
    CROSS JOIN genes b
    WHERE a.chromosome = b.chromosome
    GROUP BY a.chromosome, a.name, a.start_pos, a.end_pos
)
SELECT
    chromosome,
    AVG(nearest_distance) AS avg_distance,
    MIN(nearest_distance) AS min_distance,
    MAX(nearest_distance) AS max_distance
FROM nearest_distances
GROUP BY chromosome
ORDER BY chromosome;
```

## Performance Tips

### 1. Pre-Filter by Chromosome

**Always** include `WHERE a.chromosome = b.chromosome` to avoid unnecessary comparisons:

```sql
-- ✅ GOOD
SELECT DISTANCE(a.position, b.position)
FROM features_a a CROSS JOIN features_b b
WHERE a.chromosome = b.chromosome;

-- ❌ BAD (computes distances across all chromosomes)
SELECT DISTANCE(a.position, b.position)
FROM features_a a CROSS JOIN features_b b;
```

### 2. Use Indexes

Create indexes on commonly filtered columns:

```sql
CREATE INDEX idx_features_chrom_start ON features(chromosome, start_pos);
```

### 3. Filter Before CROSS JOIN

When possible, filter tables before joining:

```sql
-- Filter to specific chromosomes of interest
WITH filtered_a AS (
    SELECT * FROM features_a WHERE chromosome IN ('chr1', 'chr2')
),
filtered_b AS (
    SELECT * FROM features_b WHERE chromosome IN ('chr1', 'chr2')
)
SELECT DISTANCE(a.position, b.position)
FROM filtered_a a CROSS JOIN filtered_b b
WHERE a.chromosome = b.chromosome;
```

### 4. Materialize Results for Reuse

If computing distances multiple times, materialize results:

```sql
CREATE TABLE distances AS
SELECT
    a.id AS id_a,
    b.id AS id_b,
    DISTANCE(a.position, b.position) AS distance
FROM features_a a CROSS JOIN features_b b
WHERE a.chromosome = b.chromosome;

-- Then query the materialized table
SELECT * FROM distances WHERE distance < 10000 ORDER BY distance;
```

## Edge Cases

### NULL Distances

DISTANCE returns `NULL` when:
- Features are on different chromosomes
- Features are on different strands (when `stranded=true`)
- Any input column is NULL

**Handling NULLs in ORDER BY**:
```sql
-- Explicitly sort NULLs last
ORDER BY DISTANCE(a.position, b.position) NULLS LAST

-- Or filter them out
WHERE DISTANCE(a.position, b.position) IS NOT NULL
```

### Adjacent Features (Bookended)

Adjacent features where `end_a = start_b` return distance `0` (following bedtools convention):

```sql
-- Feature A: chr1:100-200
-- Feature B: chr1:200-300
-- DISTANCE: 0 (not 1)
```

### Zero-Width Features (Points)

Distance from/to point features works as expected:

```sql
-- Feature A: chr1:1000-1000 (point)
-- Feature B: chr1:1500-2000 (interval)
-- DISTANCE: 500
```

## Comparison with bedtools closest

| Feature | GIQL DISTANCE | bedtools closest |
|---------|---------------|------------------|
| Find closest feature | ✅ With window functions | ✅ Built-in |
| Basic distance (`-d`) | `DISTANCE(a.pos, b.pos)` | `-d` flag |
| Strand-specific (`-s`) | `stranded=true` | `-s` flag |
| Signed distance (`-D ref`) | `signed=true` | `-D ref` flag |
| Ignore overlaps (`-io`) | `WHERE distance > 0` | `-io` flag |
| K-nearest (`-k 3`) | `WHERE rank <= 3` | `-k 3` flag |
| Pre-sorting required | ❌ No | ✅ Yes |
| SQL composability | ✅ Full SQL power | ❌ Limited |

## Common Errors

### Error: Column not found

```
Error: Column 'a.pos' is not a genomic position column
```

**Solution**: Register the column in your schema or use default column names (`chromosome`, `start_pos`, `end_pos`)

### Error: Cannot parse genomic range

```
Error: Could not parse genomic range: 'invalid_range'
```

**Solution**: Use valid range format: `'chr1:1000-2000'` or `'chr1:1000-2000:+'`

### Performance: Query too slow

**Symptoms**: CROSS JOIN of large tables takes too long

**Solutions**:
1. Add `WHERE a.chromosome = b.chromosome` pre-filter
2. Create indexes on chromosome and start_pos columns
3. Filter to regions of interest before CROSS JOIN
4. Consider using bedtools for very large datasets (optimized for genomics scale)

## Python API Usage

```python
from giql import GIQLEngine

# Create engine
engine = GIQLEngine(target_dialect="duckdb")

# Register tables
engine.register_table_schema("peaks", genomic_column="position")
engine.register_table_schema("genes", genomic_column="position")

# Execute query
query = """
SELECT
    a.name AS peak,
    b.name AS nearest_gene,
    DISTANCE(a.position, b.position) AS distance
FROM peaks a
CROSS JOIN genes b
WHERE a.chromosome = b.chromosome
ORDER BY distance
LIMIT 10
"""

result = engine.execute(query)
for row in result:
    print(f"{row['peak']} -> {row['nearest_gene']}: {row['distance']}bp")
```

## Next Steps

- See `data-model.md` for implementation details
- See `spec.md` for complete operator specification
- See `plan.md` for development roadmap

## Quick Reference

### Syntax

```sql
DISTANCE(interval_a, interval_b [, stranded=false] [, signed=false]) -> INTEGER | NULL
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `interval_a` | Position column or range | Required | First genomic interval |
| `interval_b` | Position column or range | Required | Second genomic interval |
| `stranded` | Boolean | `false` | If true, return NULL for different strands |
| `signed` | Boolean | `false` | If true, return negative for upstream |

### Return Value

| Value | Meaning |
|-------|---------|
| `0` | Overlapping intervals |
| Positive integer | Gap distance in base pairs (unsigned) or downstream distance (signed) |
| Negative integer | Upstream distance (signed mode only) |
| `NULL` | Different chromosomes or different strands (when stranded=true) |
