# NEAREST Operator Quick Start Guide

## What is NEAREST?

The NEAREST operator finds the k-nearest genomic features to query intervals, eliminating the need to write complex window functions for k-nearest neighbor queries.

## Basic Usage

### Find 3 Nearest Genes to Each Peak

```sql
SELECT *
FROM peaks
CROSS JOIN LATERAL NEAREST(genes, k=3)
```

Result: For each peak, returns the 3 closest genes with computed distance.

### Find 5 Nearest Features to a Specific Location

```sql
SELECT *
FROM NEAREST(genes, reference='chr1:1000000-1001000', k=5)
```

Result: Returns the 5 genes closest to chr1:1000000-1001000.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `k` | integer | 1 | Maximum neighbors to return |
| `reference` | string/column | implicit | Query position (literal or column) |
| `max_distance` | integer | none | Only include neighbors within distance |
| `stranded` | boolean | false | Only consider same-strand features |
| `signed` | boolean | false | Calculate directional distances |

## Common Patterns

### Distance-Constrained Queries

Find up to 5 genes within 100kb:

```sql
SELECT *
FROM peaks
CROSS JOIN LATERAL NEAREST(genes, k=5, max_distance=100000)
```

### Strand-Specific Queries

Find 3 nearest same-strand features:

```sql
SELECT *
FROM peaks
CROSS JOIN LATERAL NEAREST(genes, k=3, stranded=true)
```

### Directional Queries

Find 3 nearest upstream genes:

```sql
SELECT *
FROM peaks
CROSS JOIN LATERAL NEAREST(genes, k=3, signed=true)
WHERE distance < 0  -- Negative = upstream
```

Find 3 nearest downstream genes:

```sql
SELECT *
FROM peaks
CROSS JOIN LATERAL NEAREST(genes, k=3, signed=true)
WHERE distance > 0  -- Positive = downstream
```

### Combined Constraints

Find up to 5 same-strand genes within 50kb:

```sql
SELECT *
FROM peaks
CROSS JOIN LATERAL NEAREST(
    genes,
    k=5,
    max_distance=50000,
    stranded=true
)
```

## Filtering Results

You can add WHERE clauses to filter results:

```sql
SELECT *
FROM peaks
CROSS JOIN LATERAL NEAREST(genes, k=5)
WHERE distance < 10000  -- Only close neighbors
  AND gene_type = 'protein_coding'  -- Only protein-coding genes
```

## Multiple Query Points

Use a CTE to query multiple positions:

```sql
WITH query_points AS (
    SELECT 'chr1:1000-2000' AS position UNION ALL
    SELECT 'chr2:5000-6000' UNION ALL
    SELECT 'chr3:10000-11000'
)
SELECT *
FROM query_points
CROSS JOIN LATERAL NEAREST(genes, reference=query_points.position, k=3)
```

## Performance Tips

1. **Create indexes** on genomic columns:
   ```sql
   CREATE INDEX idx_genes_chr_start ON genes(chromosome, start_pos);
   CREATE INDEX idx_peaks_chr_start ON peaks(chromosome, start_pos);
   ```

2. **Use max_distance** to reduce candidate set:
   - Filters candidates before ranking
   - Can improve performance 2-10x

3. **Chromosome pre-filtering** is automatic:
   - GIQL automatically adds `WHERE a.chromosome = b.chromosome`
   - ~24x reduction for human genome

4. **Small k is faster**:
   - k≤10 benefits most from optimizations
   - Larger k (50+) may be slower

## Comparison with bedtools

| bedtools | GIQL Equivalent |
|----------|-----------------|
| `bedtools closest` | `NEAREST(target, k=1)` |
| `bedtools closest -d` | Automatic (distance always included) |
| `bedtools closest -t all` | Automatic (ties included) |
| `bedtools closest -s` | `NEAREST(target, stranded=true)` |
| `bedtools closest -D ref` | `NEAREST(target, signed=true)` |
| `bedtools closest -io` | Add `WHERE distance > 0` |

## Expected Performance

With proper indexes:

| Dataset Size | k Value | Typical Query Time |
|--------------|---------|-------------------|
| <10K rows | Any | <100ms |
| 10K-1M rows | 1-10 | <5s |
| 1M+ rows | 1-10 | <30s |

## Example: Peak-Gene Association

Find the 3 closest genes to each ChIP-seq peak:

```sql
SELECT
    peaks.peak_id,
    peaks.chromosome,
    peaks.start_pos AS peak_start,
    genes.gene_name,
    genes.start_pos AS gene_start,
    distance
FROM peaks
CROSS JOIN LATERAL NEAREST(genes, k=3)
ORDER BY peaks.peak_id, distance
```

Result shows each peak with its 3 nearest genes, ordered by distance.
