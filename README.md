# GIQL - Genomic Interval Query Language

A SQL dialect for genomic range queries with multi-database support.

## Status

WIP

## Overview

GIQL extends SQL with spatial operators for genomic interval queries. It transpiles to standard SQL that works across multiple database backends.

## Installation

### From Source

Clone the repository and install locally:

```bash
# Clone the repository
git clone https://github.com/abdennlab/giql.git
cd giql

# Install in development mode
pip install -e .

# Or with development dependencies
pip install -e ".[dev]"
```

### Building Documentation

To build the documentation locally:

```bash
cd docs

# Install documentation dependencies
pip install -r requirements.txt

# Build HTML documentation
make html

# View the documentation
# The built docs will be in docs/_build/html/
# Open docs/_build/html/index.html in your browser
```

## Quick Start

```python
from giql import GIQLEngine

# Create engine with DuckDB backend
with GIQLEngine(target_dialect="duckdb") as engine:
    # Load genomic data
    engine.load_csv("variants", "variants.csv")
    engine.register_table_schema(
        "variants",
        {
            "id": "INTEGER",
            "chromosome": "VARCHAR",
            "start_pos": "BIGINT",
            "end_pos": "BIGINT",
        },
        genomic_column="position",
    )

    # Query with genomic operators (returns cursor for streaming)
    cursor = engine.execute("""
        SELECT * FROM variants
        WHERE position INTERSECTS 'chr1:1000-2000'
    """)

    # Process results lazily
    for row in cursor:
        print(row)

    # Or just transpile to SQL without executing
    sql = engine.transpile("""
        SELECT * FROM variants
        WHERE position INTERSECTS 'chr1:1000-2000'
    """)
    print(sql)  # See the generated SQL
```

## Recipes: Replicating Bedtools with GIQL

These examples show how to replicate common `bedtools intersect` operations using GIQL.

### Setup

```python
from giql import GIQLEngine

# Load two genomic datasets
engine = GIQLEngine(target_dialect="duckdb")
engine.load_csv("features_a", "file_a.bed")
engine.load_csv("features_b", "file_b.bed")

# Register schemas
for table in ["features_a", "features_b"]:
    engine.register_table_schema(
        table,
        {
            "chromosome": "VARCHAR",
            "start_pos": "BIGINT",
            "end_pos": "BIGINT",
            "name": "VARCHAR",
            "score": "FLOAT",
            "strand": "VARCHAR",
        },
        genomic_column="position",
    )
```

**Note:** All examples below use `engine.execute()` which returns a DB-API 2.0 cursor for lazy/streaming iteration. To process results:

```python
# Iterate over results
cursor = engine.execute("SELECT ...")
for row in cursor:
    print(row)

# Or materialize to pandas DataFrame
import pandas as pd
cursor = engine.execute("SELECT ...")
df = pd.DataFrame(cursor.fetchall(), columns=[desc[0] for desc in cursor.description])
```

### Basic Intersect Operations

#### Default: Report overlaps between A and B

**Bedtools:**
```bash
bedtools intersect -a file_a.bed -b file_b.bed
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT DISTINCT a.*
    FROM features_a a, features_b b
    WHERE a.position INTERSECTS b.position
""")
```

#### `-wa`: Write original A entry for each overlap

**Bedtools:**
```bash
bedtools intersect -a file_a.bed -b file_b.bed -wa
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT a.*
    FROM features_a a, features_b b
    WHERE a.position INTERSECTS b.position
""")
```

#### `-wb`: Write original B entry for each overlap

**Bedtools:**
```bash
bedtools intersect -a file_a.bed -b file_b.bed -wb
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT b.*
    FROM features_a a, features_b b
    WHERE a.position INTERSECTS b.position
""")
```

#### `-wa -wb`: Write both A and B entries

**Bedtools:**
```bash
bedtools intersect -a file_a.bed -b file_b.bed -wa -wb
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT a.*, b.*
    FROM features_a a, features_b b
    WHERE a.position INTERSECTS b.position
""")
```

### Filtering Operations

#### `-v`: Report A entries with NO overlap in B

**Bedtools:**
```bash
bedtools intersect -a file_a.bed -b file_b.bed -v
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT a.*
    FROM features_a a
    LEFT JOIN features_b b ON a.position INTERSECTS b.position
    WHERE b.chromosome IS NULL
""")
```

#### `-u`: Report A entries with ANY overlap in B (unique)

**Bedtools:**
```bash
bedtools intersect -a file_a.bed -b file_b.bed -u
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT DISTINCT a.*
    FROM features_a a
    INNER JOIN features_b b ON a.position INTERSECTS b.position
""")
```

### Counting Operations

#### `-c`: Count B overlaps for each A feature

**Bedtools:**
```bash
bedtools intersect -a file_a.bed -b file_b.bed -c
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT a.*, COUNT(b.name) as overlap_count
    FROM features_a a
    LEFT JOIN features_b b ON a.position INTERSECTS b.position
    GROUP BY a.chromosome, a.start_pos, a.end_pos, a.name, a.score, a.strand
""")
```

### Strand-Specific Operations

#### `-s`: Same strand overlaps only

**Bedtools:**
```bash
bedtools intersect -a file_a.bed -b file_b.bed -s
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT a.*
    FROM features_a a, features_b b
    WHERE a.position INTERSECTS b.position
      AND a.strand = b.strand
""")
```

#### `-S`: Opposite strand overlaps only

**Bedtools:**
```bash
bedtools intersect -a file_a.bed -b file_b.bed -S
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT a.*
    FROM features_a a, features_b b
    WHERE a.position INTERSECTS b.position
      AND a.strand != b.strand
      AND a.strand IN ('+', '-')
      AND b.strand IN ('+', '-')
""")
```

### Overlap Fraction Requirements

#### `-f`: Minimum overlap fraction of A

**Bedtools:**
```bash
bedtools intersect -a file_a.bed -b file_b.bed -f 0.5
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT a.*
    FROM features_a a, features_b b
    WHERE a.position INTERSECTS b.position
      AND (
          LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)
      ) >= 0.5 * (a.end_pos - a.start_pos)
""")
```

#### `-F`: Minimum overlap fraction of B

**Bedtools:**
```bash
bedtools intersect -a file_a.bed -b file_b.bed -F 0.5
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT a.*
    FROM features_a a, features_b b
    WHERE a.position INTERSECTS b.position
      AND (
          LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)
      ) >= 0.5 * (b.end_pos - b.start_pos)
""")
```

#### `-r`: Reciprocal overlap (both A and B must meet fraction)

**Bedtools:**
```bash
bedtools intersect -a file_a.bed -b file_b.bed -f 0.5 -r
```

**GIQL:**
```python
cursor = engine.execute("""
    WITH overlap_calcs AS (
        SELECT
            a.*,
            (LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)) as overlap_bp,
            (a.end_pos - a.start_pos) as a_length,
            (b.end_pos - b.start_pos) as b_length
        FROM features_a a, features_b b
        WHERE a.position INTERSECTS b.position
    )
    SELECT chromosome, start_pos, end_pos, name, score, strand
    FROM overlap_calcs
    WHERE overlap_bp >= 0.5 * a_length
      AND overlap_bp >= 0.5 * b_length
""")
```

### Join Operations

#### `-loj`: Left outer join (report all A, with NULL for non-overlapping)

**Bedtools:**
```bash
bedtools intersect -a file_a.bed -b file_b.bed -loj
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT a.*, b.*
    FROM features_a a
    LEFT JOIN features_b b ON a.position INTERSECTS b.position
""")
```

#### `-wo`: Write overlap amount in base pairs

**Bedtools:**
```bash
bedtools intersect -a file_a.bed -b file_b.bed -wo
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT
        a.*,
        b.*,
        (LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)) as overlap_bp
    FROM features_a a, features_b b
    WHERE a.position INTERSECTS b.position
""")
```

#### `-wao`: Write overlap amount for ALL A features

**Bedtools:**
```bash
bedtools intersect -a file_a.bed -b file_b.bed -wao
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT
        a.*,
        b.*,
        CASE
            WHEN b.chromosome IS NULL THEN 0
            ELSE LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)
        END as overlap_bp
    FROM features_a a
    LEFT JOIN features_b b ON a.position INTERSECTS b.position
""")
```

### Multiple Database Files

#### Intersect with multiple B files

**Bedtools:**
```bash
bedtools intersect -a file_a.bed -b file1.bed file2.bed file3.bed
```

**GIQL:**
```python
# Load multiple files
engine.load_csv("features_b1", "file1.bed")
engine.load_csv("features_b2", "file2.bed")
engine.load_csv("features_b3", "file3.bed")

# Register schemas for each...

# Query using UNION to combine all B features
cursor = engine.execute("""
    WITH all_b_features AS (
        SELECT * FROM features_b1
        UNION ALL
        SELECT * FROM features_b2
        UNION ALL
        SELECT * FROM features_b3
    )
    SELECT DISTINCT a.*
    FROM features_a a
    INNER JOIN all_b_features b ON a.position INTERSECTS b.position
""")
```

### Advanced Queries

#### Find features in A that overlap at least 2 features in B

**Bedtools:**
```bash
bedtools intersect -a file_a.bed -b file_b.bed -c | awk '$NF >= 2'
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT a.*
    FROM features_a a
    INNER JOIN features_b b ON a.position INTERSECTS b.position
    GROUP BY a.chromosome, a.start_pos, a.end_pos, a.name, a.score, a.strand
    HAVING COUNT(*) >= 2
""")
```

#### Complex filtering: Overlaps specific genes with quality threshold

**GIQL:**
```python
cursor = engine.execute("""
    SELECT v.*, g.name as gene_name
    FROM variants v
    INNER JOIN genes g ON v.position INTERSECTS g.position
    WHERE v.quality >= 30
      AND g.name IN ('BRCA1', 'BRCA2', 'TP53')
      AND v.chromosome = g.chromosome
    ORDER BY v.chromosome, v.start_pos
""")
```

#### Aggregate statistics per chromosome

**GIQL:**
```python
cursor = engine.execute("""
    SELECT
        a.chromosome,
        COUNT(DISTINCT a.name) as total_features,
        COUNT(b.name) as total_overlaps,
        AVG(LEAST(a.end_pos, b.end_pos) - GREATEST(a.start_pos, b.start_pos)) as avg_overlap_bp
    FROM features_a a
    LEFT JOIN features_b b ON a.position INTERSECTS b.position
    GROUP BY a.chromosome
    ORDER BY a.chromosome
""")
```

### Clustering and Merging Operations

GIQL provides dedicated operators for clustering and merging overlapping genomic intervals, replicating `bedtools cluster` and `bedtools merge` functionality.

#### Cluster: Assign cluster IDs to overlapping intervals

**Basic clustering:**

**Bedtools:**
```bash
bedtools cluster -i features.bed
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT
        *,
        CLUSTER(position) AS cluster_id
    FROM features
    ORDER BY chromosome, start_pos
""")
```

**Clustering with distance parameter** (merge intervals within 1000bp):

**Bedtools:**
```bash
bedtools cluster -i features.bed -d 1000
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT
        *,
        CLUSTER(position, 1000) AS cluster_id
    FROM features
    ORDER BY chromosome, start_pos
""")
```

**Strand-specific clustering:**

**Bedtools:**
```bash
bedtools cluster -i features.bed -s
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT
        *,
        CLUSTER(position, stranded=true) AS cluster_id
    FROM features
    ORDER BY chromosome, start_pos
""")
```

#### Merge: Combine overlapping intervals

**Basic merge:**

**Bedtools:**
```bash
bedtools merge -i features.bed
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT MERGE(position)
    FROM features
""")
```

Returns `chromosome`, `start_pos`, and `end_pos` of merged intervals.

**Merge with distance parameter:**

**Bedtools:**
```bash
bedtools merge -i features.bed -d 1000
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT MERGE(position, 1000)
    FROM features
""")
```

**Strand-specific merge:**

**Bedtools:**
```bash
bedtools merge -i features.bed -s
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT MERGE(position, stranded=true)
    FROM features
""")
```

**Merge with aggregations** (count features):

**Bedtools:**
```bash
bedtools merge -i features.bed -c 1 -o count
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT
        MERGE(position),
        COUNT(*) as feature_count
    FROM features
""")
```

**Merge with column operations** (average score):

**Bedtools:**
```bash
bedtools merge -i features.bed -c 5 -o mean
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT
        MERGE(position),
        AVG(score) as avg_score
    FROM features
""")
```

**Collect names of merged features:**

**Bedtools:**
```bash
bedtools merge -i features.bed -c 4 -o collapse
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT
        MERGE(position),
        STRING_AGG(name, ',') as feature_names
    FROM features
""")
```

### Distance Calculation

GIQL provides a `DISTANCE()` operator for calculating genomic distances between intervals, replicating `bedtools closest -d` functionality.

#### Calculate distance between intervals

**Bedtools:**
```bash
bedtools closest -a peaks.bed -b genes.bed -d
```

**GIQL:**
```python
cursor = engine.execute("""
    SELECT
        a.name AS peak,
        b.name AS gene,
        DISTANCE(a.position, b.position) AS distance
    FROM peaks a
    CROSS JOIN genes b
    WHERE a.chromosome = b.chromosome
    ORDER BY a.name, distance
""")
```

**Distance calculation rules:**
- Returns `0` for overlapping intervals
- Returns positive integer (gap in base pairs) for non-overlapping intervals
- Returns `NULL` for intervals on different chromosomes

**Performance tip:** Always include `WHERE a.chromosome = b.chromosome` to avoid unnecessary cross-chromosome comparisons.

## Transpiling GIQL to SQL

The `transpile()` method allows you to convert GIQL queries to standard SQL without executing them. This is useful for debugging, understanding the generated SQL, or integrating GIQL with external tools:

```python
from giql import GIQLEngine

with GIQLEngine(target_dialect="duckdb") as engine:
    # Transpile GIQL to target SQL dialect
    sql = engine.transpile("""
        SELECT * FROM variants
        WHERE position INTERSECTS 'chr1:1000-2000'
    """)

    print(sql)
    # Output: SELECT * FROM variants WHERE chromosome = 'chr1' AND start_pos < 2000 AND end_pos > 1000

# Different dialects generate different SQL
with GIQLEngine(target_dialect="sqlite") as engine:
    sql = engine.transpile("SELECT * FROM variants WHERE position INTERSECTS 'chr1:1000-2000'")
    # Generates SQLite-compatible SQL
```

The transpiled SQL can be executed directly on your database or used with other tools. The `transpile()` method respects the engine's `verbose` setting to print detailed transpilation information.

## Development

This project is in active development. See `giql.md` for the development plan.
