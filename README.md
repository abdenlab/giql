# GIQL - Genomic Interval Query Language

A SQL dialect for genomic range queries. Transpiles to standard SQL.


## Overview

GIQL extends SQL with spatial operators for genomic interval queries. It transpiles GIQL queries into standard SQL that can be executed on any database backend.

GIQL provides a familiar SQL syntax for bioinformatics workflows, allowing you to express complex genomic range operations without writing intricate SQL expressions. Whether you're filtering variants by genomic region, finding overlapping features, or calculating distances between intervals, GIQL makes these operations intuitive and portable.

## Features

- **SQL-based**: Familiar SQL syntax with genomic extensions
- **Spatial operators**: INTERSECTS, CONTAINS, WITHIN for range relationships
- **Distance operators**: DISTANCE, NEAREST for proximity queries
- **Aggregation operators**: CLUSTER, MERGE for combining intervals
- **Set quantifiers**: ANY, ALL for multi-range queries
- **Transpilation**: Converts GIQL to standard SQL for execution on any backend

## Installation

### From PyPI

Install the latest stable release:

```bash
pip install giql
```

Or the latest release candidate:

```bash
pip install --pre giql
```

### From Source

Clone the repository and install locally:

```bash
# Clone the repository
git clone https://github.com/abdenlab/giql.git
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
from giql import transpile

# Transpile a GIQL query to standard SQL
sql = transpile(
    "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
    tables=["peaks"],
)
print(sql)
```

With custom column mappings:

```python
from giql import Table, transpile

sql = transpile(
    "SELECT * FROM variants WHERE position INTERSECTS 'chr1:1000-2000'",
    tables=[
        Table(
            "variants",
            genomic_col="position",
            chrom_col="chromosome",
            start_col="start_pos",
            end_col="end_pos",
        )
    ],
)
```

Execution example with DuckDB:

```python
import duckdb
import oxbow as ox
from giql import transpile

conn = duckdb.connect()
peaks = ox.from_bed("peaks.bed", bed_schema="bed6+4").to_duckdb(conn)  # streaming source

sql = transpile(
    "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
    tables=["peaks"],
)
df = con.execute(sql).fetchdf()
```

## Operators at a Glance

### Spatial Relationships

| Operator | Description |
|----------|-------------|
| `INTERSECTS` | Returns true when ranges overlap by at least one base pair |
| `CONTAINS` | Returns true when one range fully contains another |
| `WITHIN` | Returns true when one range is fully within another |

### Distance and Proximity

| Operator | Description |
|----------|-------------|
| `DISTANCE` | Calculate genomic distance between two intervals |
| `NEAREST` | Find k-nearest genomic features |

### Aggregation

| Operator | Description |
|----------|-------------|
| `CLUSTER` | Assign cluster IDs to overlapping intervals |
| `MERGE` | Combine overlapping intervals into unified regions |

### Set Quantifiers

| Quantifier | Description |
|------------|-------------|
| `ANY` | Match if condition holds for any of the specified ranges |
| `ALL` | Match if condition holds for all of the specified ranges |

## Documentation

For complete documentation, build the docs locally (see above) or visit the hosted documentation.

The documentation includes:

- **Operator Reference**: Detailed documentation for each operator with examples
- **Recipes**: Common query patterns for intersections, distance calculations, and clustering
- **Bedtools Migration Guide**: How to replicate bedtools operations with GIQL
- **Guides**: Performance optimization, multi-backend configuration, and schema mapping

## Development

This project is in active development.
