# GIQL - Genomic Interval Query Language

A SQL dialect for genomic range queries with multi-database support.

## Status

WIP

## Overview

GIQL extends SQL with spatial operators for genomic interval queries. It transpiles to standard SQL that works across multiple database backends including DuckDB and SQLite.

GIQL provides a familiar SQL syntax for bioinformatics workflows, allowing you to express complex genomic range operations without writing intricate SQL expressions. Whether you're filtering variants by genomic region, finding overlapping features, or calculating distances between intervals, GIQL makes these operations intuitive and portable across databases.

## Features

- **SQL-based**: Familiar SQL syntax with genomic extensions
- **Multi-backend**: Works with DuckDB, SQLite, and more
- **Spatial operators**: INTERSECTS, CONTAINS, WITHIN for range relationships
- **Distance operators**: DISTANCE, NEAREST for proximity queries
- **Aggregation operators**: CLUSTER, MERGE for combining intervals
- **Set quantifiers**: ANY, ALL for multi-range queries
- **Transpilation**: Convert GIQL to standard SQL for debugging or external use

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
        genomic_column="interval",
    )

    # Query with genomic operators (returns cursor for streaming)
    cursor = engine.execute("""
        SELECT * FROM variants
        WHERE interval INTERSECTS 'chr1:1000-2000'
    """)

    # Process results lazily
    for row in cursor:
        print(row)

    # Or just transpile to SQL without executing
    sql = engine.transpile("""
        SELECT * FROM variants
        WHERE interval INTERSECTS 'chr1:1000-2000'
    """)
    print(sql)  # See the generated SQL
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
