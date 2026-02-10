# GIQL - Genomic Interval Query Language

<samp>
  <p align="center">
      <a href="https://giql.readthedocs.io/">docs</a> |
      <a href="https://giql.readthedocs.io/dialect">syntax</a> |
      <a href="https://giql.readthedocs.io/transpilation">transpiler</a>
  </p>
</samp>

GIQL is an extended SQL dialect that allows you to declaratively express genomic interval operations.

The `giql` Python package transpiles GIQL queries into standard SQL syntax for execution on any database or analytics engine.

> **Note:** This project is in active development — APIs, syntax, and behavior may change.

## Installation

To install the transpiler:

```bash
pip install giql
```

## Usage (transpilation)

The `giql` package transpiles GIQL queries to standard SQL.

```python
from giql import transpile

sql = transpile(
    "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
    tables=["peaks"],
)
print(sql)
```

Each table referenced in a GIQL query exposes a genomic "pseudo-column" that maps to separate logical chromosome, start, end, and strand columns. You can customize the column mappings.

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
print(sql)
```

The transpiled SQL can be executed with fast genome-unaware databases or in-memory analytic engines like DuckDB.

You can also use [oxbow](https://oxbow.readthedocs.io) to efficiently stream specialized genomics formats into DuckDB. 

```python
import duckdb
import oxbow as ox
from giql import transpile

conn = duckdb.connect()

# Load a streaming data source as a DuckDB relation
peaks = ox.from_bed("peaks.bed", bed_schema="bed6+4").to_duckdb(conn)

sql = transpile(
    "SELECT * FROM peaks WHERE interval INTERSECTS 'chr1:1000-2000'",
    tables=["peaks"],
)

# Execute and return the output as a dataframe
df = con.execute(sql).fetchdf()
```

## Development

```bash
git clone https://github.com/abdenlab/giql.git
cd giql
uv sync
```

To build the documentation locally:

```bash
uv run --group docs sphinx-build docs docs/_build
# The built docs will be in docs/_build/html/
```

For serve the docs locally with automatic rebuild:
```bash
uv run --group docs sphinx-autobuild docs docs/_build
```
