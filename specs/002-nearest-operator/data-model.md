# Data Model: NEAREST Operator

**Feature**: 002-nearest-operator
**Date**: 2025-11-07

## AST Nodes

### GIQLNearest Expression

```python
class GIQLNearest(exp.Func):
    """NEAREST function for k-nearest neighbor queries."""

    arg_types = {
        "this": True,          # Required: target table name
        "reference": False,    # Optional: position reference
        "k": False,            # Optional: number of neighbors (default=1)
        "max_distance": False, # Optional: distance threshold
        "stranded": False,     # Optional: strand-specific
        "signed": False,       # Optional: directional distance
    }
```

**Parameters**:
- `this`: Target table to search (e.g., `genes`)
- `reference`: Query position (column ref like `peaks.position` or literal like `'chr1:1000-2000'`)
- `k`: Maximum neighbors to return per query (default=1)
- `max_distance`: Only include neighbors within this distance
- `stranded`: Only consider same-strand features (default=false)
- `signed`: Calculate directional distances (default=false)

## Transpilation Context

### Mode Detection

**Standalone Mode**: Reference is required
```sql
FROM NEAREST(genes, reference='chr1:1000-2000', k=3)
```

**Correlated Mode**: Reference defaults to outer table's `.position`
```sql
FROM peaks CROSS JOIN LATERAL NEAREST(genes, k=3)
-- Implicitly uses peaks.position as reference
```

### Resolution Flow

1. **Parse NEAREST** → GIQLNearest AST node
2. **Detect mode**: standalone (explicit reference) vs correlated (LATERAL context)
3. **Resolve target table**: Extract table name, look up schema
4. **Resolve reference**:
   - Literal: Parse genomic range → (chr, start, end, strand)
   - Column ref: Resolve via `_get_column_refs()` → physical columns
   - Implicit: Use `outer_table.position`
5. **Generate SQL**: Dialect-specific (LATERAL or window functions)

## Generated SQL Structure

### PostgreSQL/DuckDB (LATERAL)

```sql
SELECT outer_table.*, nearest_results.*
FROM outer_table
CROSS JOIN LATERAL (
    SELECT target_table.*,
           DISTANCE(...) AS distance
    FROM target_table
    WHERE outer_table.chromosome = target_table.chromosome
      [AND DISTANCE(...) <= max_distance]
    ORDER BY distance
    LIMIT k
) nearest_results
```

### SQLite (Window Functions)

```sql
WITH ranked AS (
    SELECT outer_table.*, target_table.*,
           DISTANCE(...) AS distance,
           RANK() OVER (
               PARTITION BY outer_table.id
               ORDER BY DISTANCE(...)
           ) AS rank
    FROM outer_table
    CROSS JOIN target_table
    WHERE outer_table.chromosome = target_table.chromosome
      [AND DISTANCE(...) <= max_distance]
)
SELECT * FROM ranked WHERE rank <= k
```

## Schema Resolution

**Target Table** → Position Columns:
1. Look up table in `SchemaInfo`
2. Find genomic column (default: `.position`)
3. Extract physical columns: `chromosome`, `start_pos`, `end_pos`, `strand`

**Reference** → Query Position:
- **Literal**: `'chr1:1000-2000'` → Parse via `RangeParser`
- **Column**: `a.position` → Resolve via `_get_column_refs()`
- **Implicit**: Use outer table from LATERAL context

## Result Schema

**Correlated Mode** (both tables):
```
| outer_col1 | outer_col2 | target_col1 | target_col2 | distance |
```

**Standalone Mode** (target only):
```
| target_col1 | target_col2 | distance |
```

## Key Data Flows

1. **Parse**: GIQL query → AST with GIQLNearest node
2. **Analyze**: Detect mode, resolve tables/columns
3. **Transform**: Generate dialect-specific SQL
4. **Execute**: Database runs generated SQL
5. **Return**: Result rows with computed distance column
