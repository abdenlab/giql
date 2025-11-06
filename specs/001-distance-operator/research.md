# Research: DISTANCE Operator Implementation

**Date**: 2025-11-06
**Status**: Complete

## Overview

Research findings for implementing the DISTANCE() operator as a transpiled SQL expression (not an actual UDF). The operator will follow the same pattern as INTERSECTS/CONTAINS/WITHIN spatial predicates, generating pure SQL CASE expressions.

## 1. SQL Expression Pattern (Not UDF Registration)

### Decision

DISTANCE will be transpiled into a SQL CASE expression, **not** registered as a database UDF. This follows the existing pattern for spatial predicates and ensures maximum portability.

### Rationale

1. **Consistency**: IN TERSECTS, CONTAINS, and WITHIN all transpile to SQL predicates
2. **Portability**: Pure SQL works identically across all dialects without UDF registration
3. **Simplicity**: No need to manage UDF lifecycle or handle different registration APIs
4. **Performance**: Database can optimize inline SQL expressions better than UDF calls

### Implementation Pattern

From `src/giql/generators/base.py` (lines 168-243), spatial predicates transpile column-to-column operations:

```python
def _generate_column_join(self, left_col: str, right_col: str, op_type: str) -> str:
    l_chrom, l_start, l_end = self._get_column_refs(left_col, None)
    r_chrom, r_start, r_end = self._get_column_refs(right_col, None)

    if op_type == "intersects":
        return (
            f"({l_chrom} = {r_chrom} "
            f"AND {l_start} < {r_end} "
            f"AND {l_end} > {r_start})"
        )
```

DISTANCE will follow this pattern but generate a CASE expression instead of a boolean predicate.

## 2. SQL CASE Expression for Distance Calculation

### SQL Template

```sql
CASE
    -- Different chromosomes
    WHEN {chrom_a} != {chrom_b} THEN NULL

    -- Different strands (when stranded mode enabled)
    WHEN {stranded_flag} AND {strand_a} != {strand_b} THEN NULL

    -- Overlapping intervals
    WHEN {start_a} < {end_b} AND {end_a} > {start_b} THEN 0

    -- A is before B
    WHEN {end_a} <= {start_b} THEN
        CASE WHEN {signed_flag} THEN ({start_b} - {end_a})
             ELSE ({start_b} - {end_a}) END

    -- B is before A
    ELSE
        CASE WHEN {signed_flag} THEN -({start_a} - {end_b})
             ELSE ({start_a} - {end_b}) END
END
```

### Cross-Dialect Compatibility

**Tested across DuckDB, SQLite, PostgreSQL:**
- ✅ CASE expressions: Standard SQL-92, supported by all
- ✅ NULL handling: Identical across all dialects
- ✅ Integer arithmetic: Standard across all
- ✅ Boolean expressions in WHEN: All support
- ⚠️ Boolean literals: Need to handle as 0/1 for SQLite (see section 3)

## 3. Handling Named Parameters (stranded, signed)

### Current Pattern from CLUSTER/MERGE

From `src/giql/expressions.py` (lines 94-123), GIQL functions handle named parameters via `from_arg_list()`:

```python
class GIQLCluster(exp.Func):
    arg_types = {
        "this": True,  # genomic column
        "distance": False,  # optional
        "stranded": False,  # optional
    }

    @classmethod
    def from_arg_list(cls, args):
        kwargs = {}
        positional_args = []

        for arg in args:
            if isinstance(arg, exp.EQ):
                # Named parameter: extract name and value
                param_name = arg.this.name if isinstance(arg.this, exp.Column) else str(arg.this)
                kwargs[param_name.lower()] = arg.expression
            else:
                positional_args.append(arg)

        # Map positional to kwargs
        if len(positional_args) > 0:
            kwargs["this"] = positional_args[0]

        return cls(**kwargs)
```

### Implementation for DISTANCE

```python
class GIQLDistance(exp.Func):
    arg_types = {
        "this": True,          # interval_a
        "expression": True,    # interval_b
        "stranded": False,     # optional boolean
        "signed": False,       # optional boolean
    }

    @classmethod
    def from_arg_list(cls, args):
        # Same pattern as CLUSTER/MERGE
        # Extracts named parameters from EQ nodes
        # Maps to arg_types dictionary
```

### Boolean Parameter Handling in SQL Generation

**Key insight**: Named parameters are resolved at transpilation time, not runtime. The SQL generator evaluates `stranded` and `signed` flags and generates the appropriate CASE expression.

**SQLite Boolean Handling**: SQLite doesn't have a native BOOLEAN type. Store/pass as INTEGER (0=false, 1=true). This is handled transparently by Python's sqlite3 module when we evaluate the parameters at transpilation time.

Since parameters are evaluated during transpilation (not at SQL runtime), we don't need to worry about SQL boolean compatibility - we just generate different CASE expressions based on the parameter values.

## 4. NULL Handling in Window Functions (ORDER BY)

### Research Question

How do NULL distance values behave in `ORDER BY DISTANCE(...)` clauses across SQL dialects?

### Findings

**Default NULL sorting behavior:**
- **PostgreSQL**: NULLS LAST (configurable with `NULLS FIRST`/`NULLS LAST`)
- **DuckDB**: NULLS LAST by default
- **SQLite**: NULLS FIRST by default

**Impact on DISTANCE queries:**

When finding closest features:
```sql
ROW_NUMBER() OVER (PARTITION BY ... ORDER BY DISTANCE(...))
```

NULL distances (from different chromosomes/strands) will be sorted differently:
- **PostgreSQL/DuckDB**: NULLs appear after valid distances ✅ (desired behavior)
- **SQLite**: NULLs appear before valid distances ❌ (undesired - need explicit handling)

### Solution

**Explicitly specify NULL ordering** in generated SQL for consistency:

```sql
ROW_NUMBER() OVER (
    PARTITION BY a.chromosome, a.start_pos, a.end_pos
    ORDER BY DISTANCE(a.position, b.position) NULLS LAST
)
```

**Alternative**: Filter NULLs in WHERE clause before window function:
```sql
WHERE a.chromosome = b.chromosome  -- Pre-filter eliminates chromosome mismatches
```

### Recommendation

1. **Document** that users should pre-filter by chromosome for performance and correctness
2. **Examples** should always include `WHERE a.chromosome = b.chromosome`
3. **Optional**: Generator could detect DISTANCE in ORDER BY and automatically append `NULLS LAST`

## 5. Strand Information Resolution

### Research Question

How does the existing codebase resolve position columns to chromosome/start/end? Does it already handle strand columns?

### Findings from Codebase

**Current position resolution** (`src/giql/generators/base.py`, line 296-340):

```python
def _get_column_refs(self, column_ref: str, table_name: str | None = None) -> tuple[str, str, str]:
    """Get physical column names for genomic data.

    Returns: Tuple of (chromosome_col, start_col, end_col)
    """
    chrom_col = DEFAULT_CHROM_COL
    start_col = DEFAULT_START_COL
    end_col = DEFAULT_END_COL

    # Extract table from column_ref if present (e.g., "a.position" -> "a")
    # Look up schema to get actual column names
    # Return resolved column references
```

**Strand handling** (`src/giql/transformer.py`, lines 60-92):

```python
def _get_genomic_columns(self, query: exp.Select) -> tuple[str, str, str, str]:
    """Returns: (chrom_col, start_col, end_col, strand_col)"""
    # Already resolves strand_col from schema!
```

### Key Insight

✅ Strand column resolution **already exists** in transformer! The transformer's `_get_genomic_columns()` method returns a 4-tuple including `strand_col`.

**However**, the SQL generator's `_get_column_refs()` only returns 3 columns (chrom, start, end). We need to:

1. **Extend `_get_column_refs()`** to optionally return strand:
   ```python
   def _get_column_refs(..., include_strand=False):
       if include_strand:
           return (chrom_col, start_col, end_col, strand_col)
       return (chrom_col, start_col, end_col)
   ```

2. **Or** create a new method specifically for DISTANCE:
   ```python
   def _get_column_refs_with_strand(column_ref, table_name):
       return (chrom_col, start_col, end_col, strand_col)
   ```

### Default Strand Handling

**For literal ranges without strand** (e.g., `'chr1:1000-2000'`):
- Default to `strand='.'` (unstranded)
- Parser already supports strand notation: `'chr1:1000-2000:+'`

## 6. Integration with Existing Code

### Files to Modify

1. **`src/giql/expressions.py`**
   - Add `GIQLDistance` class (similar to `GIQLCluster`/`GIQLMerge`)
   - Implement `from_arg_list()` for named parameter parsing

2. **`src/giql/dialect.py`**
   - Add DISTANCE function parsing
   - Register DISTANCE as a function (if not already generic)

3. **`src/giql/generators/base.py`**
   - Add `giqldistance_sql()` method
   - Extend `_get_column_refs()` to include strand when needed
   - Generate CASE expression based on parsed parameters

4. **No transformer changes needed**
   - DISTANCE doesn't require CTE transformation like CLUSTER/MERGE

### SQL Generation Method Structure

```python
def giqldistance_sql(self, expression: GIQLDistance) -> str:
    """Generate SQL CASE expression for DISTANCE function."""

    # Extract parameters from expression
    interval_a = expression.this
    interval_b = expression.args.get("expression")
    stranded = expression.args.get("stranded")  # exp.Boolean or None
    signed = expression.args.get("signed")      # exp.Boolean or None

    # Evaluate boolean parameters (they're AST nodes)
    stranded_enabled = stranded and stranded.this  # Get actual bool value
    signed_enabled = signed and signed.this

    # Resolve column references (including strand if stranded mode)
    if stranded_enabled:
        chrom_a, start_a, end_a, strand_a = self._get_column_refs_with_strand(...)
        chrom_b, start_b, end_b, strand_b = self._get_column_refs_with_strand(...)
    else:
        chrom_a, start_a, end_a = self._get_column_refs(...)
        chrom_b, start_b, end_b = self._get_column_refs(...)
        strand_a = strand_b = None

    # Generate CASE expression
    return self._generate_distance_case(
        chrom_a, start_a, end_a, strand_a,
        chrom_b, start_b, end_b, strand_b,
        stranded_enabled, signed_enabled
    )
```

## Summary of Key Decisions

| Topic | Decision | Rationale |
|-------|----------|-----------|
| Implementation Approach | Transpile to SQL CASE expression | Follows existing pattern, maximizes portability |
| UDF Registration | Not needed | Pure SQL approach |
| Named Parameters | Parse at transpilation time | Evaluated before SQL generation |
| Boolean Handling | Resolve during transpilation | No runtime boolean compatibility issues |
| NULL Ordering | Document + recommend pre-filtering | User-controlled, consistent with best practices |
| Strand Resolution | Extend existing `_get_column_refs()` | Leverage existing infrastructure |
| Testing | Same pattern as spatial predicates | Unit + integration + property-based |

## Alternatives Considered & Rejected

### Alternative 1: Actual UDF Registration
- **Rejected**: Requires database-specific registration, less portable, more complex lifecycle management

### Alternative 2: Scalar Subquery
- **Rejected**: More complex SQL, harder to optimize by database

### Alternative 3: Stored Procedure
- **Rejected**: PostgreSQL-specific, breaks portability requirement

## Implementation Risks & Mitigations

| Risk | Mitigation |
|------|------------|
| Complex CASE expression hard to debug | Generate well-formatted SQL with clear structure |
| Performance with large CROSS JOINs | Document best practices, show pre-filtering examples |
| NULL handling differences across dialects | Explicit NULLS LAST in ORDER BY, document behavior |
| Strand column might not exist | Graceful handling, clear error messages |

## Next Steps

Proceed to Phase 1:
1. Generate `data-model.md` with DISTANCE function interface design
2. Create `quickstart.md` with usage examples
3. Update agent context with implementation decisions
