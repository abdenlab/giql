# Data Model: DISTANCE Operator

**Date**: 2025-11-06
**Status**: Design Complete

## Overview

This document defines the data model and interface design for the DISTANCE operator, including AST representation, SQL generation pipeline, and data flow through the GIQL transpiler.

## Entity Model

### 1. GIQLDistance (AST Node)

**Location**: `src/giql/expressions.py`

**Purpose**: Represents a DISTANCE function call in the abstract syntax tree

**Class Definition**:
```python
class GIQLDistance(exp.Func):
    """DISTANCE function for calculating genomic distances between intervals.

    Generates SQL CASE expression that computes distance between two genomic
    intervals, with optional strand-specific and signed (directional) modes.

    Examples:
        DISTANCE(a.position, b.position)
        DISTANCE(a.position, 'chr1:1000-2000')
        DISTANCE(a.position, b.position, stranded=true)
        DISTANCE(a.position, b.position, signed=true)
        DISTANCE(a.position, b.position, stranded=true, signed=true)
    """

    arg_types = {
        "this": True,          # Required: interval_a (column ref or literal range)
        "expression": True,    # Required: interval_b (column ref or literal range)
        "stranded": False,     # Optional: boolean for strand-specific distance
        "signed": False,       # Optional: boolean for directional distance
    }

    @classmethod
    def from_arg_list(cls, args):
        """Parse argument list, handling named parameters.

        :param args:
            List of arguments from parser
        :return:
            GIQLDistance instance with properly mapped arguments
        """
        kwargs = {}
        positional_args = []

        # Separate named (EQ) and positional arguments
        for arg in args:
            if isinstance(arg, exp.EQ):
                # Named parameter: extract name and value
                param_name = (
                    arg.this.name if isinstance(arg.this, exp.Column) else str(arg.this)
                )
                kwargs[param_name.lower()] = arg.expression
            else:
                positional_args.append(arg)

        # Map positional arguments
        if len(positional_args) >= 1:
            kwargs["this"] = positional_args[0]
        if len(positional_args) >= 2:
            kwargs["expression"] = positional_args[1]

        return cls(**kwargs)
```

**Field Specifications**:

| Field | Type | Required | Description | Example |
|-------|------|----------|-------------|---------|
| `this` | exp.Column \| exp.Literal | Yes | First genomic interval | `a.position` or `'chr1:100-200'` |
| `expression` | exp.Column \| exp.Literal | Yes | Second genomic interval | `b.position` or `'chr2:300-400'` |
| `stranded` | exp.Boolean \| None | No | Strand-specific mode flag | `true` (default: `false`) |
| `signed` | exp.Boolean \| None | No | Directional distance flag | `true` (default: `false`) |

**Validation Rules**:
- Both `this` and `expression` must be present
- `stranded` and `signed` must be boolean literals if provided
- If a literal range is provided, it must be parseable by `RangeParser.parse()`

### 2. DistanceCalculation (Conceptual Model)

**Purpose**: Represents the logical operation performed by DISTANCE

**Inputs**:
```python
@dataclass
class GenomicInterval:
    chromosome: str
    start: int
    end: int
    strand: str  # '+', '-', or '.'

@dataclass
class DistanceParameters:
    stranded: bool = False  # If true, return NULL for different strands
    signed: bool = False    # If true, return negative for upstream
```

**Output**:
```python
Optional[int]  # Distance in base pairs, or NULL
```

**Distance Calculation Logic**:
```python
def calculate_distance(
    interval_a: GenomicInterval,
    interval_b: GenomicInterval,
    params: DistanceParameters
) -> Optional[int]:
    """Calculate genomic distance between two intervals.

    Returns:
        - NULL: Different chromosomes or (different strands when stranded=true)
        - 0: Overlapping intervals
        - Positive integer: Gap distance (unsigned mode)
        - Negative/positive integer: Directional distance (signed mode)
    """

    # Different chromosomes
    if interval_a.chromosome != interval_b.chromosome:
        return None

    # Different strands (when stranded mode enabled)
    if params.stranded and interval_a.strand != interval_b.strand:
        return None

    # Check for overlap
    if interval_a.start < interval_b.end and interval_a.end > interval_b.start:
        return 0

    # Calculate gap distance
    if interval_a.end <= interval_b.start:
        # A is before B
        distance = interval_b.start - interval_a.end
        return distance if not params.signed else distance  # positive = downstream

    else:
        # B is before A
        distance = interval_a.start - interval_b.end
        return distance if not params.signed else -distance  # negative = upstream
```

## Data Flow Pipeline

### 1. Parsing Phase

**Input**: GIQL query string
```sql
SELECT a.*, b.*, DISTANCE(a.position, b.position, stranded=true) as dist
FROM features_a a
CROSS JOIN features_b b
```

**Process**:
1. Parser (`src/giql/dialect.py`) recognizes `DISTANCE` as a function
2. Calls `GIQLDistance.from_arg_list()` with parsed arguments
3. Creates `GIQLDistance` AST node with:
   - `this` = `Column(table='a', name='position')`
   - `expression` = `Column(table='b', name='position')`
   - `stranded` = `Boolean(this=True)`
   - `signed` = `None`

**Output**: AST with `GIQLDistance` nodes

### 2. SQL Generation Phase

**Input**: AST with `GIQLDistance` nodes

**Process** (`src/giql/generators/base.py`):

1. **Generator encounters GIQLDistance node**:
   ```python
   def giqldistance_sql(self, expression: GIQLDistance) -> str
   ```

2. **Extract and evaluate parameters**:
   ```python
   stranded_enabled = (expression.args.get("stranded") and
                       expression.args["stranded"].this)  # Get boolean value
   signed_enabled = (expression.args.get("signed") and
                     expression.args["signed"].this)
   ```

3. **Resolve position columns to physical columns**:
   ```python
   # If stranded mode, need strand columns
   if stranded_enabled:
       chrom_a, start_a, end_a, strand_a = self._get_column_refs_with_strand(
           interval_a_ref, table_a
       )
       chrom_b, start_b, end_b, strand_b = self._get_column_refs_with_strand(
           interval_b_ref, table_b
       )
   else:
       chrom_a, start_a, end_a = self._get_column_refs(interval_a_ref, table_a)
       chrom_b, start_b, end_b = self._get_column_refs(interval_b_ref, table_b)
       strand_a = strand_b = None
   ```

4. **Generate SQL CASE expression**:
   ```python
   sql = f"""CASE
       WHEN {chrom_a} != {chrom_b} THEN NULL
       {"WHEN " + strand_a + " != " + strand_b + " THEN NULL" if stranded_enabled else ""}
       WHEN {start_a} < {end_b} AND {end_a} > {start_b} THEN 0
       WHEN {end_a} <= {start_b} THEN
           {"(" + start_b + " - " + end_a + ")" if not signed_enabled
            else "(" + start_b + " - " + end_a + ")"}
       ELSE
           {"(" + start_a + " - " + end_b + ")" if not signed_enabled
            else "-(" + start_a + " - " + end_b + ")"}
   END"""
   ```

**Output**: SQL string with CASE expression

### 3. Example Data Flow

**Input GIQL**:
```sql
DISTANCE(a.position, b.position, signed=true)
```

**Parsed AST**:
```python
GIQLDistance(
    this=Column(table='a', name='position'),
    expression=Column(table='b', name='position'),
    stranded=None,
    signed=Boolean(this=True)
)
```

**Column Resolution** (assuming default schema):
```python
a.position → (a.chromosome, a.start_pos, a.end_pos)
b.position → (b.chromosome, b.start_pos, b.end_pos)
```

**Generated SQL**:
```sql
CASE
    WHEN a.chromosome != b.chromosome THEN NULL
    WHEN a.start_pos < b.end_pos AND a.end_pos > b.start_pos THEN 0
    WHEN a.end_pos <= b.start_pos THEN (b.start_pos - a.end_pos)
    ELSE -(a.start_pos - b.end_pos)
END
```

## Schema Integration

### Position Column Resolution

**Current Schema System** (`src/giql/schema.py`):

```python
@dataclass
class ColumnInfo:
    name: str
    is_genomic: bool = False
    chrom_col: Optional[str] = None
    start_col: Optional[str] = None
    end_col: Optional[str] = None
    strand_col: Optional[str] = None  # Already exists!
```

**Usage in DISTANCE**:

When user writes `DISTANCE(a.position, b.position)`:
1. Schema lookup finds `position` is genomic column for table `a`
2. Retrieves `chrom_col`, `start_col`, `end_col`, `strand_col` mappings
3. Generates SQL using actual column names (e.g., `a.chr`, `a.start`, `a.end`, `a.strand`)

**Default Column Names** (when no schema provided):
- Chromosome: `chromosome`
- Start: `start_pos`
- End: `end_pos`
- Strand: `strand`

## Type System

### Input Types

| GIQL Type | SQL Type | Python Type | Description |
|-----------|----------|-------------|-------------|
| Position Column | Column Reference | str | Reference to genomic position column |
| Genomic Range Literal | String Literal | str | Parseable range like `'chr1:1000-2000'` |
| Boolean Parameter | Boolean Literal | bool | `true` or `false` |

### Output Type

| GIQL Type | SQL Type | Python Type | Description |
|-----------|----------|-------------|-------------|
| Distance | INTEGER | int \| None | Base pairs between intervals |

### NULL Semantics

**DISTANCE returns NULL when**:
- Intervals on different chromosomes
- Intervals on different strands (when `stranded=true`)
- Any input column is NULL (chromosome, start, end, or strand)

**SQL NULL propagation**:
- NULL in arithmetic operations → NULL
- NULL in comparisons → NULL (treated as false in WHEN clauses)
- NULL in ORDER BY → Sorted last (with explicit `NULLS LAST`)

## Error Handling

### Parse-Time Errors

| Error | When | Message |
|-------|------|---------|
| Missing required argument | `DISTANCE(a.position)` | "DISTANCE requires 2 arguments, got 1" |
| Invalid parameter name | `DISTANCE(a.pos, b.pos, foo=true)` | "Unknown parameter 'foo' for DISTANCE" |
| Invalid parameter type | `DISTANCE(a.pos, b.pos, stranded=123)` | "Parameter 'stranded' must be boolean, got integer" |

### Runtime Errors (SQL Generation)

| Error | When | Message |
|-------|------|---------|
| Column not genomic | Position column not registered in schema | "Column 'a.pos' is not a genomic position column" |
| Missing strand column | `stranded=true` but table has no strand column | "Table 'features_a' has no strand column (required for stranded=true)" |
| Invalid range literal | `DISTANCE(a.pos, 'invalid')` | "Could not parse genomic range: 'invalid'" |

## Performance Characteristics

### Computational Complexity

| Operation | Complexity | Notes |
|-----------|------------|-------|
| Parse DISTANCE node | O(1) | Simple AST node creation |
| Resolve columns | O(1) | Schema lookup |
| Generate SQL | O(1) | String concatenation |
| Execute SQL | O(n*m) | Database-dependent (CROSS JOIN of n×m rows) |

### SQL Expression Size

- Base CASE expression: ~200-300 characters
- With stranded mode: +50 characters
- With signed mode: +30 characters per branch
- Total: ~300-400 characters per DISTANCE call

**Impact**: Negligible - modern databases handle expressions of this size efficiently

### Optimization Opportunities

1. **Pre-filtering**: `WHERE a.chromosome = b.chromosome` reduces CROSS JOIN size
2. **Indexing**: Indexes on chromosome, start_pos columns improve JOIN performance
3. **Materialization**: For repeated distance calculations, consider materializing results

## Testing Strategy

### Unit Tests

**Test AST Node Creation**:
```python
def test_giqldistance_from_arg_list_positional():
    """Test GIQLDistance parses positional arguments.

    Given:
        Argument list with two positional arguments
    When:
        from_arg_list() is called
    Then:
        GIQLDistance node created with correct this and expression
    """

def test_giqldistance_from_arg_list_named_parameters():
    """Test GIQLDistance parses named parameters.

    Given:
        Argument list with stranded=true and signed=false
    When:
        from_arg_list() is called
    Then:
        GIQLDistance node created with correct boolean values
    """
```

**Test SQL Generation**:
```python
def test_giqldistance_sql_basic():
    """Test basic DISTANCE transpilation.

    Given:
        DISTANCE(a.position, b.position)
    When:
        giqldistance_sql() generates SQL
    Then:
        CASE expression with NULL, overlap, and distance logic
    """

def test_giqldistance_sql_stranded():
    """Test stranded DISTANCE transpilation.

    Given:
        DISTANCE(a.position, b.position, stranded=true)
    When:
        giqldistance_sql() generates SQL
    Then:
        CASE expression includes strand comparison
    """
```

### Integration Tests

**End-to-End Query Tests**:
```python
def test_distance_query_execution():
    """Test DISTANCE in complete query.

    Given:
        Two tables with genomic features
    When:
        Query with DISTANCE is executed
    Then:
        Results match expected distances
    """
```

### Property-Based Tests

```python
@given(interval_a=genomic_interval(), interval_b=genomic_interval())
def test_distance_symmetric_unsigned(interval_a, interval_b):
    """DISTANCE(a, b) == DISTANCE(b, a) when signed=false."""

@given(intervals=st.lists(genomic_interval(), min_size=2, max_size=2))
def test_distance_overlap_is_zero(intervals):
    """Overlapping intervals return distance=0."""
```

## Summary

The DISTANCE operator follows GIQL's established patterns:
- ✅ AST representation similar to CLUSTER/MERGE
- ✅ SQL generation similar to INTERSECTS/WITHIN
- ✅ Schema integration via existing position column system
- ✅ Named parameter handling via from_arg_list()
- ✅ Clean separation between parse-time and generation-time logic

**Next**: See `quickstart.md` for user-facing usage examples
