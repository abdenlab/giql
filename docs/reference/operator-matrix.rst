Operator Compatibility Matrix
=============================

This reference shows which GIQL operators are supported by each database backend.

Backend Support Overview
------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 15 15 15 25

   * - Operator
     - DuckDB
     - SQLite
     - PostgreSQL
     - Notes
   * - **Spatial Operators**
     -
     -
     -
     -
   * - INTERSECTS
     - ✅
     - ✅
     - 🚧
     -
   * - CONTAINS
     - ✅
     - ✅
     - 🚧
     -
   * - WITHIN
     - ✅
     - ✅
     - 🚧
     -
   * - **Distance Operators**
     -
     -
     -
     -
   * - DISTANCE
     - ✅
     - ✅
     - 🚧
     -
   * - NEAREST
     - ✅
     - ⚠️
     - 🚧
     - SQLite: slower for large k
   * - **Aggregation Operators**
     -
     -
     -
     -
   * - CLUSTER
     - ✅
     - ✅
     - 🚧
     -
   * - MERGE
     - ✅
     - ✅
     - 🚧
     -
   * - **Set Quantifiers**
     -
     -
     -
     -
   * - ANY
     - ✅
     - ✅
     - 🚧
     -
   * - ALL
     - ✅
     - ✅
     - 🚧
     -

Legend
------

.. list-table::
   :widths: 10 90

   * - ✅
     - **Full support** - Operator works with full functionality
   * - ⚠️
     - **Partial support** - Operator works but with limitations
   * - 🚧
     - **Planned** - Support planned for future release
   * - ❌
     - **Not supported** - Operator not available for this backend

Operator Details by Backend
---------------------------

DuckDB
~~~~~~

All operators are fully supported on DuckDB. DuckDB is the recommended backend
for most use cases due to its excellent analytical query performance.

**Strengths:**

- Efficient columnar storage
- Parallel query execution
- Full LATERAL join support (used by NEAREST)
- Rich window function support (used by CLUSTER)

SQLite
~~~~~~

All operators work on SQLite, with some performance considerations:

**NEAREST operator:**

- Works correctly but may be slower for large k values
- Performance depends on table size and index availability
- Consider using ``max_distance`` to limit search space

**CLUSTER and MERGE:**

- Full functionality
- May be slower than DuckDB for very large datasets

PostgreSQL (Planned)
~~~~~~~~~~~~~~~~~~~~

PostgreSQL support is planned for a future release. Expected to have full
operator support.

SQL Feature Requirements
------------------------

GIQL operators require certain SQL features from the underlying database:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - SQL Feature
     - Used By
   * - Basic predicates (AND, OR, comparison)
     - All spatial operators
   * - CASE expressions
     - DISTANCE, conditional logic
   * - LATERAL joins
     - NEAREST
   * - Window functions
     - CLUSTER
   * - Aggregate functions
     - MERGE, COUNT, etc.
   * - Common Table Expressions (WITH)
     - Complex queries, MERGE

Version Compatibility
---------------------

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Backend
     - Minimum Version
     - Notes
   * - DuckDB
     - 0.8.0+
     - Recommended: latest stable
   * - SQLite
     - 3.25.0+
     - Requires window function support
   * - PostgreSQL
     - 12+
     - Planned

Checking Compatibility
----------------------

Verify operator support at runtime:

.. code-block:: python

   from giql import GIQLEngine

   with GIQLEngine(target_dialect="duckdb") as engine:
       # Transpile a query to verify it works
       try:
           sql = engine.transpile("""
               SELECT * FROM features
               WHERE interval INTERSECTS 'chr1:1000-2000'
           """)
           print("INTERSECTS supported")
       except Exception as e:
           print(f"INTERSECTS not supported: {e}")
