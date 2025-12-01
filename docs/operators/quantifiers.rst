Set Quantifiers
===============

Set quantifiers extend spatial operators to work with multiple ranges simultaneously.
They allow you to test whether a genomic position matches any or all of a set of
specified ranges in a single query.

.. contents::
   :local:
   :depth: 2

.. _any-quantifier:

ANY
---

Match if the condition holds for any of the specified ranges.

Description
~~~~~~~~~~~

The ``ANY`` quantifier tests whether a genomic position satisfies a spatial
relationship with at least one range from a provided set. It acts as a logical
OR across multiple range comparisons.

This is useful for:

- Filtering features that overlap any of several regions of interest
- Checking membership in a set of genomic windows
- Multi-region queries without complex OR clauses

Syntax
~~~~~~

.. code-block:: sql

   -- With INTERSECTS
   position INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000', 'chr2:1000-3000')

   -- With CONTAINS
   position CONTAINS ANY('chr1:1500', 'chr1:2500')

   -- With WITHIN
   position WITHIN ANY('chr1:0-10000', 'chr2:0-10000')

Parameters
~~~~~~~~~~

**position**
   A genomic column registered with the engine.

**ranges**
   A comma-separated list of genomic range literals.

Return Value
~~~~~~~~~~~~

Boolean: ``true`` if the spatial condition holds for at least one of the specified
ranges, ``false`` otherwise.

Examples
~~~~~~~~

**Match Multiple Regions:**

Find variants in any of several regions of interest:

.. code-block:: python

   cursor = engine.execute("""
       SELECT * FROM variants
       WHERE position INTERSECTS ANY(
           'chr1:1000-2000',
           'chr1:5000-6000',
           'chr2:1000-3000'
       )
   """)

**Check Against Gene Promoters:**

Find features overlapping any of a set of promoter regions:

.. code-block:: python

   cursor = engine.execute("""
       SELECT * FROM peaks
       WHERE position INTERSECTS ANY(
           'chr1:11869-12869',   -- Gene A promoter
           'chr1:29554-30554',   -- Gene B promoter
           'chr1:69091-70091'    -- Gene C promoter
       )
   """)

**Combine with Other Filters:**

Filter by multiple regions and additional criteria:

.. code-block:: python

   cursor = engine.execute("""
       SELECT * FROM variants
       WHERE position INTERSECTS ANY('chr1:1000-2000', 'chr2:5000-6000')
         AND quality >= 30
         AND filter = 'PASS'
   """)

**Multi-Chromosome Query:**

Query across different chromosomes efficiently:

.. code-block:: python

   cursor = engine.execute("""
       SELECT * FROM features
       WHERE position INTERSECTS ANY(
           'chr1:100000-200000',
           'chr2:100000-200000',
           'chr3:100000-200000',
           'chrX:100000-200000'
       )
   """)

Backend Compatibility
~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Backend
     - Support
     - Notes
   * - DuckDB
     - Full
     -
   * - SQLite
     - Full
     -
   * - PostgreSQL
     - Planned
     -

Performance Notes
~~~~~~~~~~~~~~~~~

- ``ANY`` expands to multiple OR conditions in the generated SQL
- For very large sets of ranges, consider using a separate table and JOIN instead
- The optimizer may benefit from indexes on chromosome and position columns

Related
~~~~~~~

- :ref:`ALL <all-quantifier>` - Match all ranges (logical AND)
- :ref:`INTERSECTS <intersects-operator>` - Base spatial operator

----

.. _all-quantifier:

ALL
---

Match if the condition holds for all of the specified ranges.

Description
~~~~~~~~~~~

The ``ALL`` quantifier tests whether a genomic position satisfies a spatial
relationship with every range in a provided set. It acts as a logical AND
across multiple range comparisons.

This is useful for:

- Finding features that span multiple specific positions
- Ensuring complete coverage of a set of points
- Strict multi-point containment queries

Syntax
~~~~~~

.. code-block:: sql

   -- With CONTAINS
   position CONTAINS ALL('chr1:1500', 'chr1:1600', 'chr1:1700')

   -- With INTERSECTS (less common, but valid)
   position INTERSECTS ALL('chr1:1000-1100', 'chr1:1050-1150')

Parameters
~~~~~~~~~~

**position**
   A genomic column registered with the engine.

**ranges**
   A comma-separated list of genomic range literals.

Return Value
~~~~~~~~~~~~

Boolean: ``true`` if the spatial condition holds for all of the specified
ranges, ``false`` otherwise.

Examples
~~~~~~~~

**Find Features Containing Multiple Points:**

Find genes that contain all specified SNP positions:

.. code-block:: python

   cursor = engine.execute("""
       SELECT * FROM genes
       WHERE position CONTAINS ALL(
           'chr1:1500',
           'chr1:1600',
           'chr1:1700'
       )
   """)

**Ensure Complete Coverage:**

Find intervals that span a set of required positions:

.. code-block:: python

   cursor = engine.execute("""
       SELECT * FROM features
       WHERE position CONTAINS ALL(
           'chr1:10000',
           'chr1:15000',
           'chr1:20000'
       )
   """)

**Find Overlapping Regions:**

Find features that overlap with all specified windows (useful for finding
features in the intersection of multiple regions):

.. code-block:: python

   cursor = engine.execute("""
       SELECT * FROM features
       WHERE position INTERSECTS ALL(
           'chr1:1000-2000',
           'chr1:1500-2500'
       )
   """)

   # This finds features that overlap BOTH ranges
   # (i.e., features in the intersection: chr1:1500-2000)

Backend Compatibility
~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Backend
     - Support
     - Notes
   * - DuckDB
     - Full
     -
   * - SQLite
     - Full
     -
   * - PostgreSQL
     - Planned
     -

Performance Notes
~~~~~~~~~~~~~~~~~

- ``ALL`` expands to multiple AND conditions in the generated SQL
- Queries with ``ALL`` may be more restrictive, potentially reducing result sets
- Consider whether ``ANY`` might be more appropriate for your use case

Related
~~~~~~~

- :ref:`ANY <any-quantifier>` - Match any range (logical OR)
- :ref:`CONTAINS <contains-operator>` - Base containment operator

Choosing Between ANY and ALL
----------------------------

Use **ANY** when you want to find features that match at least one of several criteria:

.. code-block:: python

   # Find variants in gene A OR gene B OR gene C
   WHERE position INTERSECTS ANY('gene_a_region', 'gene_b_region', 'gene_c_region')

Use **ALL** when you want to find features that satisfy all criteria simultaneously:

.. code-block:: python

   # Find features that contain ALL of these positions
   WHERE position CONTAINS ALL('pos1', 'pos2', 'pos3')

Common Patterns
---------------

**Exclusion with ANY:**

Find features that don't overlap any blacklisted region:

.. code-block:: python

   cursor = engine.execute("""
       SELECT * FROM peaks
       WHERE NOT position INTERSECTS ANY(
           'chr1:1000000-2000000',  -- Centromere
           'chr1:5000000-5500000'   -- Known artifact region
       )
   """)

**Combining ANY and ALL:**

Complex queries can combine both quantifiers:

.. code-block:: python

   cursor = engine.execute("""
       SELECT * FROM features
       WHERE position INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')
         AND position CONTAINS ALL('chr1:1100', 'chr1:1200')
   """)
