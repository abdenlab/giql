Set Quantifiers
===============

Set quantifiers extend spatial operators to work with multiple ranges simultaneously.
They allow you to test whether a genomic position matches any or all of a set of
specified ranges in a single query.

.. contents::
   :local:
   :depth: 1

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
   interval INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000', 'chr2:1000-3000')

   -- With CONTAINS
   interval CONTAINS ANY('chr1:1500', 'chr1:2500')

   -- With WITHIN
   interval WITHIN ANY('chr1:0-10000', 'chr2:0-10000')

Parameters
~~~~~~~~~~

**interval**
   A genomic column.

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

.. code-block:: sql

   SELECT * FROM variants
   WHERE interval INTERSECTS ANY(
       'chr1:1000-2000',
       'chr1:5000-6000',
       'chr2:1000-3000'
   )

**Check Against Gene Promoters:**

Find features overlapping any of a set of promoter regions:

.. code-block:: sql

   SELECT * FROM peaks
   WHERE interval INTERSECTS ANY(
       'chr1:11869-12869',   -- Gene A promoter
       'chr1:29554-30554',   -- Gene B promoter
       'chr1:69091-70091'    -- Gene C promoter
   )

**Combine with Other Filters:**

Filter by multiple regions and additional criteria:

.. code-block:: sql

   SELECT * FROM variants
   WHERE interval INTERSECTS ANY('chr1:1000-2000', 'chr2:5000-6000')
     AND quality >= 30
     AND filter = 'PASS'

**Multi-Chromosome Query:**

Query across different chromosomes efficiently:

.. code-block:: sql

   SELECT * FROM features
   WHERE interval INTERSECTS ANY(
       'chr1:100000-200000',
       'chr2:100000-200000',
       'chr3:100000-200000',
       'chrX:100000-200000'
   )

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
   interval CONTAINS ALL('chr1:1500', 'chr1:1600', 'chr1:1700')

   -- With INTERSECTS (less common, but valid)
   interval INTERSECTS ALL('chr1:1000-1100', 'chr1:1050-1150')

Parameters
~~~~~~~~~~

**interval**
   A genomic column.

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

.. code-block:: sql

   SELECT * FROM genes
   WHERE interval CONTAINS ALL(
       'chr1:1500',
       'chr1:1600',
       'chr1:1700'
   )

**Ensure Complete Coverage:**

Find intervals that span a set of required positions:

.. code-block:: sql

   SELECT * FROM features
   WHERE interval CONTAINS ALL(
       'chr1:10000',
       'chr1:15000',
       'chr1:20000'
   )

**Find Overlapping Regions:**

Find features that overlap with all specified windows (useful for finding
features in the intersection of multiple regions):

.. code-block:: sql

   SELECT * FROM features
   WHERE interval INTERSECTS ALL(
       'chr1:1000-2000',
       'chr1:1500-2500'
   )

   -- This finds features that overlap BOTH ranges
   -- (i.e., features in the intersection: chr1:1500-2000)

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


Common Patterns
---------------

**Exclusion with ANY:**

Find features that don't overlap any blacklisted region:

.. code-block:: sql

   SELECT * FROM peaks
   WHERE NOT interval INTERSECTS ANY(
       'chr1:1000000-2000000',  -- Centromere
       'chr1:5000000-5500000'   -- Known artifact region
   )

**Combining ANY and ALL:**

Complex queries can combine both quantifiers:

.. code-block:: sql

   SELECT * FROM features
   WHERE interval INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')
     AND interval CONTAINS ALL('chr1:1100', 'chr1:1200')
