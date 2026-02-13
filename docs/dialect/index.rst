The GIQL Dialect
================

GIQL extends SQL with operators specifically designed for genomic interval queries.
These operators enable powerful spatial reasoning over genomic coordinates without
requiring complex SQL expressions.

Operators are organized by functionality. All operators work across supported 
database backends (DuckDB, SQLite, with PostgreSQL planned). Each operator page 
includes a compatibility table showing backend support status.

Logical genomic range columns
-----------------------------

GIQL allows you to reference *logical* genomic range columns in queries. Such logical
columns do not need to exist explicitly or materially in any of your data sources,
and no specialized composite data types are needed. Rather, a logical genomic range 
column or "pseudo-column" can be mapped to physical columns in your source table that 
contain the required information and use conventional data types (reference sequence name, 
start and end coordinates, optional strand, etc.).

In GIQL queries, you reference a logical genomic range column using a designated name 
like ``interval``:

.. code-block:: sql

   SELECT * FROM variants WHERE interval INTERSECTS 'chr1:1000-2000'

By providing a :doc:`schema mapping <../transpilation/schema-mapping>` for the genomic 
range columns of each of the tables in a GIQL query, the GIQL transpiler can translate 
range operations into standard SQL expressions to be consumed by a general-purpose 
query engine. Alternatively, a GIQL-aware query engine could use the schema mapping
directly for optimization.

Spatial Relationship Operators
------------------------------

Test positional relationships between genomic ranges.

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - Operator
     - Description
     - Example
   * - :ref:`INTERSECTS <intersects-operator>`
     - Returns true when ranges overlap by at least one base pair
     - ``interval INTERSECTS 'chr1:1000-2000'``
   * - :ref:`CONTAINS <contains-operator>`
     - Returns true when one range fully contains another
     - ``interval CONTAINS 'chr1:1500'``
   * - :ref:`WITHIN <within-operator>`
     - Returns true when one range is fully within another
     - ``interval WITHIN 'chr1:1000-5000'``

See :doc:`spatial-operators` for detailed documentation.

Distance and Proximity Operators
--------------------------------

Calculate distances and find nearest features.

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - Operator
     - Description
     - Example
   * - :ref:`DISTANCE <distance-operator>`
     - Calculate genomic distance between two intervals
     - ``DISTANCE(a.interval, b.interval)``
   * - :ref:`NEAREST <nearest-operator>`
     - Find k-nearest genomic features
     - ``NEAREST(genes, reference=peaks.interval, k=5)``

See :doc:`distance-operators` for detailed documentation.

Aggregation Operators
---------------------

Combine and cluster genomic intervals.

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - Operator
     - Description
     - Example
   * - :ref:`CLUSTER <cluster-operator>`
     - Assign cluster IDs to overlapping intervals
     - ``CLUSTER(interval) AS cluster_id``
   * - :ref:`MERGE <merge-operator>`
     - Combine overlapping intervals into unified regions
     - ``SELECT MERGE(interval) FROM features``

See :doc:`aggregation-operators` for detailed documentation.

Set Quantifiers
---------------

Apply operators to multiple ranges simultaneously.

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - Quantifier
     - Description
     - Example
   * - :ref:`ANY <any-quantifier>`
     - Match if condition holds for any of the specified ranges
     - ``interval INTERSECTS ANY('chr1:1000-2000', 'chr2:5000-6000')``
   * - :ref:`ALL <all-quantifier>`
     - Match if condition holds for all of the specified ranges
     - ``interval CONTAINS ALL('chr1:1500', 'chr1:1600')``

See :doc:`quantifiers` for detailed documentation.


.. .. toctree::
..    :maxdepth: 2
..    :hidden:

..    spatial-operators
..    distance-operators
..    aggregation-operators
..    quantifiers
