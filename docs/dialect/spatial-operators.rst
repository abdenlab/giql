Spatial Operators
=================

Spatial relationship operators test positional relationships between genomic ranges.
These are the core operators for determining whether genomic intervals overlap,
contain, or are contained within other intervals.

.. _intersects-operator:

INTERSECTS
----------

Returns true when two genomic ranges overlap by at least one base pair.

Description
~~~~~~~~~~~

The ``INTERSECTS`` operator is the most commonly used spatial operator. It tests
whether two genomic intervals share any overlapping bases. Two intervals intersect
if they are on the same chromosome and their coordinate ranges overlap.

Mathematically, intervals ``[start_a, end_a)`` and ``[start_b, end_b)`` intersect when:

- They are on the same chromosome, AND
- ``start_a < end_b`` AND ``start_b < end_a``

Syntax
~~~~~~

.. code-block:: sql

   -- Compare against a literal range
   interval INTERSECTS 'chr1:1000-2000'

   -- Compare against another genomic column (joins)
   a.interval INTERSECTS b.interval

   -- With set quantifiers
   interval INTERSECTS ANY('chr1:1000-2000', 'chr2:5000-6000')

Parameters
~~~~~~~~~~

**interval**
   A genomic column from a registered table.

**literal_range**
   A string literal specifying a genomic range in the format ``'chromosome:start-end'``.

**other_interval**
   Another genomic column from the same or different table (for joins).

Return Value
~~~~~~~~~~~~

Boolean: ``true`` if the ranges overlap, ``false`` otherwise.

Examples
~~~~~~~~

**Basic Usage:**

Find all variants that overlap a specific genomic region:

.. code-block:: sql

   SELECT * FROM variants
   WHERE interval INTERSECTS 'chr1:1000-2000'

**Column-to-Column Joins:**

Find variants that overlap with any gene:

.. code-block:: sql

   SELECT v.*, g.name AS gene_name
   FROM variants v
   INNER JOIN genes g ON v.interval INTERSECTS g.interval

**With WHERE Clause:**

Find overlapping features with additional filtering:

.. code-block:: sql

   SELECT v.*, g.name
   FROM variants v
   INNER JOIN genes g ON v.interval INTERSECTS g.interval
   WHERE v.quality >= 30
     AND g.biotype = 'protein_coding'

**Left Outer Join:**

Find all variants, with gene information where available:

.. code-block:: sql

   SELECT v.*, g.name AS gene_name
   FROM variants v
   LEFT JOIN genes g ON v.interval INTERSECTS g.interval

Related Operators
~~~~~~~~~~~~~~~~~

- :ref:`CONTAINS <contains-operator>` - Check if one range fully contains another
- :ref:`WITHIN <within-operator>` - Check if one range is fully within another
- :ref:`DISTANCE <distance-operator>` - Calculate distance between non-overlapping ranges

----

.. _contains-operator:

CONTAINS
--------

Returns true when one genomic range fully contains another.

Description
~~~~~~~~~~~

The ``CONTAINS`` operator tests whether one genomic interval completely encompasses
another. The containing interval must cover the entire span of the contained interval.

Mathematically, interval ``[start_a, end_a)`` contains ``[start_b, end_b)`` when:

- They are on the same chromosome, AND
- ``start_a <= start_b`` AND ``end_a >= end_b``

Syntax
~~~~~~

.. code-block:: sql

   -- Check if interval contains a point
   interval CONTAINS 'chr1:1500'

   -- Check if interval contains a range
   interval CONTAINS 'chr1:1200-1800'

   -- Column-to-column comparison
   gene.interval CONTAINS variant.interval

Parameters
~~~~~~~~~~

**interval**
   A genomic column.

**literal_range**
   A string literal specifying a genomic point or range.

**other_interval**
   Another genomic column for comparisons.

Return Value
~~~~~~~~~~~~

Boolean: ``true`` if the first range fully contains the second, ``false`` otherwise.

Examples
~~~~~~~~

**Point Containment:**

Find genes that contain a specific position:

.. code-block:: sql

   SELECT * FROM genes
   WHERE interval CONTAINS 'chr1:1500'

**Range Containment:**

Find large features that fully contain smaller features:

.. code-block:: sql

   SELECT g.name AS gene_name, e.name AS exon_name
   FROM genes g
   INNER JOIN exons e ON g.interval CONTAINS e.interval

**Filtering Fully Contained Variants:**

Find variants that are completely within gene boundaries:

.. code-block:: sql

   SELECT v.*
   FROM variants v
   INNER JOIN genes g ON g.interval CONTAINS v.interval

Related Operators
~~~~~~~~~~~~~~~~~

- :ref:`WITHIN <within-operator>` - Inverse of CONTAINS
- :ref:`INTERSECTS <intersects-operator>` - Partial overlap (less strict)

----

.. _within-operator:

WITHIN
------

Returns true when one genomic range is fully contained within another.

Description
~~~~~~~~~~~

The ``WITHIN`` operator is the inverse of ``CONTAINS``. It tests whether a genomic
interval falls completely inside another interval.

Mathematically, interval ``[start_a, end_a)`` is within ``[start_b, end_b)`` when:

- They are on the same chromosome, AND
- ``start_a >= start_b`` AND ``end_a <= end_b``

Syntax
~~~~~~

.. code-block:: sql

   -- Check if interval is within a range
   interval WITHIN 'chr1:1000-5000'

   -- Column-to-column comparison
   variant.interval WITHIN gene.interval

Parameters
~~~~~~~~~~

**interval**
   A genomic column.

**literal_range**
   A string literal specifying the containing range.

**other_interval**
   Another genomic column for comparisons.

Return Value
~~~~~~~~~~~~

Boolean: ``true`` if the first range is fully within the second, ``false`` otherwise.

Examples
~~~~~~~~

**Filter to Region:**

Find all features within a specific genomic window:

.. code-block:: sql

   SELECT * FROM features
   WHERE interval WITHIN 'chr1:1000000-2000000'

**Find Nested Features:**

Find exons that are completely within their parent gene:

.. code-block:: sql

   SELECT e.*, g.name AS gene_name
   FROM exons e
   INNER JOIN genes g ON e.interval WITHIN g.interval

Related Operators
~~~~~~~~~~~~~~~~~

- :ref:`CONTAINS <contains-operator>` - Inverse of WITHIN
- :ref:`INTERSECTS <intersects-operator>` - Partial overlap (less strict)
