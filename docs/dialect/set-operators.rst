Set Operations
==============

Set operations cut a set of genomic intervals against another set of intervals,
producing finer-grained sub-intervals. Unlike the aggregation operators, which
collapse intervals into summarized regions, a set operation *multiplies* rows --
one input interval becomes one or more output rows.

.. _disjoin-operator:

DISJOIN
-------

Split genomic intervals at reference breakpoints into sub-intervals that
never partially overlap a reference interval.

Description
~~~~~~~~~~~

The ``DISJOIN`` operator is a table function: it takes a set of *target*
intervals and a *reference* set of intervals, and cuts each target interval
at every reference breakpoint (a reference ``start`` or ``end`` position)
that falls strictly inside it. Every resulting sub-interval is fully
contained by each reference interval it overlaps -- it can never *partially*
overlap one.

Unlike :ref:`MERGE <merge-operator>` and :ref:`CLUSTER <cluster-operator>`,
which aggregate intervals, ``DISJOIN`` *multiplies* rows: one target interval
becomes one or more sub-interval rows. The full target row passes through
unchanged and the sub-interval is appended as ``disjoin_chrom``,
``disjoin_start``, ``disjoin_end``.

When no ``reference`` is given it defaults to the target set, so
``DISJOIN(features)`` splits the set against its own breakpoints. Selecting
the distinct sub-intervals then yields the globally non-overlapping partition
-- the equivalent of Bioconductor's ``GenomicRanges::disjoin()``. ``DISJOIN``
deliberately departs from that set-reducing shape: as a SQL table function it
*multiplies* rows and keeps each parent row intact, so it composes with
``JOIN`` and ``SELECT``. ``SELECT DISTINCT disjoin_chrom, disjoin_start,
disjoin_end`` recovers the canonical ``disjoin()`` partition exactly.

This is useful for:

- Partitioning a set of intervals into non-overlapping segments
- Re-tiling target features onto a data-defined or external grid
- Splitting features so downstream aggregates never double-count overlaps

Syntax
~~~~~~

.. code-block:: sql

   -- Self-mode: split the set against its own breakpoints
   SELECT * FROM DISJOIN(features)

   -- Split target features against an explicit reference set
   SELECT * FROM DISJOIN(features, reference := mask)

   -- The reference may be a subquery
   SELECT * FROM DISJOIN(features, reference := (SELECT * FROM mask))

   -- The target may be a CTE assembled in the same query
   WITH x AS (SELECT chrom, "start", "end" FROM features WHERE chrom = 'chr1')
   SELECT * FROM DISJOIN(x)

Parameters
~~~~~~~~~~

**target**
   The set of intervals to split. May be a registered table (whose configured
   genomic columns are used) or a CTE defined in an enclosing ``WITH`` clause
   (assumed to expose canonical 0-based half-open ``chrom`` / ``start`` /
   ``end`` columns). A CTE shadows a registered table of the same name,
   matching SQL scoping. A bare name that is neither is rejected.

**reference** *(optional)*
   A registered table, a CTE defined in the same query, or a ``(SELECT ...)``
   subquery whose interval boundaries supply the breakpoints. Defaults to
   ``target`` when omitted. A bare name that is neither a registered table nor
   an in-query CTE is rejected.

Return Value
~~~~~~~~~~~~

Every column of the matched target row, passed through unchanged, plus the
sub-interval:

- ``disjoin_chrom`` - Chromosome of the sub-interval
- ``disjoin_start`` - Start of the sub-interval
- ``disjoin_end`` - End of the sub-interval

A sub-interval that overlaps no reference interval is dropped (the coverage
filter). In self-mode every sub-interval is covered by its own parent, so
nothing is dropped. In reference mode a target interval that overlaps no
reference interval at all yields no output rows -- ``DISJOIN`` can drop a
target entirely, not merely trim it.

.. note::

   ``disjoin_chrom`` / ``disjoin_start`` / ``disjoin_end`` are reserved output
   column names. If the target table already carries a column with one of
   those names, the output will contain two columns of that name; DuckDB
   silently renames the second, so ``SELECT disjoin_start`` then resolves to
   the passed-through parent column rather than the computed sub-interval.
   Rename the conflicting target column before the call.

.. note::

   ``DISJOIN`` reserves identifiers beginning with ``__giql_dj_`` for its
   internal CTEs. A target, reference, or enclosing CTE that starts with
   this prefix is rejected at transpile time -- rename the relation before
   the call.

Examples
~~~~~~~~

**Partition a set of intervals:**

Given two overlapping intervals ``A = [0, 20)`` and ``B = [10, 30)``,
``DISJOIN`` in self-mode cuts ``A`` at breakpoint ``10`` and ``B`` at
breakpoint ``20``:

.. code-block:: sql

   SELECT DISTINCT disjoin_chrom, disjoin_start, disjoin_end
   FROM DISJOIN(features)
   ORDER BY disjoin_start

   -- Returns the partition: [0,10), [10,20), [20,30)

**Split features against a mask:**

.. code-block:: sql

   SELECT name, disjoin_start, disjoin_end
   FROM DISJOIN(features, reference := mask)

**Re-tile against a uniform grid:**

.. code-block:: sql

   WITH bins AS (
       SELECT 'chr1' AS chrom, x AS start, x + 1000 AS "end"
       FROM range(0, 250000000, 1000) AS t(x)
   )
   SELECT * FROM DISJOIN(features, reference := bins)

.. note::

   ``range()`` is DuckDB-specific table-generating syntax; on other engines
   build the bin grid with the equivalent construct. The grid must span every
   chromosome present in ``features`` -- a feature on a chromosome the grid
   omits has no covering bin and is dropped by the coverage filter.

.. note::

   ``DISJOIN`` is a table function and appears in the ``FROM`` clause, like
   ``NEAREST``. On PostgreSQL the derived table must be given an alias
   (``FROM DISJOIN(features) AS d``).

.. note::

   The appended ``disjoin_start`` / ``disjoin_end`` columns are emitted in the
   target table's coordinate system -- the same convention as its
   passed-through ``start`` / ``end`` columns, so every column of an output row
   shares one convention. Cut positions are computed canonically inside the
   operator; only their final representation follows the target table. A CTE
   target is assumed canonical (0-based half-open), so its appended
   ``disjoin_*`` columns are emitted canonically as well.

.. note::

   ``DISJOIN`` operates on coordinates only and is strand-blind: a reference
   interval cuts a target interval regardless of either one's strand. To
   disjoin per strand, filter the target and reference to a single strand
   before the call.

.. note::

   Unlike ``DISJOIN``, ``NEAREST`` does not accept a CTE as its target -- only
   a registered table. Naming an in-query CTE in ``NEAREST(...)`` raises a
   transpilation error suggesting ``DISJOIN`` if a CTE target is desired.

See also :doc:`../recipes/disjoin` for end-to-end recipes, including the
CTE-target patterns for filtering and canonicalisation.

Related Operators
~~~~~~~~~~~~~~~~~

- :ref:`MERGE <merge-operator>` - Collapse overlapping intervals into one
- :ref:`CLUSTER <cluster-operator>` - Assign cluster IDs without splitting
