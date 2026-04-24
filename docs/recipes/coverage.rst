Coverage
========

This section covers patterns for computing genome-wide coverage and signal
summaries using GIQL's ``COVERAGE`` operator.

Basic Coverage
--------------

Binned coverage underpins most genome-wide signal summaries — read-pileup plots for ChIP-seq, exon-level depth in RNA-seq, and peak-density overviews across megabases. The recipes below start from a canonical interval-count and build toward more specialised summaries.

Count Overlapping Features
~~~~~~~~~~~~~~~~~~~~~~~~~~

Count the number of features overlapping each 1 kb bin across the genome:

.. code-block:: sql

   SELECT COVERAGE(interval, 1000) AS depth
   FROM features

**Sample output:**

.. code-block:: text

   ┌────────┬────────┬────────┬───────┐
   │ chrom  │ start  │  end   │ depth │
   ├────────┼────────┼────────┼───────┤
   │ chr1   │      0 │   1000 │     3 │
   │ chr1   │   1000 │   2000 │     1 │
   │ chr1   │   2000 │   3000 │     0 │
   │ ...    │    ... │    ... │   ... │
   └────────┴────────┴────────┴───────┘

Each row represents one genomic bin. Bins with no overlapping features appear with a count of zero. An interval that spans more than one bin is counted in each bin it overlaps (the ``bedtools coverage`` convention), so the sum of bin counts is generally greater than the number of source intervals.

**Use case:** Compute read depth or feature density at a fixed resolution.

Custom Bin Size
~~~~~~~~~~~~~~~

Use a finer resolution of 100 bp:

.. code-block:: sql

   SELECT COVERAGE(interval, 100) AS depth
   FROM reads

**Use case:** High-resolution coverage tracks for visualisation.

Named Resolution Parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~

The resolution can also be supplied by name:

.. code-block:: sql

   SELECT COVERAGE(interval, resolution := 500) AS depth
   FROM features

Both ``:=`` and ``=>`` are accepted for named parameters.

.. note::

   Weighted summary statistics (mean, sum, min, max over interval values, with bin-boundary-aware weighting) are not yet implemented. See the project tracker for the follow-up.

Filtered Coverage
-----------------

Strand-Specific Coverage
~~~~~~~~~~~~~~~~~~~~~~~~

Compute coverage for each strand separately by filtering:

.. code-block:: sql

   -- Plus strand
   SELECT COVERAGE(interval, 1000) AS depth
   FROM features
   WHERE strand = '+'

.. code-block:: sql

   -- Minus strand
   SELECT COVERAGE(interval, 1000) AS depth
   FROM features
   WHERE strand = '-'

**Use case:** Strand-specific signal tracks for RNA-seq or stranded assays.

Coverage of High-Scoring Features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Restrict coverage to features above a quality threshold:

.. code-block:: sql

   SELECT COVERAGE(interval, 1000) AS depth
   FROM features
   WHERE score > 10

5' End Counting
~~~~~~~~~~~~~~~

To count only the 5' ends of features (e.g. TSS or read starts), first
create a view or CTE that trims each interval to its 5' end, then apply
``COVERAGE``:

.. code-block:: sql

   WITH five_prime AS (
       SELECT chrom, "start", "start" + 1 AS "end"
       FROM features
       WHERE strand = '+'
       UNION ALL
       SELECT chrom, "end" - 1 AS "start", "end"
       FROM features
       WHERE strand = '-'
   )
   SELECT COVERAGE(interval, 1000) AS tss_count
   FROM five_prime

Normalised Coverage
-------------------

RPM Normalisation
~~~~~~~~~~~~~~~~~

Normalise bin counts to reads per million (RPM) by dividing by the total
number of reads:

.. code-block:: sql

   WITH bins AS (
       SELECT COVERAGE(interval, 1000) AS depth
       FROM reads
   ),
   total AS (
       SELECT COUNT(*) AS n FROM reads
   )
   SELECT
       bins.chrom,
       bins.start,
       bins.end,
       bins.depth * 1000000.0 / total.n AS rpm
   FROM bins, total
