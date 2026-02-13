Distance and Neighbors
======================

This section covers patterns for calculating genomic distances and finding
nearest features using GIQL's distance operators.

Calculating Distances
---------------------

Distance Between Feature Pairs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate the distance between features in two tables:

.. code-block:: sql

   SELECT
       a.name AS feature_a,
       b.name AS feature_b,
       DISTANCE(a.interval, b.interval) AS distance
   FROM features_a a
   CROSS JOIN features_b b
   WHERE a.chrom = b.chrom
   ORDER BY a.name, distance

**Use case:** Generate a distance matrix between regulatory elements and genes.

.. note::

   Always include ``WHERE a.chrom = b.chrom`` to avoid comparing
   features on different chromosomes (which returns NULL for distance).

Identify Overlapping vs Proximal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Classify relationships based on distance:

.. code-block:: sql

   SELECT
       p.name AS peak,
       g.name AS gene,
       DISTANCE(p.interval, g.interval) AS dist,
       CASE
           WHEN DISTANCE(p.interval, g.interval) = 0 THEN 'overlapping'
           WHEN DISTANCE(p.interval, g.interval) <= 1000 THEN 'proximal (<1kb)'
           WHEN DISTANCE(p.interval, g.interval) <= 10000 THEN 'nearby (<10kb)'
           ELSE 'distant'
       END AS relationship
   FROM peaks p
   CROSS JOIN genes g
   WHERE p.chrom = g.chrom

**Use case:** Categorize peak-gene relationships for enhancer analysis.

Filter by Maximum Distance
~~~~~~~~~~~~~~~~~~~~~~~~~~

Find feature pairs within a distance threshold:

.. code-block:: sql

   SELECT
       a.name,
       b.name,
       DISTANCE(a.interval, b.interval) AS dist
   FROM features_a a
   CROSS JOIN features_b b
   WHERE a.chrom = b.chrom
     AND DISTANCE(a.interval, b.interval) <= 50000
   ORDER BY dist

**Use case:** Find regulatory elements within 50kb of genes.

K-Nearest Neighbor Queries
--------------------------

Find K Nearest Features
~~~~~~~~~~~~~~~~~~~~~~~

For each peak, find the 3 nearest genes:

.. code-block:: sql

   SELECT
       peaks.name AS peak,
       nearest.name AS gene,
       nearest.distance
   FROM peaks
   CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=3) AS nearest
   ORDER BY peaks.name, nearest.distance

**Use case:** Annotate ChIP-seq peaks with nearby genes.

Nearest Feature to a Specific Location
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Find the 5 nearest genes to a specific genomic coordinate:

.. code-block:: sql

   SELECT name, distance
   FROM NEAREST(genes, reference='chr1:1000000-1001000', k=5)
   ORDER BY distance

**Use case:** Explore the genomic neighborhood of a position of interest.

Nearest with Distance Constraint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Find nearest features within a maximum distance:

.. code-block:: sql

   SELECT
       peaks.name AS peak,
       nearest.name AS gene,
       nearest.distance
   FROM peaks
   CROSS JOIN LATERAL NEAREST(
       genes,
       reference=peaks.interval,
       k=5,
       max_distance=100000
   ) AS nearest
   ORDER BY peaks.name, nearest.distance

**Use case:** Find regulatory targets within 100kb, ignoring distant genes.

Strand-Specific Queries
-----------------------

Same-Strand Nearest Neighbors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Find nearest features on the same strand only:

.. code-block:: sql

   SELECT
       peaks.name AS peak,
       nearest.name AS gene,
       nearest.strand,
       nearest.distance
   FROM peaks
   CROSS JOIN LATERAL NEAREST(
       genes,
       reference=peaks.interval,
       k=3,
       stranded=true
   ) AS nearest
   ORDER BY peaks.name, nearest.distance

**Use case:** Find same-strand genes for strand-specific regulatory analysis.

Directional Queries
-------------------

Upstream Features
~~~~~~~~~~~~~~~~~

Find features upstream (5') of reference positions using signed distances:

.. code-block:: sql

   SELECT
       peaks.name AS peak,
       nearest.name AS gene,
       nearest.distance
   FROM peaks
   CROSS JOIN LATERAL NEAREST(
       genes,
       reference=peaks.interval,
       k=10,
       signed=true
   ) AS nearest
   WHERE nearest.distance < 0
   ORDER BY peaks.name, nearest.distance DESC

**Use case:** Find genes upstream of regulatory elements.

.. note::

   With ``signed=true``, negative distances indicate upstream features
   and positive distances indicate downstream features.

Downstream Features
~~~~~~~~~~~~~~~~~~~

Find features downstream (3') of reference positions:

.. code-block:: sql

   SELECT
       peaks.name AS peak,
       nearest.name AS gene,
       nearest.distance
   FROM peaks
   CROSS JOIN LATERAL NEAREST(
       genes,
       reference=peaks.interval,
       k=10,
       signed=true
   ) AS nearest
   WHERE nearest.distance > 0
   ORDER BY peaks.name, nearest.distance

**Use case:** Identify downstream targets of promoter elements.

Promoter-Proximal Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~

Find features within a specific distance window around the reference:

.. code-block:: sql

   SELECT
       peaks.name AS peak,
       nearest.name AS gene,
       nearest.distance
   FROM peaks
   CROSS JOIN LATERAL NEAREST(
       genes,
       reference=peaks.interval,
       k=10,
       signed=true
   ) AS nearest
   WHERE nearest.distance BETWEEN -2000 AND 500
   ORDER BY peaks.name, ABS(nearest.distance)

**Use case:** Find genes with peaks in their promoter regions (-2kb to +500bp from TSS).

Combined Parameters
-------------------

Strand-Specific with Distance Constraint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Find nearby same-strand features:

.. code-block:: sql

   SELECT
       peaks.name AS peak,
       nearest.name AS gene,
       nearest.distance
   FROM peaks
   CROSS JOIN LATERAL NEAREST(
       genes,
       reference=peaks.interval,
       k=5,
       max_distance=50000,
       stranded=true,
       signed=true
   ) AS nearest
   WHERE nearest.distance BETWEEN -10000 AND 10000
   ORDER BY peaks.name, ABS(nearest.distance)

**Use case:** Find same-strand genes within ±10kb for promoter-enhancer analysis.

Distance Statistics
-------------------

Average Distance to Nearest Gene
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate the average distance from peaks to their nearest gene:

.. code-block:: sql

   WITH nearest_genes AS (
       SELECT
           peaks.name AS peak,
           nearest.distance
       FROM peaks
       CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=1) AS nearest
   )
   SELECT
       COUNT(*) AS peak_count,
       AVG(distance) AS avg_distance,
       MIN(distance) AS min_distance,
       MAX(distance) AS max_distance
   FROM nearest_genes

**Use case:** Characterize the genomic distribution of peaks relative to genes.

Distance Distribution by Chromosome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Analyze distance patterns per chromosome:

.. code-block:: sql

   WITH nearest_genes AS (
       SELECT
           peaks.chrom,
           peaks.name AS peak,
           nearest.distance
       FROM peaks
       CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=1) AS nearest
   )
   SELECT
       chrom,
       COUNT(*) AS peak_count,
       AVG(distance) AS avg_distance
   FROM nearest_genes
   GROUP BY chrom
   ORDER BY chrom

**Use case:** Compare regulatory element distribution across chromosomes.

Window Expansion Patterns
-------------------------

Expand Search Window
~~~~~~~~~~~~~~~~~~~~

Find features within an expanded window around each feature:

.. code-block:: sql

   WITH expanded AS (
       SELECT
           name,
           chrom,
           start - 5000 AS search_start,
           end + 5000 AS search_end
       FROM peaks
   )
   SELECT
       e.name AS peak,
       b.*
   FROM expanded e
   JOIN features_b b
       ON b.chrom = e.chrom
       AND b.start < e.search_end
       AND b.end > e.search_start

**Use case:** Find all features within 5kb flanking regions.

.. note::

   This pattern uses raw coordinate manipulation rather than the NEAREST
   operator, which is useful when you need custom window shapes.
