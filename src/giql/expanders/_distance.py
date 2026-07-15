"""Shared distance-CASE string builder for the NEAREST expander.

The private ``_``-prefix keeps this module out of the ``expanders`` package
auto-discovery (see :mod:`giql.expanders`), so it is a plain helper module rather
than a registered expander.

.. note::

   KEEP IN SYNC: :func:`generate_distance_case` here and the AST builder in
   :mod:`giql.expanders.distance` (``expand_distance``) produce the *same*
   distance CASE by two routes. DISTANCE migrated to the AST expander (epic
   #137, issue #140); this string form survives because NEAREST assembles its
   ORDER BY / filter math as SQL fragments (string-first) and splices the CASE
   into them. Any change to the distance math here must be mirrored in the
   expander (and vice versa); the parity test in ``tests/test_distance_udf.py``
   guards drift.
"""

from __future__ import annotations


def generate_distance_case(
    chrom_a: str,
    start_a: str,
    end_a: str,
    strand_a: str | None,
    chrom_b: str,
    start_b: str,
    end_b: str,
    strand_b: str | None,
    stranded: bool = False,
    signed: bool = False,
) -> str:
    """Generate a SQL CASE expression for distance calculation.

    Distances follow bedtools ``closest -d`` semantics: overlapping
    intervals report ``0``, book-ended (adjacent) intervals where
    ``A.end == B.start`` in half-open coordinates report ``1``, and a raw
    half-open gap of N bases reports ``N + 1``. The ``+ 1`` is applied to
    the absolute gap magnitude before any directional sign, so a downstream
    book-ended pair reports ``+1`` and an upstream one ``-1`` in signed mode.

    :param chrom_a:
        Chromosome column for interval A
    :param start_a:
        Start column for interval A
    :param end_a:
        End column for interval A
    :param strand_a:
        Strand column for interval A (None if not stranded)
    :param chrom_b:
        Chromosome column for interval B
    :param start_b:
        Start column for interval B
    :param end_b:
        End column for interval B
    :param strand_b:
        Strand column for interval B (None if not stranded)
    :param stranded:
        Whether to use strand-aware distance calculation
    :param signed:
        Whether to return signed distance (negative for upstream, positive for
        downstream)
    :return:
        SQL CASE expression
    """
    if not stranded or strand_a is None or strand_b is None:
        # Basic distance calculation without strand awareness
        if signed:
            # Signed distance: negative for upstream (B before A),
            # positive for downstream (B after A)
            return (
                f"CASE WHEN {chrom_a} != {chrom_b} THEN NULL "
                f"WHEN {start_a} < {end_b} AND {end_a} > {start_b} THEN 0 "
                f"WHEN {end_a} <= {start_b} THEN ({start_b} - {end_a} + 1) "
                f"ELSE -({start_a} - {end_b} + 1) END"
            )
        # Unsigned (absolute) distance
        return (
            f"CASE WHEN {chrom_a} != {chrom_b} THEN NULL "
            f"WHEN {start_a} < {end_b} AND {end_a} > {start_b} THEN 0 "
            f"WHEN {end_a} <= {start_b} THEN ({start_b} - {end_a} + 1) "
            f"ELSE ({start_a} - {end_b} + 1) END"
        )

    # Stranded distance calculation
    # Return NULL if either strand is '.', '?', or NULL
    # Calculate distance and multiply by -1 if first interval is on '-' strand
    if signed:
        # Stranded + signed: apply strand flip AND directional sign
        return (
            f"CASE WHEN {chrom_a} != {chrom_b} THEN NULL "
            f"WHEN {strand_a} IS NULL OR {strand_b} IS NULL THEN NULL "
            f"WHEN {strand_a} = '.' OR {strand_a} = '?' THEN NULL "
            f"WHEN {strand_b} = '.' OR {strand_b} = '?' THEN NULL "
            f"WHEN {start_a} < {end_b} AND {end_a} > {start_b} THEN 0 "
            f"WHEN {end_a} <= {start_b} THEN "
            f"CASE WHEN {strand_a} = '-' THEN -({start_b} - {end_a} + 1) "
            f"ELSE ({start_b} - {end_a} + 1) END "
            f"ELSE "
            f"CASE WHEN {strand_a} = '-' THEN ({start_a} - {end_b} + 1) "
            f"ELSE -({start_a} - {end_b} + 1) END END"
        )
    # Stranded but not signed: apply strand flip only
    return (
        f"CASE WHEN {chrom_a} != {chrom_b} THEN NULL "
        f"WHEN {strand_a} IS NULL OR {strand_b} IS NULL THEN NULL "
        f"WHEN {strand_a} = '.' OR {strand_a} = '?' THEN NULL "
        f"WHEN {strand_b} = '.' OR {strand_b} = '?' THEN NULL "
        f"WHEN {start_a} < {end_b} AND {end_a} > {start_b} THEN 0 "
        f"WHEN {end_a} <= {start_b} THEN "
        f"CASE WHEN {strand_a} = '-' THEN -({start_b} - {end_a} + 1) "
        f"ELSE ({start_b} - {end_a} + 1) END "
        f"ELSE "
        f"CASE WHEN {strand_a} = '-' THEN -({start_a} - {end_b} + 1) "
        f"ELSE ({start_a} - {end_b} + 1) END END"
    )
