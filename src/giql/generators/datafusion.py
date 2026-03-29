"""DataFusion SQL generator for GIQL transpilation.

Emits ``giql_intersects()`` function calls for column-to-column
INTERSECTS joins instead of expanding to raw overlap predicates.
A DataFusion logical optimizer rule matches on that function call
and rewrites it to a binned equi-join with adaptive bin sizing.
"""

from __future__ import annotations

from giql.generators.base import BaseGIQLGenerator


class DataFusionGIQLGenerator(BaseGIQLGenerator):
    """Generator that preserves INTERSECTS semantics for DataFusion.

    For column-to-column INTERSECTS joins, emits::

        (l.chrom = r.chrom AND giql_intersects(l.start, l.end, r.start, r.end))

    instead of the standard overlap predicates. The chrom equi-key is
    preserved as plain SQL so DataFusion can use it for hash
    partitioning. All other operations (literal range queries,
    CONTAINS, WITHIN) fall through to the base generator.
    """

    def _generate_column_join(self, left_col: str, right_col: str, op_type: str) -> str:
        if op_type == "intersects":
            l_chrom, l_start, l_end = self._get_column_refs(left_col, None)
            r_chrom, r_start, r_end = self._get_column_refs(right_col, None)
            return (
                f"({l_chrom} = {r_chrom} "
                f"AND giql_intersects({l_start}, {l_end}, {r_start}, {r_end}))"
            )
        return super()._generate_column_join(left_col, right_col, op_type)
