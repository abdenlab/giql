"""Query transformers for GIQL operations.

This module contains transformers that rewrite queries containing GIQL-specific
operations (like CLUSTER, MERGE, and binned INTERSECTS joins) into equivalent
SQL with CTEs.
"""

import itertools

from sqlglot import exp

from giql.constants import DEFAULT_BIN_SIZE
from giql.constants import DEFAULT_CHROM_COL
from giql.constants import DEFAULT_END_COL
from giql.constants import DEFAULT_START_COL
from giql.constants import DEFAULT_STRAND_COL
from giql.expressions import GIQLCluster
from giql.expressions import GIQLCoverage
from giql.expressions import GIQLMerge
from giql.expressions import Intersects
from giql.table import Tables

# Mapping from COVERAGE stat parameter to SQL aggregate function
COVERAGE_STAT_MAP = {
    "count": "COUNT",
    "mean": "AVG",
    "sum": "SUM",
    "min": "MIN",
    "max": "MAX",
}


class ClusterTransformer:
    """Transforms queries containing CLUSTER into CTE-based queries.

    CLUSTER cannot be a simple window function because it requires nested
    window functions (LAG inside SUM). Instead, we transform:

        SELECT *, CLUSTER(interval) AS cluster_id FROM features

    Into:

        WITH lag_calc AS (
            SELECT *, LAG(end_pos) OVER (...) AS prev_end FROM features
        )
        SELECT *, SUM(CASE WHEN prev_end >= start_pos ...) AS cluster_id
        FROM lag_calc
    """

    def __init__(self, tables: Tables):
        """Initialize transformer.

        :param tables:
            Table configurations for column mapping
        """
        self.tables = tables

    def _get_table_name(self, query: exp.Select) -> str | None:
        """Extract table name from query's FROM clause.

        :param query:
            Query to extract table name from
        :return:
            Table name if FROM contains a simple table, None otherwise
        """
        from_clause = query.args.get("from_")
        if not from_clause:
            return None

        if isinstance(from_clause.this, exp.Table):
            return from_clause.this.name

        return None

    def _get_genomic_columns(self, query: exp.Select) -> tuple[str, str, str, str]:
        """Get genomic column names from table config or defaults.

        :param query:
            Query to extract table and column info from
        :return:
            Tuple of (chrom_col, start_col, end_col, strand_col)
        """
        table_name = self._get_table_name(query)

        # Default column names
        chrom_col = DEFAULT_CHROM_COL
        start_col = DEFAULT_START_COL
        end_col = DEFAULT_END_COL
        strand_col = DEFAULT_STRAND_COL

        if table_name:
            table = self.tables.get(table_name)
            if table:
                chrom_col = table.chrom_col
                start_col = table.start_col
                end_col = table.end_col
                if table.strand_col:
                    strand_col = table.strand_col

        return chrom_col, start_col, end_col, strand_col

    def transform(self, query: exp.Expression) -> exp.Expression:
        """Transform query if it contains CLUSTER expressions.

        :param query:
            Parsed query AST
        :return:
            Transformed query AST
        """
        if not isinstance(query, exp.Select):
            return query

        # First, recursively transform any CTEs that might contain CLUSTER
        if query.args.get("with_"):
            cte = query.args["with_"]
            for cte_expr in cte.expressions:
                if isinstance(cte_expr, exp.CTE):
                    # Transform the CTE's subquery
                    cte_expr.set("this", self.transform(cte_expr.this))

        # Recursively transform subqueries in FROM clause
        if query.args.get("from_"):
            from_clause = query.args["from_"]
            self._transform_subqueries_in_node(from_clause)

        # Recursively transform subqueries in JOIN clauses
        if query.args.get("joins"):
            for join in query.args["joins"]:
                self._transform_subqueries_in_node(join)

        # Recursively transform subqueries in WHERE clause
        if query.args.get("where"):
            self._transform_subqueries_in_node(query.args["where"])

        # Find all CLUSTER expressions in the SELECT clause
        cluster_exprs = self._find_cluster_expressions(query)

        if not cluster_exprs:
            return query

        # Transform query for each CLUSTER expression
        for cluster_expr in cluster_exprs:
            query = self._transform_for_cluster(query, cluster_expr)

        return query

    def _transform_subqueries_in_node(self, node: exp.Expression):
        """Recursively transform subqueries within an expression node.

        :param node:
            Expression node to search for subqueries
        """
        # Find and transform any Subquery nodes
        for subquery in node.find_all(exp.Subquery):
            if isinstance(subquery.this, exp.Select):
                transformed = self.transform(subquery.this)
                subquery.set("this", transformed)

    def _find_cluster_expressions(self, query: exp.Select) -> list[GIQLCluster]:
        """Find all CLUSTER expressions in query.

        :param query:
            Query to search
        :return:
            List of CLUSTER expressions
        """
        cluster_exprs = []

        for expression in query.expressions:
            # Check if this is a CLUSTER expression or an alias containing one
            if isinstance(expression, GIQLCluster):
                cluster_exprs.append(expression)
            elif isinstance(expression, exp.Alias):
                if isinstance(expression.this, GIQLCluster):
                    cluster_exprs.append(expression.this)

        return cluster_exprs

    def _transform_for_cluster(
        self, query: exp.Select, cluster_expr: GIQLCluster
    ) -> exp.Select:
        """Transform query to compute CLUSTER using CTEs.

        :param query:
            Original query
        :param cluster_expr:
            CLUSTER expression to transform
        :return:
            Transformed query with CTEs
        """
        # Extract CLUSTER parameters
        distance_expr = cluster_expr.args.get("distance")

        # Handle distance parameter - could be int literal or None
        if distance_expr:
            if isinstance(distance_expr, exp.Literal):
                distance = int(distance_expr.this)
            else:
                # Try to extract as string and convert
                try:
                    distance = int(str(distance_expr.this))
                except (ValueError, AttributeError):
                    distance = 0
        else:
            distance = 0

        stranded_expr = cluster_expr.args.get("stranded")
        if stranded_expr:
            # Handle different types of boolean expressions
            if isinstance(stranded_expr, exp.Boolean):
                stranded = stranded_expr.this
            elif isinstance(stranded_expr, exp.Literal):
                stranded = str(stranded_expr.this).upper() == "TRUE"
            else:
                # Try to extract the value as a string
                stranded = str(stranded_expr).upper() in ("TRUE", "1", "YES")
        else:
            stranded = False

        # Get column names from table config or use defaults
        chrom_col, start_col, end_col, strand_col = self._get_genomic_columns(query)

        # Build partition clause
        partition_cols = [exp.column(chrom_col, quoted=True)]
        if stranded:
            partition_cols.append(exp.column(strand_col, quoted=True))

        # Build ORDER BY for window
        order_by = [exp.Ordered(this=exp.column(start_col, quoted=True))]

        # Create LAG window spec
        lag_window = exp.Window(
            this=exp.Anonymous(
                this="LAG", expressions=[exp.column(end_col, quoted=True)]
            ),
            partition_by=partition_cols,
            order=exp.Order(expressions=order_by),
        )

        # Add distance offset if specified
        if distance > 0:
            lag_with_distance = exp.Add(
                this=lag_window, expression=exp.Literal.number(distance)
            )
        else:
            lag_with_distance = lag_window

        # Create CASE expression for is_new_cluster
        case_expr = exp.Case(
            ifs=[
                exp.If(
                    this=exp.GTE(
                        this=lag_with_distance,
                        expression=exp.column(start_col, quoted=True),
                    ),
                    true=exp.Literal.number(0),
                )
            ],
            default=exp.Literal.number(1),
        )

        # Build CTE SELECT expressions (all original except CLUSTER, plus is_new_cluster)
        cte_expressions = []
        for expression in query.expressions:
            # Skip CLUSTER expressions
            if isinstance(expression, GIQLCluster):
                continue
            elif isinstance(expression, exp.Alias) and isinstance(
                expression.this, GIQLCluster
            ):
                continue
            else:
                cte_expressions.append(expression)

        # Ensure required columns for window functions are included
        required_cols = {chrom_col, start_col, end_col}
        if stranded:
            required_cols.add(strand_col)

        # Check if required columns are already in the select list
        selected_cols = set()
        for expr in cte_expressions:
            if isinstance(expr, exp.Column):
                selected_cols.add(expr.name)
            elif isinstance(expr, exp.Alias):
                # Don't count aliases as the source column
                pass
            elif isinstance(expr, exp.Star):
                # SELECT * includes all columns
                selected_cols = required_cols  # Assume all are covered
                break

        # Add missing required columns
        for col in required_cols - selected_cols:
            cte_expressions.append(exp.column(col, quoted=True))

        # Add is_new_cluster calculation
        cte_expressions.append(exp.alias_(case_expr, "is_new_cluster", quoted=False))

        # Build CTE query
        cte_select = exp.Select()
        cte_select.select(*cte_expressions, copy=False)

        # Copy FROM, WHERE, GROUP BY, HAVING from original (but not ORDER BY)
        # Use copy() to avoid sharing references between queries
        if query.args.get("from_"):
            from_clause = query.args["from_"].copy()
            cte_select.set("from_", from_clause)
        if query.args.get("where"):
            cte_select.set("where", query.args["where"].copy())
        if query.args.get("group"):
            cte_select.set("group", query.args["group"].copy())
        if query.args.get("having"):
            cte_select.set("having", query.args["having"].copy())

        # Create outer query with SUM over is_new_cluster
        sum_window = exp.Window(
            this=exp.Sum(this=exp.column("is_new_cluster")),
            partition_by=partition_cols,
            order=exp.Order(expressions=order_by),
        )

        # Build outer SELECT expressions (replace CLUSTER with SUM)
        new_expressions = []
        for expression in query.expressions:
            if isinstance(expression, GIQLCluster):
                new_expressions.append(sum_window)
            elif isinstance(expression, exp.Alias) and isinstance(
                expression.this, GIQLCluster
            ):
                # Keep the alias but replace the expression
                new_expressions.append(
                    exp.alias_(sum_window, expression.alias, quoted=False)
                )
            else:
                new_expressions.append(expression)

        # Build new query
        new_query = exp.Select()
        new_query.select(*new_expressions, copy=False)

        # Wrap CTE in subquery and set as FROM clause
        subquery = exp.Subquery(
            this=cte_select,
            alias=exp.TableAlias(this=exp.Identifier(this="lag_calc")),
        )
        new_query.from_(subquery, copy=False)

        # Copy ORDER BY from original to outer query
        if query.args.get("order"):
            new_query.order_by(*query.args["order"].expressions, copy=False)

        return new_query


class MergeTransformer:
    """Transforms queries containing MERGE into GROUP BY queries.

    MERGE combines overlapping intervals using CLUSTER + aggregation:

        SELECT MERGE(interval) FROM features

    Into:

        WITH clustered AS (
            SELECT *, CLUSTER(interval) AS __giql_cluster_id FROM features
        )
        SELECT
            chromosome,
            MIN(start_pos) AS start_pos,
            MAX(end_pos) AS end_pos
        FROM clustered
        GROUP BY chromosome, __giql_cluster_id
        ORDER BY chromosome, start_pos
    """

    def __init__(self, tables: Tables):
        """Initialize transformer.

        :param tables:
            Table configurations for column mapping
        """
        self.tables = tables
        self.cluster_transformer = ClusterTransformer(tables)

    def transform(self, query: exp.Expression) -> exp.Expression:
        """Transform query if it contains MERGE expressions.

        :param query:
            Parsed query AST
        :return:
            Transformed query AST
        """
        if not isinstance(query, exp.Select):
            return query

        # First, recursively transform any CTEs that might contain MERGE
        if query.args.get("with_"):
            cte = query.args["with_"]
            for cte_expr in cte.expressions:
                if isinstance(cte_expr, exp.CTE):
                    # Transform the CTE's subquery
                    cte_expr.set("this", self.transform(cte_expr.this))

        # Recursively transform subqueries in FROM clause
        if query.args.get("from_"):
            from_clause = query.args["from_"]
            self._transform_subqueries_in_node(from_clause)

        # Recursively transform subqueries in JOIN clauses
        if query.args.get("joins"):
            for join in query.args["joins"]:
                self._transform_subqueries_in_node(join)

        # Recursively transform subqueries in WHERE clause
        if query.args.get("where"):
            self._transform_subqueries_in_node(query.args["where"])

        # Find all MERGE expressions in the SELECT clause
        merge_exprs = self._find_merge_expressions(query)

        if not merge_exprs:
            return query

        # For now, support only one MERGE expression
        if len(merge_exprs) > 1:
            raise ValueError("Multiple MERGE expressions not yet supported")

        merge_expr = merge_exprs[0]
        return self._transform_for_merge(query, merge_expr)

    def _transform_subqueries_in_node(self, node: exp.Expression):
        """Recursively transform subqueries within an expression node.

        :param node:
            Expression node to search for subqueries
        """
        # Find and transform any Subquery nodes
        for subquery in node.find_all(exp.Subquery):
            if isinstance(subquery.this, exp.Select):
                transformed = self.transform(subquery.this)
                subquery.set("this", transformed)

    def _find_merge_expressions(self, query: exp.Select) -> list[GIQLMerge]:
        """Find all MERGE expressions in query.

        :param query:
            Query to search
        :return:
            List of MERGE expressions
        """
        merge_exprs = []

        for expression in query.expressions:
            if isinstance(expression, GIQLMerge):
                merge_exprs.append(expression)
            elif isinstance(expression, exp.Alias):
                if isinstance(expression.this, GIQLMerge):
                    merge_exprs.append(expression.this)

        return merge_exprs

    def _transform_for_merge(
        self, query: exp.Select, merge_expr: GIQLMerge
    ) -> exp.Select:
        """Transform query to compute MERGE using CLUSTER + GROUP BY.

        :param query:
            Original query
        :param merge_expr:
            MERGE expression to transform
        :return:
            Transformed query with clustering and aggregation
        """
        # Extract MERGE parameters (same as CLUSTER)
        distance_expr = merge_expr.args.get("distance")
        stranded_expr = merge_expr.args.get("stranded")

        # Get column names from table config or use defaults
        (
            chrom_col,
            start_col,
            end_col,
            strand_col,
        ) = self.cluster_transformer._get_genomic_columns(query)

        # Build CLUSTER expression with same parameters
        cluster_kwargs = {"this": merge_expr.this}
        if distance_expr:
            cluster_kwargs["distance"] = distance_expr
        if stranded_expr:
            cluster_kwargs["stranded"] = stranded_expr

        cluster_expr = GIQLCluster(**cluster_kwargs)

        # Create intermediate query with CLUSTER
        # Start with original query's FROM/WHERE/etc
        cluster_query = exp.Select()
        cluster_query.select(exp.Star(), copy=False)
        cluster_query.select(
            exp.alias_(cluster_expr, "__giql_cluster_id", quoted=False),
            append=True,
            copy=False,
        )

        # Copy FROM, WHERE from original
        # Use copy() to avoid sharing references between queries
        if query.args.get("from_"):
            cluster_query.set("from_", query.args["from_"].copy())
        if query.args.get("where"):
            cluster_query.set("where", query.args["where"].copy())

        # Apply CLUSTER transformation to get the CTE-based query
        cluster_query = self.cluster_transformer.transform(cluster_query)

        # Build GROUP BY columns
        group_by_cols = [exp.column(chrom_col)]

        # Handle stranded parameter
        if stranded_expr:
            if isinstance(stranded_expr, exp.Boolean):
                stranded = stranded_expr.this
            elif isinstance(stranded_expr, exp.Literal):
                stranded = str(stranded_expr.this).upper() == "TRUE"
            else:
                stranded = str(stranded_expr).upper() in ("TRUE", "1", "YES")
        else:
            stranded = False

        if stranded:
            group_by_cols.append(exp.column(strand_col, quoted=True))

        group_by_cols.append(exp.column("__giql_cluster_id"))

        # Build SELECT expressions for merged intervals
        select_exprs = []

        # Add group-by columns (non-aggregated)
        select_exprs.append(exp.column(chrom_col, quoted=True))
        if stranded:
            select_exprs.append(exp.column(strand_col, quoted=True))

        # Add merged interval bounds
        select_exprs.append(
            exp.alias_(
                exp.Min(this=exp.column(start_col, quoted=True)), start_col, quoted=False
            )
        )
        select_exprs.append(
            exp.alias_(
                exp.Max(this=exp.column(end_col, quoted=True)), end_col, quoted=False
            )
        )

        # Process other columns from original SELECT
        for expression in query.expressions:
            # Skip the MERGE expression itself
            if isinstance(expression, GIQLMerge):
                continue
            elif isinstance(expression, exp.Alias) and isinstance(
                expression.this, GIQLMerge
            ):
                continue
            # Include other columns (they should be aggregates or in GROUP BY)
            else:
                select_exprs.append(expression)

        # Build final query
        final_query = exp.Select()
        final_query.select(*select_exprs, copy=False)

        # FROM the clustered subquery
        subquery = exp.Subquery(
            this=cluster_query,
            alias=exp.TableAlias(this=exp.Identifier(this="clustered")),
        )
        final_query.from_(subquery, copy=False)

        # Add GROUP BY
        final_query.group_by(*group_by_cols, copy=False)

        # Add ORDER BY (chromosome, start)
        final_query.order_by(
            exp.Ordered(this=exp.column(chrom_col, quoted=True)), copy=False
        )
        final_query.order_by(
            exp.Ordered(this=exp.column(start_col, quoted=True)), append=True, copy=False
        )

        # Preserve any existing CTEs from the original query
        if query.args.get("with_"):
            final_query.set("with_", query.args["with_"].copy())

        return final_query


class IntersectsBinnedJoinTransformer:
    """Transforms column-to-column INTERSECTS into binned equi-joins.

    Handles both explicit JOIN ON and implicit cross-join (WHERE) patterns.
    Two rewrite strategies are selected based on the SELECT list:

    **No wildcards** (``SELECT a.chrom, b.start, ...``) — ``__giql_bin``
    cannot appear in the output regardless of CTE content, so the simpler
    1-join full-CTE approach is used:

        WITH __giql_a_binned AS (
            SELECT *, UNNEST(range(
                CAST("start" / B AS BIGINT),
                CAST(("end" - 1) / B + 1 AS BIGINT)
            )) AS __giql_bin FROM peaks
        ),
        __giql_b_binned AS (...)
        SELECT DISTINCT a.chrom, b.start, ...
        FROM __giql_a_binned AS a
        JOIN __giql_b_binned AS b
          ON a."chrom" = b."chrom" AND a.__giql_bin = b.__giql_bin
          AND a."start" < b."end" AND a."end" > b."start"

    **Wildcards present** (``SELECT a.*, b.*``) — ``__giql_bin`` would leak
    into ``a.*`` expansion if ``a`` aliases a full-select CTE. A key-only
    bridge CTE pattern is used instead, keeping original table references:

        WITH __giql_peaks_bins AS (
            SELECT "chrom", "start", "end",
                   UNNEST(range(...)) AS __giql_bin FROM peaks
        ),
        __giql_genes_bins AS (...)
        SELECT DISTINCT a.*, b.*
        FROM peaks a
        JOIN __giql_peaks_bins __giql_c0
          ON a."chrom" = __giql_c0."chrom" AND a."start" = __giql_c0."start"
          AND a."end" = __giql_c0."end"
        JOIN __giql_genes_bins __giql_c1
          ON __giql_c0."chrom" = __giql_c1."chrom"
          AND __giql_c0.__giql_bin = __giql_c1.__giql_bin
        JOIN genes b
          ON b."chrom" = __giql_c1."chrom" AND b."start" = __giql_c1."start"
          AND b."end" = __giql_c1."end"
          AND a."start" < b."end" AND a."end" > b."start"

    Literal-range INTERSECTS (e.g., ``WHERE interval INTERSECTS 'chr1:...'``)
    are left untouched.

    SELECT DISTINCT is added to deduplicate rows produced by multi-bin
    matches.  This means rows that are identical across every selected
    column will be collapsed — include a distinguishing column (e.g., an
    id or score) to preserve duplicates that differ only in unselected
    columns.  The bridge path's key-match joins on ``(chrom, start,
    end)`` and may fan out if multiple source rows share those values;
    DISTINCT corrects for this.
    """

    def __init__(self, tables: Tables, bin_size: int | None = None):
        """Initialize transformer.

        :param tables:
            Table configurations for column mapping
        :param bin_size:
            Bin width for the equi-join rewrite. Defaults to
            DEFAULT_BIN_SIZE if not specified.
        """
        self.tables = tables
        resolved = bin_size if bin_size is not None else DEFAULT_BIN_SIZE
        if not isinstance(resolved, int) or resolved <= 0:
            raise ValueError(f"bin_size must be a positive integer, got {resolved!r}")
        self.bin_size = resolved

    def transform(self, query: exp.Expression) -> exp.Expression:
        if not isinstance(query, exp.Select):
            return query

        # The pairs-CTE approach computes matching (left_key, right_key)
        # pairs via an INNER binned join with DISTINCT on key columns,
        # then joins the original tables through the pairs CTE.  This
        # avoids adding SELECT DISTINCT to the output query, which would
        # collapse legitimately different source rows that happen to
        # have identical selected columns.  It also handles outer joins
        # correctly by preventing bin fan-out from creating spurious
        # NULL rows.
        return self._transform_with_pairs(query)

    def _select_has_wildcards(self, query: exp.Select) -> bool:
        """Return True if any SELECT item is a wildcard (* or table.*)."""
        for expr in query.expressions:
            if isinstance(expr, exp.Star):
                return True
            if isinstance(expr, exp.Column) and isinstance(expr.this, exp.Star):
                return True
        return False

    def _has_outer_join_intersects(self, query: exp.Select) -> bool:
        """Return True if any outer JOIN has an INTERSECTS predicate."""
        for join in query.args.get("joins") or []:
            if join.args.get("side") and join.args.get("on"):
                if self._find_column_intersects_in(join.args["on"]):
                    return True
        return False

    def _transform_with_pairs(self, query: exp.Select) -> exp.Select:
        """Transform INTERSECTS joins using the pairs-CTE approach.

        Computes matching (left_key, right_key) pairs via an INNER
        binned join with DISTINCT on key columns, then joins the
        original tables through the pairs CTE.  Unlike the full-CTE
        and bridge paths, this does NOT add SELECT DISTINCT to the
        output — deduplication happens inside the pairs CTE on
        (chrom, start, end) keys, preserving all source rows.
        """
        joins = query.args.get("joins") or []
        key_binned: dict[str, str] = {}
        pairs_idx = 0
        new_joins: list[exp.Join] = []
        rewrote_any = False

        for join in joins:
            on = join.args.get("on")
            if on:
                intersects = self._find_column_intersects_in(on)
                if intersects:
                    extra = self._extract_non_intersects(on, intersects)
                    replacement = self._build_pairs_replacement_joins(
                        query, join, intersects, extra, key_binned, pairs_idx
                    )
                    new_joins.extend(replacement)
                    pairs_idx += 1
                    rewrote_any = True
                    continue
            new_joins.append(join)

        where = query.args.get("where")
        if where:
            intersects = self._find_column_intersects_in(where.this)
            if intersects:
                cross_join = self._find_cross_join_for_intersects(
                    query, intersects, new_joins
                )
                if cross_join is not None:
                    new_joins.remove(cross_join)
                    replacement = self._build_pairs_replacement_joins(
                        query,
                        cross_join,
                        intersects,
                        None,
                        key_binned,
                        pairs_idx,
                    )
                    new_joins.extend(replacement)
                    self._remove_intersects_from_where(query, intersects)
                    pairs_idx += 1
                    rewrote_any = True

        if rewrote_any:
            query.set("joins", new_joins)

        return query

    def _build_pairs_cte(
        self,
        name: str,
        l_cte: str,
        r_cte: str,
        l_cols: tuple[str, str, str],
        r_cols: tuple[str, str, str],
    ) -> exp.CTE:
        """Build a DISTINCT inner-join pairs CTE.

        Returns a CTE named *name* that selects the six key columns
        (__giql_l_chrom, __giql_l_start, __giql_l_end, __giql_r_chrom,
        __giql_r_start, __giql_r_end) from an INNER join of the two bin
        CTEs on chrom, __giql_bin, and the overlap predicate.
        """
        l_alias = "__giql_l"
        r_alias = "__giql_r"

        select = exp.Select()
        select.set("distinct", exp.Distinct())

        first = True
        for tbl_alias, cols, prefix in [
            (l_alias, l_cols, "__giql_l"),
            (r_alias, r_cols, "__giql_r"),
        ]:
            for col, suffix in zip(cols, ["_chrom", "_start", "_end"]):
                col_expr = exp.Alias(
                    this=exp.column(col, table=tbl_alias, quoted=True),
                    alias=exp.Identifier(this=f"{prefix}{suffix}"),
                )
                select.select(col_expr, append=not first, copy=False)
                first = False

        select.from_(
            exp.Table(
                this=exp.Identifier(this=l_cte),
                alias=exp.TableAlias(this=exp.Identifier(this=l_alias)),
            ),
            copy=False,
        )

        join_on = exp.And(
            this=exp.And(
                this=exp.EQ(
                    this=exp.column(l_cols[0], table=l_alias, quoted=True),
                    expression=exp.column(r_cols[0], table=r_alias, quoted=True),
                ),
                expression=exp.EQ(
                    this=exp.column("__giql_bin", table=l_alias),
                    expression=exp.column("__giql_bin", table=r_alias),
                ),
            ),
            expression=self._build_overlap(l_alias, r_alias, l_cols, r_cols),
        )

        select.join(
            exp.Table(
                this=exp.Identifier(this=r_cte),
                alias=exp.TableAlias(this=exp.Identifier(this=r_alias)),
            ),
            on=join_on,
            copy=False,
        )

        return exp.CTE(
            this=select,
            alias=exp.TableAlias(this=exp.Identifier(this=name)),
        )

    def _build_pairs_replacement_joins(
        self,
        query: exp.Select,
        join: exp.Join,
        intersects: Intersects,
        extra: exp.Expression | None,
        key_binned: dict[str, str],
        pairs_idx: int,
    ) -> list[exp.Join]:
        """Build a pairs CTE and two replacement joins for one INTERSECTS.

        Returns two joins:
        - join1: from_alias [SIDE] JOIN __giql_pairs ON from.key = pairs.from_key
        - join2: [SIDE] JOIN join_table ON join.key = pairs.join_key [AND extra]
        """
        from_table = query.args["from_"].this
        join_table = join.this
        if not isinstance(from_table, exp.Table) or not isinstance(
            join_table, exp.Table
        ):
            return [join]

        from_alias = from_table.alias or from_table.name
        join_alias = join_table.alias or join_table.name
        from_table_name = from_table.name
        join_table_name = join_table.name

        left_alias = intersects.this.table
        from_cols = self._get_columns(from_table_name)
        join_cols = self._get_columns(join_table_name)

        # Determine which INTERSECTS side maps to FROM vs JOIN table
        if left_alias == from_alias:
            l_table_name, r_table_name = from_table_name, join_table_name
            l_cols, r_cols = from_cols, join_cols
            from_prefix, join_prefix = "__giql_l", "__giql_r"
        else:
            l_table_name, r_table_name = join_table_name, from_table_name
            l_cols, r_cols = join_cols, from_cols
            from_prefix, join_prefix = "__giql_r", "__giql_l"

        # Ensure key-only bin CTEs exist
        l_cte = self._ensure_key_binned(query, l_table_name, key_binned)
        r_cte = self._ensure_key_binned(query, r_table_name, key_binned)

        # Build and attach the pairs CTE
        pairs_name = f"__giql_pairs_{pairs_idx}"
        pairs_cte = self._build_pairs_cte(pairs_name, l_cte, r_cte, l_cols, r_cols)
        existing_with = query.args.get("with_")
        if existing_with:
            existing_with.append("expressions", pairs_cte)
        else:
            query.set("with_", exp.With(expressions=[pairs_cte]))

        side = join.args.get("side")
        p_alias = f"__giql_p{pairs_idx}"

        # join1: [SIDE] JOIN pairs ON from.key = pairs.from_key
        join1_on = self._build_key_match(from_alias, from_cols, p_alias, from_prefix)
        join1_kwargs: dict = {
            "this": exp.Table(
                this=exp.Identifier(this=pairs_name),
                alias=exp.TableAlias(this=exp.Identifier(this=p_alias)),
            ),
            "on": join1_on,
        }
        if side:
            join1_kwargs["side"] = side
        join1 = exp.Join(**join1_kwargs)

        # join2: [SIDE] JOIN join_table ON join.key = pairs.join_key
        join2_on = self._build_key_match(join_alias, join_cols, p_alias, join_prefix)
        if extra:
            join2_on = exp.And(this=join2_on, expression=extra)
        join2_kwargs: dict = {
            "this": exp.Table(
                this=exp.Identifier(this=join_table_name),
                alias=exp.TableAlias(this=exp.Identifier(this=join_alias)),
            ),
            "on": join2_on,
        }
        if side:
            join2_kwargs["side"] = side
        join2 = exp.Join(**join2_kwargs)

        return [join1, join2]

    def _build_key_match(
        self,
        table_alias: str,
        cols: tuple[str, str, str],
        pairs_alias: str,
        prefix: str,
    ) -> exp.And:
        """Build ``table.chrom = pairs.prefix_chrom AND ...`` for all three keys."""
        return exp.And(
            this=exp.And(
                this=exp.EQ(
                    this=exp.column(cols[0], table=table_alias, quoted=True),
                    expression=exp.column(f"{prefix}_chrom", table=pairs_alias),
                ),
                expression=exp.EQ(
                    this=exp.column(cols[1], table=table_alias, quoted=True),
                    expression=exp.column(f"{prefix}_start", table=pairs_alias),
                ),
            ),
            expression=exp.EQ(
                this=exp.column(cols[2], table=table_alias, quoted=True),
                expression=exp.column(f"{prefix}_end", table=pairs_alias),
            ),
        )

    def _transform_full_cte(self, query: exp.Select) -> exp.Select:
        joins = query.args.get("joins") or []
        binned: dict[str, tuple[str, str, str]] = {}
        rewrote_any = False

        for join in joins:
            on = join.args.get("on")
            if on:
                intersects = self._find_column_intersects_in(on)
                if intersects:
                    self._rewrite_join_on_full_cte(query, join, intersects, binned)
                    rewrote_any = True

        # Implicit cross-join: FROM a, b WHERE a.interval INTERSECTS b.interval
        where = query.args.get("where")
        if where:
            intersects = self._find_column_intersects_in(where.this)
            if intersects:
                cross_join = self._find_cross_join_for_intersects(
                    query, intersects, joins
                )
                if cross_join is not None:
                    self._rewrite_cross_join_full_cte(
                        query, cross_join, intersects, binned
                    )
                    rewrote_any = True

        if rewrote_any:
            query.set("distinct", exp.Distinct())

        return query

    def _rewrite_join_on_full_cte(
        self,
        query: exp.Select,
        join: exp.Join,
        intersects: Intersects,
        binned: dict[str, tuple[str, str, str]],
    ) -> None:
        from_table = query.args["from_"].this
        join_table = join.this
        if not isinstance(from_table, exp.Table) or not isinstance(
            join_table, exp.Table
        ):
            return

        from_alias, from_cols = self._ensure_table_binned_full(
            query, from_table, query.args["from_"], binned
        )
        join_alias, join_cols = self._ensure_table_binned_full(
            query, join_table, join, binned
        )

        extra = self._extract_non_intersects(join.args.get("on"), intersects)

        equi_join = exp.And(
            this=exp.EQ(
                this=exp.column(from_cols[0], table=from_alias, quoted=True),
                expression=exp.column(join_cols[0], table=join_alias, quoted=True),
            ),
            expression=exp.EQ(
                this=exp.column("__giql_bin", table=from_alias),
                expression=exp.column("__giql_bin", table=join_alias),
            ),
        )
        # Place both equi-join and overlap in ON so LEFT/RIGHT/FULL semantics hold.
        new_on = exp.And(
            this=equi_join,
            expression=self._build_overlap(from_alias, join_alias, from_cols, join_cols),
        )
        if extra:
            new_on = exp.And(this=new_on, expression=extra)
        join.set("on", new_on)

    def _rewrite_cross_join_full_cte(
        self,
        query: exp.Select,
        cross_join: exp.Join,
        intersects: Intersects,
        binned: dict[str, tuple[str, str, str]],
    ) -> None:
        from_table = query.args["from_"].this
        join_table = cross_join.this
        if not isinstance(from_table, exp.Table) or not isinstance(
            join_table, exp.Table
        ):
            return

        from_alias, from_cols = self._ensure_table_binned_full(
            query, from_table, query.args["from_"], binned
        )
        join_alias, join_cols = self._ensure_table_binned_full(
            query, join_table, cross_join, binned
        )

        equi_join = exp.And(
            this=exp.EQ(
                this=exp.column(from_cols[0], table=from_alias, quoted=True),
                expression=exp.column(join_cols[0], table=join_alias, quoted=True),
            ),
            expression=exp.EQ(
                this=exp.column("__giql_bin", table=from_alias),
                expression=exp.column("__giql_bin", table=join_alias),
            ),
        )
        cross_join.set(
            "on",
            exp.And(
                this=equi_join,
                expression=self._build_overlap(
                    from_alias, join_alias, from_cols, join_cols
                ),
            ),
        )
        self._remove_intersects_from_where(query, intersects)

    def _ensure_table_binned_full(
        self,
        query: exp.Select,
        table: exp.Table,
        parent: exp.Expression,
        binned: dict[str, tuple[str, str, str]],
    ) -> tuple[str, tuple[str, str, str]]:
        """Create a full SELECT * CTE for *table* if needed; replace ref in *parent*."""
        alias = table.alias or table.name
        if alias in binned:
            return alias, binned[alias]

        table_name = table.name
        cols = self._get_columns(table_name)
        cte_name = f"__giql_{alias}_binned"

        cte = exp.CTE(
            this=self._build_full_binned_select(table_name, cols),
            alias=exp.TableAlias(this=exp.Identifier(this=cte_name)),
        )
        existing_with = query.args.get("with_")
        if existing_with:
            existing_with.append("expressions", cte)
        else:
            query.set("with_", exp.With(expressions=[cte]))

        parent.set(
            "this",
            exp.Table(
                this=exp.Identifier(this=cte_name),
                alias=exp.TableAlias(this=exp.Identifier(this=alias)),
            ),
        )
        binned[alias] = cols
        return alias, cols

    def _build_bin_range(
        self, start: str, end: str
    ) -> tuple[exp.Expression, exp.Expression]:
        """Build the (low, high) bin-index expressions for UNNEST(range(...)).

        Returns ``start // bin_size`` and ``(end - 1) // bin_size + 1``.
        Uses integer floor division to avoid rounding errors from
        float division + CAST.
        """
        bs = self.bin_size

        low = exp.IntDiv(
            this=exp.column(start, quoted=True),
            expression=exp.Literal.number(bs),
        )
        high = exp.Add(
            this=exp.IntDiv(
                this=exp.Paren(
                    this=exp.Sub(
                        this=exp.column(end, quoted=True),
                        expression=exp.Literal.number(1),
                    ),
                ),
                expression=exp.Literal.number(bs),
            ),
            expression=exp.Literal.number(1),
        )
        return low, high

    def _build_full_binned_select(
        self, table_name: str, cols: tuple[str, str, str]
    ) -> exp.Select:
        """Build ``SELECT *, UNNEST(range(...)) AS __giql_bin FROM <table>``."""
        _chrom, start, end = cols
        low, high = self._build_bin_range(start, end)

        range_fn = exp.Anonymous(this="range", expressions=[low, high])
        unnest_fn = exp.Anonymous(this="UNNEST", expressions=[range_fn])
        bin_alias = exp.Alias(
            this=unnest_fn,
            alias=exp.Identifier(this="__giql_bin"),
        )

        select = exp.Select()
        select.select(exp.Star(), copy=False)
        select.select(bin_alias, append=True, copy=False)
        select.from_(exp.Table(this=exp.Identifier(this=table_name)), copy=False)
        return select

    def _transform_bridge(self, query: exp.Select) -> exp.Select:
        joins = query.args.get("joins") or []
        key_binned: dict[str, str] = {}
        connector_counter = itertools.count()
        new_joins: list[exp.Join] = []
        rewrote_any = False

        for join in joins:
            on = join.args.get("on")
            if on:
                intersects = self._find_column_intersects_in(on)
                if intersects:
                    extra = self._build_join_back_joins(
                        query,
                        join,
                        intersects,
                        key_binned,
                        connector_counter,
                        preserve_kind=True,
                    )
                    new_joins.extend(extra)
                    rewrote_any = True
                    continue
            new_joins.append(join)

        where = query.args.get("where")
        if where:
            intersects = self._find_column_intersects_in(where.this)
            if intersects:
                cross_join = self._find_cross_join_for_intersects(
                    query, intersects, new_joins
                )
                if cross_join is not None:
                    new_joins.remove(cross_join)
                    extra = self._build_join_back_joins(
                        query,
                        cross_join,
                        intersects,
                        key_binned,
                        connector_counter,
                        preserve_kind=False,
                    )
                    new_joins.extend(extra)
                    self._remove_intersects_from_where(query, intersects)
                    rewrote_any = True

        if rewrote_any:
            query.set("joins", new_joins)
            query.set("distinct", exp.Distinct())

        return query

    def _find_column_intersects_in(self, expr: exp.Expression) -> Intersects | None:
        """Return the first column-to-column Intersects node in *expr*, or None.

        Only the first match is returned.  A single JOIN with multiple
        INTERSECTS conditions in its ON clause is not supported; only the
        first will be rewritten.
        """
        for node in expr.find_all(Intersects):
            if (
                isinstance(node.this, exp.Column)
                and node.this.table
                and isinstance(node.expression, exp.Column)
                and node.expression.table
            ):
                return node
        return None

    def _find_cross_join_for_intersects(
        self,
        query: exp.Select,
        intersects: Intersects,
        current_joins: list[exp.Join],
    ) -> exp.Join | None:
        """Find the implicit cross-join entry for the table in a WHERE INTERSECTS."""
        from_table = query.args["from_"].this
        if not isinstance(from_table, exp.Table):
            return None
        from_alias = from_table.alias or from_table.name

        left_alias = intersects.this.table
        right_alias = intersects.expression.table
        if left_alias == from_alias:
            target_alias = right_alias
        elif right_alias == from_alias:
            target_alias = left_alias
        else:
            return None

        for join in current_joins:
            if isinstance(join.this, exp.Table):
                alias = join.this.alias or join.this.name
                if alias == target_alias:
                    return join
        return None

    def _remove_intersects_from_where(
        self, query: exp.Select, intersects: Intersects
    ) -> None:
        """Remove the INTERSECTS predicate from the WHERE clause."""
        where = query.args.get("where")
        if not where:
            return
        remainder = self._extract_non_intersects(where.this, intersects)
        if remainder is None:
            query.set("where", None)
        else:
            query.set("where", exp.Where(this=remainder))

    def _extract_non_intersects(
        self, expr: exp.Expression | None, intersects: Intersects
    ) -> exp.Expression | None:
        """Return the parts of an AND tree that are not the INTERSECTS node."""
        if expr is None or expr is intersects:
            return None
        if isinstance(expr, exp.And):
            if expr.this is intersects:
                return expr.expression
            if expr.expression is intersects:
                return expr.this
            left = self._extract_non_intersects(expr.this, intersects)
            right = self._extract_non_intersects(expr.expression, intersects)
            if left is None:
                return right
            if right is None:
                return left
            return exp.And(this=left, expression=right)
        return expr

    def _get_columns(self, table_name: str) -> tuple[str, str, str]:
        """Return (chrom, start, end) column names for a table."""
        table = self.tables.get(table_name)
        if table:
            return (table.chrom_col, table.start_col, table.end_col)
        return (DEFAULT_CHROM_COL, DEFAULT_START_COL, DEFAULT_END_COL)

    def _build_overlap(
        self,
        from_alias: str,
        join_alias: str,
        from_cols: tuple[str, str, str],
        join_cols: tuple[str, str, str],
    ) -> exp.And:
        """Build ``from.start < join.end AND from.end > join.start``."""
        return exp.And(
            this=exp.LT(
                this=exp.column(from_cols[1], table=from_alias, quoted=True),
                expression=exp.column(join_cols[2], table=join_alias, quoted=True),
            ),
            expression=exp.GT(
                this=exp.column(from_cols[2], table=from_alias, quoted=True),
                expression=exp.column(join_cols[1], table=join_alias, quoted=True),
            ),
        )

    def _find_table_name_for_alias(self, query: exp.Select, alias: str) -> str:
        """Resolve an alias to its underlying table name."""
        from_table = query.args["from_"].this
        if isinstance(from_table, exp.Table):
            if (from_table.alias or from_table.name) == alias:
                return from_table.name
        for join in query.args.get("joins") or []:
            if isinstance(join.this, exp.Table):
                t = join.this
                if (t.alias or t.name) == alias:
                    return t.name
        return alias  # fallback: alias == table name

    def _build_key_only_bins_select(
        self, table_name: str, cols: tuple[str, str, str]
    ) -> exp.Select:
        """Build ``SELECT chrom, start, end, UNNEST(range(...)) AS __giql_bin FROM table``."""
        chrom, start, end = cols
        low, high = self._build_bin_range(start, end)

        range_fn = exp.Anonymous(this="range", expressions=[low, high])
        unnest_fn = exp.Anonymous(this="UNNEST", expressions=[range_fn])
        bin_alias = exp.Alias(
            this=unnest_fn,
            alias=exp.Identifier(this="__giql_bin"),
        )

        select = exp.Select()
        select.select(exp.column(chrom, quoted=True), copy=False)
        select.select(exp.column(start, quoted=True), append=True, copy=False)
        select.select(exp.column(end, quoted=True), append=True, copy=False)
        select.select(bin_alias, append=True, copy=False)
        select.from_(exp.Table(this=exp.Identifier(this=table_name)), copy=False)
        return select

    def _ensure_key_binned(
        self,
        query: exp.Select,
        table_name: str,
        key_binned: dict[str, str],
    ) -> str:
        """Ensure a key-only bins CTE exists for *table_name*; return its name."""
        if table_name in key_binned:
            return key_binned[table_name]

        cte_name = f"__giql_{table_name}_bins"
        cols = self._get_columns(table_name)
        cte = exp.CTE(
            this=self._build_key_only_bins_select(table_name, cols),
            alias=exp.TableAlias(this=exp.Identifier(this=cte_name)),
        )

        existing_with = query.args.get("with_")
        if existing_with:
            existing_with.append("expressions", cte)
        else:
            query.set("with_", exp.With(expressions=[cte]))

        key_binned[table_name] = cte_name
        return cte_name

    def _build_join_back_joins(
        self,
        query: exp.Select,
        join: exp.Join,
        intersects: Intersects,
        key_binned: dict[str, str],
        connector_counter: itertools.count,
        *,
        preserve_kind: bool,
    ) -> list[exp.Join]:
        """Build three replacement JOINs for one INTERSECTS using the join-back pattern.

        join1 is always INNER because it key-matches a table against its
        own bin CTE — every row has a corresponding bin entry by
        construction, so the join side has no effect.

        join2 and join3 inherit the original join's side (LEFT, RIGHT)
        when *preserve_kind* is True.
        """
        join_table = join.this
        if not isinstance(join_table, exp.Table):
            return [join]

        join_alias = join_table.alias or join_table.name
        join_table_name = join_table.name

        left_alias = intersects.this.table
        right_alias = intersects.expression.table
        other_alias = left_alias if right_alias == join_alias else right_alias
        if other_alias == join_alias:
            return [join]  # can't determine structure

        extra = self._extract_non_intersects(join.args.get("on"), intersects)

        other_table_name = self._find_table_name_for_alias(query, other_alias)
        other_cols = self._get_columns(other_table_name)
        join_cols = self._get_columns(join_table_name)

        other_cte = self._ensure_key_binned(query, other_table_name, key_binned)
        join_cte = self._ensure_key_binned(query, join_table_name, key_binned)

        c0 = f"__giql_c{next(connector_counter)}"
        c1 = f"__giql_c{next(connector_counter)}"

        join_side = None
        if preserve_kind:
            join_side = join.args.get("side")

        # join1: key-match from other_alias to its bin CTE
        join1 = exp.Join(
            this=exp.Table(
                this=exp.Identifier(this=other_cte),
                alias=exp.TableAlias(this=exp.Identifier(this=c0)),
            ),
            on=exp.And(
                this=exp.And(
                    this=exp.EQ(
                        this=exp.column(other_cols[0], table=other_alias, quoted=True),
                        expression=exp.column(other_cols[0], table=c0, quoted=True),
                    ),
                    expression=exp.EQ(
                        this=exp.column(other_cols[1], table=other_alias, quoted=True),
                        expression=exp.column(other_cols[1], table=c0, quoted=True),
                    ),
                ),
                expression=exp.EQ(
                    this=exp.column(other_cols[2], table=other_alias, quoted=True),
                    expression=exp.column(other_cols[2], table=c0, quoted=True),
                ),
            ),
        )

        # join2: bin equi-join (chrom + __giql_bin match)
        join2_kwargs: dict = {
            "this": exp.Table(
                this=exp.Identifier(this=join_cte),
                alias=exp.TableAlias(this=exp.Identifier(this=c1)),
            ),
            "on": exp.And(
                this=exp.EQ(
                    this=exp.column(other_cols[0], table=c0, quoted=True),
                    expression=exp.column(join_cols[0], table=c1, quoted=True),
                ),
                expression=exp.EQ(
                    this=exp.column("__giql_bin", table=c0),
                    expression=exp.column("__giql_bin", table=c1),
                ),
            ),
        }
        if join_side:
            join2_kwargs["side"] = join_side
        join2 = exp.Join(**join2_kwargs)

        # join3: key-match from join CTE back to actual join table + overlap
        key_match = exp.And(
            this=exp.And(
                this=exp.EQ(
                    this=exp.column(join_cols[0], table=join_alias, quoted=True),
                    expression=exp.column(join_cols[0], table=c1, quoted=True),
                ),
                expression=exp.EQ(
                    this=exp.column(join_cols[1], table=join_alias, quoted=True),
                    expression=exp.column(join_cols[1], table=c1, quoted=True),
                ),
            ),
            expression=exp.EQ(
                this=exp.column(join_cols[2], table=join_alias, quoted=True),
                expression=exp.column(join_cols[2], table=c1, quoted=True),
            ),
        )
        join3_on = exp.And(
            this=key_match,
            expression=self._build_overlap(
                other_alias, join_alias, other_cols, join_cols
            ),
        )
        if extra:
            join3_on = exp.And(this=join3_on, expression=extra)
        join3_kwargs: dict = {
            "this": exp.Table(
                this=exp.Identifier(this=join_table_name),
                alias=exp.TableAlias(this=exp.Identifier(this=join_alias)),
            ),
            "on": join3_on,
        }
        if join_side:
            join3_kwargs["side"] = join_side

        join3 = exp.Join(**join3_kwargs)

        return [join1, join2, join3]


class CoverageTransformer:
    """Transforms queries containing COVERAGE into binned coverage queries.

    COVERAGE tiles the genome into fixed-width bins and aggregates overlapping
    intervals per bin:

        SELECT COVERAGE(interval, 1000) FROM features

    Into:

        WITH __giql_bins AS (
            SELECT chrom, bin_start AS start, bin_start + 1000 AS "end"
            FROM (
                SELECT DISTINCT chrom, MAX("end") AS __max_end
                FROM features GROUP BY chrom
            ) AS __giql_chroms,
            LATERAL generate_series(0, __max_end, 1000) AS t(bin_start)
        )
        SELECT bins.chrom, bins.start, bins."end", COUNT(source.*)
        FROM __giql_bins AS bins
        LEFT JOIN features AS source
          ON source.start < bins."end"
          AND source."end" > bins.start
          AND source.chrom = bins.chrom
        GROUP BY bins.chrom, bins.start, bins."end"
        ORDER BY bins.chrom, bins.start
    """

    def __init__(self, tables: Tables):
        """Initialize transformer.

        :param tables:
            Table configurations for column mapping
        """
        self.tables = tables
        self.cluster_transformer = ClusterTransformer(tables)

    def _get_table_alias(self, query: exp.Select) -> str | None:
        """Extract table alias from query's FROM clause.

        :param query:
            Query to extract alias from
        :return:
            Table alias if present, None otherwise
        """
        from_clause = query.args.get("from_")
        if not from_clause:
            return None
        if isinstance(from_clause.this, exp.Table):
            return from_clause.this.alias
        return None

    def transform(self, query: exp.Expression) -> exp.Expression:
        """Transform query if it contains COVERAGE expressions.

        :param query:
            Parsed query AST
        :return:
            Transformed query AST
        """
        if not isinstance(query, exp.Select):
            return query

        # Recursively transform CTEs
        if query.args.get("with_"):
            cte = query.args["with_"]
            for cte_expr in cte.expressions:
                if isinstance(cte_expr, exp.CTE):
                    cte_expr.set("this", self.transform(cte_expr.this))

        # Recursively transform subqueries in FROM/JOIN/WHERE
        for key in ("from_", "where"):
            if query.args.get(key):
                self._transform_subqueries_in_node(query.args[key])
        if query.args.get("joins"):
            for join in query.args["joins"]:
                self._transform_subqueries_in_node(join)

        # Find COVERAGE expressions in SELECT
        coverage_exprs = self._find_coverage_expressions(query)
        if not coverage_exprs:
            return query

        if len(coverage_exprs) > 1:
            raise ValueError("Multiple COVERAGE expressions not yet supported")

        return self._transform_for_coverage(query, coverage_exprs[0])

    def _transform_subqueries_in_node(self, node: exp.Expression):
        """Recursively transform subqueries within an expression node.

        :param node:
            Expression node to search for subqueries
        """
        for subquery in node.find_all(exp.Subquery):
            if isinstance(subquery.this, exp.Select):
                transformed = self.transform(subquery.this)
                subquery.set("this", transformed)

    def _find_coverage_expressions(self, query: exp.Select) -> list[GIQLCoverage]:
        """Find all COVERAGE expressions in query.

        :param query:
            Query to search
        :return:
            List of COVERAGE expressions
        """
        coverage_exprs = []
        for expression in query.expressions:
            if isinstance(expression, GIQLCoverage):
                coverage_exprs.append(expression)
            elif isinstance(expression, exp.Alias):
                if isinstance(expression.this, GIQLCoverage):
                    coverage_exprs.append(expression.this)
        return coverage_exprs

    def _transform_for_coverage(
        self, query: exp.Select, coverage_expr: GIQLCoverage
    ) -> exp.Select:
        """Transform query to compute COVERAGE using bins CTE + JOIN + GROUP BY.

        :param query:
            Original query
        :param coverage_expr:
            COVERAGE expression to transform
        :return:
            Transformed query
        """
        # Extract parameters
        resolution_expr = coverage_expr.args.get("resolution")
        if isinstance(resolution_expr, exp.Literal):
            resolution = int(resolution_expr.this)
        elif (
            isinstance(resolution_expr, exp.Neg)
            and isinstance(resolution_expr.this, exp.Literal)
        ):
            resolution = -int(resolution_expr.this.this)
        else:
            raise ValueError("COVERAGE resolution must be an integer literal")

        if resolution <= 0:
            raise ValueError(
                f"COVERAGE resolution must be positive, got {resolution}"
            )

        stat_expr = coverage_expr.args.get("stat")
        if stat_expr:
            if isinstance(stat_expr, exp.Literal):
                stat = stat_expr.this.strip("'\"").lower()
            else:
                stat = str(stat_expr).strip("'\"").lower()
        else:
            stat = "count"

        if stat not in COVERAGE_STAT_MAP:
            raise ValueError(
                f"Unknown COVERAGE stat '{stat}'. "
                f"Must be one of: {', '.join(COVERAGE_STAT_MAP)}"
            )

        sql_agg = COVERAGE_STAT_MAP[stat]

        # Extract target parameter
        target_expr = coverage_expr.args.get("target")
        if target_expr:
            if isinstance(target_expr, exp.Literal):
                target_col = target_expr.this.strip("'\"")
            else:
                target_col = str(target_expr).strip("'\"")
        else:
            target_col = None

        # Get column names and table info
        chrom_col, start_col, end_col, _ = (
            self.cluster_transformer._get_genomic_columns(query)
        )
        table_name = self.cluster_transformer._get_table_name(query)
        table_alias = self._get_table_alias(query)
        source_ref = table_alias or table_name or "source"

        # Build __giql_chroms subquery:
        #   SELECT DISTINCT chrom, MAX("end") AS __max_end FROM <table> GROUP BY chrom
        chroms_select = exp.Select()
        chroms_select.select(
            exp.column(chrom_col, quoted=True),
            copy=False,
        )
        chroms_select.select(
            exp.alias_(
                exp.Max(this=exp.column(end_col, quoted=True)),
                "__max_end",
                quoted=False,
            ),
            append=True,
            copy=False,
        )

        if table_name:
            if table_alias:
                chroms_select.from_(
                    exp.alias_(exp.to_table(table_name), table_alias, table=True),
                    copy=False,
                )
            else:
                chroms_select.from_(exp.to_table(table_name), copy=False)

        # Apply WHERE from original query to the chroms subquery too,
        # qualifying unqualified column references with the table name
        if query.args.get("where"):
            chroms_where = query.args["where"].copy()
            if table_name:
                for col in chroms_where.find_all(exp.Column):
                    if not col.table:
                        col.set("table", exp.Identifier(this=table_name))
            chroms_select.set("where", chroms_where)

        chroms_select.group_by(exp.column(chrom_col, quoted=True), copy=False)

        chroms_subquery = exp.Subquery(
            this=chroms_select,
            alias=exp.TableAlias(this=exp.Identifier(this="__giql_chroms")),
        )

        # Build bins CTE using raw SQL for generate_series + LATERAL
        # since SQLGlot doesn't natively support generate_series
        bins_select = exp.Select()
        bins_select.select(
            exp.column(chrom_col, table="__giql_chroms", quoted=True),
            copy=False,
        )
        bins_select.select(
            exp.alias_(
                exp.column("bin_start"),
                start_col,
                quoted=True,
            ),
            append=True,
            copy=False,
        )
        bins_select.select(
            exp.alias_(
                exp.Add(
                    this=exp.column("bin_start"),
                    expression=exp.Literal.number(resolution),
                ),
                end_col,
                quoted=True,
            ),
            append=True,
            copy=False,
        )

        # FROM __giql_chroms subquery
        bins_select.from_(chroms_subquery, copy=False)

        # CROSS JOIN LATERAL generate_series(0, __max_end, resolution) AS t(bin_start)
        lateral_join = exp.Join(
            this=exp.Lateral(
                this=exp.Anonymous(
                    this="generate_series",
                    expressions=[
                        exp.Literal.number(0),
                        exp.column("__max_end"),
                        exp.Literal.number(resolution),
                    ],
                ),
                alias=exp.TableAlias(
                    this=exp.Identifier(this="t"),
                    columns=[exp.Identifier(this="bin_start")],
                ),
            ),
            kind="CROSS",
        )
        bins_select.append("joins", lateral_join)

        # Wrap bins_select as a CTE named __giql_bins
        bins_cte = exp.CTE(
            this=bins_select,
            alias=exp.TableAlias(this=exp.Identifier(this="__giql_bins")),
        )
        with_clause = exp.With(expressions=[bins_cte])

        # Build the aggregate expression
        if stat == "count":
            if target_col:
                agg_expr = exp.Anonymous(
                    this="COUNT",
                    expressions=[
                        exp.column(target_col, table=source_ref, quoted=True),
                    ],
                )
            else:
                agg_expr = exp.Anonymous(
                    this="COUNT",
                    expressions=[
                        exp.column(chrom_col, table=source_ref, quoted=True),
                    ],
                )
        else:
            if target_col:
                agg_expr = exp.Anonymous(
                    this=sql_agg,
                    expressions=[
                        exp.column(target_col, table=source_ref, quoted=True),
                    ],
                )
            else:
                agg_expr = exp.Anonymous(
                    this=sql_agg,
                    expressions=[
                        exp.Sub(
                            this=exp.column(end_col, table=source_ref, quoted=True),
                            expression=exp.column(
                                start_col, table=source_ref, quoted=True
                            ),
                        )
                    ],
                )

        # Build main SELECT
        final_query = exp.Select()

        # Add bin coordinate columns
        final_query.select(
            exp.column(chrom_col, table="bins", quoted=True),
            copy=False,
        )
        final_query.select(
            exp.column(start_col, table="bins", quoted=True),
            append=True,
            copy=False,
        )
        final_query.select(
            exp.column(end_col, table="bins", quoted=True),
            append=True,
            copy=False,
        )

        # Replace COVERAGE(...) in select list with aggregate, and add other columns
        for expression in query.expressions:
            if isinstance(expression, GIQLCoverage):
                final_query.select(
                    exp.alias_(agg_expr, "value", quoted=False),
                    append=True,
                    copy=False,
                )
            elif isinstance(expression, exp.Alias) and isinstance(
                expression.this, GIQLCoverage
            ):
                final_query.select(
                    exp.alias_(agg_expr, expression.alias, quoted=False),
                    append=True,
                    copy=False,
                )
            else:
                final_query.select(expression, append=True, copy=False)

        # FROM __giql_bins AS bins
        final_query.from_(
            exp.Table(
                this=exp.Identifier(this="__giql_bins"),
                alias=exp.TableAlias(this=exp.Identifier(this="bins")),
            ),
            copy=False,
        )

        # LEFT JOIN source ON overlap conditions
        source_table = exp.to_table(table_name) if table_name else exp.to_table("source")
        source_table.set("alias", exp.TableAlias(this=exp.Identifier(this=source_ref)))

        join_condition = exp.And(
            this=exp.And(
                this=exp.LT(
                    this=exp.column(start_col, table=source_ref, quoted=True),
                    expression=exp.column(end_col, table="bins", quoted=True),
                ),
                expression=exp.GT(
                    this=exp.column(end_col, table=source_ref, quoted=True),
                    expression=exp.column(start_col, table="bins", quoted=True),
                ),
            ),
            expression=exp.EQ(
                this=exp.column(chrom_col, table=source_ref, quoted=True),
                expression=exp.column(chrom_col, table="bins", quoted=True),
            ),
        )

        # Merge original WHERE into the JOIN ON condition so that
        # LEFT JOIN still produces zero-coverage bins (WHERE would filter
        # them out because source columns are NULL for non-matching bins)
        if query.args.get("where"):
            where_condition = query.args["where"].this.copy()
            # Qualify unqualified column references with source_ref
            for col in where_condition.find_all(exp.Column):
                if not col.table:
                    col.set("table", exp.Identifier(this=source_ref))
            join_condition = exp.And(
                this=join_condition,
                expression=where_condition,
            )

        left_join = exp.Join(
            this=source_table,
            on=join_condition,
            kind="LEFT",
        )
        final_query.append("joins", left_join)

        # GROUP BY bins.chrom, bins.start, bins.end
        final_query.group_by(
            exp.column(chrom_col, table="bins", quoted=True),
            copy=False,
        )
        final_query.group_by(
            exp.column(start_col, table="bins", quoted=True),
            append=True,
            copy=False,
        )
        final_query.group_by(
            exp.column(end_col, table="bins", quoted=True),
            append=True,
            copy=False,
        )

        # ORDER BY bins.chrom, bins.start
        final_query.order_by(
            exp.Ordered(this=exp.column(chrom_col, table="bins", quoted=True)),
            copy=False,
        )
        final_query.order_by(
            exp.Ordered(this=exp.column(start_col, table="bins", quoted=True)),
            append=True,
            copy=False,
        )

        # Attach the WITH clause, preserving any user CTEs from the input query
        existing_with = query.args.get("with_")
        if existing_with:
            merged_ctes = [cte.copy() for cte in existing_with.expressions] + [bins_cte]
            final_query.set("with_", exp.With(expressions=merged_ctes))
        else:
            final_query.set("with_", with_clause)

        return final_query
