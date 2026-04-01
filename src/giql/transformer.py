"""Query transformers for GIQL operations.

This module contains transformers that rewrite queries containing GIQL-specific
operations (like CLUSTER, MERGE, and binned INTERSECTS joins) into equivalent
SQL with CTEs.
"""

from sqlglot import exp

from giql.constants import DEFAULT_CHROM_COL
from giql.constants import DEFAULT_END_COL
from giql.constants import DEFAULT_START_COL
from giql.constants import DEFAULT_STRAND_COL
from giql.expressions import GIQLCluster
from giql.expressions import GIQLMerge
from giql.expressions import Intersects
from giql.table import Tables

DEFAULT_BIN_SIZE = 10_000


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
    """

    def __init__(self, tables: Tables, bin_size: int | None = None):
        self.tables = tables
        self.bin_size = bin_size if bin_size is not None else DEFAULT_BIN_SIZE
        if self.bin_size <= 0:
            raise ValueError(f"bin_size must be a positive integer, got {self.bin_size}")

    def transform(self, query: exp.Expression) -> exp.Expression:
        if not isinstance(query, exp.Select):
            return query

        # The bridge path can't faithfully represent FULL OUTER JOIN
        # because the three-join chain's bin fan-out creates spurious
        # unmatched rows.  Fall back to full-CTE for those queries.
        if self._select_has_wildcards(query) and not self._has_full_outer_join(query):
            return self._transform_bridge(query)
        return self._transform_full_cte(query)

    def _select_has_wildcards(self, query: exp.Select) -> bool:
        """Return True if any SELECT item is a wildcard (* or table.*)."""
        for expr in query.expressions:
            if isinstance(expr, exp.Star):
                return True
            if isinstance(expr, exp.Column) and isinstance(expr.this, exp.Star):
                return True
        return False

    def _has_full_outer_join(self, query: exp.Select) -> bool:
        """Return True if any JOIN in the query is a FULL OUTER JOIN."""
        for join in query.args.get("joins") or []:
            if join.args.get("side") == "FULL":
                return True
        return False

    def _transform_full_cte(self, query: exp.Select) -> exp.Select:
        joins = query.args.get("joins") or []
        # binned: alias -> (chrom_col, start_col, end_col)
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

    def _build_full_binned_select(
        self, table_name: str, cols: tuple[str, str, str]
    ) -> exp.Select:
        """Build ``SELECT *, UNNEST(range(...)) AS __giql_bin FROM <table>``."""
        _chrom, start, end = cols
        B = self.bin_size

        low = exp.Cast(
            this=exp.Div(
                this=exp.column(start, quoted=True),
                expression=exp.Literal.number(B),
            ),
            to=exp.DataType(this=exp.DataType.Type.BIGINT),
        )
        high = exp.Cast(
            this=exp.Add(
                this=exp.Div(
                    this=exp.Paren(
                        this=exp.Sub(
                            this=exp.column(end, quoted=True),
                            expression=exp.Literal.number(1),
                        ),
                    ),
                    expression=exp.Literal.number(B),
                ),
                expression=exp.Literal.number(1),
            ),
            to=exp.DataType(this=exp.DataType.Type.BIGINT),
        )

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
        connector_idx = [0]
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
                        connector_idx,
                        preserve_kind=True,
                    )
                    new_joins.extend(extra)
                    rewrote_any = True
                    continue
            new_joins.append(join)

        # Implicit cross-join: FROM a, b WHERE a.interval INTERSECTS b.interval
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
                        connector_idx,
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
        """Find an Intersects node where both sides are table-qualified columns."""
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
        where_expr = where.this
        if where_expr is intersects:
            query.set("where", None)
        elif isinstance(where_expr, exp.And):
            if where_expr.this is intersects:
                query.set("where", exp.Where(this=where_expr.expression))
            elif where_expr.expression is intersects:
                query.set("where", exp.Where(this=where_expr.this))
            else:
                intersects.replace(exp.true())
        else:
            intersects.replace(exp.true())

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
        B = self.bin_size

        low = exp.Cast(
            this=exp.Div(
                this=exp.column(start, quoted=True),
                expression=exp.Literal.number(B),
            ),
            to=exp.DataType(this=exp.DataType.Type.BIGINT),
        )
        high = exp.Cast(
            this=exp.Add(
                this=exp.Div(
                    this=exp.Paren(
                        this=exp.Sub(
                            this=exp.column(end, quoted=True),
                            expression=exp.Literal.number(1),
                        ),
                    ),
                    expression=exp.Literal.number(B),
                ),
                expression=exp.Literal.number(1),
            ),
            to=exp.DataType(this=exp.DataType.Type.BIGINT),
        )

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
        connector_idx: list[int],
        *,
        preserve_kind: bool,
    ) -> list[exp.Join]:
        """Build three replacement JOINs for one INTERSECTS using the join-back pattern.

        join1: JOIN key_cte_for_other  connector_a  ON other_alias key-matches connector_a
        join2: JOIN key_cte_for_join   connector_b  ON connector_a equi-joins  connector_b
        join3: JOIN original_join_table join_alias  ON join_alias key-matches  connector_b
                                                    AND overlap predicate
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

        c0 = f"__giql_c{connector_idx[0]}"
        c1 = f"__giql_c{connector_idx[0] + 1}"
        connector_idx[0] += 2

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
