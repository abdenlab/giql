"""Query transformers for GIQL operations.

This module contains transformers that rewrite queries containing GIQL-specific
operations (like CLUSTER and MERGE) into equivalent SQL with CTEs.
"""

from sqlglot import exp

from giql.constants import DEFAULT_CHROM_COL
from giql.constants import DEFAULT_END_COL
from giql.constants import DEFAULT_START_COL
from giql.constants import DEFAULT_STRAND_COL
from giql.expressions import GIQLCluster
from giql.expressions import GIQLMerge
from giql.table import Tables


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
