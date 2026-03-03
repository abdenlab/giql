"""Tests for BaseGIQLGenerator.

Test specification: specs/test_base.md
"""

import pytest
from hypothesis import HealthCheck
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st
from sqlglot import exp
from sqlglot import parse_one

from giql import Table
from giql.dialect import GIQLDialect
from giql.expressions import GIQLNearest
from giql.generators import BaseGIQLGenerator
from giql.table import Tables


@pytest.fixture
def tables_info():
    """Basic Tables with a single table containing genomic columns."""
    tables = Tables()
    tables.register("variants", Table("variants"))
    return tables


@pytest.fixture
def tables_with_two_tables():
    """Tables with two tables for column-to-column tests."""
    tables = Tables()
    tables.register("features_a", Table("features_a"))
    tables.register("features_b", Table("features_b"))
    return tables


@pytest.fixture
def tables_with_closed_intervals():
    """Tables with CLOSED interval type for bedtools compatibility tests."""
    tables = Tables()
    tables.register("bed_features", Table("bed_features", interval_type="closed"))
    return tables


@pytest.fixture
def tables_with_peaks_and_genes():
    """Tables with peaks and genes tables for NEAREST tests."""
    tables = Tables()
    tables.register("peaks", Table("peaks"))
    tables.register("genes", Table("genes"))
    return tables


class TestBaseGIQLGenerator:
    """Tests for BaseGIQLGenerator class."""

    def test_instantiation_defaults(self):
        """
        GIVEN no tables provided
        WHEN Generator is instantiated with defaults
        THEN Generator has empty Tables and SUPPORTS_LATERAL is True.
        """
        generator = BaseGIQLGenerator()

        assert generator.tables is not None
        assert "variants" not in generator.tables
        assert generator.SUPPORTS_LATERAL is True

    def test_instantiation_with_tables(self, tables_info):
        """
        GIVEN a valid Tables object with table definitions
        WHEN Generator is instantiated with tables
        THEN Generator stores tables and can resolve column references.
        """
        generator = BaseGIQLGenerator(tables=tables_info)

        assert generator.tables is tables_info
        assert "variants" in generator.tables

    def test_instantiation_kwargs_forwarding(self):
        """
        GIVEN Generator with custom kwargs
        WHEN Generator is instantiated with **kwargs
        THEN Generator passes kwargs to parent class.
        """
        # The parent Generator class accepts various kwargs like 'pretty'
        generator = BaseGIQLGenerator(pretty=True)

        # If kwargs forwarding works, generator should have pretty attribute
        assert generator.pretty is True

    def test_select_sql_basic(self, tables_info):
        """
        GIVEN a SELECT expression with FROM clause containing a table
        WHEN select_sql is called
        THEN Table context is tracked and alias mapping is built.
        """
        sql = "SELECT * FROM variants"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_info)
        output = generator.generate(ast)

        expected = "SELECT * FROM variants"
        assert output == expected

    def test_select_sql_with_alias(self, tables_info):
        """
        GIVEN a SELECT with aliased table (e.g., FROM table AS t)
        WHEN select_sql is called
        THEN Alias-to-table mapping includes the alias.
        """
        sql = "SELECT * FROM variants AS v WHERE v.interval INTERSECTS 'chr1:1000-2000'"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_info)
        output = generator.generate(ast)

        expected = (
            "SELECT * FROM variants AS v WHERE "
            '(v."chrom" = \'chr1\' AND v."start" < 2000 '
            'AND v."end" > 1000)'
        )
        assert output == expected

    def test_select_sql_with_joins(self, tables_with_two_tables):
        """
        GIVEN a SELECT with JOINs
        WHEN select_sql is called
        THEN All joined tables and aliases are tracked.
        """
        sql = "SELECT * FROM features_a AS a JOIN features_b AS b ON a.id = b.id"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_two_tables)
        output = generator.generate(ast)

        expected = "SELECT * FROM features_a AS a JOIN features_b AS b ON a.id = b.id"
        assert output == expected

    def test_intersects_sql_with_literal(self):
        """
        GIVEN an Intersects expression with literal range 'chr1:1000-2000'
        WHEN intersects_sql is called
        THEN SQL with chrom = 'chr1' AND start < 2000 AND end > 1000 is generated.
        """
        sql = "SELECT * FROM variants WHERE interval INTERSECTS 'chr1:1000-2000'"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator()
        output = generator.generate(ast)

        expected = (
            "SELECT * FROM variants WHERE "
            '("chrom" = \'chr1\' AND "start" < 2000 AND "end" > 1000)'
        )
        assert output == expected

    def test_intersects_sql_column_join(self, tables_with_two_tables):
        """
        GIVEN an Intersects expression with column-to-column
            (a.interval INTERSECTS b.interval)
        WHEN intersects_sql is called
        THEN SQL with column-to-column comparison is generated.
        """
        sql = (
            "SELECT * FROM features_a AS a CROSS JOIN features_b AS b "
            "WHERE a.interval INTERSECTS b.interval"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_two_tables)
        output = generator.generate(ast)

        expected = (
            "SELECT * FROM features_a AS a CROSS JOIN features_b AS b WHERE "
            '(a."chrom" = b."chrom" AND a."start" < b."end" '
            'AND a."end" > b."start")'
        )
        assert output == expected

    @given(
        chrom_num=st.sampled_from(["1", "2", "3", "X", "Y", "M"]),
        start=st.integers(min_value=0, max_value=1_000_000_000),
        length=st.integers(min_value=1, max_value=1_000_000),
    )
    def test_intersects_sql_validity_property(self, chrom_num, start, length):
        """
        GIVEN any valid chromosome and coordinates (start < end)
        WHEN Intersects SQL is generated and parsed
        THEN Generated SQL is syntactically valid.
        """
        chrom = f"chr{chrom_num}"
        end = start + length
        sql = f"SELECT * FROM variants WHERE interval INTERSECTS '{chrom}:{start}-{end}'"

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator()
        output = generator.generate(ast)

        # Verify we can parse the output SQL (proves it's syntactically valid)
        parsed = parse_one(output)
        assert parsed is not None

    def test_contains_sql_point_query(self):
        """
        GIVEN a Contains expression with point query 'chr1:1000'
        WHEN contains_sql is called
        THEN SQL with start <= 1000 AND end > 1000 is generated.
        """
        sql = "SELECT * FROM variants WHERE interval CONTAINS 'chr1:1500'"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator()
        output = generator.generate(ast)

        expected = (
            "SELECT * FROM variants WHERE "
            '("chrom" = \'chr1\' AND "start" <= 1500 AND "end" > 1500)'
        )
        assert output == expected

    def test_contains_sql_range_query(self):
        """
        GIVEN a Contains expression with range query 'chr1:1000-2000'
        WHEN contains_sql is called
        THEN SQL with start <= 1000 AND end >= 2000 is generated.
        """
        sql = "SELECT * FROM variants WHERE interval CONTAINS 'chr1:1500-2000'"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator()
        output = generator.generate(ast)

        expected = (
            "SELECT * FROM variants WHERE "
            '("chrom" = \'chr1\' AND "start" <= 1500 '
            'AND "end" >= 2000)'
        )
        assert output == expected

    def test_contains_sql_column_join(self, tables_with_two_tables):
        """
        GIVEN a Contains expression with column-to-column join
        WHEN contains_sql is called
        THEN SQL with left contains right comparison is generated.
        """
        sql = (
            "SELECT * FROM features_a AS a CROSS JOIN features_b AS b "
            "WHERE a.interval CONTAINS b.interval"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_two_tables)
        output = generator.generate(ast)

        expected = (
            "SELECT * FROM features_a AS a CROSS JOIN features_b AS b WHERE "
            '(a."chrom" = b."chrom" '
            'AND a."start" <= b."start" AND a."end" >= b."end")'
        )
        assert output == expected

    @given(
        chrom_num=st.sampled_from(["1", "2", "3", "X", "Y", "M"]),
        start=st.integers(min_value=0, max_value=1_000_000_000),
        length=st.integers(min_value=2, max_value=1_000_000),  # Min length 2 for range
    )
    def test_contains_sql_coordinate_validity_property(self, chrom_num, start, length):
        """
        GIVEN any valid genomic range coordinates
        WHEN Contains SQL is generated
        THEN the generated SQL is syntactically valid and references the range.
        """
        chrom = f"chr{chrom_num}"
        end = start + length
        sql = f"SELECT * FROM variants WHERE interval CONTAINS '{chrom}:{start}-{end}'"

        ast = parse_one(sql, dialect=GIQLDialect)
        generator = BaseGIQLGenerator()
        output = generator.generate(ast)

        # Verify we can parse the output SQL (proves it's syntactically valid)
        parsed = parse_one(output)
        assert parsed is not None

        # The output should reference the start value
        assert str(start) in output

    def test_within_sql_with_literal(self):
        """
        GIVEN a Within expression with literal range 'chr1:1000-2000'
        WHEN within_sql is called
        THEN SQL with start >= 1000 AND end <= 2000 is generated.
        """
        sql = "SELECT * FROM variants WHERE interval WITHIN 'chr1:1000-5000'"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator()
        output = generator.generate(ast)

        expected = (
            "SELECT * FROM variants WHERE "
            '("chrom" = \'chr1\' AND "start" >= 1000 AND "end" <= 5000)'
        )
        assert output == expected

    def test_within_sql_column_join(self, tables_with_two_tables):
        """
        GIVEN a Within expression with column-to-column join
        WHEN within_sql is called
        THEN SQL with left within right comparison is generated.
        """
        sql = (
            "SELECT * FROM features_a AS a CROSS JOIN features_b AS b "
            "WHERE a.interval WITHIN b.interval"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_two_tables)
        output = generator.generate(ast)

        expected = (
            "SELECT * FROM features_a AS a CROSS JOIN features_b AS b WHERE "
            '(a."chrom" = b."chrom" '
            'AND a."start" >= b."start" AND a."end" <= b."end")'
        )
        assert output == expected

    def test_spatialsetpredicate_sql_any(self):
        """
        GIVEN a SpatialSetPredicate with ANY quantifier and multiple ranges
        WHEN spatialsetpredicate_sql is called
        THEN SQL with OR-combined conditions is generated.
        """
        sql = (
            "SELECT * FROM variants "
            "WHERE interval INTERSECTS ANY('chr1:1000-2000', 'chr1:5000-6000')"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator()
        output = generator.generate(ast)

        expected = (
            "SELECT * FROM variants WHERE "
            '(("chrom" = \'chr1\' AND "start" < 2000 AND "end" > 1000) '
            'OR ("chrom" = \'chr1\' AND "start" < 6000 AND "end" > 5000))'
        )
        assert output == expected

    def test_spatialsetpredicate_sql_all(self):
        """
        GIVEN a SpatialSetPredicate with ALL quantifier and multiple ranges
        WHEN spatialsetpredicate_sql is called
        THEN SQL with AND-combined conditions is generated.
        """
        sql = (
            "SELECT * FROM variants "
            "WHERE interval INTERSECTS ALL('chr1:1000-2000', 'chr1:1500-1800')"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator()
        output = generator.generate(ast)

        expected = (
            "SELECT * FROM variants WHERE "
            '(("chrom" = \'chr1\' AND "start" < 2000 AND "end" > 1000) '
            'AND ("chrom" = \'chr1\' AND "start" < 1800 AND "end" > 1500))'
        )
        assert output == expected

    def test_giqlnearest_sql_standalone(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQLNearest in standalone mode with literal reference
        WHEN giqlnearest_sql is called
        THEN Subquery with ORDER BY distance LIMIT k is generated.
        """
        sql = "SELECT * FROM NEAREST(genes, reference='chr1:1000-2000', k=3)"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.generate(ast)

        expected = (
            "SELECT * FROM (\n"
            "                SELECT genes.*, "
            "CASE WHEN 'chr1' != genes.\"chrom\" THEN NULL "
            'WHEN 1000 < genes."end" AND 2000 > genes."start" THEN 0 '
            'WHEN 2000 <= genes."start" '
            'THEN (genes."start" - 2000) '
            'ELSE (1000 - genes."end") END AS distance\n'
            "                FROM genes\n"
            "                WHERE 'chr1' = genes.\"chrom\"\n"
            "                ORDER BY ABS("
            "CASE WHEN 'chr1' != genes.\"chrom\" THEN NULL "
            'WHEN 1000 < genes."end" AND 2000 > genes."start" THEN 0 '
            'WHEN 2000 <= genes."start" '
            'THEN (genes."start" - 2000) '
            'ELSE (1000 - genes."end") END)\n'
            "                LIMIT 3\n"
            "            )"
        )
        assert output == expected

    def test_giqlnearest_sql_correlated(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQLNearest in correlated mode (LATERAL join context)
        WHEN giqlnearest_sql is called
        THEN LATERAL-compatible subquery is generated.
        """
        sql = (
            "SELECT * FROM peaks "
            "CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=3)"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.generate(ast)

        expected = (
            "SELECT * FROM peaks CROSS JOIN LATERAL (\n"
            "                SELECT genes.*, "
            'CASE WHEN peaks."chrom" != genes."chrom" THEN NULL '
            'WHEN peaks."start" < genes."end" '
            'AND peaks."end" > genes."start" THEN 0 '
            'WHEN peaks."end" <= genes."start" '
            'THEN (genes."start" - peaks."end") '
            'ELSE (peaks."start" - genes."end") END AS distance\n'
            "                FROM genes\n"
            '                WHERE peaks."chrom" = genes."chrom"\n'
            "                ORDER BY ABS("
            'CASE WHEN peaks."chrom" != genes."chrom" THEN NULL '
            'WHEN peaks."start" < genes."end" '
            'AND peaks."end" > genes."start" THEN 0 '
            'WHEN peaks."end" <= genes."start" '
            'THEN (genes."start" - peaks."end") '
            'ELSE (peaks."start" - genes."end") END)\n'
            "                LIMIT 3\n"
            "            )"
        )
        assert output == expected

    def test_giqlnearest_sql_with_max_distance(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQLNearest with max_distance parameter
        WHEN giqlnearest_sql is called
        THEN WHERE clause includes distance filter.
        """
        sql = (
            "SELECT * FROM peaks "
            "CROSS JOIN LATERAL NEAREST("
            "genes, reference=peaks.interval, k=5, max_distance=100000)"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.generate(ast)

        expected = (
            "SELECT * FROM peaks CROSS JOIN LATERAL (\n"
            "                SELECT genes.*, "
            'CASE WHEN peaks."chrom" != genes."chrom" THEN NULL '
            'WHEN peaks."start" < genes."end" '
            'AND peaks."end" > genes."start" THEN 0 '
            'WHEN peaks."end" <= genes."start" '
            'THEN (genes."start" - peaks."end") '
            'ELSE (peaks."start" - genes."end") END AS distance\n'
            "                FROM genes\n"
            '                WHERE peaks."chrom" = genes."chrom" '
            "AND (ABS("
            'CASE WHEN peaks."chrom" != genes."chrom" THEN NULL '
            'WHEN peaks."start" < genes."end" '
            'AND peaks."end" > genes."start" THEN 0 '
            'WHEN peaks."end" <= genes."start" '
            'THEN (genes."start" - peaks."end") '
            'ELSE (peaks."start" - genes."end") END)) <= 100000\n'
            "                ORDER BY ABS("
            'CASE WHEN peaks."chrom" != genes."chrom" THEN NULL '
            'WHEN peaks."start" < genes."end" '
            'AND peaks."end" > genes."start" THEN 0 '
            'WHEN peaks."end" <= genes."start" '
            'THEN (genes."start" - peaks."end") '
            'ELSE (peaks."start" - genes."end") END)\n'
            "                LIMIT 5\n"
            "            )"
        )
        assert output == expected

    def test_giqlnearest_sql_stranded(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQLNearest with stranded=True
        WHEN giqlnearest_sql is called
        THEN Strand matching is included in WHERE clause.
        """
        sql = (
            "SELECT * FROM peaks "
            "CROSS JOIN LATERAL NEAREST("
            "genes, reference=peaks.interval, k=3, stranded=true)"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.generate(ast)

        expected = (
            "SELECT * FROM peaks CROSS JOIN LATERAL (\n"
            "                SELECT genes.*, "
            'CASE WHEN peaks."chrom" != genes."chrom" THEN NULL '
            'WHEN peaks."strand" IS NULL OR genes."strand" IS NULL THEN NULL '
            "WHEN peaks.\"strand\" = '.' OR peaks.\"strand\" = '?' THEN NULL "
            "WHEN genes.\"strand\" = '.' OR genes.\"strand\" = '?' THEN NULL "
            'WHEN peaks."start" < genes."end" '
            'AND peaks."end" > genes."start" THEN 0 '
            'WHEN peaks."end" <= genes."start" '
            "THEN CASE WHEN peaks.\"strand\" = '-' "
            'THEN -(genes."start" - peaks."end") '
            'ELSE (genes."start" - peaks."end") END '
            "ELSE CASE WHEN peaks.\"strand\" = '-' "
            'THEN -(peaks."start" - genes."end") '
            'ELSE (peaks."start" - genes."end") END END AS distance\n'
            "                FROM genes\n"
            '                WHERE peaks."chrom" = genes."chrom" '
            'AND peaks."strand" = genes."strand"\n'
            "                ORDER BY ABS("
            'CASE WHEN peaks."chrom" != genes."chrom" THEN NULL '
            'WHEN peaks."strand" IS NULL OR genes."strand" IS NULL THEN NULL '
            "WHEN peaks.\"strand\" = '.' OR peaks.\"strand\" = '?' THEN NULL "
            "WHEN genes.\"strand\" = '.' OR genes.\"strand\" = '?' THEN NULL "
            'WHEN peaks."start" < genes."end" '
            'AND peaks."end" > genes."start" THEN 0 '
            'WHEN peaks."end" <= genes."start" '
            "THEN CASE WHEN peaks.\"strand\" = '-' "
            'THEN -(genes."start" - peaks."end") '
            'ELSE (genes."start" - peaks."end") END '
            "ELSE CASE WHEN peaks.\"strand\" = '-' "
            'THEN -(peaks."start" - genes."end") '
            'ELSE (peaks."start" - genes."end") END END)\n'
            "                LIMIT 3\n"
            "            )"
        )
        assert output == expected

    def test_giqlnearest_sql_signed(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQLNearest with signed=True
        WHEN giqlnearest_sql is called
        THEN Distance expression includes signed calculation.
        """
        sql = (
            "SELECT * FROM peaks "
            "CROSS JOIN LATERAL NEAREST("
            "genes, reference=peaks.interval, k=3, signed=true)"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.generate(ast)

        expected = (
            "SELECT * FROM peaks CROSS JOIN LATERAL (\n"
            "                SELECT genes.*, "
            'CASE WHEN peaks."chrom" != genes."chrom" THEN NULL '
            'WHEN peaks."start" < genes."end" '
            'AND peaks."end" > genes."start" THEN 0 '
            'WHEN peaks."end" <= genes."start" '
            'THEN (genes."start" - peaks."end") '
            'ELSE -(peaks."start" - genes."end") END AS distance\n'
            "                FROM genes\n"
            '                WHERE peaks."chrom" = genes."chrom"\n'
            "                ORDER BY ABS("
            'CASE WHEN peaks."chrom" != genes."chrom" THEN NULL '
            'WHEN peaks."start" < genes."end" '
            'AND peaks."end" > genes."start" THEN 0 '
            'WHEN peaks."end" <= genes."start" '
            'THEN (genes."start" - peaks."end") '
            'ELSE -(peaks."start" - genes."end") END)\n'
            "                LIMIT 3\n"
            "            )"
        )
        assert output == expected

    def test_giqlnearest_sql_no_lateral_support(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQLNearest on a generator with SUPPORTS_LATERAL=False
        WHEN giqlnearest_sql is called in correlated mode
        THEN ValueError is raised with helpful message.
        """

        # Create a generator subclass without LATERAL support
        class NoLateralGenerator(BaseGIQLGenerator):
            SUPPORTS_LATERAL = False

        # Use query without explicit reference to trigger correlated mode
        sql = "SELECT * FROM peaks CROSS JOIN LATERAL NEAREST(genes, k=3)"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = NoLateralGenerator(tables=tables_with_peaks_and_genes)

        with pytest.raises(ValueError, match="LATERAL"):
            generator.generate(ast)

    @settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
    @given(
        k=st.integers(min_value=1, max_value=100),
        max_distance=st.integers(min_value=1, max_value=10_000_000),
    )
    def test_giqlnearest_sql_parameter_handling_property(
        self, tables_with_peaks_and_genes, k, max_distance
    ):
        """
        GIVEN any valid k value (positive integer) and max_distance
        WHEN Nearest SQL is generated
        THEN k appears in LIMIT clause, max_distance in WHERE.
        """
        sql = (
            f"SELECT * FROM NEAREST("
            f"genes, reference='chr1:1000-2000', k={k}, max_distance={max_distance})"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.generate(ast)

        # k should appear in LIMIT
        assert f"LIMIT {k}" in output
        # max_distance should appear in WHERE
        assert str(max_distance) in output

    def test_giqldistance_sql_basic(self, tables_with_two_tables):
        """
        GIVEN a GIQLDistance with two column references
        WHEN giqldistance_sql is called
        THEN CASE expression for distance calculation is generated.
        """
        sql = (
            "SELECT DISTANCE(a.interval, b.interval) as dist "
            "FROM features_a a CROSS JOIN features_b b"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_two_tables)
        output = generator.generate(ast)

        expected = (
            'SELECT CASE WHEN a."chrom" != b."chrom" THEN NULL '
            'WHEN a."start" < b."end" AND a."end" > b."start" '
            'THEN 0 WHEN a."end" <= b."start" '
            'THEN (b."start" - a."end") '
            'ELSE (a."start" - b."end") END AS dist '
            "FROM features_a AS a CROSS JOIN features_b AS b"
        )
        assert output == expected

    def test_giqldistance_sql_stranded(self, tables_with_two_tables):
        """
        GIVEN a GIQLDistance with stranded=True
        WHEN giqldistance_sql is called
        THEN Strand-aware distance CASE expression is generated.
        """
        sql = (
            "SELECT DISTANCE(a.interval, b.interval, stranded=true) as dist "
            "FROM features_a a CROSS JOIN features_b b"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_two_tables)
        output = generator.generate(ast)

        expected = (
            'SELECT CASE WHEN a."chrom" != b."chrom" THEN NULL '
            'WHEN a."strand" IS NULL OR b."strand" IS NULL THEN NULL '
            "WHEN a.\"strand\" = '.' OR a.\"strand\" = '?' THEN NULL "
            "WHEN b.\"strand\" = '.' OR b.\"strand\" = '?' THEN NULL "
            'WHEN a."start" < b."end" '
            'AND a."end" > b."start" THEN 0 '
            'WHEN a."end" <= b."start" '
            "THEN CASE WHEN a.\"strand\" = '-' "
            'THEN -(b."start" - a."end") '
            'ELSE (b."start" - a."end") END '
            "ELSE CASE WHEN a.\"strand\" = '-' "
            'THEN -(a."start" - b."end") '
            'ELSE (a."start" - b."end") END END AS dist '
            "FROM features_a AS a CROSS JOIN features_b AS b"
        )
        assert output == expected

    def test_giqldistance_sql_signed(self, tables_with_two_tables):
        """
        GIVEN a GIQLDistance with signed=True
        WHEN giqldistance_sql is called
        THEN Signed distance CASE expression is generated.
        """
        sql = (
            "SELECT DISTANCE(a.interval, b.interval, signed=true) as dist "
            "FROM features_a a CROSS JOIN features_b b"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_two_tables)
        output = generator.generate(ast)

        expected = (
            'SELECT CASE WHEN a."chrom" != b."chrom" THEN NULL '
            'WHEN a."start" < b."end" AND a."end" > b."start" '
            'THEN 0 WHEN a."end" <= b."start" '
            'THEN (b."start" - a."end") '
            'ELSE -(a."start" - b."end") END AS dist '
            "FROM features_a AS a CROSS JOIN features_b AS b"
        )
        assert output == expected

    def test_giqldistance_sql_stranded_and_signed(self, tables_with_two_tables):
        """
        GIVEN a GIQLDistance with both stranded and signed=True
        WHEN giqldistance_sql is called
        THEN Combined stranded+signed distance expression is generated.
        """
        sql = (
            "SELECT "
            "DISTANCE(a.interval, b.interval, stranded=true, signed=true) as dist "
            "FROM features_a a CROSS JOIN features_b b"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_two_tables)
        output = generator.generate(ast)

        expected = (
            'SELECT CASE WHEN a."chrom" != b."chrom" THEN NULL '
            'WHEN a."strand" IS NULL OR b."strand" IS NULL THEN NULL '
            "WHEN a.\"strand\" = '.' OR a.\"strand\" = '?' THEN NULL "
            "WHEN b.\"strand\" = '.' OR b.\"strand\" = '?' THEN NULL "
            'WHEN a."start" < b."end" '
            'AND a."end" > b."start" THEN 0 '
            'WHEN a."end" <= b."start" '
            "THEN CASE WHEN a.\"strand\" = '-' "
            'THEN -(b."start" - a."end") '
            'ELSE (b."start" - a."end") END '
            "ELSE CASE WHEN a.\"strand\" = '-' "
            'THEN (a."start" - b."end") '
            'ELSE -(a."start" - b."end") END END AS dist '
            "FROM features_a AS a CROSS JOIN features_b AS b"
        )
        assert output == expected

    def test_giqldistance_with_closed_intervals(self, tables_with_closed_intervals):
        """
        GIVEN intervals from table with CLOSED interval type
        WHEN Distance calculation is performed
        THEN Distance includes +1 adjustment (bedtools compatibility).
        """
        # Add a second table with closed intervals for distance calculation
        tables_with_closed_intervals.register(
            "bed_features_b", Table("bed_features_b", interval_type="closed")
        )

        sql = (
            "SELECT DISTANCE(a.interval, b.interval) as dist "
            "FROM bed_features a CROSS JOIN bed_features_b b"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_closed_intervals)
        output = generator.generate(ast)

        expected = (
            'SELECT CASE WHEN a."chrom" != b."chrom" THEN NULL '
            'WHEN a."start" < b."end" '
            'AND a."end" > b."start" THEN 0 '
            'WHEN a."end" <= b."start" '
            'THEN (b."start" - a."end" + 1) '
            'ELSE (a."start" - b."end" + 1) END AS dist '
            "FROM bed_features AS a CROSS JOIN bed_features_b AS b"
        )
        assert output == expected

    def test_error_handling_invalid_range(self):
        """
        GIVEN invalid genomic range string in Intersects
        WHEN intersects_sql is called
        THEN ValueError with descriptive message is raised.
        """
        sql = "SELECT * FROM variants WHERE interval INTERSECTS 'invalid'"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator()

        with pytest.raises(ValueError, match="Could not parse genomic range"):
            generator.generate(ast)

    def test_error_handling_unknown_operation(self):
        """
        GIVEN unknown operation type in spatial operations
        WHEN a spatial operation with unknown op_type is attempted
        THEN ValueError is raised.

        Note: This test verifies internal error handling by directly calling
        a method with invalid input, which would only occur through code errors.
        """
        # This is an indirect test - we verify the generator raises ValueError
        # when given malformed range strings as that's how errors surface
        sql = "SELECT * FROM variants WHERE interval INTERSECTS 'chr:a-b'"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator()

        with pytest.raises(ValueError):
            generator.generate(ast)

    def test_select_sql_join_without_alias(self, tables_with_two_tables):
        """
        GIVEN a SELECT with JOIN where joined table has no alias
        WHEN select_sql is called
        THEN Table name is used directly in alias mapping.
        """
        sql = "SELECT * FROM features_a JOIN features_b ON features_a.id = features_b.id"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_two_tables)
        output = generator.generate(ast)

        expected = (
            "SELECT * FROM features_a JOIN features_b ON features_a.id = features_b.id"
        )
        assert output == expected

    def test_giqlnearest_sql_stranded_literal_with_strand(
        self, tables_with_peaks_and_genes
    ):
        """
        GIVEN a GIQLNearest with stranded=True and literal reference containing strand
        WHEN giqlnearest_sql is called
        THEN Strand from literal range is parsed and used in filtering.
        """
        sql = (
            "SELECT * FROM NEAREST("
            "genes, reference='chr1:1000-2000:+', k=3, stranded=true)"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.generate(ast)

        # Should contain strand literal '+' and strand filtering
        assert "'+'" in output
        assert 'genes."strand"' in output

    def test_giqlnearest_sql_stranded_implicit_reference(
        self, tables_with_peaks_and_genes
    ):
        """
        GIVEN a GIQLNearest in correlated mode with implicit reference and stranded=True
        WHEN giqlnearest_sql is called
        THEN Strand column is resolved from outer table and used.
        """
        sql = "SELECT * FROM peaks CROSS JOIN LATERAL NEAREST(genes, k=3, stranded=true)"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.generate(ast)

        # Should have strand columns from both tables
        assert 'peaks."strand"' in output
        assert 'genes."strand"' in output

    def test_giqlnearest_sql_closed_intervals(self):
        """
        GIVEN a GIQLNearest with target table using CLOSED interval type
        WHEN giqlnearest_sql is called
        THEN Distance calculation includes +1 adjustment for bedtools compatibility.
        """
        tables = Tables()
        tables.register("genes_closed", Table("genes_closed", interval_type="closed"))

        sql = "SELECT * FROM NEAREST(genes_closed, reference='chr1:1000-2000', k=3)"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables)
        output = generator.generate(ast)

        # Should have +1 adjustment for closed intervals
        assert "+ 1)" in output

    def test_giqldistance_sql_literal_first_arg_error(self, tables_with_two_tables):
        """
        GIVEN a GIQLDistance with literal range as first argument
        WHEN giqldistance_sql is called
        THEN ValueError is raised indicating literals not supported.
        """
        sql = "SELECT DISTANCE('chr1:1000-2000', b.interval) as dist FROM features_b b"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_two_tables)

        with pytest.raises(ValueError, match="Literal range as first argument"):
            generator.generate(ast)

    def test_giqldistance_sql_literal_second_arg_error(self, tables_with_two_tables):
        """
        GIVEN a GIQLDistance with literal range as second argument
        WHEN giqldistance_sql is called
        THEN ValueError is raised indicating literals not supported.
        """
        sql = "SELECT DISTANCE(a.interval, 'chr1:1000-2000') as dist FROM features_a a"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_two_tables)

        with pytest.raises(ValueError, match="Literal range as second argument"):
            generator.generate(ast)

    def test_giqlnearest_sql_missing_outer_table_error(
        self, tables_with_peaks_and_genes
    ):
        """
        GIVEN a GIQLNearest in correlated mode without reference where outer table
            cannot be found
        WHEN giqlnearest_sql is called
        THEN ValueError is raised with helpful message about specifying reference.
        """

        nearest = GIQLNearest(
            this=exp.Table(this=exp.Identifier(this="genes")),
            k=exp.Literal.number(3),
        )

        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)

        with pytest.raises(ValueError, match="Could not find outer table"):
            generator.giqlnearest_sql(nearest)

    def test_giqlnearest_sql_outer_table_not_in_tables(self):
        """
        GIVEN a GIQLNearest in correlated mode where outer table is not registered
        WHEN giqlnearest_sql is called
        THEN ValueError is raised listing the issue.
        """
        tables = Tables()
        tables.register("genes", Table("genes"))

        nearest = GIQLNearest(
            this=exp.Table(this=exp.Identifier(this="genes")),
            k=exp.Literal.number(3),
        )

        generator = BaseGIQLGenerator(tables=tables)
        generator._alias_to_table = {"unknown_table": "unknown_table"}
        generator._find_outer_table_in_lateral_join = lambda x: "unknown_table"

        with pytest.raises(ValueError, match="not found in tables"):
            generator.giqlnearest_sql(nearest)

    def test_giqlnearest_sql_invalid_reference_range(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQLNearest with invalid/unparseable reference range string
        WHEN giqlnearest_sql is called
        THEN ValueError is raised with parse error details.
        """
        sql = "SELECT * FROM NEAREST(genes, reference='invalid_range', k=3)"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)

        with pytest.raises(ValueError, match="Could not parse reference genomic range"):
            generator.generate(ast)

    def test_giqlnearest_sql_no_tables_error(self):
        """
        GIVEN a GIQLNearest without tables registered
        WHEN giqlnearest_sql is called
        THEN ValueError is raised because target table cannot be resolved.
        """
        sql = "SELECT * FROM NEAREST(genes, reference='chr1:1000-2000', k=3)"
        ast = parse_one(sql, dialect=GIQLDialect)

        # Generator with empty tables - table won't be found
        generator = BaseGIQLGenerator()

        with pytest.raises(ValueError, match="not found in tables"):
            generator.generate(ast)

    def test_giqlnearest_sql_target_not_in_tables(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQLNearest with target table not registered
        WHEN giqlnearest_sql is called
        THEN ValueError is raised listing available tables.
        """
        sql = "SELECT * FROM NEAREST(unknown_table, reference='chr1:1000-2000', k=3)"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)

        with pytest.raises(ValueError, match="not found in tables"):
            generator.generate(ast)

    def test_intersects_sql_unqualified_column(self):
        """
        GIVEN an unqualified column reference (no table prefix) in spatial operation
        WHEN intersects_sql is called
        THEN Default column names are used without table qualifier.
        """
        sql = "SELECT * FROM variants WHERE interval INTERSECTS 'chr1:1000-2000'"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator()
        output = generator.generate(ast)

        expected = (
            "SELECT * FROM variants WHERE "
            '("chrom" = \'chr1\' AND "start" < 2000 AND "end" > 1000)'
        )
        assert output == expected

    def test_giqlnearest_sql_stranded_unqualified_reference(
        self, tables_with_peaks_and_genes
    ):
        """
        GIVEN a GIQLNearest with stranded=True and unqualified column reference
        WHEN giqlnearest_sql is called
        THEN Strand column is resolved without table prefix.
        """

        # Create NEAREST with stranded=True and an unqualified column reference
        # The reference is an unqualified column (no table prefix)
        nearest = GIQLNearest(
            this=exp.Table(this=exp.Identifier(this="genes")),
            reference=exp.Column(this=exp.Identifier(this="interval")),
            k=exp.Literal.number(3),
            stranded=exp.Boolean(this=True),
        )

        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.giqlnearest_sql(nearest)

        # Should produce valid output with unqualified strand column
        assert "LIMIT 3" in output
        # The strand column should be unqualified (no table prefix)
        assert '"strand"' in output

    def test_giqlnearest_sql_identifier_target(self, tables_with_peaks_and_genes):
        """
        GIVEN a GIQLNearest where target is an Identifier (not Table or Column)
        WHEN giqlnearest_sql is called
        THEN Target is converted to string and lookup proceeds.
        """

        # Use exp.Identifier directly - not Table or Column
        # This triggers the else branch at line 830 where str(target) is called
        nearest = GIQLNearest(
            this=exp.Identifier(this="genes"),
            reference=exp.Literal.string("chr1:1000-2000"),
            k=exp.Literal.number(3),
        )

        generator = BaseGIQLGenerator(tables=tables_with_peaks_and_genes)
        output = generator.giqlnearest_sql(nearest)

        # Should succeed and produce valid SQL
        assert "genes" in output
        assert "LIMIT 3" in output

    @given(
        bool_repr=st.sampled_from(["true", "TRUE", "True", "1", "yes", "YES"]),
    )
    @settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_giqldistance_stranded_param_truthy_values_property(
        self, tables_with_two_tables, bool_repr
    ):
        """
        GIVEN a GIQLDistance with stranded parameter in various truthy representations
        WHEN giqldistance_sql is called
        THEN The parameter is parsed as True and strand-aware distance is calculated.
        """
        sql = (
            f"SELECT DISTANCE(a.interval, b.interval, stranded={bool_repr}) as dist "
            "FROM features_a a, features_b b"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_two_tables)
        output = generator.generate(ast)

        # Should include strand handling (NULL checks for strand columns)
        assert "strand" in output.lower()
        assert "NULL" in output

    @given(
        bool_repr=st.sampled_from(["false", "FALSE", "False", "0", "no", "NO"]),
    )
    @settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_giqldistance_stranded_param_falsy_values_property(
        self, tables_with_two_tables, bool_repr
    ):
        """
        GIVEN a GIQLDistance with stranded parameter in various falsy representations
        WHEN giqldistance_sql is called
        THEN The parameter is parsed as False and basic distance is calculated.
        """
        sql = (
            f"SELECT DISTANCE(a.interval, b.interval, stranded={bool_repr}) as dist "
            "FROM features_a a, features_b b"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_two_tables)
        output = generator.generate(ast)

        # Should NOT include strand NULL checks (basic distance)
        assert "strand" not in output.lower()

    @given(
        bool_repr=st.sampled_from(["true", "TRUE", "True", "1", "yes", "YES"]),
    )
    @settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_giqldistance_signed_param_truthy_values_property(
        self, tables_with_two_tables, bool_repr
    ):
        """
        GIVEN a GIQLDistance with signed parameter in various truthy representations
        WHEN giqldistance_sql is called
        THEN The parameter is parsed as True and signed distance is calculated.
        """
        sql = (
            f"SELECT DISTANCE(a.interval, b.interval, signed={bool_repr}) as dist "
            "FROM features_a a, features_b b"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_two_tables)
        output = generator.generate(ast)

        # Signed distance has negative sign for upstream intervals
        assert "-(" in output

    @given(
        bool_repr=st.sampled_from(["false", "FALSE", "False", "0", "no", "NO"]),
    )
    @settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_giqldistance_signed_param_falsy_values_property(
        self, tables_with_two_tables, bool_repr
    ):
        """
        GIVEN a GIQLDistance with signed parameter in various falsy representations
        WHEN giqldistance_sql is called
        THEN The parameter is parsed as False and unsigned distance is calculated.
        """
        sql = (
            f"SELECT DISTANCE(a.interval, b.interval, signed={bool_repr}) as dist "
            "FROM features_a a, features_b b"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator(tables=tables_with_two_tables)
        output = generator.generate(ast)

        # Unsigned distance has no negative sign (both ELSE branches are positive)
        # Count occurrences of "-(" - signed has 1, unsigned has 0
        assert output.count("-(") == 0
