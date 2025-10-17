from sqlglot import parse_one

from giql.dialect import GIQLDialect
from giql.generators import BaseGIQLGenerator
from giql.generators import GIQLDuckDBGenerator


class TestBaseGenerator:
    def test_generate_simple_intersects(self):
        """
        GIVEN a SQL query with INTERSECTS operator
        WHEN generating SQL code
        THEN should produce standard SQL with range conditions
        """
        sql = "SELECT * FROM variants WHERE position INTERSECTS 'chr1:1000-2000'"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator()
        output = generator.generate(ast)

        # Should expand to chromosome/position checks
        assert "\"chromosome\" = 'chr1'" in output
        assert '"start_pos" < 2000' in output
        assert '"end_pos" > 1000' in output

    def test_generate_contains(self):
        """
        GIVEN a SQL query with CONTAINS operator
        WHEN generating SQL code
        THEN should produce containment conditions
        """
        sql = "SELECT * FROM variants WHERE position CONTAINS 'chr1:1500'"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator()
        output = generator.generate(ast)

        # Point query: start <= point < end
        assert "\"chromosome\" = 'chr1'" in output
        assert '"start_pos" <= 1500' in output
        assert '"end_pos" > 1500' in output

    def test_generate_within(self):
        """
        GIVEN a SQL query with WITHIN operator
        WHEN generating SQL code
        THEN should produce within conditions
        """
        sql = "SELECT * FROM variants WHERE position WITHIN 'chr1:1000-5000'"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator()
        output = generator.generate(ast)

        # Left within right: start1 >= start2 AND end1 <= end2
        assert "\"chromosome\" = 'chr1'" in output
        assert '"start_pos" >= 1000' in output
        assert '"end_pos" <= 5000' in output

    def test_generate_intersects_any(self):
        """
        GIVEN a SQL query with INTERSECTS ANY operator
        WHEN generating SQL code
        THEN should produce OR conditions
        """
        sql = (
            "SELECT * FROM v WHERE position INTERSECTS ANY("
            "'chr1:1000-2000', 'chr1:5000-6000')"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator()
        output = generator.generate(ast)

        # Should have two conditions combined with OR
        assert " OR " in output
        assert output.count("\"chromosome\" = 'chr1'") == 2

    def test_generate_intersects_all(self):
        """
        GIVEN a SQL query with INTERSECTS ALL operator
        WHEN generating SQL code
        THEN should produce AND conditions
        """
        sql = (
            "SELECT * FROM v WHERE position INTERSECTS ALL("
            "'chr1:1000-2000', 'chr1:1500-1800')"
        )
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator()
        output = generator.generate(ast)

        # Should have two conditions combined with AND
        assert " AND " in output
        assert output.count("\"chromosome\" = 'chr1'") == 2

    def test_generate_with_table_alias(self):
        """
        GIVEN a SQL query with table alias
        WHEN generating SQL code
        THEN should properly qualify column names
        """
        sql = "SELECT * FROM variants v WHERE v.position INTERSECTS 'chr1:1000-2000'"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator()
        output = generator.generate(ast)

        # Should use table alias in generated conditions
        assert 'v."chromosome"' in output
        assert 'v."start_pos"' in output
        assert 'v."end_pos"' in output

    def test_contains_range_query(self):
        """
        GIVEN a SQL query with CONTAINS on a range (not a point)
        WHEN generating SQL code
        THEN should use range containment logic
        """
        sql = "SELECT * FROM variants WHERE position CONTAINS 'chr1:1500-2000'"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator()
        output = generator.generate(ast)

        # Range containment: start1 <= start2 AND end1 >= end2
        assert "\"chromosome\" = 'chr1'" in output
        assert '"start_pos" <= 1500' in output
        assert '"end_pos" >= 2000' in output

    def test_invalid_range_string(self):
        """
        GIVEN a SQL query with invalid range format
        WHEN generating SQL code
        THEN should raise ValueError
        """
        sql = "SELECT * FROM variants WHERE position INTERSECTS 'invalid'"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = BaseGIQLGenerator()
        try:
            _ = generator.generate(ast)
            assert False, "Should have raised ValueError"
        except ValueError as e:
            assert "Could not parse genomic range" in str(e)


class TestDuckDBGenerator:
    def test_duckdb_generator_basic(self):
        """
        GIVEN a SQL query with INTERSECTS operator
        WHEN using DuckDB generator
        THEN should produce valid DuckDB SQL
        """
        sql = "SELECT * FROM variants WHERE position INTERSECTS 'chr1:1000-2000'"
        ast = parse_one(sql, dialect=GIQLDialect)

        generator = GIQLDuckDBGenerator()
        output = generator.generate(ast)

        # Should still have the basic range conditions
        assert "\"chromosome\" = 'chr1'" in output
        assert '"start_pos" < 2000' in output
        assert '"end_pos" > 1000' in output
