"""Functional suite for the GIQL operator usage-pattern catalogue.

Executes every catalogued DISJOIN usage pattern's transpiled SQL against a
matrix of database engines and snapshots the result rows with pytest-manifest,
verifying each pattern is valid SQL and functionally correct in its query
context. The functional matrix crosses two dimensions: the usage patterns of
``usage_patterns.py`` and the :data:`~tests.usage_patterns.DISJOIN_PROFILES`
semantic profiles -- one per DISJOIN edge case. The committed snapshots under
``tests/manifests/`` are the oracle: regenerate with ``pytest -O`` and review
the diff before committing.

Patterns GIQL cannot transpile (a subquery, CTE, or nested operator result in a
table-function's target position) are catalogued but strict-``xfail``.

The module also carries unit-tier guard tests for the catalogue's descriptor
machinery: :class:`TestTableFixture`, :class:`TestRenderTableFunction`, and
:class:`TestDisjoinProfiles`.
"""

import pytest

from giql import transpile
from giql.table import Table

from .usage_patterns import CANONICAL_COLUMNS
from .usage_patterns import DISJOIN_PROFILES
from .usage_patterns import DISJOIN_USAGE
from .usage_patterns import OperatorType
from .usage_patterns import OperatorUsage
from .usage_patterns import TableFixture
from .usage_patterns import UsagePattern
from .usage_patterns import render
from .usage_patterns import templated_patterns


class _TemplatelessUsage(OperatorUsage):
    """An OperatorUsage subclass with an empty template table, for guard tests."""

    name = "NO_OP"
    operator_type = OperatorType.TABLE_FUNCTION

    @property
    def tables(self):
        """Return no tables -- this descriptor is never executed."""
        return ()

    def template_table(self):
        """Return an empty template table so render() always rejects."""
        return {}

    def slots(self):
        """Return no slot values -- no template is ever rendered."""
        return {}


def _run_duckdb(fixtures: tuple[TableFixture, ...], sql: str) -> list:
    """Load `fixtures` into in-memory DuckDB, execute `sql`, return the rows.

    Each fixture's physical schema (`fixture.columns`) drives the `CREATE
    TABLE`, so a fixture carrying a custom `Table` config is materialized with
    its own physical column names and types.
    """
    import duckdb

    conn = duckdb.connect(":memory:")
    try:
        for fixture in fixtures:
            column_ddl = ", ".join(
                f'"{column}" {sql_type}' for column, sql_type in fixture.columns
            )
            conn.execute(f'CREATE TABLE "{fixture.name}" ({column_ddl})')
            if fixture.rows:
                placeholders = ", ".join("?" for _ in fixture.columns)
                conn.executemany(
                    f'INSERT INTO "{fixture.name}" VALUES ({placeholders})',
                    [list(row) for row in fixture.rows],
                )
        return conn.execute(sql).fetchall()
    finally:
        conn.close()


# Database engines the functional matrix runs against. DuckDB only for now;
# SQLite and PostgreSQL are future entries keyed in here.
ENGINES = {"duckdb": _run_duckdb}


def _execute(usage: OperatorUsage, query: str, engine: str) -> list:
    """Transpile `query`, run it on `engine`, return rows sorted for stability.

    The descriptor's `tables` yields a `Table` config for any fixture with a
    custom physical layout, so GIQL resolves custom column names and
    coordinate systems during transpilation.
    """
    sql = transpile(query, tables=list(usage.tables))
    rows = ENGINES[engine](usage.fixtures, sql)
    return sorted(([list(row) for row in rows]), key=repr)


def _profile_cases() -> list:
    """Return (profile, mode, usage, pattern, query, engine) params for the matrix.

    One param is produced per (profile, mode, pattern, engine). Each param
    carries a stable id, the pattern's `usage` mark, and -- for an
    untranspilable pattern -- a strict `xfail`. Reference-mode and self-mode
    are flattened into one parametrization so a single test method covers the
    whole matrix; `mode` is threaded into the manifest key so a profile's two
    modes never collide on the same snapshot row.
    """
    cases: list = []
    for profile in DISJOIN_PROFILES:
        modes = (("self", profile.self_usage), ("reference", profile.reference_usage))
        for mode, usage in modes:
            for pattern in templated_patterns(usage):
                marks = [pytest.mark.usage(pattern), pytest.mark.integration]
                if pattern in usage.xfails:
                    marks.append(
                        pytest.mark.xfail(
                            reason=usage.xfails[pattern], strict=True
                        )
                    )
                for engine in ENGINES:
                    cases.append(
                        pytest.param(
                            profile.profile_id,
                            mode,
                            usage,
                            pattern,
                            render(pattern, usage),
                            engine,
                            id=(
                                f"{profile.profile_id}-{mode}-"
                                f"{engine}-{pattern.name}"
                            ),
                            marks=marks,
                        )
                    )
    return cases


class TestDisjoinUsageMatrix:
    """Every DISJOIN usage pattern executes to a snapshotted result per profile."""

    @pytest.mark.parametrize(
        ("profile_id", "mode", "usage", "pattern", "query", "engine"),
        _profile_cases(),
    )
    def test_disjoin_pattern_executes_to_snapshot_across_profiles(
        self, profile_id, mode, usage, pattern, query, engine, manifest
    ):
        """Test every DISJOIN usage pattern executes to its per-profile snapshot.

        Given:
            A canonical usage pattern rendered for one DISJOIN semantic profile
            in either self-mode or reference-mode.
        When:
            Its transpiled SQL is executed against the engine.
        Then:
            It should produce result rows matching the captured manifest
            snapshot keyed by profile, mode, engine, and pattern.
        """
        # Arrange
        key = f"{profile_id}-{mode}-{engine}-{pattern.name}"

        # Act
        rows = _execute(usage, query, engine)

        # Assert
        assert manifest[key] == rows

    def test_templated_patterns_should_cover_every_canonical_pattern(self):
        """Test the DISJOIN catalogue templates every canonical pattern.

        Given:
            The canonical table-function usage-pattern list, including the
            OUTER_COUNT, OUTER_LEFT_JOIN, SELF_JOIN_OVERLAP, NESTED_AGGREGATE,
            and ORDERED_CTE members.
        When:
            Collecting the patterns the DISJOIN catalogue can render.
        Then:
            It should render every canonical pattern, the new members
            included.
        """
        # Act
        renderable = set(templated_patterns(DISJOIN_USAGE))

        # Assert
        assert renderable == set(DISJOIN_USAGE.canonical_patterns)
        new_members = {
            UsagePattern.OUTER_COUNT,
            UsagePattern.OUTER_LEFT_JOIN,
            UsagePattern.SELF_JOIN_OVERLAP,
            UsagePattern.NESTED_AGGREGATE,
            UsagePattern.ORDERED_CTE,
        }
        assert new_members <= renderable


class TestTableFixture:
    """Tests for the TableFixture catalogue fixture descriptor."""

    def test___init___with_name_and_rows_only(self):
        """Test a TableFixture defaults to the canonical schema and no Table.

        Given:
            A TableFixture constructed with only a name and rows.
        When:
            The fixture is inspected.
        Then:
            It should expose the canonical column schema and a None Table
            config.
        """
        # Act
        fixture = TableFixture("peaks", (("chr1", 0, 10, "p"),))

        # Assert
        assert fixture.columns == CANONICAL_COLUMNS
        assert fixture.table is None

    def test___init___with_custom_table_and_schema(self):
        """Test a TableFixture exposes a supplied Table config and column schema.

        Given:
            A TableFixture constructed with a custom Table config and an
            explicit physical column schema.
        When:
            The fixture is inspected.
        Then:
            It should expose both the custom schema and the Table config.
        """
        # Arrange
        schema = (("seqid", "VARCHAR"), ("lo", "INTEGER"), ("hi", "INTEGER"))
        config = Table("v", chrom_col="seqid", start_col="lo", end_col="hi")

        # Act
        fixture = TableFixture("v", (("chr1", 1, 9),), columns=schema, table=config)

        # Assert
        assert fixture.columns == schema
        assert fixture.table is config

    def test___init___with_zero_rows(self):
        """Test a TableFixture with no rows is a valid empty fixture.

        Given:
            A TableFixture constructed with an empty rows tuple.
        When:
            The fixture is inspected.
        Then:
            It should be a valid fixture exposing zero rows and the canonical
            schema.
        """
        # Act
        fixture = TableFixture("empty", ())

        # Assert
        assert fixture.rows == ()
        assert fixture.columns == CANONICAL_COLUMNS


class TestRenderTableFunction:
    """Tests for render() applied to table-function usage descriptors."""

    def test_render_with_new_member_pattern(self):
        """Test render() produces a complete query for a new-member pattern.

        Given:
            A DISJOIN_PROFILES descriptor and the NESTED_AGGREGATE usage
            pattern (a member added for issue #87).
        When:
            The pattern is rendered through render().
        Then:
            It should return a complete query string with every slot filled
            and no unsubstituted braces.
        """
        # Act
        query = render(UsagePattern.NESTED_AGGREGATE, DISJOIN_USAGE)

        # Assert
        assert "{" not in query and "}" not in query
        assert "DISJOIN(features)" in query
        assert "disjoin_end" in query and "disjoin_start" in query

    def test_render_with_pattern_having_no_template(self):
        """Test render() rejects a pattern its operator class cannot template.

        Given:
            An operator usage descriptor whose class has an empty template
            table, and a usage pattern absent from it.
        When:
            render() is called for that pattern.
        Then:
            It should raise ValueError naming the operator and the pattern.
        """
        # Arrange
        templateless = _TemplatelessUsage()

        # Act & assert
        with pytest.raises(ValueError, match=r"NO_OP.*NESTED_AGGREGATE"):
            render(UsagePattern.NESTED_AGGREGATE, templateless)


class TestDisjoinProfiles:
    """Tests for the DISJOIN_PROFILES semantic-profile collection."""

    @pytest.mark.parametrize("profile", DISJOIN_PROFILES, ids=lambda p: p.profile_id)
    def test_profile_patterns_reference_only_registered_fixtures(self, profile):
        """Test every table a profile's patterns reference is a known fixture.

        Given:
            A DISJOIN profile and its self-mode and reference-mode usages.
        When:
            Collecting every table name the rendered patterns reference.
        Then:
            It should reference only tables present in the usage's fixtures.
        """
        # Arrange
        for usage in (profile.self_usage, profile.reference_usage):
            fixture_names = {fixture.name for fixture in usage.fixtures}
            referenced = {
                usage.target,
                usage.join_table,
                usage.ref_operand,
            }

            # Act & assert
            assert referenced <= fixture_names

    @pytest.mark.parametrize("profile", DISJOIN_PROFILES, ids=lambda p: p.profile_id)
    def test_profile_names_target_and_reference_in_fixtures(self, profile):
        """Test every profile names a target and reference present in fixtures.

        Given:
            A DISJOIN profile's self-mode and reference-mode usages.
        When:
            Resolving the declared target and reference table names.
        Then:
            It should find the target in self-mode fixtures and both the
            target and the reference in reference-mode fixtures.
        """
        # Arrange
        self_names = {f.name for f in profile.self_usage.fixtures}
        ref_names = {f.name for f in profile.reference_usage.fixtures}

        # Act & assert
        assert profile.self_usage.target in self_names
        assert profile.reference_usage.target in ref_names
        assert profile.reference_usage.reference in ref_names

    def test_profile_ids_are_unique(self):
        """Test every DISJOIN profile id is unique across the collection.

        Given:
            The DISJOIN_PROFILES collection, whose ids prefix manifest keys.
        When:
            Collecting every profile id.
        Then:
            It should contain no duplicate id, so manifest keys never collide.
        """
        # Act
        ids = [profile.profile_id for profile in DISJOIN_PROFILES]

        # Assert
        assert len(ids) == len(set(ids))

    def test_profile_ids_match_nested_usage_profile_ids(self):
        """Test each profile's id matches the profile_id of its bound usages.

        Given:
            A DISJOIN profile and the two TableFunctionUsage descriptors it
            binds.
        When:
            Comparing the profile id to each usage's profile_id field.
        Then:
            It should match on both the self-mode and reference-mode usage.
        """
        # Act & assert
        for profile in DISJOIN_PROFILES:
            assert profile.self_usage.profile_id == profile.profile_id
            assert profile.reference_usage.profile_id == profile.profile_id
