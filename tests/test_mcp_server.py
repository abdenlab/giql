"""Tests for the GIQL MCP server."""

from __future__ import annotations

import pytest

from giql.mcp import server
from giql.mcp.server import (
    DOC_PATHS,
    OPERATORS,
    clean_rst,
    explain_operator,
    find_docs_root,
    get_documentation,
    get_syntax_reference,
    list_operators,
    search_docs,
)


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture()
def docs_root(tmp_path, monkeypatch):
    """Set up a temporary docs tree with synthetic RST files."""
    for rel_path in DOC_PATHS.values():
        doc_file = tmp_path / rel_path
        doc_file.parent.mkdir(parents=True, exist_ok=True)
        stem = doc_file.stem.replace("-", " ").title()
        doc_file.write_text(
            f"{stem}\n{'=' * len(stem)}\n\nContent about {stem}.\n",
            encoding="utf-8",
        )

    # Write a richer spatial-operators file so explain_operator can extract
    spatial = tmp_path / "dialect" / "spatial-operators.rst"
    spatial.write_text(
        "Spatial Operators\n"
        "=================\n\n"
        "INTERSECTS\n"
        "----------\n"
        "Overlap check details.\n\n"
        "CONTAINS\n"
        "--------\n"
        "Containment check details.\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(server, "DOCS_ROOT", tmp_path)
    return tmp_path


@pytest.fixture()
def no_docs(monkeypatch):
    """Simulate missing documentation directory."""
    monkeypatch.setattr(server, "DOCS_ROOT", None)


# =============================================================================
# TestCleanRst
# =============================================================================


class TestCleanRst:
    """Tests for clean_rst utility."""

    def test_removes_toctree(self):
        """
        GIVEN RST content with a toctree directive
        WHEN cleaned
        THEN the toctree block is removed
        """
        rst = "Title\n=====\n\n.. toctree::\n   :maxdepth: 2\n\n   page1\n   page2\n\nBody text."
        result = clean_rst(rst)
        assert "toctree" not in result
        assert "Body text." in result

    def test_converts_code_blocks(self):
        """
        GIVEN RST content with a code-block directive
        WHEN cleaned
        THEN the directive is converted to markdown fenced code
        """
        rst = "Example\n=======\n\n.. code-block:: sql\n\n   SELECT * FROM t\n"
        result = clean_rst(rst)
        assert "```sql" in result
        assert "SELECT * FROM t" in result

    def test_cleans_cross_references(self):
        """
        GIVEN RST content with :ref: and :doc: cross-references
        WHEN cleaned
        THEN they are replaced with plain text labels
        """
        rst = "See :ref:`label <target>` and :doc:`guide`."
        result = clean_rst(rst)
        assert ":ref:" not in result
        assert ":doc:" not in result
        assert "label" in result
        assert "guide" in result


# =============================================================================
# TestFindDocsRoot
# =============================================================================


class TestFindDocsRoot:
    """Tests for find_docs_root path discovery."""

    def test_env_var_strategy(self, tmp_path, monkeypatch):
        """
        GIVEN the GIQL_DOCS environment variable points to an existing directory
        WHEN find_docs_root is called
        THEN it returns that directory
        """
        monkeypatch.setenv("GIQL_DOCS", str(tmp_path))
        result = find_docs_root()
        assert result == tmp_path

    def test_nonexistent_env_var_fallthrough(self, tmp_path, monkeypatch):
        """
        GIVEN GIQL_DOCS points to a nonexistent directory
        WHEN find_docs_root is called
        THEN it does not return that path
        """
        monkeypatch.setenv("GIQL_DOCS", str(tmp_path / "nonexistent"))
        result = find_docs_root()
        assert result != tmp_path / "nonexistent"


# =============================================================================
# TestDocPaths
# =============================================================================


class TestDocPaths:
    """Tests for the DOC_PATHS mapping."""

    def test_all_values_end_in_rst(self):
        """
        GIVEN the DOC_PATHS dict
        WHEN inspecting all values
        THEN every value ends with .rst
        """
        for key, value in DOC_PATHS.items():
            assert value.endswith(".rst"), f"{key} -> {value} does not end in .rst"

    def test_valid_directory_prefixes(self):
        """
        GIVEN the DOC_PATHS dict
        WHEN inspecting directory prefixes of values
        THEN only known subdirectories are used
        """
        allowed = {"dialect", "transpilation", "recipes"}
        for key, value in DOC_PATHS.items():
            parts = value.split("/")
            if len(parts) > 1:
                assert parts[0] in allowed, f"Unexpected directory: {parts[0]}"

    def test_operator_doc_files_use_dialect(self):
        """
        GIVEN the OPERATORS dict
        WHEN inspecting doc_file values
        THEN all start with 'dialect/'
        """
        for name, info in OPERATORS.items():
            doc_file = info["doc_file"]
            assert doc_file.startswith("dialect/"), (
                f"{name} doc_file {doc_file!r} does not start with 'dialect/'"
            )


# =============================================================================
# TestListOperators
# =============================================================================


class TestListOperators:
    """Tests for the list_operators tool."""

    def test_returns_all_operators(self):
        """
        GIVEN the GIQL MCP server
        WHEN list_operators is called
        THEN it returns all 9 operators
        """
        result = list_operators()
        assert len(result) == 9

    def test_each_has_required_fields(self):
        """
        GIVEN the GIQL MCP server
        WHEN list_operators is called
        THEN each entry has name, category, description, and syntax
        """
        required = {"name", "category", "description", "syntax"}
        for entry in list_operators():
            assert required.issubset(entry.keys()), (
                f"Missing fields in {entry.get('name', '???')}"
            )


# =============================================================================
# TestExplainOperator
# =============================================================================


class TestExplainOperator:
    """Tests for the explain_operator tool."""

    def test_known_operator_with_docs(self, docs_root):
        """
        GIVEN a docs root with spatial-operators.rst containing an INTERSECTS section
        WHEN explain_operator is called for INTERSECTS
        THEN the result includes full_documentation
        """
        result = explain_operator("INTERSECTS")
        assert result["name"] == "INTERSECTS"
        assert result["full_documentation"] is not None
        assert "Overlap check" in result["full_documentation"]

    def test_case_insensitive(self):
        """
        GIVEN a valid operator name in lowercase
        WHEN explain_operator is called
        THEN it resolves the operator correctly
        """
        result = explain_operator("distance")
        assert result["name"] == "DISTANCE"
        assert "error" not in result

    def test_unknown_operator_error(self):
        """
        GIVEN an invalid operator name
        WHEN explain_operator is called
        THEN it returns an error dict with available operators
        """
        result = explain_operator("BOGUS")
        assert "error" in result
        assert "available" in result

    def test_no_docs_returns_none_full_documentation(self, no_docs):
        """
        GIVEN DOCS_ROOT is None
        WHEN explain_operator is called for a valid operator
        THEN full_documentation is None
        """
        result = explain_operator("MERGE")
        assert result["name"] == "MERGE"
        assert result["full_documentation"] is None


# =============================================================================
# TestGetSyntaxReference
# =============================================================================


class TestGetSyntaxReference:
    """Tests for the get_syntax_reference tool."""

    def test_returns_embedded_reference(self):
        """
        GIVEN the MCP server
        WHEN get_syntax_reference is called
        THEN it returns a string containing all operator categories
        """
        result = get_syntax_reference()
        assert "GIQL Syntax Quick Reference" in result
        assert "INTERSECTS" in result
        assert "NEAREST" in result
        assert "CLUSTER" in result
        assert "ANY" in result


# =============================================================================
# TestGetDocumentation
# =============================================================================


class TestGetDocumentation:
    """Tests for the get_documentation resource."""

    def test_valid_path(self, docs_root):
        """
        GIVEN a docs root with synthetic RST files
        WHEN get_documentation is called with a valid path
        THEN it returns the cleaned file content
        """
        result = get_documentation("quickstart")
        assert "Quickstart" in result

    def test_invalid_path(self, docs_root):
        """
        GIVEN a docs root
        WHEN get_documentation is called with an invalid path
        THEN it returns an error listing available paths
        """
        result = get_documentation("nonexistent/page")
        assert "Documentation not found" in result
        assert "quickstart" in result

    def test_docs_unavailable(self, no_docs):
        """
        GIVEN DOCS_ROOT is None
        WHEN get_documentation is called with a valid path
        THEN it returns a message about inaccessible docs
        """
        result = get_documentation("quickstart")
        assert "not accessible" in result


# =============================================================================
# TestSearchDocs
# =============================================================================


class TestSearchDocs:
    """Tests for the search_docs tool."""

    def test_matching_term(self, docs_root):
        """
        GIVEN a docs root with synthetic RST files
        WHEN search_docs is called with a term present in a file
        THEN it returns results with path, title, and excerpt
        """
        results = search_docs("Quickstart")
        assert any(r.get("path") == "quickstart" for r in results)
        match = next(r for r in results if r.get("path") == "quickstart")
        assert "title" in match
        assert "excerpt" in match

    def test_no_results(self, docs_root):
        """
        GIVEN a docs root with synthetic RST files
        WHEN search_docs is called with a term not in any file
        THEN it returns a no-results message
        """
        results = search_docs("xyzzy_nonexistent_term_12345")
        assert any("No results" in r.get("message", "") for r in results)

    def test_docs_unavailable(self, no_docs):
        """
        GIVEN DOCS_ROOT is None
        WHEN search_docs is called
        THEN it returns an error about inaccessible docs
        """
        results = search_docs("anything")
        assert any("error" in r for r in results)
