"""GIQL MCP Server implementation.

Provides documentation access via the Model Context Protocol.
"""

from __future__ import annotations

import importlib.resources
import os
import re
from pathlib import Path
from typing import Any

from fastmcp import FastMCP

from giql.table import Table
from giql.transpile import transpile as giql_transpile

# =============================================================================
# Operator Metadata
# =============================================================================

OPERATORS: dict[str, dict[str, Any]] = {
    "INTERSECTS": {
        "category": "spatial",
        "description": "Returns true when two genomic ranges overlap by at least one base pair",
        "syntax": "interval INTERSECTS 'chr1:1000-2000'",
        "parameters": [
            {"name": "interval", "description": "A genomic column or range"},
            {"name": "target", "description": "A genomic range literal or column"},
        ],
        "returns": "Boolean: true if ranges overlap, false otherwise",
        "example": "SELECT * FROM variants WHERE interval INTERSECTS 'chr1:1000-2000'",
        "doc_file": "dialect/spatial-operators.rst",
    },
    "CONTAINS": {
        "category": "spatial",
        "description": "Returns true when one genomic range fully contains another",
        "syntax": "interval CONTAINS 'chr1:1500'",
        "parameters": [
            {"name": "interval", "description": "A genomic column (the container)"},
            {"name": "target", "description": "Range or point to check containment of"},
        ],
        "returns": "Boolean: true if first range fully contains second",
        "example": "SELECT * FROM genes WHERE interval CONTAINS 'chr1:1500'",
        "doc_file": "dialect/spatial-operators.rst",
    },
    "WITHIN": {
        "category": "spatial",
        "description": "Returns true when one genomic range is fully contained within another",
        "syntax": "interval WITHIN 'chr1:1000-5000'",
        "parameters": [
            {"name": "interval", "description": "A genomic column (must be within)"},
            {"name": "container", "description": "The containing range"},
        ],
        "returns": "Boolean: true if first range is fully within second",
        "example": "SELECT * FROM variants WHERE interval WITHIN 'chr1:1000-5000'",
        "doc_file": "dialect/spatial-operators.rst",
    },
    "DISTANCE": {
        "category": "distance",
        "description": "Calculate the genomic distance between two intervals",
        "syntax": "DISTANCE(interval_a, interval_b)",
        "parameters": [
            {"name": "interval_a", "description": "First genomic interval"},
            {"name": "interval_b", "description": "Second genomic interval"},
        ],
        "returns": "0 for overlapping, positive integer for gap, NULL for different chromosomes",
        "example": "SELECT DISTANCE(a.interval, b.interval) AS dist FROM a CROSS JOIN b",
        "doc_file": "dialect/distance-operators.rst",
    },
    "NEAREST": {
        "category": "distance",
        "description": "Find the k-nearest genomic features to a reference position",
        "syntax": "CROSS JOIN LATERAL NEAREST(table, reference=interval, k=5)",
        "parameters": [
            {"name": "target_table", "description": "Table to search for nearest features"},
            {"name": "reference", "description": "Reference position or column"},
            {"name": "k", "description": "Number of nearest neighbors (default: 1)"},
            {"name": "max_distance", "description": "Maximum distance threshold (optional)"},
            {"name": "stranded", "description": "Same-strand only (default: false)"},
            {"name": "signed", "description": "Return signed distances (default: false)"},
        ],
        "returns": "Rows from target table with distance column",
        "example": "SELECT * FROM peaks CROSS JOIN LATERAL NEAREST(genes, reference=peaks.interval, k=3) AS nearest",
        "doc_file": "dialect/distance-operators.rst",
    },
    "CLUSTER": {
        "category": "aggregation",
        "description": "Assign cluster IDs to overlapping or nearby genomic intervals",
        "syntax": "CLUSTER(interval) AS cluster_id",
        "parameters": [
            {"name": "interval", "description": "Genomic column to cluster"},
            {"name": "distance", "description": "Max gap to consider same cluster (default: 0)"},
            {"name": "stranded", "description": "Cluster by strand (default: false)"},
        ],
        "returns": "Integer cluster ID",
        "example": "SELECT *, CLUSTER(interval) AS cluster_id FROM features",
        "doc_file": "dialect/aggregation-operators.rst",
    },
    "MERGE": {
        "category": "aggregation",
        "description": "Combine overlapping genomic intervals into unified regions",
        "syntax": "SELECT MERGE(interval) FROM table",
        "parameters": [
            {"name": "interval", "description": "Genomic column to merge"},
            {"name": "distance", "description": "Max gap to merge (default: 0)"},
            {"name": "stranded", "description": "Merge by strand (default: false)"},
        ],
        "returns": "Merged interval coordinates (chromosome, start_pos, end_pos)",
        "example": "SELECT MERGE(interval), COUNT(*) FROM features",
        "doc_file": "dialect/aggregation-operators.rst",
    },
    "ANY": {
        "category": "quantifier",
        "description": "Match if condition holds for any of the specified ranges",
        "syntax": "interval INTERSECTS ANY('chr1:1000-2000', 'chr2:5000-6000')",
        "parameters": [
            {"name": "ranges", "description": "Comma-separated list of genomic ranges"},
        ],
        "returns": "Boolean: true if condition holds for at least one range",
        "example": "SELECT * FROM variants WHERE interval INTERSECTS ANY('chr1:1000-2000', 'chr2:5000-6000')",
        "doc_file": "dialect/quantifiers.rst",
    },
    "ALL": {
        "category": "quantifier",
        "description": "Match if condition holds for all of the specified ranges",
        "syntax": "interval CONTAINS ALL('chr1:1500', 'chr1:1600')",
        "parameters": [
            {"name": "ranges", "description": "Comma-separated list of genomic ranges"},
        ],
        "returns": "Boolean: true if condition holds for all ranges",
        "example": "SELECT * FROM genes WHERE interval CONTAINS ALL('chr1:1500', 'chr1:1600', 'chr1:1700')",
        "doc_file": "dialect/quantifiers.rst",
    },
}

# =============================================================================
# Path Discovery
# =============================================================================


def find_docs_root() -> Path | None:
    """Locate GIQL documentation directory.

    Searches in order:
    1. GIQL_DOCS environment variable
    2. Package data (bundled docs)
    3. Development repo structure
    """
    # Strategy 1: Environment variable
    if "GIQL_DOCS" in os.environ:
        docs = Path(os.environ["GIQL_DOCS"])
        if docs.exists():
            return docs

    # Strategy 2: Package data via importlib.resources
    try:
        # For Python 3.11+, use files() API
        files = importlib.resources.files("giql.mcp")
        docs_path = files / "docs"
        # Check if it's a real path (not in a zip)
        if hasattr(docs_path, "_path"):
            path = Path(docs_path._path)
            if path.exists():
                return path
        # Try to get as a traversable
        if docs_path.is_dir():
            # For installed packages, extract path
            return Path(str(docs_path))
    except (TypeError, AttributeError, FileNotFoundError):
        pass

    # Strategy 3: Development structure (src/giql/mcp -> repo/docs)
    module_path = Path(__file__).resolve()
    repo_root = module_path.parent.parent.parent.parent
    docs = repo_root / "docs"
    if docs.exists():
        return docs

    return None


DOCS_ROOT = find_docs_root()

# =============================================================================
# Documentation Paths
# =============================================================================

DOC_PATHS: dict[str, str] = {
    "index": "index.rst",
    "quickstart": "quickstart.rst",
    "dialect/index": "dialect/index.rst",
    "dialect/spatial-operators": "dialect/spatial-operators.rst",
    "dialect/distance-operators": "dialect/distance-operators.rst",
    "dialect/aggregation-operators": "dialect/aggregation-operators.rst",
    "dialect/quantifiers": "dialect/quantifiers.rst",
    "transpilation/index": "transpilation/index.rst",
    "transpilation/api-reference": "transpilation/api-reference.rst",
    "transpilation/execution": "transpilation/execution.rst",
    "transpilation/performance": "transpilation/performance.rst",
    "transpilation/schema-mapping": "transpilation/schema-mapping.rst",
    "recipes/index": "recipes/index.rst",
    "recipes/intersect": "recipes/intersect.rst",
    "recipes/distance": "recipes/distance.rst",
    "recipes/clustering": "recipes/clustering.rst",
    "recipes/bedtools-migration": "recipes/bedtools-migration.rst",
    "recipes/advanced": "recipes/advanced.rst",
}

# =============================================================================
# RST Utilities
# =============================================================================


def clean_rst(content: str) -> str:
    """Convert RST content to clean text for LLM consumption."""
    # Remove directive markers but keep content
    content = re.sub(
        r"\.\. (toctree|contents)::.*?(?=\n[^\s\n]|\Z)", "", content, flags=re.DOTALL
    )

    # Simplify code blocks
    content = re.sub(r"\.\. code-block:: (\w+)\n\n", r"```\1\n", content)

    # Clean up RST cross-references
    content = re.sub(r":ref:`([^`<]+)(?:\s*<[^>]+>)?`", r"\1", content)
    content = re.sub(r":doc:`([^`]+)`", r"\1", content)
    content = re.sub(r"``([^`]+)``", r"`\1`", content)

    # Clean up list-table directives
    content = re.sub(
        r"\.\. list-table::.*?(?=\n[^\s*]|\Z)", "", content, flags=re.DOTALL
    )

    # Remove excess blank lines
    content = re.sub(r"\n{3,}", "\n\n", content)

    return content.strip()


def read_doc_file(relative_path: str) -> str | None:
    """Read and clean a documentation file."""
    if DOCS_ROOT is None:
        return None

    if not relative_path.endswith(".rst"):
        relative_path = relative_path + ".rst"

    file_path = DOCS_ROOT / relative_path
    if not file_path.exists():
        return None

    content = file_path.read_text(encoding="utf-8")
    return clean_rst(content)


# =============================================================================
# Server Setup
# =============================================================================

mcp = FastMCP(
    name="giql",
)

# =============================================================================
# Resources
# =============================================================================


@mcp.resource("giql://docs/{path}")
def get_documentation(path: str) -> str:
    """Get GIQL documentation by path.

    Available paths:
    - index, quickstart
    - dialect/index, dialect/spatial-operators, dialect/distance-operators,
      dialect/aggregation-operators, dialect/quantifiers
    - transpilation/index, transpilation/api-reference,
      transpilation/execution, transpilation/performance,
      transpilation/schema-mapping
    - recipes/index, recipes/intersect, recipes/distance,
      recipes/clustering, recipes/bedtools-migration, recipes/advanced
    """
    if path not in DOC_PATHS:
        available = ", ".join(sorted(DOC_PATHS.keys()))
        return f"Documentation not found for path: {path}. Available: {available}"

    content = read_doc_file(DOC_PATHS[path])
    if content is None:
        return "Documentation file not accessible. Ensure GIQL_DOCS is set correctly."

    return content


# =============================================================================
# Tools
# =============================================================================


@mcp.tool()
def list_operators() -> list[dict[str, str]]:
    """List all available GIQL operators with brief descriptions.

    Returns a list of operators organized by category:
    - Spatial: INTERSECTS, CONTAINS, WITHIN
    - Distance: DISTANCE, NEAREST
    - Aggregation: CLUSTER, MERGE
    - Quantifiers: ANY, ALL
    """
    result = []
    for name, info in OPERATORS.items():
        result.append(
            {
                "name": name,
                "category": info["category"],
                "description": info["description"],
                "syntax": info["syntax"],
            }
        )
    return result


@mcp.tool()
def explain_operator(name: str) -> dict[str, Any]:
    """Get detailed documentation for a specific GIQL operator.

    Args:
        name: Operator name (case-insensitive). One of:
              INTERSECTS, CONTAINS, WITHIN, DISTANCE, NEAREST,
              CLUSTER, MERGE, ANY, ALL
    """
    name_upper = name.upper().strip()

    if name_upper not in OPERATORS:
        return {"error": f"Unknown operator: {name}", "available": list(OPERATORS.keys())}

    op = OPERATORS[name_upper]

    # Try to get full docs from RST file
    full_docs = None
    if DOCS_ROOT and "doc_file" in op:
        content = read_doc_file(op["doc_file"])
        if content:
            # Extract section for this operator
            pattern = rf"^{name_upper}\n[~=\-]+\n(.*?)(?=\n[A-Z]+\n[~=\-]+|\Z)"
            match = re.search(pattern, content, re.MULTILINE | re.DOTALL)
            if match:
                full_docs = f"{name_upper}\n{'=' * len(name_upper)}\n{match.group(1).strip()}"

    return {
        "name": name_upper,
        "category": op["category"],
        "description": op["description"],
        "syntax": op["syntax"],
        "parameters": op["parameters"],
        "returns": op["returns"],
        "example": op["example"],
        "full_documentation": full_docs,
    }


@mcp.tool()
def get_syntax_reference() -> str:
    """Get the complete GIQL syntax quick reference.

    Returns a formatted reference covering:
    - Genomic range literal format
    - All operators with syntax
    - Common query patterns
    """
    return """
GIQL Syntax Quick Reference
===========================

Genomic Range Literals
----------------------
Format: 'chromosome:start-end'
Examples:
  'chr1:1000-2000'     -- Range from 1000 to 2000 on chr1
  'chr1:1000'          -- Point at position 1000
  'chrX:50000-100000'  -- Range on chrX

Coordinates: 0-based, half-open [start, end)

Spatial Operators
-----------------
INTERSECTS - Test if ranges overlap
  interval INTERSECTS 'chr1:1000-2000'
  a.interval INTERSECTS b.interval

CONTAINS - Test if range contains another
  interval CONTAINS 'chr1:1500'
  gene.interval CONTAINS variant.interval

WITHIN - Test if range is within another
  interval WITHIN 'chr1:1000-5000'

Distance Operators
------------------
DISTANCE - Calculate distance between intervals
  DISTANCE(a.interval, b.interval)

NEAREST - Find k-nearest neighbors
  CROSS JOIN LATERAL NEAREST(table, reference=interval, k=5) AS nearest

Aggregation Operators
---------------------
CLUSTER - Assign cluster IDs
  CLUSTER(interval) AS cluster_id
  CLUSTER(interval, 1000) AS cluster_id  -- with distance

MERGE - Combine overlapping intervals
  SELECT MERGE(interval) FROM table

Set Quantifiers
---------------
ANY - Match any of multiple ranges
  interval INTERSECTS ANY('chr1:1000-2000', 'chr2:5000-6000')

ALL - Match all of multiple ranges
  interval CONTAINS ALL('chr1:1500', 'chr1:1600')
"""


@mcp.tool()
def search_docs(query: str) -> list[dict[str, str]]:
    """Search across GIQL documentation.

    Args:
        query: Search term or phrase (case-insensitive)

    Returns:
        List of matching documents with path, title, and excerpt
    """
    if DOCS_ROOT is None:
        return [{"error": "Documentation not accessible. Set GIQL_DOCS."}]

    results = []
    query_lower = query.lower()

    for path, filename in DOC_PATHS.items():
        content = read_doc_file(filename)
        if content is None:
            continue

        if query_lower in content.lower():
            # Extract title (first heading)
            title_match = re.search(r"^([^\n]+)\n[=\-~]+", content)
            title = title_match.group(1) if title_match else path

            # Find excerpt around match
            idx = content.lower().find(query_lower)
            start = max(0, idx - 100)
            end = min(len(content), idx + len(query) + 100)
            excerpt = content[start:end].strip()
            if start > 0:
                excerpt = "..." + excerpt
            if end < len(content):
                excerpt = excerpt + "..."

            results.append(
                {
                    "path": path,
                    "title": title,
                    "excerpt": excerpt,
                }
            )

    if not results:
        return [{"message": f"No results found for '{query}'"}]

    return results


@mcp.tool()
def transpile(
    query: str, tables: list[str | dict[str, Any]] | None = None
) -> dict[str, str]:
    """Transpile a GIQL query to SQL.

    Converts a GIQL query string (with genomic operators like INTERSECTS,
    CONTAINS, WITHIN, etc.) into standard SQL.

    Args:
        query: The GIQL query string.
        tables: Optional list of table specifications. Each entry can be:
            - A string table name (uses default column mappings:
              chrom, start, end, strand)
            - A dict with keys: name (required), plus optional chrom_col,
              start_col, end_col, strand_col, genomic_col,
              coordinate_system, interval_type
    """
    try:
        converted: list[str | Table] | None = None
        if tables is not None:
            converted = []
            for entry in tables:
                if isinstance(entry, str):
                    converted.append(entry)
                else:
                    converted.append(Table(**entry))
        result = giql_transpile(query, tables=converted)
        return {"sql": result}
    except (ValueError, Exception) as e:
        return {"error": str(e)}


# =============================================================================
# Entry Point
# =============================================================================


def main() -> None:
    """Run the MCP server."""
    mcp.run()


if __name__ == "__main__":
    main()
