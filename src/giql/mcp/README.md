# GIQL MCP Server

An [MCP](https://modelcontextprotocol.io/) server that exposes GIQL documentation, operator references, and syntax guides to LLM-powered tools.

## Installation

```bash
pip install giql[mcp]
```

## Running the server

The `giql-mcp` console script starts the server using stdio transport:

```bash
giql-mcp
```

## Client configuration

Add the server to your MCP client config (e.g. Claude Desktop `claude_desktop_config.json`):

```json
{
  "mcpServers": {
    "giql": {
      "command": "giql-mcp"
    }
  }
}
```

## Interactive testing

Use the [FastMCP Inspector](https://gofastmcp.com/getting-started/inspector) to explore the server in a web UI:

```bash
fastmcp dev inspector src/giql/mcp/server.py --with-editable .
```

Then open `http://localhost:6274` in your browser.

## Available tools and resources

### Tools

| Tool | Description |
|---|---|
| `list_operators` | List all GIQL operators with brief descriptions |
| `explain_operator` | Get detailed documentation for a specific operator |
| `get_syntax_reference` | Get the complete GIQL syntax quick reference |
| `search_docs` | Search across GIQL documentation |
| `transpile` | Transpile a GIQL query to SQL |

### Resources

| URI Template | Description |
|---|---|
| `giql://docs/{path}` | Access documentation by path (e.g. `dialect/spatial-operators`, `recipes/intersect`) |
