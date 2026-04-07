---
name: understand-chat
description: >
  Ask questions about the codebase using the knowledge graph. Use this skill
  whenever the user says "/understand-chat <query>", "ask about the code",
  "what does X do", or similar. Requires a knowledge graph at
  .understand-anything/knowledge-graph.json (run /understand to generate one).
argument-hint: "[query]"
---

The key words MUST, MUST NOT, SHALL, SHALL NOT, SHOULD, SHOULD NOT, REQUIRED, RECOMMENDED, MAY, and OPTIONAL in this document are to be interpreted as described in RFC 2119.

# Understand-Chat Skill

Answer questions about this codebase using the knowledge graph at `.understand-anything/knowledge-graph.json`.

## Pipeline Context

This skill is a **standalone utility** — it is not part of the SDLC pipeline sequence (`/issue` → `/implement` → `/test` → `/commit` → `/pr`). It can be invoked at any time to query the knowledge graph for architectural context, component relationships, or dependency chains.

## Arguments

A query string MUST be provided as the argument (e.g., `/understand-chat how does transpilation work`).

## Graph Structure Reference

The knowledge graph JSON has this structure:

- `project` — `{name, description, languages, frameworks, analyzedAt, gitCommitHash}`
- `nodes[]` — each has `{id, type, name, filePath, summary, tags[], complexity, languageNotes?}`
  - Node types: `file`, `function`, `class`, `module`, `concept`, `config`, `resource`, `pipeline`
  - ID conventions: `file:path`, `function:path:name`, `class:path:name`
- `edges[]` — each has `{source, target, type, direction, weight}`
  - Key edge types: `imports`, `contains`, `calls`, `depends_on`, `tested_by`, `configures`, `extends`
- `layers[]` — each has `{id, name, description, nodeIds[]}`
- `tour[]` — each has `{order, title, description, nodeIds[]}`

## How to Read Efficiently

1. Use Grep to search within the JSON for relevant entries BEFORE reading the full file.
2. Only read sections you need — do NOT load the entire graph into context.
3. Node `name` and `summary` fields are the most useful for understanding.
4. Edges tell you how components connect — follow `imports` and `calls` for dependency chains.

## Workflow

### TL;DR

1. Verify graph exists
2. Read project metadata
3. Search for relevant nodes
4. Find connected edges
5. Read layer context
6. Answer the query

### 1. Verify graph exists

```bash
test -f .understand-anything/knowledge-graph.json && echo "exists" || echo "missing"
```

If the graph does not exist, inform the user and suggest running `/understand` to generate one. MUST NOT proceed without a graph.

### 2. Read project metadata

Read the `"project"` section (first ~25 lines) from the top of the knowledge graph file for high-level context — project name, description, languages, and frameworks.

### 3. Search for relevant nodes

Search the knowledge graph file for nodes matching the user's query keywords:

- Grep `"name"` fields for keyword matches.
- Grep `"summary"` fields for semantic matches.
- Grep `"tags"` arrays for topic matches.
- Note the `id` values of all matching nodes.

### 4. Find connected edges

For each matched node ID, grep for that ID in the `edges` section to find:

- What it imports or depends on (downstream dependencies).
- What calls or imports it (upstream dependents).

This gives the 1-hop subgraph around the query.

### 5. Read layer context

Grep for `"layers"` to understand which architectural layers the matched nodes belong to. Layer assignments reveal whether matched components are in the core domain, infrastructure, testing, or other architectural tiers.

### 6. Answer the query

Using only the relevant subgraph, answer the user's query:

- Reference specific files, functions, and relationships from the graph.
- Explain which layer(s) are relevant and why.
- Be concise but thorough — link concepts to actual code locations.
- If the query does not match any nodes, say so and suggest related terms from the graph.
