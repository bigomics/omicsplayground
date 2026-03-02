## Diagram Rules

Classify nodes into:
{{node_names}}

Classify links into:
{{link_names}}

Required edge semantics:
- `phenotype -> process` with relation `association` or `risk`
- `process -> module` with relation `supports` or `association`

Prioritization:
- Include phenotype/context node(s) and major process nodes first.
- Include top hub genes and mark Steiner/connector roles where available.
- Keep the graph readable: prefer fewer, evidence-supported nodes.

Constraints:
- Use only entities present in the report text.
- Do not invent mechanisms, targets, or causal claims.
- Keep node labels short and legible.

## Detail Fields

- `detail_1`: Hub or Steiner genes, comma-separated (e.g. "TP53, AKT1")
- `detail_2`: Biological process or pathway name (e.g. "apoptosis signaling")
