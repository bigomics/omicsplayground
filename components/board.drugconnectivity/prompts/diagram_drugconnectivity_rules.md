## Diagram Rules

> **Node types for this board** (override base instructions):
> Do NOT use `module`, `phenotype`, or `process`.

Classify nodes into:
{{node_names}}

Classify links into:
{{link_names}}

Required edge semantics:
- `contrast -> drug` with relation `opposes` or `mimics`
- `drug -> moa` with relation `has_moa`
- `drug -> target` with relation `hits`

Prioritization:
- Include the strongest opposing drugs first.
- Include mimicking drugs only when they add mechanistic context.
- Keep the graph readable: prefer fewer, evidence-supported nodes.

Constraints:
- Use only entities present in the report text.
- Do not invent mechanisms, targets, or causal claims.
