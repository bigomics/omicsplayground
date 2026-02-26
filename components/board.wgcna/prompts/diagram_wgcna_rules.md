## Diagram Rules

Classify nodes into:
{{node_names}}

Classify links into:
{{link_names}}

Required edge semantics:
- `module -> phenotype` with relation `positive` or `negative`
- `module -> process` with relation `positive` or `negative`

Node labels:
- Format: `[ModuleID]: biology` (2-5 words describing function)
- Example: `[turquoise]: cell cycle regulation`
- Do NOT use full gene names in labels -- genes go in the `genes` field only

Edge strength calibration:
- Strong co-expression or known direct interaction -> strength 0.8-1.0
- Moderate relationship or inferred -> strength 0.4-0.7
- Weak or speculative -> strength 0.1-0.3

Prioritization:
- Hub nodes (high connectivity) should be phenotype or key process nodes.
- Prefer 6-12 edges total for readability.

Constraints:
- Use only entities present in the report text.
- Do not invent mechanisms, targets, or causal claims.
- All text fields must use plain ASCII characters only.
