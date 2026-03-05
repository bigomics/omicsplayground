## Diagram Rules

Classify nodes into:
{{node_names}}

Classify links into:
{{link_names}}

Required edge semantics:
- `contrast -> drug` with relation `opposes` or `mimics`
- `drug -> moa` with relation `has_moa`
- `drug -> target` with relation `hits`

Primary extraction rule for sample-to-drug edges:
- Prioritize explicit polarity phrases in the report text such as
  `<sample/contrast> is opposed by <drug>` and `<sample/contrast> is mimicked by
  <drug>`.
- Create `contrast -> drug` edges only from explicit polarity statements in the
  report text.
- Do not infer `contrast -> drug` polarity from MOA links, target links, or
  general narrative context alone.

Prioritization:
- Include the strongest opposing drugs first.
- Include mimicking drugs only when they add mechanistic context.
- Keep the graph readable: prefer fewer, evidence-supported nodes.

Constraints:
- Use only entities present in the report text.
- Do not invent mechanisms, targets, or causal claims.

## Detail Fields

- `detail_1`: Target genes or key biomarkers, comma-separated. Leave empty for non-gene nodes.
- `detail_2`: Mechanism of action or drug class (e.g. "HDAC inhibitor", "kinase pathway")
