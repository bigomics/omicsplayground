## Diagram Rules

Classify nodes into:
{{node_names}}

Classify links into:
{{link_names}}

Required edge semantics:
- `contrast -> drug` with relation `opposes` or `mimics`
- `drug -> moa` with relation `has_moa`
- `drug -> target` with relation `hits`

### Contrast coverage

Every contrast discussed in the **Main Findings** section of the report MUST
appear as a `contrast` node in the diagram. If the report's cross-contrast
section identifies additional contrasts with distinct pharmacological themes,
include those too. A diagram that covers only a subset of the report's
contrasts is incomplete.

### Primary extraction rule for sample-to-drug edges
- Prioritize explicit polarity phrases in the report text such as
  `<sample/contrast> is opposed by <drug>` and `<sample/contrast> is mimicked by
  <drug>`.
- Create `contrast -> drug` edges only from explicit polarity statements in the
  report text.
- Do not infer `contrast -> drug` polarity from MOA links, target links, or
  general narrative context alone.

### Prioritization
- Include the strongest opposing drugs first.
- Include mimicking drugs only when they add mechanistic context.
- Keep the graph readable: prefer fewer, evidence-supported nodes.

### Constraints
- Use only entities present in the report text or supporting data tables.
- Do not invent mechanisms, targets, or causal claims.

## Detail Fields

Every node type should have informative detail fields that help a biologist
interpret the diagram at a glance:

- **contrast nodes**:
  - `detail_1`: Key upregulated or dysregulated pathways (e.g. "mTOR, PI3K-Akt")
  - `detail_2`: Experimental context in 2-5 words (e.g. "48h T-cell activation")
- **drug nodes**:
  - `detail_1`: Target genes or key biomarkers from the report or data tables,
    comma-separated. If no specific targets are mentioned, leave empty.
  - `detail_2`: Mechanism of action or drug class (e.g. "HDAC inhibitor")
- **moa nodes**:
  - `detail_1`: Representative drugs sharing this MOA (e.g. "vorinostat, panobinostat")
  - `detail_2`: Biological effect in 2-5 words (e.g. "chromatin remodeling")
- **target nodes**:
  - `detail_1`: Gene symbol(s) if mentioned in the report (e.g. "CDK4, CDK6")
  - `detail_2`: Functional role (e.g. "cell cycle checkpoint")

All detail field content MUST come from the report text or supporting data
tables. Do not fill detail fields with general pharmacological knowledge — if
a gene, target, or MOA is not mentioned in the report or data, leave the
field empty rather than infer it.
