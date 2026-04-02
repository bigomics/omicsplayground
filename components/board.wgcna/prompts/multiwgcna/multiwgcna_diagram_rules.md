## Diagram Rules

### The Story the Diagram Must Tell

A MultiWGCNA diagram should show how layer-specific modules converge into a coherent multi-omics program.

1. Cross-layer module relationships are the backbone. Every strong cross-layer correlation in the report should appear.
2. Trait anchoring provides directionality. Modules should connect to phenotype nodes with positive or negative relationships when supported by the report.
3. Local biological interpretation still matters. Each module should connect to at least one process node grounded in its enrichment results.
4. Layer context should be visible so the reader can distinguish transcript, protein, metabolite, or other modality-specific modules.

Classify nodes into:
{{node_names}}

Classify links into:
{{link_names}}

Required edge semantics:
- `module -> phenotype` with relation `positive` or `negative`
- `module -> process` with relation `positive` or `negative`
- `module -> module` with relation `cross_layer` when modules from different layers are strongly correlated
- `layer -> module` with relation `association`

Node labels:
- Module labels: `[layer] module: biology`
- Layer labels: short datatype names such as `gx`, `px`, `mx`
- Keep labels short and readable

Constraints:
- Use only entities present in the report text or supporting data
- Do not invent mechanisms or causal claims
- Prefer 15-25 edges total for readability

## Detail Fields

- **layer nodes**:
  - `detail_1`: short datatype description
  - `detail_2`: number of modules represented
- **module nodes**:
  - `detail_1`: top 3-5 hub features
  - `detail_2`: module biology in 2-5 words
- **phenotype nodes**:
  - `detail_1`: key conditions or sample groups
  - `detail_2`: experimental context
- **process nodes**:
  - `detail_1`: top 2-3 features driving the process
  - `detail_2`: pathway or functional category
