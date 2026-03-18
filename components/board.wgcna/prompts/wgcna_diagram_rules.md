## Diagram Rules

### The Story the Diagram Must Tell

A WGCNA diagram's value is in showing the INTERPLAY between co-expression programs, not in listing them. Follow this hierarchy:

1. **The opposing programs** (anticorrelated module pairs) are the backbone. Every strong anti-correlation in the report MUST appear. This is THE insight: what goes up when something else goes down, and what that trade-off means biologically.

2. **The trait anchoring** shows directionality. Each module's green/red arrow to the trait node tells the reader which programs are activated vs. suppressed in the experimental condition. Together with the anti-correlation edges, this creates a readable narrative: "X activates Y (green) while suppressing Z (red); these programs are mutually exclusive."

3. **Hub gene validation.** detail_1 should list hub genes that a domain expert would recognize as relevant to the module's function. For human data: known oncogenes, tumor suppressors, druggable kinases. For plants: known stress-response or hormonal regulators. For C. elegans: characterized gene families (nhr, str, col, etc.).

1. **Cross-module convergence** (pathway nodes). Only worth showing when the same pathway appears enriched in two modules with opposite trait directions, this reveals a pathway that is being activated in one arm and suppressed in another.

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
- Prefer 15-25 edges total for readability.
- ALWAYS INCLUDE SOME NODES FOR BIOLOGICAL PROCESSES FOR EACH MODULE

Constraints:
- Use only entities present in the report text.
- Do not invent mechanisms, targets, or causal claims.

## Detail Fields

Every node type should have informative detail fields that help a biologist
interpret the diagram at a glance:

- **module nodes**:
  - `detail_1`: Top 3-5 hub genes, comma-separated (e.g. "TP53, BRCA1, MYC")
  - `detail_2`: Module biological function in 2-5 words (e.g. "cell cycle regulation")
- **phenotype nodes**:
  - `detail_1`: Key sample groups or conditions (e.g. "treated, control, 48h")
  - `detail_2`: Experimental context in 2-5 words (e.g. "T-cell activation time-course")
- **process nodes**:
  - `detail_1`: Top 2-3 genes driving this process, comma-separated
  - `detail_2`: Pathway or GO category (e.g. "oxidative phosphorylation")

All detail field content MUST come from the report text or supporting data
tables. Do not fill detail fields with general biological knowledge — if a
gene or condition is not mentioned in the report or data, leave the field
empty rather than infer it.
