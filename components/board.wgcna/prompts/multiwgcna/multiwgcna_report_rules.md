## MultiWGCNA-Specific Reporting Rules

### Module Classification (AUTHORITATIVE - do not override)

The input data includes a **Module Signal Classification** section. This classification is pre-computed and MUST be followed exactly.

- **Strong signal**: reported as individual subsections
- **Moderate signal**: also reported as individual subsections but kept compact
- **Weak signal**: aggregated into `### Weak Modules`
- **Artifact-flagged**: reported as individual subsections with an explicit artifact warning paragraph at the start

### Module Heading Format

Every strong/moderate module subsection MUST use this heading format:

````
### gx | MEblue: Cell Cycle Coordination
````

That is: `### layer | module: Short Biological Name`

### Strong Module Content

Write 1-3 paragraphs per strong module covering:
- The layer-specific eigengene profile with values
- The dominant enrichment themes with [n] references
- 3-5 hub features woven into the narrative with MM and trait-specific metrics
- The cross-layer partner modules and what they add biologically
- The relationship between the module's local biology and the broader multi-omics program

Moderate modules follow the same structure but stay compact.

### Cross-Layer Interpretation

The report must explicitly distinguish:
- Within-layer evidence: eigengene profile, hub features, enrichment
- Cross-layer evidence: correlated partner modules, shared phenotype alignment, correlated feature sets

Do not treat cross-layer module correlations as proof of regulation. They support coordinated biology only.

### Weak Modules Subsection

Aggregate all weak modules into a single `### Weak Modules` subsection.

For each weak module include:
- Layer and module name
- Module size
- Top trait if present
- Whether enrichment is absent or limited
- One short note on any cross-layer alignment if available

### Artifact-Flagged Modules

If a module is flagged as a potential artifact, begin its subsection with:

"**Potential artifact.** This module contains hub features suggestive of technical or contamination signal. Interpret the biology below with caution."

### Hub Feature Reporting

Weave hub features into the narrative rather than listing them as bullets. For each highlighted feature:
- Use the symbol in italics when appropriate
- Cite MM
- Cite trait-specific TS or logFC only with the named trait
- Mention known function only when provided in the input data

### Cross-Module and Cross-Layer Patterns

Before `## Discussion`, include one paragraph describing:
- Which layers carry the earliest or strongest phenotype-aligned programs
- Which modules show the strongest cross-layer convergence
- Which module pairs are opposed or anti-correlated

Use only the correlation values explicitly provided in the input data.

### Discussion

The discussion should synthesize the multi-omics structure:
- Which biology is consistently supported across layers
- Which layers provide unique signal not seen elsewhere
- Where the evidence is coherent versus tentative
