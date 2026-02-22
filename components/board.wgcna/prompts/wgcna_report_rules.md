## WGCNA-Specific Reporting Rules

### Module Classification (AUTHORITATIVE — do not override)

The input data includes a **Module Signal Classification** section. This
classification is pre-computed and **MUST be followed exactly**. Do NOT
reclassify modules based on your own judgment of their enrichment or
correlation strength.

- **Strong signal**: reported as individual `### MEcolor: Name` subsections
- **Moderate signal**: also reported as individual subsections but kept
  compact (2-3 sentences: eigengene pattern, enrichment status, 1-2 hub
  genes). A module listed as "Moderate" MUST get its own subsection even
  if it has zero significant enrichment terms.
- **Weak signal**: aggregated into `### Weak Modules`
- **Artifact-flagged**: reported as individual subsections BUT with an
  explicit artifact warning paragraph at the start

If a module is classified as Moderate but has no enrichment, write its
subsection focusing on eigengene pattern, hub genes, and the absence of
enrichment — do NOT move it to Weak Modules.

### Module Heading Format

Every strong/moderate module subsection MUST use this heading format:

```
### MEblue: Cell Cycle Machinery
```

That is: `### MEcolor: Short Biological Name` where the biological name
is 2-4 words derived from the dominant enrichment theme. Do NOT use long
compound names. Examples:
- GOOD: "Translational Scaling", "Quiescence Program", "Early Transient Response"
- BAD: "Ribosome Biogenesis and Amino Acid Transport Module"

### Weak Modules Subsection

All modules classified as **Weak** in the Module Signal Classification
are aggregated into a single `### Weak Modules` subsection. If there are
no Weak modules, **omit this subsection entirely** — do not write it.

Format as a compact paragraph per module (2-3 sentences max):
- Name, gene count, and eigengene pattern
- State absence of enrichment explicitly
- Note hub genes briefly and any possible interpretation
- Keep speculative language minimal

Example: "**MEred** (69 genes) shows a sharp transient spike at 12h
(eigengene: -0.12 → +0.52 → -0.13, r = +0.99 with 12h contrast) but
lacks significant pathway enrichment (0 of 1000 terms at q < 0.05).
Hub genes *SERAC1* and *HIP1* suggest membrane remodeling."

### Artifact-Flagged Modules

If a module is classified as strong/moderate but flagged as a potential
artifact (see the ⚠ flag in Module Signal Classification), it still
gets its own subsection but MUST begin with a warning paragraph:

"**⚠ Potential artifact.** This module's hub genes include [list 2-3
specific genes] which are [reason, e.g., classical plasma proteins not
expected in sorted cell preparations]. The enrichment themes [list 1-2]
further support this interpretation. The biological findings below should
be treated with caution."

Then proceed with the normal narrative.

### Strong Module Content (per subsection)

Write 1-3 paragraphs per strong module covering:
- The eigengene temporal/condition profile (with values)
- The dominant enrichment themes with [n] references
- 3-5 hub genes woven into the narrative (with MM, logFC)
- How the enrichment and eigengene patterns connect
- Enrichment overlap counts when available (see Enrichment Overlap below)

Moderate modules follow the same structure but kept compact (2-3 sentences).

### Eigengene Interpretation

Eigengene = first principal component of module expression.
- Positive value = module genes tend to be upregulated
- Negative value = module genes tend to be downregulated

Describe patterns narratively with values:
  GOOD: "The MEblue eigengene rises steadily from resting (-0.22) through
  48h (+0.27) to peak at 72h (+0.33)."
  BAD: "MEblue is an activation module."

### Hub Gene Reporting

Weave 3-4 hub genes per strong module into the narrative (fewer for
moderate modules). Do NOT list all hub genes from the input — select the
most representative ones. For each:
- Gene symbol in *italics*
- Module membership (MM) — how representative of the module
- logFC vs the strongest trait contrast
- Brief functional context (from the "Known function" column)

Integrate naturally: "*CDK1* (MM = 0.97, logFC = +7.1), a cyclin-dependent
kinase central to mitotic entry, leads this module alongside replication
factors *MCM2* (MM = 0.95) and *RRM1* (MM = 0.92)."

Do NOT create bullet lists of hub genes. Weave them into the narrative.

### Grey Module

The grey module contains unassigned genes. Mention its size once in the
introduction or summary (e.g., "149 genes were not assigned to any
module"). Do not interpret it as a biological program.

### Enrichment References

- Reference enrichment terms via [n] brackets in the narrative
- Use plain-English descriptions in prose: "DNA replication [1,2]"
  not "REACTOME_DNA_REPLICATION [1]"
- The Data References section at the end maps [n] to exact term names
- When >100 terms are significant, summarize by theme — do not attempt
  to reference them all
- When 0 terms are significant, state it explicitly

### Cross-Module Patterns

After all module subsections (before ## Discussion), include one paragraph
describing relationships between modules:
- Temporal ordering (which module eigengenes peak first, second, etc.)
- Opposing pairs (anti-correlated eigengenes)
- Sequential handoffs (e.g., translation peaks before proliferation)

Support every cross-module claim with eigengene values from both modules.
When MODULE-MODULE EIGENGENE CORRELATIONS data is provided, use ONLY
the r values listed there. Do NOT estimate, infer, or fabricate
correlation values that are not in the provided table. If a module pair
is not listed (meaning |r| < 0.70), do not claim they are correlated.

### Enrichment Overlap (REQUIRED when present)

When the input includes enrichment overlap counts (e.g., "85/147"),
you MUST cite at least one overlap count per strong module:
"85 of 147 DNA metabolic process genes are in this module [1]".
This quantifies how much of a pathway the module captures and is
more informative than the q-value alone.

### Gene Family Enrichment (REQUIRED when present)

When "Gene family enrichment:" lines appear in a module's input data,
you MUST mention the top 1-2 families for strong modules in the
narrative. Example: "the module contains 23 of 242 kinases and 17 of
119 CD molecules, consistent with a signaling-heavy composition."
Gene families provide an independent validation layer beyond pathway
enrichment — they confirm what types of proteins dominate the module.
For the Discussion section, compare gene family representation across
modules when the same family appears in multiple modules.
