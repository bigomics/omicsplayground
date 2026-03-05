## What This Analysis Reveals

The Prize-Collecting Steiner Forest (PCSF) algorithm maps differentially
expressed genes onto a protein–protein interaction backbone (STRING / GRAPHITE)
and finds the optimal subnetwork balancing inclusion of high-prize nodes against
connection cost. The result is a **network portrait of the biological state**:
hub genes with high centrality are topological drivers; Steiner nodes — genes
NOT in the original differential set but recruited as connectors — are **hidden
regulators** that the expression data alone would miss.

The report articulates this portrait in biological language. The network
topology is evidence. The biological programmes are the subject.

When `## Dataset Context` is present in the input data, treat it as the
authoritative framing for the experiment and prefer it over generic
species-level defaults.

## Analytical Unit: the Functional Module

Each `### heading` in Main Findings names a **biological programme or functional
module** revealed by the network topology — not a contrast ID, not a gene list.

Derive the heading from the dominant biological theme of the module and the
experimental context:

```
GOOD: ### late_vs_ctrl: Mitotic Drive via Replication Hub
GOOD: ### treated_vs_ctrl: Coagulation–Complement Axis with Metabolic Shutdown
BAD:  ### contrast_A_vs_B
BAD:  ### Top Hub Genes
```

**Multi-contrast data with temporal or ordered conditions**: you MUST group
contrasts into 2–4 functional phases by which hub modules and pathway themes
dominate — do NOT write one subsection per contrast. Inspect the hub overlap
and pathway enrichment across contrasts first: find which modules are stable
(core programme), which emerge at specific time points, and which fade. Build
phases from these patterns. Name each phase by its biology, not by contrast IDs.

**Multi-contrast data without clear ordering**: identify the 2–4 dominant
functional modules that emerge across contrasts; each becomes one subsection.

**Single-contrast data**: one subsection per dominant functional module (usually
1–2). Use `### Minor Signals` only if a weak secondary module exists.

Use the **Contrast Ranking** tier (strong / moderate / weak / data-limited) to
calibrate depth:
- **Strong / Moderate**: full narrative treatment
- **Weak / Data-limited**: aggregate into `### Minor Signals` with explicit
  caveats about network sparsity or missing pathway support

Omit `### Minor Signals` entirely if all contrasts are strong or moderate.

## Per-Unit Content: Three Narrative Beats (always prose)

**Beat 1 — The biological programme (lead sentence)**
State what functional module the network topology reveals. Lead with the
biology, not the gene names.

> GOOD: "The network centres on a replication-licensing module whose hub
> architecture links nucleotide synthesis to checkpoint release."
> BAD: "12 hub genes were identified with centrality scores above 4.2."

**Beat 2 — The network evidence**
Name the hub genes that anchor the module. For each, weave centrality, logFC
direction, and functional annotation into prose. When the same hub appears
across multiple contrasts, say so explicitly — this is stronger evidence than
a single-contrast signal.

**Steiner nodes**: If the data reports Steiner nodes (steiner_fraction > 0),
highlight them prominently — these are the PCSF's unique discovery: regulators
invisible to standard differential expression. If steiner_fraction = 0 for a
contrast, state this once and do NOT describe any gene in that contrast as a
"hidden regulator" or "not in the original differential set" — every node in a
zero-Steiner network IS a differentially expressed terminal by definition.

> GOOD (Steiner present): "The algorithm recruited *GeneX* as a Steiner
> connector (centrality = 5.8), linking the proliferative hubs to a
> transcription factor not itself differentially expressed."
> GOOD (Steiner absent): "The network contains no Steiner nodes, so the
> topology is fully explained by the measured expression changes."
> BAD: "No Steiner nodes were recruited... *GeneY* provides a hidden link
> not differentially expressed" — this is a self-contradiction.

**Beat 3 — Pathway convergence and confidence**
Does pathway enrichment corroborate the module? Name the enriched pathway(s)
and their overlap size. State confidence explicitly when:
- Pathway signal is FALSE (no enrichment detected)
- Network contains < 100 nodes or zero Steiner nodes
- Hub logFC values are all < 0.5 (modest effect sizes)

## Cross-Unit Synthesis Paragraph (REQUIRED when > 1 unit)

After all unit subsections, one paragraph with no heading, describing:
- How the network architecture evolves across conditions or time
- Which hub modules are stable (core programme) vs. condition-specific
- Whether Steiner nodes bridge otherwise disconnected modules
- The overall biological narrative that emerges from integrating all contrasts

Every cross-unit claim must be grounded in values from the data.

## Discussion

1–2 paragraphs interpreting the overall network portrait:
- What biological programme does the network most strongly represent?
- What do the Steiner nodes reveal that differential expression alone could not?
- Translational caveats once, briefly: PPI backbone bias (well-studied proteins
  over-represented), organism-specific annotation gaps, network parameter
  sensitivity

## Word Limit

800–1200 words total across Main Findings + Discussion + Conclusion.

## What NOT to Write

- Do NOT create section headings named "Top Hub Drivers", "Contrast Ranking
  and Confidence Tiers", "Decision Snapshot", "Validation Next Steps",
  "Limitations and Risk Notes" — these describe tables and checklists, not
  biological narratives
- Do NOT reproduce hub tables or pathway tables — the researcher already has
  them in the interface
- Do NOT write bullet lists of gene names with centrality scores; weave them
  into prose
- Do NOT enumerate every contrast separately when a thematic or temporal
  grouping captures the pattern more clearly
- Do NOT repeat logFC or centrality values for every gene — select the 3–5
  most telling values per module

## Hard Constraints

- Use only facts explicitly present in the input data.
- Use only `[n]` references supplied under `## Referenceable Evidence`.
- Do NOT invent new reference numbers, renumber existing ones, or repurpose a
  reference number for a different claim.
- If a contrast reports `steiner = 0`, do NOT describe any node in that
  contrast as a hidden regulator, recruited connector, or non-differentially
  expressed bridge.
- If a contrast reports `pathway signal = no`, state that explicitly and keep
  confidence modest.
- Do NOT synthesize a numeric range across contrasts unless the underlying
  values are explicit in the input data. Cross-contrast summaries must be
  anchored by at least one supplied reference.

## Data References

Use `[n]` citations to trace quantitative claims back to the input data.
References should cite **pathway enrichment terms, cross-contrast patterns,
and network-level metrics** — NOT individual gene centrality/logFC values
(those are already woven into the prose and do not need a separate citation).

The input data provides `## Referenceable Evidence` as a pre-numbered catalog.
The final `## Data References` section MUST be a subset copy of the cited
entries from that catalog, preserving wording and numbering exactly.

Up to 15 references. Every `[n]` in the narrative MUST have a matching entry
in `Referenceable Evidence` and in the final `## Data References` section.
Do NOT cite a number that is absent from `Referenceable Evidence`.
