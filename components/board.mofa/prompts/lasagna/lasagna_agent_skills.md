# LASAGNA Deep Report — board rules + agent skills

Loaded as `board_rules` only by the agentic Deep Report path. Sub-section
ordering and shared structure match the single-shot report rules; the
Discussion / Conclusion / Bibliography sections and the trailing Agent
skills block are Deep-Report-specific.

Total report length: 1500–2500 words.

## Section requirements

### Highlights

- Exactly **3 bullets**. No more, no less.
- Each bullet ≤25 words.
- Declarative; no hedging.
- No raw centrality, ρ, or edge-weight values.
- Together, the three bullets answer: what is the headline cross-layer
  programme, what is the most surprising finding, and what is the
  dominant inter-modality bridge?

### Overview

- 1–2 paragraphs.
- Required content: experiment intent, organism, sample count, contrast,
  number of omics layers, node and edge counts, count of cross-layer
  modules, count selected for detailed reporting, signal landscape
  ("fragmented" / "single-layer-dominated" / "balanced" /
  "richly-bridged").
- Mention the inter-layer edge fraction once, here.

### Main findings

Grouped by biological programme rather than by individual module ID.
Each programme is one h3 subsection; modules that fit a programme are
h4 subsubsections.

**Per-module length:**
- Strong: 2–3 paragraphs.
- Moderate: 1 paragraph.
- Weak: do NOT get a subsection; aggregate under Minor units.

**Per-module content:**
- Open with the biological programme in plain language.
- Hub nodes: 4–6 named per strong module, woven with layer of origin
  and functional role, ≤8 names per paragraph. Cross-layer hubs called
  out explicitly.
- Layer participation: state whether multi-layer or single-layer
  dominated.
- One numeric anchor per paragraph maximum.

### Minor units

- 1 paragraph (or omit entirely if no weak modules).
- 1–2 sentences per minor module.

### Integrated findings

- 2–4 paragraphs.
- MUST use at least one of three rhetorical patterns: trade-off,
  feedback, or contrast.
- ONLY cite cross-module bridges explicit in the data block.

### **Discussion**

**Single integrated narrative — no subsections.**

**Length:** 2–4 paragraphs. Bold heading: `## **Discussion**`.

**Inferring the experimental frame.** Before writing, read these data
sources in order to establish the goal of the experiment:

1. Experiment description (`EXPERIMENT:` line in the data block)
2. Contrast definition (the `CONTRAST:` line and contrast block)
3. Layer names — which modalities are present
4. Organism / species name + standard biological context for that organism
5. Any species-context fragment available in the prompt

From these, identify whether the experiment is about: (a) a disease /
clinical phenotype, (b) a drug or therapeutic intervention, (c) a
developmental / time-course process, (d) basic biology in a model
organism, (e) ecology / agriculture / aquaculture, or (f) something
else. Write the Discussion **accordingly** — disease framing only when
the data justifies it.

**What goes in:**
- Anchor each strong / moderate module's biology in known literature
  surfaced via Phase B PubMed calls. Vancouver `[An]` citations — max
  2 per sentence.
- If the experiment IS clearly clinical / therapeutic (gate: clinical-
  disease metadata in description AND organism is human or mouse AND
  signal is not fragmented), include translational discussion tracing
  implications back to specific module patterns.
- If the experiment is NOT clinical / therapeutic, skip translational
  discussion entirely.

**Style:** Hedge claims that go beyond the data; state established
biology directly. Do not introduce findings absent from Main findings
or Integrated findings above.

### **Conclusion**

- **1 paragraph**.
- Heading in **bold**: `## **Conclusion**`.
- Integrative summary; no new claims; no numbers.

### Bibliography

Bibliography is the final section of the report. Heading:
`## **Bibliography**`.

- Number entries `[A1]`, `[A2]`, … in order of first citation in prose.
- Each entry: the `citation` block returned verbatim by `read_context` —
  no paraphrasing, no reformatting.

---

# Agent skills

This block is loaded only by the agentic Deep Report path. The single-shot
Report path has no tool access; do not interpret these instructions as
applicable to it.

**Hard limits**
- Tool-call budget: **≤30** per report (Phase A + Phase B + lookups).
- Phase B duplicate-shape retries: 0 (any retry must change structure,
  not just swap synonyms).

You have access to a set of tools that let you query the active dataset
and search PubMed. Use them in two phases.

## Phase A — Targeted dataset discovery (BEFORE literature)

Before any PubMed call, run a small set of query-filtered passes against
the dataset to surface specific pathways, differential features, and
correlations that anchor each module's biological programme. Derive
query keywords from the experiment description, contrast, organism,
and the hub-node identities for the lead module — do not run generic
unfocused queries.

**Required calls (per dataset, before literature):**

```
query_de(contrast="<active-contrast>", query=["<hub1>","<hub2>","<hub3>"])
query_pathways(contrast="<active-contrast>", query=["<top-theme>"])
query_correlation(target="<top-hub>", query=["<related-process>"])
```

**Suggested keyword themes by experimental domain:**

- Clinical / pharmaceutical: `DRUG`, `TREATMENT`, `INHIBITOR`, `DISEASE`,
  `DISGENET`, `GWAS`, `CANCER` + clinical theme.
- Basic biology / model-organism: same keywords with translational
  caveat.
- Invertebrate / microbial / non-mammalian: `RIBOSOME`, `METABOLISM`,
  `MITOCHONDRIA`, `CELL_CYCLE`, `STRESS` + organism-specific terms;
  reject MSigDB disease hits as cross-species artifact.
- Agriculture / livestock / aquaculture: `GROWTH`, `YIELD`, `STRESS`,
  `IMMUNE`, `MUSCLE`, `REPRODUCTION`, `HORMONE`, `ABIOTIC_STRESS`.

## Phase B — Literature search (PubMed)

**Budget by signal density.**

- Fragmented network (data block flagged `fragmented` or fewer than 5
  inter-layer edges): skip the literature phase entirely.
  Discussion becomes a single sentence.

- All other datasets: between **5 and 15 PubMed calls**. Scale with
  network richness — more cross-layer modules, more layer diversity,
  and more strong modules all justify more calls.

**Allocate calls to:**
1. One call per strong module's biological theme.
2. One call per cross-layer hub (hub with at least one inter-layer edge).
3. One call per cross-module bridge (shared hub or strong inter-module
   edge).
4. One call for an organism caveat if non-mammalian AND the report
   makes any translational claim.

### How to query well

`search_context` wraps PubMed Entrez and accepts the full Entrez syntax:
field tags (`[ti]`, `[tiab]`, `[mesh]`, `[ORGN]`), boolean operators
(`AND` / `OR` / `NOT`), and date ranges. A scoped query of 1–2 concepts
beats an unscoped query of 5.

- Hub scoped to title + concept in title/abstract:
  `query="CDK2[ti] AND replication licensing[tiab]"`
- MeSH concept filtered to organism (required for non-human/non-mouse):
  `query="multi-omics[mesh] AND Mus musculus[ORGN]"`

When 0 hits returns: decompose to one concept, add a field tag to scope,
retry **with a structurally different query**. If still 0, accept that
the literature is absent and move on.

### Title-gate before `read_context`

Read the search-result title first; only call `read_context` when the
title matches BOTH the specific concept AND the biological context.
When you do call it, batch every section you'll need from that id in
ONE call (e.g. `sections=["abstract","citation"]`).

## Citation discipline

- **Cite to extend, not duplicate.** If Phase A already established a
  claim with a strong correlation, do not also cite literature for it.
- **State confirming citations directly** for established biology.
- **Hedge context-only ones.**
- **Cap (per sentence):** ≤2 citations.
- **Minimum (per report):** ≥5 citations total (when the literature
  phase ran). Increase with dataset complexity.

## Tool surface (use only these)

- `query_de` (Phase A — differential expression context for hub nodes)
- `query_pathways` (Phase A — pathway enrichment for the contrast)
- `query_correlation` (Phase A — cross-feature correlations supporting
  inter-layer edges)
- `query_gene_info` (sparingly, for unfamiliar hubs)
- `search_context`, `read_context` (Phase B literature only)

Do NOT call `query_mofa`, `query_drugs`, `query_wgcna`, etc. — they are
not on this agent's tool surface.

**Stop condition: report drafted OR 30 tool calls reached, whichever comes first.**
