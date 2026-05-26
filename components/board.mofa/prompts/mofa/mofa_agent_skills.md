# MOFA Deep Report — board rules + agent skills

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
- No raw weights, NES, or padj values.
- Together, the three bullets answer: what is the headline multi-omics
  programme, what is the most surprising finding, and what is the dominant
  cross-view pattern?

### Overview

- 1–2 paragraphs.
- Required content: experiment intent, organism, sample count, number of
  omics views, number of factors extracted, count selected for detailed
  reporting, signal landscape (use one of: "low-signal", "single-axis",
  "richly-resolved", "multi-axis").
- Mention the variance landscape **once**, here. Do not return to it.

### Main findings

Grouped by biological theme rather than by individual factor. Each theme
is one h3 subsection; each theme contains the factors that fit it as h4
subsubsections.

```
### Theme 1 — [biological framing in 2–4 words]
1–3 prose paragraphs introducing the theme.

#### Factor 1: [short theme label]
1–3 prose paragraphs.

#### Factor 3: [short theme label]
1–3 prose paragraphs.

### Theme 2 — [biological framing]
#### Factor 2: [short theme label]
...
```

**Per-factor length:**
- Strong: 2–3 paragraphs.
- Moderate: 1 paragraph.
- Weak: do NOT get a subsection; aggregate under Minor units.

**Per-factor content:**
- Open with the dominant pathway theme and trait association in plain
  biological terms.
- Top features: 4–6 named per strong factor, woven with view of origin
  and functional roles, ≤8 names per paragraph. Cross-view features
  (high weight in more than one view) called out explicitly.
- Pathways: themed plain-English clusters (no IDs).
- One numeric anchor per paragraph maximum.

### Minor units

- 1 paragraph (or omit entirely if no weak factors).
- 1–2 sentences per minor factor: name it, give the dominant view, note
  enrichment absence or weakness, optionally name 1–2 top features.
- No themes, no h3 subdivisions.

### Integrated findings

- 2–4 paragraphs.
- MUST use at least one of three rhetorical patterns: trade-off,
  feedback, or contrast.
- ONLY cite cross-factor relationships explicit in the data block.

### **Discussion**

**Single integrated narrative — no subsections.**

**Length:** 2–4 paragraphs. Bold heading: `## **Discussion**`.

**Inferring the experimental frame.** Before writing, read these data
sources in order to establish the goal of the experiment:

1. Experiment description (`EXPERIMENT:` line in the data block)
2. Phenotype / trait names (the trait columns referenced in factor
   correlations)
3. Sample group names from contrast definitions
4. Organism / species name + standard biological context for that organism
5. Any species-context fragment available in the prompt

From these, identify whether the experiment is about: (a) a disease /
clinical phenotype, (b) a drug or therapeutic intervention, (c) a
developmental / time-course process, (d) basic biology in a model
organism, (e) ecology / agriculture / aquaculture, or (f) something else.
Write the Discussion **accordingly** — disease framing only when the
data justifies it.

**What goes in:**
- Anchor each strong / moderate factor's biology in known literature
  surfaced via Phase B PubMed calls. Vancouver `[An]` citations — max
  2 per sentence.
- If the experiment IS clearly clinical / therapeutic (gate: clinical-
  disease metadata present in the description or traits AND organism is
  human or mouse AND signal is not low-signal), include translational
  discussion tracing implications back to specific factor patterns.
  Cite literature for any drug / pathway / clinical claim.
- If the experiment is NOT clinical / therapeutic, skip translational
  discussion entirely — do not narrate the absence. **Do NOT write
  parenthetical "reasons" like `(reason: no clinical metadata)`** —
  those are template artefacts, never acceptable in prose.

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
- Every `[An]` in prose appears as an entry here; every entry here is
  cited at least once in prose.

Example of one well-formed entry:

> [A1] Argelaguet R, Velten B, Arnol D, et al. Multi-Omics Factor
> Analysis—a framework for unsupervised integration of multi-omics data
> sets. Mol Syst Biol. 2018 Jun 20;14(6):e8124.
> doi:10.15252/msb.20178124. PMID: 29925568.

MOFA methods literature (`[M1]`, `[M2]`) is appended deterministically by
the pipeline under a SEPARATE `## **Methods literature**` heading —
do NOT generate it yourself, and do NOT include `[M*]` in prose.

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
the dataset to surface specific factors, pathway enrichments, and top
features that anchor the biological theme. Derive query keywords from
the experiment description, phenotypes, organism, and the dominant
pathway theme of each strong factor — do not run generic unfocused
queries.

**Required calls (per dataset, before literature):**

Always start with `query_mofa(what="summary")` to see which factors carry
the strongest signal and which views drive them. Then drill on the
top-ranked factors with `query_mofa_factor`. `query=` accepts a JSON
array of short keywords with OR-semantics — prefer 2–3 targeted
small-keyword arrays over one long natural-language string.

```
query_mofa(what="summary")
query_mofa_factor(factor="<top-factor>", what="features",
                  query=["JAK","STAT","cytokine"], top_n=15)
query_mofa_factor(factor="<top-factor>", what="pathways",
                  query=["inflammation","cytokine","NF-kB"])
```

**Suggested keyword themes by experimental domain:**

- Clinical / pharmaceutical / drug discovery: `DRUG`, `TREATMENT`,
  `INHIBITOR`, `DISEASE`, `DISGENET`, `GWAS`, `CANCER` + clinical theme.
- Basic biology / model-organism (mouse/rat): same keywords with
  translational caveat.
- Invertebrate / microbial / non-mammalian: `RIBOSOME`, `METABOLISM`,
  `MITOCHONDRIA`, `CELL_CYCLE`, `STRESS` + organism-specific terms;
  reject MSigDB disease hits as cross-species artifact.
- Agriculture / livestock / aquaculture: `GROWTH`, `YIELD`, `STRESS`,
  `IMMUNE`, `MUSCLE`, `REPRODUCTION`, `HORMONE`, `ABIOTIC_STRESS`.
- Ecology: `ADAPTATION`, `STRESS`, `IMMUNE_DIVERSITY`,
  `HOST_PATHOGEN` + ortholog-aware caveat.

## Phase B — Literature search (PubMed)

**Budget by signal density.**

- Low-signal dataset (no factor clears 5 significant pathways AND
  max |weight| < 0.5): skip the literature phase entirely.
  Discussion becomes a single sentence.

- All other datasets: between **5 and 15 PubMed calls**. Scale with
  dataset richness — more multi-view factors, more biological themes,
  and more cross-factor structure all justify more calls.

**Allocate calls to:**
1. One call per strong factor's biological theme.
2. One call per cross-view feature (high weight in ≥2 views).
3. One call per strongly correlated factor-trait pair.
4. One call for an organism caveat if non-mammalian AND the report
   makes any translational claim.

### How to query well

`search_context` wraps PubMed Entrez and accepts the full Entrez syntax:
field tags (`[ti]`, `[tiab]`, `[mesh]`, `[ORGN]`), boolean operators
(`AND` / `OR` / `NOT`), and date ranges. A scoped query of 1–2 concepts
beats an unscoped query of 5.

- Feature scoped to title + concept in title/abstract:
  `query="IL6[ti] AND JAK-STAT[tiab]"`
- MeSH concept filtered to organism (required for non-human/non-mouse):
  `query="multi-omics[mesh] AND Salmo salar[ORGN]"`

When 0 hits returns: decompose to one concept, add a field tag to scope,
retry **with a structurally different query**. If still 0, accept that
the literature is absent and move on.

### Query derivation, in priority order

1. Cross-view feature symbols (verbatim, not aliases).
2. The 2-word biological core of a top enrichment term.
3. Organism + tissue/condition phrase from sample metadata, scoped via
   `[ORGN]` when non-human/non-mouse.

Stop at the first level that produces a focused query; do not pile them.

### Title-gate before `read_context`

Read the search-result title first; only call `read_context` when the
title matches BOTH the specific concept AND the biological context.
When you do call it, batch every section you'll need from that id in
ONE call (e.g. `sections=["abstract","citation"]`) — do not split.

## Citation discipline

- **Cite to extend, not duplicate.** If Phase A already established a
  claim with `padj < 1e-10`, do not also cite literature for it.
- **State confirming citations directly:**
  *"JAK-STAT drives inflammatory cytokine response [A1]."*
- **Hedge context-only ones:**
  *"This may extend to the chronic phase [A2]."*
- **Cap (per sentence):** ≤2 citations.
- **Minimum (per report):** ≥5 citations total (when the literature
  phase ran). Increase with dataset complexity.
- **Bibliography format:** cite verbatim from the `citation` block
  returned by a single `read_context(sections=["abstract","citation"],
  id=...)` call. Do not paraphrase. Numbering follows order of first
  citation in prose, starting `[A1]`.

## Tool surface (use only these)

- `query_mofa` (Phase A only — do not free-explore)
- `query_mofa_factor` (Phase A only — drill on top factors)
- `query_gene_info` (sparingly, for unfamiliar features)
- `search_context`, `read_context` (Phase B literature only)

Do NOT call `query_de`, `query_drugs`, `query_pathways`, `query_wgcna`,
etc. — they are not on this agent's tool surface.

**Stop condition: report drafted OR 30 tool calls reached, whichever comes first.**
