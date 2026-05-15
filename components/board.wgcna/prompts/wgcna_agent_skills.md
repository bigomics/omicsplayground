# WGCNA Deep Report — board rules + agent skills

Loaded as `board_rules` only by the agentic Deep Report path. Sub-section
ordering and shared structure match the single-shot report rules; the
Discussion / Conclusion / Bibliography sections and the trailing Agent
skills block are Deep-Report-specific.

Total report length: 2000–3000 words.

## Section requirements

### Highlights

- Exactly **3 bullets**. No more, no less.
- Each bullet ≤25 words.
- Declarative; no hedging.
- No `r`, `q`, `MM` values.
- Together, the three bullets answer: what is the dataset's headline
  biology, what is the most surprising finding, and what is the dominant
  cross-module pattern?

### Overview

- 1–2 paragraphs.
- Required content: experiment intent, organism, sample count, primary
  contrasts, total module count, count selected for detailed reporting,
  signal landscape (use one of: "low-signal", "single-axis", "richly-
  resolved", "multi-axis").
- Mention the grey module size **once**, here. Do not return to it.
- If the data block carries a `low-signal` flag (`r_q75 < 0.4`), open
  the section with that flag.

### Main findings

The lead findings, GROUPED by biological theme rather than by individual
module. Each theme is one h3 subsection; each theme contains the modules
that fit it as h4 subsubsections.

```
### Theme 1 — [biological framing in 2–4 words]
1–3 prose paragraphs introducing the theme.

#### MEturquoise: [short theme label]
1–3 prose paragraphs.

#### MEblue: [short theme label]
1–3 prose paragraphs.

### Theme 2 — [biological framing]
#### MEgreen: [short theme label]
...
```

**Theme heading format:** `### Theme N — [2–4 words]`. Sentence case.
**Module heading format:** `#### MEcolor: [short theme phrase]`. Theme is
2–6 words, sentence case, descriptive (e.g. *"Cell cycle machinery"*),
not technical (e.g. *"GO_0007049 enriched"*).

**Per-module length:**
- Strong: 2–3 paragraphs.
- Moderate: 1 paragraph.
- Weak: do NOT get a subsection; aggregate under Minor units.

**Per-module content:**
- Open with the eigengene direction and trait association in plain
  biological terms (not just the verbal label).
- Hub genes: 4–6 named per strong module, woven with functional roles,
  ≤8 names per paragraph.
- Pathways: themed plain-English clusters (no IDs).
- One numeric anchor per paragraph maximum (typically the trait `r` for
  the lead module's lead trait).

**Theme grouping.** Group modules into themes when biology suggests
coherent groups; otherwise list modules directly. Use your judgement.

### Minor units

- 1 paragraph (or omit entirely if no weak modules).
- 1–2 sentences per minor module: name it, give the eigengene shape,
  note enrichment absence, optionally name 1–2 hub genes.
- No themes, no h3 subdivisions.

### Integrated findings

- 2–4 paragraphs.
- MUST use at least one of the three rhetorical patterns from the
  writing style (trade-off / feedback / contrast).
- ONLY cite cross-module correlations with verbal label `correlated`
  (`|r| ≥ 0.7`) or stronger. Do not infer correlations not present in
  the data block.
- This is where opposing programs and sequential handoffs are explicitly
  named.

### **Discussion**

**Single integrated narrative — no subsections.**

**Length:** 2–4 paragraphs. Bold heading: `## **Discussion**`.

**Inferring the experimental frame.** Before writing, read these data sources
in order to establish the goal of the experiment:

1. Experiment description (`EXPERIMENT:` line in the data block)
2. Phenotype / trait names (the columns in the modules summary table —
   `condition=Basal`, `Viability percent`, `notact`, `cluster=c2`, etc.)
3. Sample group names (in the `SAMPLES` line and contrast definitions)
4. Organism / species name + standard biological context for that organism
5. Any species-context fragment available in the prompt

From these, identify whether the experiment is about: (a) a disease / clinical
phenotype, (b) a drug or therapeutic intervention, (c) a developmental /
time-course process, (d) basic biology in a model organism, (e) ecology /
agriculture / aquaculture, or (f) something else. Write the Discussion
**accordingly** — disease framing only when the data justifies it.

**What goes in:**
- Anchor each strong / moderate module's biology in known literature surfaced
  via Phase B PubMed calls. Vancouver `[An]` citations — max 2 per sentence.
- If the experiment IS clearly clinical / therapeutic (gate: clinical-disease
  metadata present in the description or traits AND organism is human or
  mouse AND signal is not low-signal), include translational discussion
  tracing implications back to specific module patterns. Cite literature for
  any drug / pathway / clinical claim.
- If the experiment is NOT clinical / therapeutic, skip translational
  discussion entirely — do not narrate the absence in any form. **Do NOT
  write parenthetical "reasons" like `(reason: no clinical metadata)` or
  `[reason: non-mammalian organism]` — those are template artefacts, never
  acceptable in prose.** If translational implications are not warranted,
  the Discussion simply does not contain them; that's it.

**Style:** Hedge claims that go beyond the data; state established biology
directly. Do not introduce findings absent from Main findings or Integrated
findings above.

### **Conclusion**

- **1 paragraph**.
- Heading in **bold**: `## **Conclusion**`.
- Integrative summary; no new claims; no numbers.
- Do not introduce findings that were not in Main findings or Integrated
  findings above.

### Bibliography

Bibliography is the final section of the report. Heading: `## **Bibliography**`.

- Number entries `[A1]`, `[A2]`, … in order of first citation in prose.
- Each entry: the `citation` block returned verbatim by `read_context` —
  no paraphrasing, no reformatting.
- Every `[An]` in prose appears as an entry here; every entry here is
  cited at least once in prose.

Example of one well-formed entry:

> [A1] Pearce EL, Poffenberger MC, Chang CH, Jones RG. Fueling immunity:
> insights into metabolism and lymphocyte function. Science. 2013 Oct
> 11;342(6155):1242454. doi:10.1126/science.1242454. PMID: 24115444.

WGCNA methods literature (`[W1]`, `[W2]`, `[W3]`) is appended
deterministically by the pipeline under a SEPARATE `## **Methods literature**`
heading — do NOT generate it yourself, and do NOT include `[W*]` in prose.

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
the dataset to surface specific genesets, enrichments, and hubs that anchor
the biological theme. Derive query keywords from the experiment description,
phenotypes, organism, and the Top trait of each strong module — do not run
generic unfocused queries.

**Required calls (per dataset, before literature):**

Always start with `what="summary"` to see which module carries the theme —
modules may be color-named (`MEturquoise`) or integer-named (`ME0`–`ME4`)
depending on the dataset; never guess. Then drill on the top-ranked module
returned. `query=` accepts a JSON array of short keywords with OR-semantics
— prefer 2–3 targeted small-keyword arrays over one long natural-language
string.

```
query_wgcna(what="summary", rank_by="query_match", query=["MCM","CDK1","E2F"])
query_wgcna(what="enrichment", query=["E2F","replication","cell cycle"])
query_wgcna(what="module:<top-ranked-module-from-summary>",
            query=["MCM","CDK1","E2F"], top_n=15)
```

**Suggested keyword themes by experimental domain** (these are starting
points — derive your own from the experiment description):

- Clinical / pharmaceutical / drug discovery: `DRUG`, `TREATMENT`, `INHIBITOR`, `DISEASE`, `DISGENET`, `GWAS`, `CANCER` + clinical theme; use `query_drugs` aggressively.
- Basic biology / model-organism (mouse/rat): same keywords as clinical with translational caveat (e.g. *"mouse model of …"*); `query_drugs` OK, flag immune divergence.
- Invertebrate / microbial / non-mammalian: `RIBOSOME`, `METABOLISM`, `MITOCHONDRIA`, `CELL_CYCLE`, `STRESS`, `AUTOPHAGY`, `PROTEOSTASIS` + organism-specific (DAF for worm; TOLL/IMD/NOTCH for fly); reject MSigDB disease hits as cross-species artifact.
- Agriculture / livestock / aquaculture: `GROWTH`, `YIELD`, `STRESS`, `IMMUNE`, `MUSCLE`, `REPRODUCTION`, `QUALITY`, `DROUGHT`, `SALINITY`, `LIGHT`, `PHOTOSYNTHESIS`, `CIRCADIAN`, `HORMONE`, `ABIOTIC_STRESS`; suppress drug/disease tools.
- Ecology: `ADAPTATION`, `STRESS`, `SECONDARY_METABOLITE`, `IMMUNE_DIVERSITY`, `HOST_PATHOGEN`, `POPULATION` + ortholog-aware caveat; suppress drug/disease tools.

## Phase B — Literature search (PubMed)

**Budget by signal density.**

- Low-signal dataset (the data block carries `r_q75 < 0.4` OR no
  enrichment passes `q < 0.05`): skip the literature phase entirely.
  Discussion becomes a single sentence.

- All other datasets: between **5 and 15 PubMed calls**. Scale with
  dataset richness — more correlated modules, more biological themes,
  and more cross-module structure all justify more calls. A sparse
  dataset sits near 5; a rich, multi-axis one near 15.

**Allocate calls to:**
1. One call per strong module's biological theme.
2. One call per drug / mechanism-of-action surfaced by Phase A's
   `query_drugs`.
3. One call per `|r_ME×ME| > 0.85` module-module handoff (the strongest
   anti-correlated pair, etc.).
4. One call for an organism caveat if non-mammalian AND the report makes
   any translational claim.

### How to query well

`search_context` wraps PubMed Entrez and accepts the full Entrez syntax:
field tags (`[ti]` for title, `[tiab]` for title+abstract, `[mesh]` for
MeSH controlled vocabulary, `[ORGN]` for organism), boolean operators
(`AND` / `OR` / `NOT`), and date ranges. A scoped query of 1–2 concepts
beats an unscoped query of 5. Two illustrative shapes (do not copy
verbatim — adapt to your dataset):

- Hub-gene scoped to title + concept in title/abstract:
  `query="MCM2[ti] AND replication licensing[tiab]"`
- MeSH concept filtered to organism (required for non-human/non-mouse):
  `query="oligodendrocyte[mesh] AND Salmo trutta[ORGN]"`

When 0 hits returns: decompose to one concept, add a field tag to scope
it, retry **with a structurally different query**. If still 0, accept
that the literature is absent and move on.

### Query derivation, in priority order

1. Hub gene symbols from strong modules (verbatim, not aliases or
   families).
2. The 2-word biological core of a top enrichment term (e.g.
   `replication licensing`, not `HALLMARK_E2F_TARGETS`).
3. Organism + tissue/condition phrase from sample metadata, scoped via
   `[ORGN]` when non-human/non-mouse.
4. Phase A drug name, if `query_drugs` was warranted.

Stop at the first level that produces a focused query; do not pile them.

### Title-gate before `read_context`

Read the search-result title first; only call `read_context` when the
title matches BOTH the specific term AND the biological context. When
you do call it, batch every section you'll need from that id in ONE
call (e.g. `sections=["abstract","citation"]`) — do not split into two
round-trips.

- Title *"GSK3 inhibition sustains naive T-cell stemness"* → accept.
- Title *"CHIR-99021 induced osteogenesis in MSCs"* → reject (wrong tissue).
- Title *"WGCNA reveals modules in T cells"* → reject (methods paper).

## Citation discipline

- **Cite to extend, not duplicate.** If Phase A already established a
  claim with `q < 1e-10`, do not also cite literature for it. Use the data.
- **State confirming citations directly:**
  *"MCMs license replication origins [A1]."*
- **Hedge context-only ones:**
  *"A recent study suggests this may extend to the resting state [A2]."*
- **Cap (per sentence):** ≤2 citations per sentence — avoid citation salad.
- **Minimum (per report):** ≥5 citations total in the report (when the literature phase ran). Increase citation count with dataset complexity — more strong modules and richer biology warrant more citations.
- **Bibliography format:** cite verbatim from the `citation` block returned
  by a single `read_context(sections=["abstract","citation"], id=...)` call —
  the abstract grounds the claim, the citation formats the bibliography.
  Do not split into two calls. Do not paraphrase or reformat the citation.
  Numbering follows order of first citation in prose, starting `[A1]`.

## Tool surface (use only these)

- `manage_pgx`, `query_wgcna` (Phase A only — do not free-explore)
- `query_gene_info` (sparingly, for unfamiliar hubs)
- `search_context`, `read_context` (Phase B literature only)

Do NOT call `query_de`, `query_drugs`, `query_pathways`, etc. unless the
domain table above says to. Stay disciplined: tool calls are billable and
the budget targets <$0.10 per report.

**Stop condition: report drafted OR 30 tool calls reached, whichever comes first.**
