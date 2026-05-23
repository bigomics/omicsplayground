# Drug Connectivity Deep Report — board rules + agent skills

Loaded as `board_rules` only by the agentic Deep Report path. Section
ordering and shared structure match the single-shot report rules; the
Discussion / Conclusion / Bibliography sections and the trailing Agent
skills block are Deep-Report-specific.

Total report length: 2000–3000 words.

Drugs are evidence; the biological programmes they implicate are the
subject of the report. Never let the writing devolve into a drug
catalogue.

## Section requirements

### Highlights

- Exactly **3 bullets**, ≤25 words each, declarative.
- Together answer: headline pharmacological programme; dominant MOA
  class; what (if anything) converges across contrasts.
- No raw NES / q-values.

### Overview

- 1–2 paragraphs.
- Required content: experiment intent, organism, sample/contrast count,
  active backend (L1000 vs sensitivity), how many contrasts fall into
  each tier (strong / moderate / weak / data-limited).
- If all contrasts are data-limited, open with that flag and suppress
  translational claims for the remainder of the report.

### Main findings

Each `### heading` names a **biological programme** revealed by the
connectivity pattern — not a contrast ID, not a drug class.

For multi-contrast data, derive headings from the Cross-Contrast MOA
Convergence matrix: identify 2–4 dominant pharmacological themes
(stable across contrasts, emergent late, contrast-specific). Group
contrasts under whichever programme they fit.

For single-contrast data: one h3 programme. Use `### Minor signals`
only if a weak secondary theme exists.

```
GOOD: ### Proteostasis-driven stress response
GOOD: ### Suppressed inflammatory signalling
BAD:  ### treated_48h_vs_ctrl
BAD:  ### Top opposing drugs
```

Use the **tier classification** from the data block for depth:

- Strong / Moderate: full prose treatment (2–3 paragraphs).
- Weak / Data-limited: aggregate into `### Minor units`.

**Per-programme content — three narrative beats, never as labels:**

1. *State the biology.* Lead sentence names the programme the
   connectivity pattern reveals; name the contrast(s) explicitly.
2. *Name the pharmacological evidence.* MOA class(es) that oppose or
   mimic the signature; 2–3 exemplar compounds in parentheses with NES.
   Use Preferred-exemplar list first. Never cite an unannotated
   compound as a mechanistic exemplar. State polarity explicitly:
   `<contrast> is opposed by <drug>` or `<contrast> is mimicked by <drug>`.
3. *Target convergence and confidence.* State whether target enrichment
   corroborates the MOA. Flag annotation coverage when < 30 % of top
   drugs are annotated — do not suppress the flag.

### Minor units

- 1 paragraph, omitted entirely if no weak/data-limited contrasts.
- 1–2 sentences per minor contrast: name it, dominant direction,
  limiting factor (data-limited / weak / low coverage). No exemplar drugs.

### Integrated findings

- 2–4 paragraphs synthesising across programmes.
- Ground every cross-contrast claim in the convergence matrix: which
  MOA classes are **stable**, which **emerge** late, which **fade**,
  which are **contrast-specific**.
- Whether opposing and mimicking programmes together tell a coherent
  biological story.
- Cross-programme claims belong here only — do not pre-empt in
  per-programme prose.

### **Discussion**

**Single integrated narrative — no subsections.**

**Length:** 2–4 paragraphs. Bold heading: `## **Discussion**`.

**Inferring the experimental frame.** Before writing, read these data
sources in order to establish the goal of the experiment:

1. Experiment description (`EXPERIMENT:` line in the data block)
2. Contrast names and phenotype variable names
3. Sample group names and organism
4. Phase A DE / pathway calls that surfaced the biological context

From these, identify whether the experiment is about: (a) a disease /
clinical phenotype, (b) a drug or therapeutic intervention,
(c) developmental / time-course biology, (d) basic biology in a model
organism, (e) agriculture / ecology, or (f) something else. Write the
Discussion accordingly.

**What goes in:**

- Anchor each strong / moderate MOA programme in known literature
  surfaced via Phase B PubMed calls. Vancouver `[An]` citations —
  max 2 per sentence.
- If the experiment IS clearly clinical / therapeutic (gate: clinical
  or disease metadata present AND organism is human or mouse AND signal
  is not data-limited), include translational discussion tracing
  implications back to specific MOA programmes. Cite literature for
  any mechanism or drug class claim.
- If the experiment is NOT clinical / therapeutic, skip translational
  discussion entirely. **Do NOT write parenthetical reasons like
  `(reason: non-clinical)` — just omit the section.**
- Restate translational caveats **once**, briefly: L1000 cell-line
  context, dose / exposure / cell-type mismatch, annotation gaps.

**Style:** Hedge claims that go beyond the data; state established
biology directly.

### **Conclusion**

- **1 paragraph**. Bold heading: `## **Conclusion**`.
- Integrative summary; no new claims; no NES values.

### Bibliography

Bibliography is the final section of the report. Heading:
`## **Bibliography**`.

- Number entries `[A1]`, `[A2]`, … in order of first citation in prose.
- Each entry: the `citation` block returned verbatim by `read_context`
  — no paraphrasing, no reformatting.
- Every `[An]` in prose appears as an entry here; every entry is cited
  at least once.

Example of one well-formed entry:

> [A1] Subramanian A, Narayan R, Corsello SM, et al. A next generation
> Connectivity Map: L1000 platform and the first 1,000,000 profiles.
> Cell. 2017;171(6):1437-1452. doi:10.1016/j.cell.2017.10.049.

Drug connectivity methods literature is appended deterministically by
the pipeline under a SEPARATE `## **Methods literature**` heading —
do NOT generate it yourself.

---

## What NOT to write

- Section headings named "Top opposing drugs", "MOA class convergence",
  "Decision snapshot", "Candidate prioritization", "Validation next
  steps", "Limitations" — these describe tables and tasks, not biology.
- Reproductions of drug tables — the researcher already has them.
- Bullet lists of drug names with NES values; weave them into prose.
- The labels "Beat 1 / Beat 2 / Beat 3".
- Unsupported or nominal-only MOA terms framed as dominant programmes.
- Unannotated probe compounds (`BG-STK33-*`, `MW-STK33-*`, `TUL-*`,
  `TUL-XXI039`, `BRD-K*`) cited as mechanistic exemplars.

---

# Agent skills

This block is loaded only by the agentic Deep Report path. The
single-shot Report path has no tool access.

**Hard limits**
- Tool-call budget: **≤30** per report (Phase A + Phase B + lookups).
- Phase B duplicate-shape retries: 0 (any retry must change structure).

You have access to tools that let you query the active dataset and
search PubMed. Use them in two phases.

## Phase A — Drug data triage and biological context (BEFORE literature)

Before any PubMed call, run a focused set of passes to:

1. **Establish what data you have** — list contrasts, active backend, tiers.
2. **Triage annotation coverage** — pull the rank table (top 30) for
   each strong/moderate contrast and count what fraction of top drugs
   carry MOA + target metadata. If < 30 % annotated, that contrast is
   data-limited.
3. **Surface MOA convergence** — use `target:GENE` queries to confirm
   at the gene level that the MOA families visible in the rank table
   are supported.
4. **Anchor the biological context** — call `query_de` and/or
   `query_pathways` on the strongest contrast to understand what
   biological signal is driving the signature; this anchors the
   Discussion.

### Required call sequence

Always list contrasts before any rank calls — contrast names vary
across datasets and must not be guessed:

```
query_drugs(include_primer=false)                        # lists contrasts + backends
```

Then for **each strong/moderate contrast**, pull the full rank table:

```
query_drugs(what="rank", contrast="HC.IL17A_vs_HC.Control", top_n=30, include_primer=false)
query_drugs(what="rank", contrast="SF.IL17A_vs_SF.Control", top_n=30, include_primer=false)
query_drugs(what="rank", contrast="act48h_vs_notact",        top_n=30, include_primer=false)
```

From each result, **count annotated rows** (rows where `moa` ≠ NA).
Apply the annotation gate below before proceeding.

### Annotation coverage gate (enforce before every MOA claim)

| Coverage (annotated / top 30) | Tier       | Action |
|-------------------------------|------------|--------|
| ≥ 50 %                        | Grain-rich | Full MOA narrative; target enrichment expected |
| 30 – 50 %                     | Mixed      | Partial narrative; hedge MOA claims; target validation required |
| < 30 %                        | Data-limited | Suppress MOA claims; name only individual exemplars that carry MOA |

### Annotation triage — chaff vs grain

Unannotated probe series carry no pharmacological information and must
be silently skipped when building the narrative. They do not count
toward the grain side of the coverage gate:

```
CHAFF (skip silently):   BG-STK33-03, BG-STK33-55, BG-STK33-59
                         MW-STK33-100, MW-STK33-2B, MW-STK33-1C
                         TUL, TUL-XX023, TUL-XXI039
                         BRD-K95985487, BRD-K83509924, BRD-A62200266
```

Annotated drugs with coherent MOA + target are the grain — these are
the only entries eligible as mechanistic exemplars:

```
GRAIN (eligible as exemplars):
  AZ-628              RAF inhibitor        BRAF|RAF1
  PD-98059            MEK inhibitor        MAP2K1|MAPK1|MAPK14|…
  alvocidib           CDK inhibitor        CDK1|CDK2|CDK4|CDK7|CDK9
  fostamatinib        SYK inhibitor        SYK
  SB-239063           p38 MAPK inhibitor   MAPK11|MAPK14
  dactinomycin        RNA pol inhibitor    POLR2A
  palbociclib         CDK inhibitor        CDK4|CDK6
  doxorubicin         topoisomerase inh.   TOP2A
  CHIR-99021          GSK3 inhibitor       GSK3A|GSK3B|MAPK1
```

When 2–3 grain entries share a target family (e.g. MAPK14 appears in
PD-98059, SB-239063, VX-745, doramapimod), that family is a **supported
MOA class** and warrants a target query.

### Target validation calls

Derive gene symbols from the grain drugs you saw in the rank table —
not from generic biological intuition. Run one call per MOA family:

```
# If p38/MEK family seen:
query_drugs(what="target:MAPK14", contrast="HC.IL17A_vs_HC.Control", top_n=15, include_primer=false)
query_drugs(what="target:MAP2K1", contrast="HC.IL17A_vs_HC.Control", top_n=15, include_primer=false)

# If CDK family seen:
query_drugs(what="target:CDK4",   contrast="Resistant_vs_Sensitive",  top_n=15, include_primer=false)
query_drugs(what="target:CDK2",   contrast="act48h_vs_notact",         top_n=15, include_primer=false)

# If SYK / src family seen:
query_drugs(what="target:SYK",    contrast="HC.IL17A_vs_HC.Control", top_n=15, include_primer=false)
query_drugs(what="target:LCK",    contrast="HC.IL17A_vs_HC.Control", top_n=15, include_primer=false)
```

A target query that returns ≥ 3 grain drugs at q < 0.05 corroborates
the MOA; fewer than 3 is nominal — hedge.

### Biological context calls

Pull DE genes and/or pathways for the strongest contrast to understand
what biological programme is driving the signature. Use keywords
derived from the experiment description, phenotype names, and organism:

```
# IL-17 cytokine experiment, human synovial fibroblasts:
query_de(contrast="HC.IL17A_vs_HC.Control", top_n=20, include_primer=false)
query_pathways(contrast="HC.IL17A_vs_HC.Control", top_n=15, include_primer=false)

# Glucocorticoid resistance, human cancer cell lines:
query_de(contrast="Resistant_vs_Sensitive", top_n=20, include_primer=false)
query_pathways(contrast="Resistant_vs_Sensitive", top_n=15, include_primer=false)

# T-cell activation time course, human proteomics:
query_de(contrast="act48h_vs_notact", top_n=20, include_primer=false)
```

**Suggested keyword themes by experimental domain** (starting points —
derive your own from the experiment description):

- Clinical / disease / drug response: `MAPK`, `MTOR`, `CDK`, `HDAC`,
  `INFLAMMATION`, `CANCER`, `RESISTANCE`, `APOPTOSIS`; query target
  families aggressively; translational framing appropriate.
- Immunology / cytokine biology (human): `CYTOKINE`, `NF-KB`, `JAK`,
  `STAT`, `MAPK`, `TLR`, `INFLAMMASOME`; `query_drugs` on kinase and
  RNA-pol families; L1000 cell-line mismatch caveat required.
- Time-course / activation (human or mouse): `CELL_CYCLE`, `CDK`,
  `E2F`, `MTOR`, `PROTEASOME`, `STRESS`; scan both opposing and
  mimicking tails — mimicking drugs often reveal the active programme.
- Basic biology / model organism (mouse): same as clinical; flag that
  L1000 profiles are human cell lines; immune-gene ortholog divergence.
- Non-mammalian / non-human (C. elegans, yeast, fish, plant): suppress
  translational claims entirely; L1000 connectivity is cross-species
  artifact; use drug data only to note convergence on conserved
  pathways (e.g. proteasome, ribosome, MAPK core) with explicit caveat.

### Backend awareness (check before rank calls)

`query_drugs()` without a contrast returns the active backend. Read it:

```
L1000/activity or L1000/gene  →  negative NES = drug opposes (reversal candidate)
                                  positive NES = drug mimics (mechanistic similarity)

CTRPv2/sensitivity             →  positive NES = predicted sensitivity / vulnerability
GDSC/sensitivity               →  negative NES = predicted resistance
```

Never mix these semantics. If the dataset exposes multiple backends,
the data block's `analysis_type` field is authoritative.

### Phase A call budget by dataset complexity

| Contrasts | Recommended Phase A calls |
|-----------|--------------------------|
| 1         | 5–7 (status + list + 1 rank + 1–2 targets + 1 DE) |
| 2–3       | 9–13 (rank × contrasts + targets × top 2 MOA + DE + pathways × 1) |
| 4+        | 13–17 (rank × top 3 contrasts + targets × top 2 MOA + DE + pathways × 1–2) |

Do not query rank tables for weak/data-limited contrasts — spend the
budget on target and DE queries for strong ones.

## Phase B — Literature search (PubMed)

**Budget by signal density.**

- All contrasts data-limited (annotation coverage < 30 % throughout,
  or all contrasts weak): skip Phase B entirely. Discussion becomes a
  single hedged paragraph naming the annotation gap.
- Otherwise: **5–12 PubMed calls**, scaled with signal richness — more
  grain-rich contrasts and more distinct MOA families justify more calls.

**Allocate calls to:**

1. One call per dominant MOA programme identified in Phase A
   (e.g., `"MEK inhibitor AND IL-17"`, `"CDK inhibitor AND glucocorticoid resistance"`).
2. One call per cross-contrast convergence (a MOA class appearing at
   q < 0.05 in ≥ 2 contrasts).
3. One call for the biological context of the signature — anchor via
   top DE gene symbols from Phase A.
4. One call for organism caveat if non-human/non-mouse AND the report
   makes any translational claim.

### How to query well

`search_context` accepts full Entrez syntax: field tags (`[ti]` title,
`[tiab]` title+abstract, `[mesh]` MeSH, `[ORGN]` organism), boolean
operators, date ranges. A scoped query of 1–2 concepts beats an
unscoped query of 5. Three illustrative shapes (do not copy verbatim —
adapt to your dataset):

- MOA class anchored to the biological context of the experiment:
  `query="p38 MAPK inhibitor[tiab] AND synovial fibroblast[tiab]"`
  `query="MEK inhibitor[tiab] AND T cell activation[tiab]"`
  `query="CDK inhibitor[tiab] AND glucocorticoid resistance[tiab]"`

- Target gene anchored to disease / phenotype:
  `query="MAPK14[ti] AND inflammation[tiab]"`
  `query="SYK[ti] AND rheumatoid arthritis[tiab]"`
  `query="CDK4[tiab] AND drug resistance[tiab]"`

- Organism-scoped (required for non-human/non-mouse datasets):
  `query="MAPK[mesh] AND Caenorhabditis elegans[ORGN]"`
  `query="proteasome inhibitor[tiab] AND Saccharomyces cerevisiae[ORGN]"`
  `query="cell cycle[mesh] AND Danio rerio[ORGN]"`

When 0 hits return: decompose to one concept, add a field tag, retry
with a **structurally different query** (swap term, change field tag,
try MeSH). If still 0, accept the absence and move on.

### Query derivation, in priority order

1. Named MOA class from a supported evidence term in Phase A
   (verbatim: `p38 MAPK inhibitor`, `CDK inhibitor`, `SYK inhibitor`).
2. Key target gene symbols from Phase A target queries
   (verbatim, not family aliases: `MAPK14`, `MAP2K1`, `SYK`, `CDK4`).
3. The biological context phrase derived from contrast / phenotype names
   (`glucocorticoid resistance`, `IL-17 stimulation`, `T cell activation`).
4. Organism + tissue, scoped via `[ORGN]` when non-human/non-mouse.

Stop at the first level that produces a focused query; do not pile them.

### Title-gate before `read_context`

Read the search-result title first; only call `read_context` when the
title matches BOTH the specific MOA/target AND the biological context.
Batch every section you need in ONE call:
`sections=["abstract","citation"]` — never split into two round-trips.

```
ACCEPT:
  "p38 MAPK inhibitors suppress IL-17-driven synovial fibroblast activation"
    → matches MOA (p38 MAPK) + context (IL-17, fibroblast)

  "SYK inhibitor fostamatinib reduces joint inflammation in rheumatoid arthritis"
    → matches MOA (SYK) + context (inflammation, arthritis)

  "CDK4/6 inhibition reverses glucocorticoid resistance in ALL"
    → matches MOA (CDK4/6) + context (glucocorticoid resistance)

REJECT:
  "WGCNA identifies p38 modules in synovial tissue"
    → methods paper, not MOA mechanism

  "SB-202190 promotes osteogenic differentiation in MSCs"
    → correct drug, wrong tissue / context

  "CDK inhibitors in breast cancer: a meta-analysis"
    → correct MOA, wrong disease context (experiment is not breast cancer)
```

## Citation discipline

- **Cite to extend, not duplicate.** If Phase A established a MOA claim
  via annotation convergence (e.g., ≥4 p38 inhibitors at q < 0.05),
  do not also cite literature for it — use the data.
- **State confirming citations directly:**
  *"p38 MAPK drives IL-17-induced cytokine release in synovial
  fibroblasts [A1]."*
  *"SYK is a central node in B-cell receptor and Fc-receptor
  inflammatory signalling [A2]."*
- **Hedge context-only citations:**
  *"MEK inhibition has been proposed as a reversal strategy in this
  setting [A3]."*
  *"Whether this extends to primary patient fibroblasts remains
  unexplored [A4]."*
- **Cap:** ≤2 citations per sentence — avoid citation salad.
- **Minimum:** ≥5 total when Phase B ran; scale up with dataset richness.
- **Format:** cite verbatim from the `citation` block returned by a
  single `read_context(sections=["abstract","citation"], id=...)` call.
  Do not paraphrase or reformat. Number in order of first citation,
  starting `[A1]`.

## Tool surface (use only these)

- `query_drugs` — rank, target, moa views (Phase A)
- `query_de`, `query_pathways` — biological context for Discussion (Phase A)
- `query_gene_info` — sparingly, for unfamiliar target gene symbols
- `search_context`, `read_context` — Phase B literature only

Do NOT call `query_wgcna`, `query_deconv`, `query_mofa`, or
`query_expression` unless they surface information about the
experimental biology that drug connectivity alone cannot provide.

**Stop condition: report drafted OR 30 tool calls reached, whichever
comes first.**
