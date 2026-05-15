# WGCNA report — board rules

Total report length: 2000–3000 words.

## Evidence scope

Your sole evidence is the data block provided with this prompt. You have no
access to PubMed, external databases, or prior experimental results. All
biological interpretations must be traceable to the data block (module-trait
associations, enrichment themes, hub gene functions). Do not introduce claims
from parametric knowledge that are not anchored in the data block.

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

**Theme grouping.** Group modules into themes when the enrichment patterns
and trait associations suggest coherent biological programs; otherwise list
modules directly. Use your judgement.

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

### Discussion

Brief synthesis layer for the candidate pipeline.

- 1–2 paragraphs.
- Heading in **bold**: `## **Discussion**`.
- Synthesize the biology revealed by the modules into a coherent narrative — what does the dataset *as a whole* suggest about the underlying biology?
- **No literature references.** This pipeline has no access to PubMed; do not cite papers, authors, or prior work.
- **No translational implications.** Do not propose drug targets, clinical implications, or therapeutic strategies — those require literature grounding the agentic pipeline can do but this one cannot.
- This is purely intra-dataset synthesis. Stay grounded in what the data block contains.
- **Low-signal datasets** (data block flagged `low-signal` or `r_q75 < 0.4`):
  limit the Discussion to **1 paragraph of no more than 3 sentences**.

### Conclusion

- **1 paragraph**.
- Heading in **bold**: `## **Conclusion**`.
- Integrative summary; no new claims; no numbers.
- Do not introduce findings that were not in Main findings or Integrated
  findings above.
