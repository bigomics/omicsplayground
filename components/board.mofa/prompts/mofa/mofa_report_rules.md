# MOFA report — board rules

Total report length: 800–1200 words.

## Evidence scope

Your sole evidence is the data block provided with this prompt. You have no
access to PubMed, external databases, or prior experimental results. All
biological interpretations must be traceable to the data block (factor
loadings, variance explained, pathway enrichments, trait correlations).
Do not introduce claims from parametric knowledge that are not anchored
in the data block.

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
- Do not mention methods details — refer to methods as "see Methods
  section below" or similar.

### Main findings

The lead findings, GROUPED by biological theme rather than by individual
factor. Each theme is one h3 subsection; each theme contains the factors
that fit it as h4 subsubsections.

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

**Theme heading format:** `### Theme N — [2–4 words]`. Sentence case.
**Factor heading format:** `#### Factor N: [short theme phrase]`. Theme is
2–6 words, sentence case, descriptive (e.g. *"Inflammatory signalling
axis"*), not technical (e.g. *"NES > 2 in MSIGDB"*).

**Per-factor length:**
- Strong: 2–3 paragraphs.
- Moderate: 1 paragraph.
- Weak: do NOT get a subsection; aggregate under Minor units.

**Per-factor content:**
- Open with the dominant pathway theme and trait association in plain
  biological terms (not just the verbal label).
- Top features: 4–6 named per strong factor, woven with view of origin
  and functional roles, ≤8 names per paragraph. When a feature has high
  weight in multiple views, say so explicitly — this is stronger evidence.
- Pathways: themed plain-English clusters (no IDs).
- One numeric anchor per paragraph maximum (typically the top NES for
  the lead enrichment, or top weight for the lead feature).

**Theme grouping.** Group factors into themes when the enrichment patterns
and view participation suggest coherent biological programs; otherwise
list factors directly. Use your judgement.

### Minor units

- 1 paragraph (or omit entirely if no weak factors).
- 1–2 sentences per minor factor: name it, give the dominant view, note
  enrichment absence or weakness, optionally name 1–2 top features.
- No themes, no h3 subdivisions.

### Integrated findings

- 2–4 paragraphs.
- MUST use at least one of three rhetorical patterns: trade-off,
  feedback, or contrast.
- ONLY cite cross-factor relationships that are explicit in the data
  block (factor correlations, opposing trait associations, shared
  features across factors). Do not infer relationships not present
  in the data.
- This is where opposing programs and view-spanning convergence are
  explicitly named.

### **Discussion**

- 1–2 paragraphs.
- Heading in **bold**: `## **Discussion**`.
- Synthesize the biology revealed by the factors into a coherent
  narrative — what does the dataset *as a whole* suggest about the
  underlying multi-omics biology?
- **No literature references.** This pipeline has no access to PubMed;
  do not cite papers, authors, or prior work.
- **No translational implications.** Do not propose drug targets,
  clinical implications, or therapeutic strategies — those require
  literature grounding that this pipeline cannot provide.
- This is purely intra-dataset synthesis. Stay grounded in what the
  data block contains.
- **Low-signal datasets** (no factor clears 5 significant pathways AND
  max |weight| < 0.5): limit the Discussion to **1 paragraph of no more
  than 3 sentences**.

### **Conclusion**

- **1 paragraph**.
- Heading in **bold**: `## **Conclusion**`.
- Integrative summary; no new claims; no numbers.
- Do not introduce findings that were not in Main findings or Integrated
  findings above.

## Hard Constraints

- Use only facts explicitly present in the input data.
- Do NOT infer causality from factor weights or correlations alone.
- If a factor has no significant enrichment, state that explicitly and
  keep confidence modest.
- Do NOT reproduce feature tables or pathway tables — the researcher
  already has them in the interface.
- Do NOT create bullet lists of feature names with weights; weave them
  into prose.
