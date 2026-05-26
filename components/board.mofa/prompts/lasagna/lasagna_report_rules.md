# LASAGNA report — board rules

Total report length: 800–1200 words.

## Evidence scope

Your sole evidence is the data block provided with this prompt. You have no
access to PubMed, external databases, or prior experimental results. All
biological interpretations must be traceable to the data block (community
structure, layer participation, cross-layer edges, hub-node identities).
Do not introduce claims from parametric knowledge that are not anchored
in the data block.

## Section requirements

### Highlights

- Exactly **3 bullets**. No more, no less.
- Each bullet ≤25 words.
- Declarative; no hedging.
- No raw centrality, ρ, or edge-weight values.
- Together, the three bullets answer: what is the dataset's headline
  cross-layer programme, what is the most surprising finding, and what
  is the dominant inter-modality bridge?

### Overview

- 1–2 paragraphs.
- Required content: experiment intent, organism, sample count, contrast
  name, number of omics layers, node and edge counts, count of
  cross-layer modules detected, count selected for detailed reporting,
  signal landscape (use one of: "fragmented", "single-layer-dominated",
  "balanced", "richly-bridged").
- Mention the inter-layer edge fraction **once**, here. Do not return
  to it.
- Do not mention methods details — refer to methods as "see Methods
  section below" or similar.

### Main findings

The lead findings, GROUPED by biological programme rather than by
individual module ID. Each programme is one h3 subsection; each
programme contains the modules that fit it as h4 subsubsections.

```
### Theme 1 — [biological framing in 2–4 words]
1–3 prose paragraphs introducing the theme.

#### M1: [short theme label]
1–3 prose paragraphs.

#### M3: [short theme label]
1–3 prose paragraphs.

### Theme 2 — [biological framing]
#### M2: [short theme label]
...
```

**Theme heading format:** `### Theme N — [2–4 words]`. Sentence case.
**Module heading format:** `#### M<N>: [short theme phrase]`. Theme is
2–6 words, sentence case, descriptive (e.g. *"Kinase–TF relay"*),
not technical (e.g. *"5 cross-layer edges"*).

**Per-module length:**
- Strong: 2–3 paragraphs.
- Moderate: 1 paragraph.
- Weak: do NOT get a subsection; aggregate under Minor units.

**Per-module content:**
- Open with the biological programme the module captures in plain
  language; lead with biology, not with node names.
- Top nodes: 4–6 named per strong module, woven with layer of origin
  and functional role, ≤8 names per paragraph. Cross-layer bridges
  (hubs connecting two or more layers) called out explicitly — this is
  the strongest evidence LASAGNA provides.
- Layer participation: state whether the module is multi-layer or
  dominated by one layer.
- One numeric anchor per paragraph maximum (typically the node count
  for the lead module).

**Theme grouping.** Group modules into themes when biology suggests
coherent groups; otherwise list modules directly. Use your judgement.

### Minor units

- 1 paragraph (or omit entirely if no weak modules).
- 1–2 sentences per minor module: name it, give the dominant layer,
  note absence of cross-layer evidence, optionally name 1–2 top nodes.
- No themes, no h3 subdivisions.

### Integrated findings

- 2–4 paragraphs.
- MUST use at least one of three rhetorical patterns: trade-off,
  feedback, or contrast.
- ONLY cite cross-module bridges, opposing programmes, or shared hub
  nodes that are explicit in the data block. Do not infer topology
  not present.
- This is where layer-spanning convergence and inter-module hand-offs
  are explicitly named.

### **Discussion**

- 1–2 paragraphs.
- Heading in **bold**: `## **Discussion**`.
- Synthesise the biology revealed by the multi-layer network into a
  coherent narrative — what does the dataset *as a whole* suggest
  about the underlying cross-omics biology?
- **No literature references.** This pipeline has no access to PubMed.
- **No translational implications.** Do not propose drug targets,
  clinical implications, or therapeutic strategies.
- Stay grounded in what the data block contains.
- **Sparse networks** (data block flagged `fragmented` or fewer than
  5 inter-layer edges): limit the Discussion to **1 paragraph of no
  more than 3 sentences**.

### **Conclusion**

- **1 paragraph**.
- Heading in **bold**: `## **Conclusion**`.
- Integrative summary; no new claims; no numbers.
- Do not introduce findings that were not in Main findings or Integrated
  findings above.

## Hard Constraints

- Use only facts explicitly present in the input data.
- Do NOT infer causality from network topology alone.
- If the network is sparse or lacks inter-layer edges, state that
  explicitly and keep confidence modest.
- Do NOT reproduce node tables or edge tables — the researcher already
  has them in the interface.
- Do NOT create bullet lists of node names with centrality scores; weave
  them into prose.
