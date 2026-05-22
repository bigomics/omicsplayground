# Drug connectivity report — board rules

Total report length: 2000–3000 words.

Drugs are evidence; the biological programmes they implicate are the
subject of the report. Never let the writing devolve into a drug
catalogue.

## Evidence scope

Your sole evidence is the data block. No PubMed, no external databases.
For each contrast, prefer the "Interpretation evidence summary" lines
over the raw drug tables — supported MOA terms first, corroborating
targets second, exemplar drugs third, raw top drugs only when no MOA
signal exists.

Read the data block's `analysis_type_description` before naming
therapeutic candidates — NES direction semantics differ between L1000
and sensitivity backends.

## Section requirements

### Highlights

- Exactly 3 bullets, ≤25 words each, declarative.
- Together answer: headline pharmacological programme; which MOA class
  is the dominant evidence; what (if anything) converges across contrasts.
- No raw NES / q-values.

### Overview

- 1–2 paragraphs.
- State: experiment intent, organism, sample/contrast count, active
  backend (L1000 vs sensitivity), how many contrasts fall into each
  tier (strong / moderate / weak / data-limited).
- If all contrasts are data-limited, open with that flag and suppress
  translational claims for the remainder of the report.

### Main findings

Each `### heading` names a **biological programme** revealed by the
connectivity pattern — not a contrast ID, not a drug class.

For multi-contrast data, derive headings by inspecting the
Cross-Contrast MOA Convergence matrix: identify 2–4 dominant
pharmacological themes (stable across contrasts, emergent late,
contrast-specific). Group contrasts under whichever programme they fit.

For single-contrast data: one h3 programme. Use `### Minor signals`
only if a weak secondary theme exists.

```
GOOD: ### Proteostasis-driven stress response
GOOD: ### Suppressed inflammatory signalling
BAD:  ### treated_48h_vs_ctrl
BAD:  ### Top opposing drugs
```

Use the **tier classification** from the data block to decide which
contrasts merit individual treatment:

- Strong / Moderate: full prose treatment within the programme's section.
- Weak / Data-limited: aggregate into `### Minor units` (below).

**Per-programme content** — three narrative beats as prose, never as
labels:

1. *State the biology.* Lead sentence names the programme the
   connectivity pattern reveals; name the contrast(s) explicitly so the
   diagram pipeline can anchor the programme to the sample.

   > GOOD: "The treated signature is consistent with an activated
   > ubiquitin–proteasome axis engaged in resolving misfolded-protein
   > load."
   > BAD: "32 drugs showed significant connectivity at q < 0.05."

2. *Name the pharmacological evidence.* MOA class(es) that oppose or
   mimic the signature; 2–3 exemplar compounds in parentheses with NES.
   When the convergence matrix shows a class is consistent across
   multiple contrasts, say so — this is stronger evidence than a
   single-contrast hit. For every exemplar named in prose, state the
   polarity explicitly:
   `<contrast> is opposed by <drug>` or `<contrast> is mimicked by <drug>`.

   Use Preferred-exemplar list first. Never cite an unannotated compound
   as a mechanistic exemplar. Never cite a compound whose annotation
   contradicts the selected mechanism.

3. *Target convergence and confidence.* State whether target enrichment
   corroborates the MOA pattern. State confidence explicitly when
   annotation coverage is low (< 30 %) or effect sizes are weak
   (|NES| < 1.5 throughout). If target support is absent, nominal, or
   mixed, say so.

### Minor units

- 1 paragraph, omitted entirely if no weak/data-limited contrasts.
- 1–2 sentences per minor contrast: name it, the dominant direction,
  flag the limitation (data-limited / weak signal / low annotation
  coverage). No exemplar drugs here.

### Integrated findings

- 2–4 paragraphs synthesising across programmes.
- Ground every cross-contrast claim in the convergence matrix:
  which MOA classes are **stable** across all contrasts, which **emerge**
  only at later time points, which **fade**, which are
  **contrast-specific**.
- Whether the opposing and mimicking programmes together tell a
  coherent biological story.
- This is the only section where cross-programme claims belong; do not
  pre-empt it in the per-programme prose.

### Discussion

- 1–2 paragraphs. Heading in bold: `## **Discussion**`.
- What programme does the signature most strongly represent?
- Balance of opposing vs. mimicking programmes — what does it say about
  the experimental state?
- Translational caveats stated **once**, briefly: L1000 cell-line
  context, dose / exposure / cell-type mismatch, annotation gaps.

### Conclusion

- 1 paragraph. Heading in bold: `## **Conclusion**`.
- Integrative summary; no new claims; no NES values.

## What NOT to write

- Section headings named "Top opposing drugs", "MOA class convergence",
  "Decision snapshot", "Candidate prioritization", "Validation next
  steps", "Limitations" — these describe tables and tasks, not biology.
- Reproductions of drug tables — the researcher already has them in the
  interface.
- Bullet lists of drug names with NES values; weave them into prose.
- The labels "Beat 1 / Beat 2 / Beat 3" — the beats structure your
  thinking, not the output.
- Unsupported or nominal-only MOA terms framed as dominant programmes.
- Unannotated compounds cited as mechanistic exemplars.

## Data references

Use `[n]` citations for specific drug names and MOA class terms as they
appear in the data block. Collect them in `## Data references` after
Conclusion:

- `[1] drug or MOA term — NES = …, q = …`

Up to 25 references. Every `[n]` in the narrative must appear in Data
references; every entry must be cited at least once.
