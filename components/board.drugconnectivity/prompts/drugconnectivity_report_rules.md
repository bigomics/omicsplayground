## What This Analysis Reveals

Drug connectivity maps the experimental gene expression signature against thousands
of known drug perturbation profiles (L1000 / CMap). The result is a
**pharmacological portrait of the biological state**: drugs that oppose the
signature (negative NES) identify pathways active in the experimental condition;
drugs that mimic it (positive NES) reveal which biological programme the
experimental state resembles.

The report articulates this portrait in biological language. The drugs are
evidence. The biological programmes are the subject.

## Evidence Priority

When the input includes a deterministic **Interpretation evidence summary** for a
contrast, use it as the primary evidence anchor. Apply this priority order:

1. Supported MOA summary terms
2. Corroborating target summary terms
3. Raw drug tables for specific cited examples

If the summary and a raw table appear to disagree, prefer the summary wording
and use the raw table only to illustrate a claim that is already supported by
the summary.

## Analytical Unit: the Biological Programme

Each `### heading` in Main Findings names a **biological programme or condition
phase** revealed by the connectivity pattern — not a contrast ID, not a drug class.

Derive the heading from the dominant opposing MOA theme and the experimental
condition:

```
GOOD: ### treated_early–mid: Proteostasis-Driven Stress Response
GOOD: ### KO_vs_WT: Suppressed Inflammatory Signalling
BAD:  ### treated_48h_vs_ctrl
BAD:  ### Top Opposing Drugs
```

**Multi-contrast data with temporal or ordered conditions**: you MUST group contrasts
into 2–4 phases by which pharmacological themes dominate — do NOT write one
subsection per contrast. Inspect the **Cross-Contrast MOA Convergence** matrix first:
find which MOA classes are consistent across all time points (stable programme), which
emerge only at later time points, and which fade. Build phases from these patterns.
A single subsection covering all five time points is acceptable if the pharmacological
theme is identical throughout. Name each phase by its biology, not by the contrast IDs.

**Multi-contrast data without clear ordering**: identify the 2–4 dominant
pharmacological themes that emerge across contrasts; each becomes one subsection.

**Single-contrast data**: one subsection for the dominant biological programme;
use `### Minor Signals` only if a weak secondary theme exists.

Use the **Contrast Signal Classification** tier (strong / moderate / weak /
data-limited) to decide which contrasts merit individual subsections:
- **Strong / Moderate**: individual `### <programme>` subsection
- **Weak / Data-limited**: aggregate into `### Minor Signals`

Omit `### Minor Signals` entirely if all contrasts are strong or moderate.

## Per-Unit Content: Three Narrative Beats (always prose)

**Beat 1 — The biological state (lead sentence)**
State what programme the connectivity pattern reveals. Lead with the biology.
Name the experimental unit explicitly in the opening sentence using the
experiment or contrast label so downstream diagram generation can anchor the
programme to the sample.

> GOOD: "The treated signature is consistent with an activated ubiquitin–proteasome
> axis engaged in resolving misfolded-protein load."
> BAD: "32 drugs showed significant connectivity at q < 0.05."

**Beat 2 — The pharmacological evidence**
Name the MOA class(es) that oppose or mimic the signature. Cite 2–3 exemplar
compounds in parentheses with their NES. When the Cross-Contrast MOA Convergence
matrix shows a class is consistent across multiple contrasts, say so explicitly —
this is stronger evidence than a single-contrast signal.

For every exemplar drug named in the prose, state its relation to the
experimental sample explicitly in plain syntax: `<sample/contrast> is opposed by
<drug>` or `<sample/contrast> is mimicked by <drug>`.

Use summary-supported MOA terms first. If no significant MOA term exists for a
direction, say so explicitly and downgrade that direction to nominal, secondary,
or data-limited rather than promoting a raw-drug mechanism to the main story.

When the input includes **Preferred ... exemplars matching supported MOA
terms**, use those compounds first. Only fall back to the broader drug tables if
no preferred exemplar is available for that mechanism.

Do not infer a mechanism from an unannotated compound. A drug can be cited as an
example only when its MOA or target annotation supports the mechanism being
described.
Do not use a compound as an exemplar when its annotation contradicts the
selected mechanism or direction of effect.

> GOOD: "Proteasome inhibitors oppose the signature consistently across all four
> conditions (MOA NES = −2.1 to −2.5, q < 0.05 in three of four contrasts),
> exemplified by *compound-X* (NES = −2.23) and *compound-Y* (NES = −2.11) [1,2]."

**Beat 3 — Target convergence and confidence**
Does the molecular target enrichment corroborate the MOA pattern? Name the
molecular targets where they converge with the MOA signal. State confidence
explicitly when annotation coverage is low (< 30 %) or effect sizes are weak
(|NES| < 1.5 throughout).

Only describe target enrichment as corroborating when the input summary presents
those targets as supporting the selected MOA theme. If target support is absent,
nominal, mixed, or data-limited, say that explicitly.

## Cross-Unit Synthesis Paragraph (REQUIRED when > 1 unit)

After all unit subsections, one paragraph with no heading, describing:
- How the pharmacological landscape evolves across conditions or time
- Which themes are stable vs. which emerge late or are contrast-specific
- Whether the opposing and mimicking programmes together tell a coherent story
  about the biological state

Every cross-unit claim must be grounded in values from the data.

## Discussion

1–2 paragraphs interpreting the overall connectivity portrait:
- What biological programme does the signature most strongly represent?
- What does the balance of opposing vs. mimicking drugs suggest about the system?
- Translational caveats once, briefly: L1000 cell-line context, annotation gaps,
  dose–time mismatch

## Word Limit

800–1200 words total across Main Findings + Discussion + Conclusion.

## What NOT to Write

- Do NOT create section headings named "Top Opposing Drugs", "Top Mimicking
  Drugs", "MOA Class Convergence", "Decision Snapshot", "Candidate
  Prioritization", "Validation Next Steps", "Limitations" — these describe
  tables and tasks, not biological narratives
- Do NOT reproduce drug tables — the researcher already has them in the interface
- Do NOT write bullet lists of drug names with NES values; weave them into prose
- Do NOT enumerate every contrast separately when a thematic or temporal grouping
  captures the pattern more clearly
- Do NOT reproduce prompt labels such as "Beat 1", "Beat 2", or "Beat 3"
- Do NOT treat unsupported or nominal-only MOA terms as dominant programmes
- Do NOT cite unannotated compounds as mechanistic exemplars

## Data References

Use `[n]` citations for specific drug names and MOA class terms as they appear in
the input tables. Collect them in `## Data References` as:

- [1] drug or MOA term — NES, q-value (from data)

Up to 25 references. Every `[n]` in the narrative must appear in Data References,
and every entry must be cited at least once.
