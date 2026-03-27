## What This Analysis Reveals

Multi-Omics Factor Analysis (MOFA) decomposes multi-omics datasets into
latent factors — axes of shared or modality-specific variation. Each factor
captures a biological programme that may span transcriptomics, proteomics,
metabolomics, or other views. The report translates these factors into a
biological narrative: what programmes dominate the data, which features
drive them, and how they relate to experimental conditions.

When `## Factor Ranking` is present in the input data, use the tier
(strong / moderate / weak) to calibrate the depth of each factor's
treatment.

## Analytical Unit: the Latent Factor

Each `### heading` in Main Findings names a **biological programme**
revealed by the factor — not a factor number alone.

Derive the heading from the dominant enrichment theme and the factor
identity:

```
GOOD: ### Factor 1: Inflammatory Signaling Axis
GOOD: ### Factor 3: Epigenetic Remodeling
BAD:  ### Factor 1
BAD:  ### Top Weighted Features
```

Use the **Factor Ranking** tier to calibrate depth:
- **Strong / Moderate**: full narrative treatment as individual subsections
- **Weak**: aggregate into `### Minor Factors` with explicit caveats
  about low variance explained or absent enrichment

Omit `### Minor Factors` entirely if all factors are strong or moderate.

## Per-Unit Content: Three Narrative Beats (always prose)

**Beat 1 — The biological programme (lead sentence)**
State what biological process the factor captures. Lead with the biology,
not the feature names.

> GOOD: "Factor 1 captures an inflammatory signalling axis spanning
> cytokine receptors and NF-κB pathway components across both
> transcriptomic and proteomic views."
> BAD: "Factor 1 has 12 features with high weights."

**Beat 2 — The driving features**
Name the top weighted features that anchor the factor. For each, weave
the loading weight, view of origin, and functional annotation into prose.
When a feature has high weight in multiple views, say so explicitly —
this is stronger evidence than single-view signal.

> GOOD: "*IL6* (weight = 0.91, transcriptome) and its receptor *IL6R*
> (weight = 0.78, proteome) jointly anchor this factor, consistent with
> an active JAK-STAT signalling loop."
> BAD: "IL6 weight 0.91. IL6R weight 0.78."

**Beat 3 — Pathway convergence and confidence**
Does pathway enrichment corroborate the factor? Name the enriched
pathway(s) with NES and q-values. State confidence explicitly when:
- No significant enrichment is detected (padj > 0.05 for all terms)
- Only a single view contributes to the factor
- Loading weights are modest (all |weight| < 0.5)

## Cross-Factor Synthesis Paragraph (REQUIRED when > 1 factor)

After all factor subsections, one paragraph with no heading, describing:
- Whether factors capture shared vs complementary biology
- Opposing or correlated factor pairs and what they suggest
- Whether certain views dominate specific factors while others are
  multi-view
- The overall biological narrative that emerges from the factor landscape

Every cross-factor claim must be grounded in values from the data.

## Discussion

1-2 paragraphs interpreting the overall multi-omics integration:
- What biological programmes does MOFA reveal that single-omics would miss?
- Which factors are multi-view (shared biology) vs single-view (modality-specific)?
- Caveats once, briefly: factor count sensitivity, view imbalance
  (unequal feature counts across modalities), correlation ≠ causation

## Word Limit

800-1200 words total across Main Findings + Discussion + Conclusion.

## Hard Constraints

- Use only facts explicitly present in the input data.
- Do NOT invent reference numbers or repurpose a reference for a
  different claim.
- If a factor has no significant enrichment, state that explicitly and
  keep confidence modest.
- Do NOT synthesize a numeric range across factors unless the underlying
  values are explicit in the input data.
- Do NOT reproduce feature tables or pathway tables — the researcher
  already has them in the interface.
- Do NOT create bullet lists of feature names with weights; weave them
  into prose.
