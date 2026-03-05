## EXPERIMENT CONTEXT

This analysis uses drug connectivity mapping and drug set enrichment analysis (DSEA) to identify drugs whose perturbation signatures match or oppose the experimental gene expression signature.

Experiment: {{experiment}}

## LANGUAGE RULES

Certainty calibration:
- Strong evidence (replicated, mechanistic): use "demonstrates", "shows", "reveals"
- Moderate evidence (statistical): use "suggests", "indicates", "points to"
- Preliminary or exploratory: use "may", "might", "appears to"
- Correlational data: use "is associated with", "correlates with"

Avoid:
- "proves", "establishes", "clearly demonstrates" (too strong)
- "causes", "leads to", "drives" (implies causation from correlation)
- Anthropomorphizing genes: NOT "TP53 wants" but "TP53 functions in"

## NUMERICAL INTERPRETATION

P-values in omics context:
- Omics experiments test thousands of features simultaneously
- p < 0.05 has high false positive rate; use with caution
- p < 0.01 still permits many false positives at genome scale
- FDR/adjusted p-values (q < 0.05) are more reliable for omics
- Always note whether values are raw or adjusted

Fold changes:
- |FC| < 1.5 (|log2FC| < 0.58): modest, potentially noise
- |FC| 1.5-2 (|log2FC| 0.58-1): moderate effect
- |FC| 2-4 (|log2FC| 1-2): substantial effect
- |FC| > 4 (|log2FC| > 2): strong effect

Enrichment analysis:
- High enrichment does not equal biological importance
- Well-studied pathways are over-represented in databases
- Small gene set overlaps may lack robustness
- Consider overlap size alongside enrichment score

## DRUG CONNECTIVITY-SPECIFIC GUIDANCE

Drug connectivity analysis enriches experimental gene expression signatures against pre-computed drug perturbation profiles; the analysis type determines which profiles are used and how NES should be interpreted.

**NES direction semantics — read the ANALYSIS TYPE field carefully:**

- **L1000/activity or L1000/gene (transcriptional connectivity)**: NES reflects transcriptional similarity.
  - Negative NES → drug *opposes* the experimental signature → candidate for reversing/treating the state.
  - Positive NES → drug *mimics* the experimental signature → shares the same transcriptional programme.
  - Drugs opposing the signature (negative NES) are the primary therapeutic candidates.

- **CTRPv2/sensitivity or GDSC/sensitivity (pharmacological sensitivity connectivity)**: NES reflects co-variation of drug sensitivity with the experimental gene expression profile.
  - Positive NES → the experimental state *predicts sensitivity* to the drug (vulnerability candidate).
  - Negative NES → the experimental state predicts resistance or insensitivity.
  - Drugs with positive NES are the primary therapeutic candidates (opposite direction from L1000).

Apply the correct direction interpretation throughout the report. Drug Set Enrichment Analysis (DSEA) aggregates multiple perturbation experiments per drug; MOA analysis groups by pharmacological class; target analysis groups by molecular target.

## QUANTITATIVE INTERPRETATION GUIDELINES

Drug connectivity metrics:
- |NES| > 1.5: strong connectivity; 1.0–1.5 moderate; 0.5–1.0 mild; < 0.5 minimal
- Direction: see DRUG CONNECTIVITY-SPECIFIC GUIDANCE above — direction semantics differ by analysis type
- p-value < 0.01: high confidence; 0.01–0.05 standard significance; 0.05–0.10 marginal; > 0.10 not significant
- q-value (adjusted p-value): accounts for multiple testing across all drugs tested

MOA and target analysis:
- MOA NES reflects whether drugs sharing a mechanism collectively match or oppose the signature
- Target gene NES indicates whether drugs acting on the same molecular target show concordant connectivity
- Higher absolute NES with lower q-values indicates stronger and more reliable MOA/target enrichment
- MOA classes or targets with few drugs should be interpreted cautiously

Integration rules:
- Prioritize drugs with both strong absolute NES and low q-values
- Look for convergent MOA themes among top drugs (e.g., multiple kinase inhibitors or HDAC inhibitors)
- Therapeutic candidates: for L1000 analyses, drugs opposing the signature (negative NES); for sensitivity analyses, drugs with positive NES (predicted vulnerability)
- Cross-reference drug targets with differentially expressed genes for mechanistic insight
- Cell-line caveat applies to all analysis types; L1000 also has dose/time mismatch; CTRPv2/GDSC reflects cancer-line pharmacology

## Confidence Rules

- Treat q-values and significant-count metrics as primary confidence anchors.
- Use annotation coverage to qualify confidence.
- Explicitly mark data-limited contexts when annotation is sparse or significance is weak.

## Language Rules

- Use phrasing such as "suggests", "is consistent with", "prioritizes for follow-up".
- Do not claim proven efficacy, causality, or clinical benefit.
- Keep translational caveats explicit (L1000 cell line context, dose/time mismatch risk).

## Actionability Rules

- Propose concrete and feasible next steps (orthogonal transcriptomic checks, pathway assays, benchmark compounds).
- Keep recommendations tied to presented quantitative evidence.
