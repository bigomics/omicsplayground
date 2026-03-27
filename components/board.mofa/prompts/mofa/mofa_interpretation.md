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

## EXPERIMENT CONTEXT

This analysis uses Multi-Omics Factor Analysis (MOFA) to identify latent
factors that capture the major sources of variation across multiple omics
data types.

Experiment: {{experiment}}

## MOFA-SPECIFIC GUIDANCE

- MOFA factors are latent variables that capture shared and unique sources
  of variation across multiple data modalities (e.g., transcriptomics,
  proteomics, metabolomics)
- Factor loadings (weights) indicate the contribution of each feature to a
  given factor; high absolute weight means the feature is strongly
  associated with that factor
- Factor scores represent sample-level values for each factor and can be
  correlated with phenotypes/traits
- Enriched pathways for a factor reveal the biological processes captured
  by that factor
- Variance explained per factor and per view indicates which factors and
  data types drive the main heterogeneity

### Multi-view variance decomposition

- A factor with high variance explained across multiple views captures
  shared biology; a factor dominated by a single view captures
  view-specific variation
- When interpreting a multi-view factor, look for convergent biological
  themes across the contributing data types
- A factor driven by a single view is not necessarily less important —
  it may reflect biology visible only in that modality

## QUANTITATIVE INTERPRETATION GUIDELINES

Definitions and relationships:
- Factor loadings (weights) quantify the relationship between each feature
  and a latent factor; features with high absolute loadings are the most
  important contributors to that factor
- Factor scores represent each sample's position along a factor axis;
  differences in scores between conditions suggest the factor captures
  condition-specific variation
- Variance explained (R2) measures how much of the total variance in each
  view is captured by each factor
- Factor-trait correlation quantifies the association between a factor and
  a phenotype/trait variable

Important note on thresholds:
- MOFA does not prescribe universal cutoffs; interpretation depends on the
  dataset and number of factors
- Treat any numeric threshold as a guideline, not a rule

Enrichment metrics:
- NES (Normalized Enrichment Score): |NES| > 2.0 very strong; 1.5-2.0
  strong; 1.0-1.5 moderate; < 1.0 weak
- padj < 0.01: high confidence; 0.01-0.05 standard; 0.05-0.10 marginal;
  > 0.10 not significant
- Pathway size: larger gene sets provide more robust enrichment estimates;
  very small sets (< 10) should be interpreted cautiously

Gene/feature metrics:
- Loading weight: features with the highest absolute weights are key
  drivers of the factor
- Centrality: high centrality features are well-connected in the factor's
  feature network
- Consider both the magnitude and sign of loadings when interpreting
  factor biology

### Metric interpretation rules (IMPORTANT)

- **Weights and centrality are factor-specific but view-independent.** A
  feature's weight describes its role within a given factor regardless of
  trait. Cite without trait qualifier: "*PLCG2* has a strong loading
  (weight = 0.85, centrality = 1.12) on Factor 1."
- **Trait correlations are trait-specific.** Factor-trait correlations are
  computed relative to one particular phenotype and ALWAYS require stating
  which trait. Write "Factor 1 correlates with treatment response
  (r = 0.72)" — never just "r = 0.72" without the trait.
- **When multiple traits are present**, compare how factors associate
  across them (e.g., "Factor 1 associates more strongly with survival
  (r = 0.72) than with age (r = 0.31)").

Integration rules:
- Prioritize factors with high variance explained and significant trait
  correlations
- Favor pathways with strong NES and low padj that are consistent with
  the top loading features
- When multiple omics layers contribute to a factor, look for convergent
  biological themes across data types
- Treat marginal enrichments or low-weight features as tentative and avoid
  over-interpretation

Use quantitative values from the tables and avoid causal overclaims.

## Writing Style: Pathway and Gene-Set References

Write the narrative using only natural biological language:

1. **No raw identifiers in prose.** Never embed database accessions or raw
   pathway names (e.g., GO_CHROMOSOME_SEGREGATION, R-HSA-68886) into the
   narrative text.

2. **Thematic grouping with inline references.** Cluster enriched pathways
   into 2-3 biological themes per factor. Pair every theme mention with
   bracketed [n] references. Never group more than 3 references in a
   single bracket.

3. **Plain-English descriptions.** Refer to pathways by their biological
   meaning (e.g., "mitotic prometaphase," "lipid metabolism") rather than
   database labels.

### Enrichment References

- Reference enrichment terms via [n] brackets in the narrative
- Use plain-English descriptions in prose: "DNA replication [1,2]"
  not "REACTOME_DNA_REPLICATION [1]"
- The Data References section at the end maps [n] to exact term names
- When >100 terms are significant, summarize by theme
- When 0 terms are significant, state it explicitly
