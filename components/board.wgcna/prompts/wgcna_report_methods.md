## Methods

Weighted gene co-expression network analysis (WGCNA) was performed on
{{n_features_wgcna}} {{feature_type}} across {{n_samples}} samples using a
{{network_type}} network with soft-thresholding power {{power}} and minimum
module size {{min_mod_size}} (merge cut height: {{merge_cut_height}},
minimum kME: {{min_kme}}). The analysis identified {{n_modules}} co-expression
modules (excluding the grey/unassigned bin of {{grey_size}} {{feature_type}}).

Module-trait correlations were computed using Pearson correlation between
module eigengenes and binary trait indicators. Pathway enrichment was
tested against {{n_genesets_tested}} gene sets. Hub genes were ranked by
module membership (MM), defined as the correlation between a gene's
expression profile and the module eigengene.

Module signal strength was classified as strong (≥10 significant enrichments
at q < 0.05, |trait r| > 0.7, ≥50 genes), moderate (1-9 significant
enrichments at q < 0.05), or weak (no significant enrichment).

### Compute Settings

| Parameter | Value |
|-----------|-------|
| Network type | {{network_type}} |
| Soft-thresholding power | {{power}} |
| Minimum module size | {{min_mod_size}} |
| Merge cut height | {{merge_cut_height}} |
| Minimum kME | {{min_kme}} |
| Features used | {{n_features_wgcna}} {{feature_type}} |
| Samples | {{n_samples}} |
| Gene sets tested | {{n_genesets_tested}} |
