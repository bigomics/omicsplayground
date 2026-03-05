EXPERIMENT: {{experiment}}
ORGANISM: {{organism}}
SAMPLES: {{n_samples}} ({{sample_groups}})
FEATURES: {{n_features}} ({{n_wgcna_features}} used for WGCNA)
WGCNA PARAMETERS: {{wgcna_params}}

## Module Overview
(Summary of all modules; sorted by significant enrichment count descending)

{{overview_table}}

## Per-Module Detail
(Eigengene profile, top trait, enrichment, hub genes, overlap, and gene families per module; use as evidence for the narrative, do not reproduce as tables)

{{module_detail}}

## Experimental Contrasts
(Contrast definitions from the experimental design)

{{contrasts}}

## Module-Module Eigengene Correlations
(Pairwise Pearson correlations between module eigengenes where |r| >= 0.70)

{{module_cors}}
