EXPERIMENT: {{experiment}}
ORGANISM: {{organism}}
SAMPLES: {{n_samples}} ({{sample_groups}})
LAYERS: {{layers}}
FEATURES: {{n_features}}
MULTIWGCNA PARAMETERS: {{wgcna_params}}

## Module Overview
(Summary of selected modules across layers; sorted by signal strength descending)

{{overview_table}}

## Per-Module Detail
(Layer, eigengene profile, top trait, enrichment, hub features, and cross-layer partners per module; use as evidence for the narrative, do not reproduce as tables)

{{module_detail}}

## Experimental Contrasts
(Contrast definitions from the experimental design)

{{contrasts}}

## Cross-Layer Module Correlations
(Pairwise Pearson correlations between module eigengenes across layers where |r| >= 0.60)

{{module_cors}}
