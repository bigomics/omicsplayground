# Experiment metadata

## Overview

EXPERIMENT: {{experiment}}
ORGANISM: {{organism}}
SAMPLES: {{n_samples}}
FEATURES: {{n_features_total}} ({{n_features_used}} used for WGCNA)
WGCNA PARAMETERS: signed, power={{power}}, minModSize={{min_mod_size}}, mergeCutHeight={{merge_cut_height}}

## Experimental contrasts

{{contrasts_block}}

---

# Cross-module overview

## Modules summary

(Sorted by significant enrichment count, descending. Only the top
`n_top_in_report` modules are shown in detail; the rest are aggregated under
"Minor units" in the report.)

| Module | Genes | Top correlated trait | Top anti-correlated trait | Enrichment hits |
| ------ | ----- | -------------------- | ------------------------- | --------------- |
{{modules_summary_table}}

{{modules_summary_footnote}}

## Module-module eigengene correlations

(Pairs with `|r| ≥ 0.70` only; weaker pairs omitted.)

{{module_module_correlations}}

---

# Per-module detail

(For each module in the top-modules selection. Used as evidence for the
narrative; do not reproduce as tables in the prose.)

{{module_detail}}
