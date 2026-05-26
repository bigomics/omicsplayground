# Experiment metadata

## Overview

EXPERIMENT: {{experiment}}
ORGANISM: {{organism}}
CONTRAST: {{contrast}}
SAMPLES: {{n_samples}}
LAYERS: {{n_layers}} ({{layer_names}})
NETWORK: {{n_nodes}} nodes, {{n_edges}} edges
INTER-LAYER FRACTION: {{inter_layer_pct}}
MODULES: {{n_modules_total}} ({{n_modules_used}} used for the report)

## Experimental contrast definition

{{contrast_block}}

---

# Cross-module overview

## Layer participation

(Number of nodes contributed by each omics layer, verbalised.)

{{layer_participation}}

## Modules summary

(Sorted by signal score, descending. Only the top `n_modules_used`
modules are shown in detail; the rest are aggregated under
"Minor units" in the report.)

| Module | Tier | Nodes | Layers spanned | Cross-layer edges | Top hub |
| ------ | ---- | ----- | -------------- | ----------------- | ------- |
{{modules_summary_table}}

_Lead module: {{lead_module}} (highest signal score)._

## Cross-module bridges

(Module pairs sharing one or more hub nodes or strong cross-layer edges.
Lists only the pairs with shared evidence; omits isolated modules.)

{{module_bridges}}

---

# Per-module detail

(For each module in the top-modules selection. Used as evidence for the
narrative; do not reproduce as tables in the prose.)

{{module_detail}}
