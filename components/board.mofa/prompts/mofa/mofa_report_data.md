# Experiment metadata

## Overview

EXPERIMENT: {{experiment}}
ORGANISM: {{organism}}
SAMPLES: {{n_samples}}
VIEWS: {{n_views}} ({{view_names}})
FACTORS: {{n_factors_total}} ({{n_factors_used}} used for the report)

## Experimental contrasts

{{contrasts_block}}

---

# Cross-factor overview

## Variance explained per view

(Per-factor variance explained across views; multi-view factors are
strongest evidence of integrated biology.)

{{variance_block}}

## Factors summary

(Sorted by signal score, descending. Only the top `n_factors_used`
factors are shown in detail; the rest are aggregated under
"Minor units" in the report.)

| Factor | Tier | Top trait | Top pathway theme | Sig. pathways |
| ------ | ---- | --------- | ----------------- | ------------- |
{{factors_summary_table}}

_Lead factor: {{lead_factor}} (highest signal score)._

## Factor-factor correlations

(Pairs with `|r| ≥ 0.30` only; weaker pairs omitted.)

{{factor_correlations}}

---

# Per-factor detail

(For each factor in the top-factors selection. Used as evidence for the
narrative; do not reproduce as tables in the prose.)

{{factor_detail}}
