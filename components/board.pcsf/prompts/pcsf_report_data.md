EXPERIMENT: {{experiment}}
ORGANISM: {{organism}}
DATA TYPE: {{datatype}}
CONTRASTS: {{n_contrasts}} ({{contrast_list}})
SAMPLES: {{sample_info}}

## Per-Contrast Factual Constraints
(Reliability metrics per contrast: steiner node fraction, pathway signal strength, hub gene effects — use to calibrate confidence in each contrast)

{{fact_constraints}}

## Referenceable Evidence
(Pre-numbered evidence entries drawn from network metrics, pathway data, and recurrent pathways — cite by [n] in the narrative)

{{reference_catalog}}

## Contrast Ranking
(Quality score, tier, node/hub counts, and pathway signal per contrast — prioritise strong/moderate contrasts in the narrative)

{{contrast_ranking}}

## Per-Contrast Network Details
(Hub genes with centrality/logFC, pathway overlap, network metrics, and caveats for each contrast; use as evidence for the narrative, do not reproduce as tables)

{{contrast_details}}
