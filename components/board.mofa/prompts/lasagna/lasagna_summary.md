## EXPERIMENT CONTEXT

This analysis uses LASAGNA (Layered Stacked Analysis of Gene Network
Architecture) to build a stacked multi-layer network from multiple
omics modalities for a single experimental contrast and identify
cross-layer modules of coordinated features.

Experiment: {{experiment}}

## Task

Interpret the biological function and significance of LASAGNA module
**{{module}}** in contrast **{{contrast}}** based on the layer
participation, hub nodes, and cross-layer connectivity provided below.

## Module data

{{module_detail}}

## Output Instructions

Always begin your summary with this heading: **Module {{module}}**

Then synthesise the data block above to explain:
- The dominant biological programme this module captures
- Which omics layers contribute and whether the module is cross-layer
  (multi-layer) or single-layer
- How the hub nodes anchor the programme — especially nodes with
  inter-layer edges, which carry the strongest evidence
- The biological relevance of the cross-layer edges to the contrast

Treat every verbal label in the data block (`hub` / `central` /
`well-bridged` / `strong` / `correlated` …) as authoritative — do not
re-derive from raw numbers; raw numbers are omitted by design.
