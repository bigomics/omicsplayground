### Module {{module}} ({{n_nodes}} nodes) — {{tier}}

Layer participation: {{layer_mix}}
[e.g. "transcriptome=12, proteome=4" — which omics layers contribute
nodes to this module]

Cross-layer connectivity: {{connectivity_label}}
[Verbal bin: {{n_inter_edges}} inter-layer edges across {{n_layers_spanned}}
layers. Bins: rich-bridged / well-bridged / lightly-bridged /
single-layer.]

Trait coordination: {{trait_summary}}
[Verbal correlation between the module's aggregate phenotypic signal
and the active contrast. Omit interpretation if labelled "—".]

Top hub nodes (top {{n_top}} by network role):
{{top_nodes_table}}
[Columns: Symbol, Layer, Network role (verbal centrality bucket),
Cross-layer (yes/no — does this hub have at least one inter-layer
edge), Function (one-line role from gene annotation). Raw centrality
values are omitted; they remain available in the underlying R
structures.]

Top cross-layer edges (top {{n_top_edges}} by edge strength):
{{top_edges_table}}
[Columns: From layer / From node, To layer / To node, Edge strength
(verbal ρ bucket). Raw ρ values are omitted.]

[Optional sections, only if non-empty:]
Bridging features: {{cross_layer_hubs}}
[Hub nodes with at least one inter-layer edge — the strongest evidence
of integrated biology this module provides.]

Partner modules (sharing hubs or members):
{{partner_modules}}
[Modules whose top hubs or membership overlap this one — use to
identify cross-module hand-offs and biological themes that span
several modules.]
