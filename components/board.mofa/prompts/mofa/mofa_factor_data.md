### Factor {{factor}} ({{n_features}} features) — {{tier}}

View participation: {{view_mix}}
[e.g. "transcriptome=8, proteome=4" — which omics views drive this factor]

Per-view variance explained: {{variance_line}}
[Verbal magnitude per view (negligible / minor / moderate / major).
Raw percentage shown as parenthetical anchor on the lead view only.]

Trait coordination: {{trait_summary}}
[Top correlated traits with verbal correlation label. Omit the line
entirely if no trait clears `weakly associated`.]

Pathway enrichment ({{n_sig_pathways}} significant of {{n_total_pathways}}):
{{pathway_themes_table}}
[Columns: Rank, Pathway, Direction & strength (verbal, sign-aware NES
bucket), Significance (verbal padj bucket). Themed groupings of pathway
names, plain English; raw MSigDB/Reactome/GO IDs in a parenthetical at
the end of each theme line.]

Top weighted features (top {{n_top}} by |weight|):
{{top_features_table}}
[Columns: Symbol, Contribution (verbal |w| bucket), Network role (verbal
centrality bucket). Raw weight and centrality values are omitted; they
remain available in the underlying R structures for verification.]

[Optional sections, only if non-empty:]
Cross-view features: {{cross_view_features}}
[Features with `strong` or `dominant` contribution in more than one view —
strongest evidence of integrated biology.]
