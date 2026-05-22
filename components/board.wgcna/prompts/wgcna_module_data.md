### {{ME_color}} ({{n_genes}} genes) — {{tier}}

Eigengene profile: {{eigengene_profile_qualitative}}
[e.g. "rises sharply at 72h, peaks at 96h"]
Trait coordination: ↑ in {{top_pos_trait}} ({{top_pos_verbal}}), ↓ in {{top_neg_trait}} ({{top_neg_verbal}})
[Omit a side if no trait clears |r| >= 0.5. If neither side has a qualifying
trait, omit the Trait coordination line entirely.]

Enrichment ({{n_sig_terms}} significant of {{n_total_terms}}):
{{enrichment_themes_table}}
[Themed groupings of pathway names, plain English, with verbal q label per theme.
Raw GO/Reactome/etc IDs in a parenthetical at the end of each theme line.]

Hub genes (top {{n_hub}} by MM):
{{hub_genes_table}}
[Two columns: gene symbol, known function (1 short clause). MM value omitted
from this table; available raw in the underlying R structures for any
verification needed.]

[Optional sections, only if non-empty:]
Gene families: {{gene_families_summary}}
Enrichment overlaps: {{enrichment_overlaps}}
