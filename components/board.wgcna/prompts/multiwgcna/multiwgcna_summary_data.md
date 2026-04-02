## EXPERIMENT CONTEXT

This analysis uses Multi-Omics WGCNA to identify co-expression modules within each omics layer and to connect related modules across layers.

Experiment: {{experiment}}

## Task

Interpret the biological function and significance of Multi-Omics WGCNA module **{{module}}** from layer **{{layer}}** based on the within-layer enrichment, trait association, and cross-layer integration evidence provided below.

## Module Overview

**Layer:** {{layer}}

**Module:** {{module}}

**Associated Phenotypes:** {{phenotypes}}

{{module_stats}}

Treat the eigengene profile and the cross-layer partner modules as complementary evidence. Explain both the local biology of this module and how it participates in the broader multi-omics program.

## Enrichment Analysis Results

The following pathways and gene sets show significant enrichment in this module:

{{genesets}}

Focus on identifying:
1. The dominant biological theme or process in this layer
2. How multiple enriched pathways relate to each other
3. How the module connects to phenotype-associated biology across layers

## Hub Features
{{keygenes_section}}

## Cross-Layer Evidence
{{cross_layer_section}}

## Output Instructions

Always begin your summary with a heading in this exact format (maintain strict adherence):
**{{layer}} | {{module}} module** | Correlated phenotypes: {{phenotypes}}

Then synthesize the above information to explain:
- The primary biological function of this module in its own layer
- How the enriched pathways collectively support this interpretation
- The biological relevance to the experimental phenotypes
- How hub features might organize this module's function
- How this module aligns with or opposes related modules from other omics layers
- Use quantitative metrics to support your interpretation (cite specific scores, q-values, MM, TS, LogFC, centrality, and cross-layer correlation values when available)
