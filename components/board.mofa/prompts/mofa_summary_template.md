## Task

Interpret the biological function and significance of MOFA factor **{{factor}}** based on the enrichment analysis and feature loading results provided below.

## Factor Overview

**Factor:** {{factor}}

**Associated Traits:** {{traits}}

The factor scores show correlation with the traits listed above. Consider how the factor's biological function might relate to these experimental conditions.

## Enrichment Analysis Results

The following pathways and gene sets show significant enrichment for this factor:

{{genesets}}

Focus on identifying:
1. The dominant biological theme or process
2. How multiple enriched pathways relate to each other
3. Potential mechanistic connections to the associated traits

## Top Loading Features
{{top_genes}}

## Output Instructions

Always begin your summary with a heading in this exact format (maintain strict adherence):
**{{factor}}** | Associated traits: {{traits}}

Then synthesize the above information to explain:
- The primary biological function captured by this factor
- How the enriched pathways collectively support this interpretation
- The biological relevance to the experimental traits
- How top loading features might drive or characterize this factor's biology
- Whether the factor captures shared variation across multiple omics layers or is specific to one data type
- Use quantitative metrics to support your interpretation (cite specific NES, padj, loading weights, and centrality values when available)
