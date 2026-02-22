## Task

Interpret the biological function and temporal dynamics of time series cluster **{{module}}** based on the enrichment analysis and gene expression patterns provided below.

## Cluster Overview

**Cluster:** {{module}}

**Contrast:** {{contrast}}

{{cluster_info}}

This cluster contains genes that share a similar temporal expression profile. Consider how the cluster's biological function relates to the observed temporal pattern and experimental conditions.

## Enrichment Analysis Results

The following gene sets show significant correlation with this cluster's temporal expression pattern:

{{genesets}}

Focus on identifying:
1. The dominant biological theme or process
2. How multiple enriched gene sets relate to each other temporally
3. Potential mechanistic connections between the temporal pattern and the biological function

## Output Instructions

Always begin your summary with a heading in this exact format (maintain strict adherence):
**{{module}}** | Contrast: {{contrast}}

Then synthesize the above information to explain:
- The primary biological function of this gene cluster
- The temporal dynamics: when are these genes activated or repressed, and what does this timing suggest?
- How the enriched gene sets collectively support this interpretation
- The biological relevance to the experimental contrast and time course design
- Use quantitative metrics to support your interpretation (cite specific correlation values, p-values, and fold changes when available)
