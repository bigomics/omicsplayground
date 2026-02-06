## Task

Interpret the biological function and significance of WGCNA module **{module}** based on the enrichment analysis and correlation results provided below.

## Module Overview

**Module:** {module}

**Associated Phenotypes:** {phenotypes}

{module_stats}

The module eigengene shows correlation with the phenotypes listed above. Consider how the module's biological function might relate to these experimental conditions.

## Enrichment Analysis Results

The following pathways and gene sets show significant enrichment in this module:

{genesets}

Focus on identifying:
1. The dominant biological theme or process
2. How multiple enriched pathways relate to each other
3. Potential mechanistic connections to the associated phenotypes

## Hub Genes
{keygenes_section}

## Output Instructions

Always begin your summary with a heading in this exact format (mantain stric adherence):
**{module} module** | Correlated phenotypes: {phenotypes}

Then synthesize the above information to explain:
- The primary biological function of this module
- How the enriched pathways collectively support this interpretation
- The biological relevance to the experimental phenotypes
- How hub genes might drive or regulate this module's function
- Use quantitative metrics to support your interpretation (cite specific scores, q-values, MM, TS, LogFC, and centrality values when available)
