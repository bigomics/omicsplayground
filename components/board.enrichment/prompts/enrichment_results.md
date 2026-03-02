## Task

Interpret the biological significance of the gene set enrichment analysis results for the contrast **{{contrast}}** based on the enrichment data provided below.

## Analysis Overview

**Contrast:** {{contrast}}

**Phenotype:** {{phenotype}}

The enrichment analysis compares gene set activity between the experimental conditions defined by this contrast. Positive logFC indicates upregulation in the first group relative to the second.

## Top Enriched Gene Sets

The following gene sets show the strongest differential enrichment for this contrast:

{{genesets}}

Focus on identifying:
1. The dominant biological themes among the top enriched pathways
2. How upregulated and downregulated pathways relate to each other
3. Potential mechanistic connections to the experimental contrast

## Leading Edge Genes

{{top_genes}}

## Output Instructions

Always begin your summary with a heading in this exact format (maintain strict adherence):
**{{contrast}}** | Enrichment summary

Then synthesize the above information to explain:
- The primary biological processes activated or suppressed in this contrast
- How the enriched pathways collectively support a coherent biological interpretation
- The biological relevance to the experimental phenotype
- How leading edge genes might drive the observed enrichment patterns
- Use quantitative metrics to support your interpretation (cite specific logFC, q-values, and gene-level statistics when available)
