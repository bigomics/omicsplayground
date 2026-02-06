## Task

Interpret the biological significance of the correlation analysis results for the gene **{gene}** based on the correlation data provided below.

## Analysis Overview

**Selected Gene:** {gene} ({gene_symbol})

**Correlation Method:** Pearson correlation with Glasso-based partial correlation

**Gene Filter:** {gene_filter}

The correlation analysis identifies genes that are co-expressed with the selected gene across all samples in the dataset.

## Top Correlated Features

The following features show the strongest correlation with **{gene_symbol}**:

{correlated_genes}

Focus on identifying:
1. The dominant biological themes among the top positively and negatively correlated genes
2. Whether the correlated genes suggest a specific pathway or functional module
3. Genes with strong partial correlation that may indicate direct regulatory relationships

## Partial Correlation Network

{partial_correlation}

## Output Instructions

Always begin your summary with a heading in this exact format (maintain strict adherence):
**{gene_symbol}** | Correlation summary

Then synthesize the above information to explain:
- The primary biological processes and pathways associated with genes co-expressed with {gene_symbol}
- How the correlation pattern (positive vs. negative) reveals functional relationships
- Which genes show evidence of direct interaction (strong partial correlation) versus indirect co-regulation
- The potential biological role of {gene_symbol} based on its correlation neighborhood
- Use quantitative metrics to support your interpretation (cite specific correlation and partial correlation values when available)
