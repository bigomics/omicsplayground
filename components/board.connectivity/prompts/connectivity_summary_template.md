## Task

Interpret the biological significance of the connectivity analysis results for the contrast **{contrast}** based on the similarity data provided below.

## Analysis Overview

**Contrast:** {contrast}

The connectivity analysis compares the gene expression signature of this contrast against a reference database of public dataset signatures. Positive correlation (rho) indicates concordant regulation; negative correlation indicates discordant (opposite) regulation.

## Top Similar Signatures

The following reference signatures show the strongest connectivity to the query contrast:

{top_signatures}

Focus on identifying:
1. The dominant biological themes among the most similar experiments
2. Whether the top matches suggest concordant or discordant regulation patterns
3. Potential shared biological mechanisms or drug reversal opportunities

## Leading Edge Genes

{leading_edge_genes}

## Enriched Pathways Across Similar Signatures

{enriched_pathways}

## Output Instructions

Always begin your summary with a heading in this exact format (maintain strict adherence):
**{contrast}** | Connectivity summary

Then synthesize the above information to explain:
- The primary biological connections between the query experiment and the most similar public signatures
- Whether the top matches are concordant (same direction) or discordant (opposite direction) and what this implies
- Key driver genes shared across multiple similar signatures and their biological roles
- How the enriched pathways across similar signatures support a coherent biological interpretation
- Use quantitative metrics to support your interpretation (cite specific scores, rho values, and gene-level statistics when available)
