## Task

Interpret the biological significance of the pathway enrichment results for the contrast **{contrast}** using {pathway_type} pathway analysis.

## Contrast

**Comparison:** {contrast}

## Enriched Pathways

The following pathways show significant enrichment for the selected contrast:

{genesets}

Focus on identifying:
1. The dominant biological themes or processes
2. How multiple enriched pathways relate to each other mechanistically
3. The direction of pathway regulation (activated vs. repressed)
4. Potential biological implications for the experimental comparison

## Output Instructions

Always begin your summary with a heading in this exact format (maintain strict adherence):
**{contrast}** | Pathway database: {pathway_type}

Then synthesize the above information to explain:
- The primary biological processes affected in this comparison
- How the enriched pathways collectively support a coherent biological interpretation
- The balance between upregulated and downregulated pathways and what this implies
- Potential mechanistic connections between the top enriched pathways
- Use quantitative metrics to support your interpretation (cite specific logFC values, q-values when available)
