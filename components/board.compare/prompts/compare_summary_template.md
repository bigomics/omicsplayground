## Task

Interpret the biological significance of comparing fold-change signatures between the contrasts **{contrast1}** and **{contrast2}** based on the comparison data provided below.

## Analysis Overview

**Contrast 1:** {contrast1}

**Contrast 2:** {contrast2}

This analysis compares the transcriptional responses of two experimental contrasts by correlating their gene-level fold-change profiles. Concordant changes suggest shared biological mechanisms while discordant changes highlight condition-specific responses.

## Correlation Statistics

{correlation_stats}

Focus on identifying:
1. The degree of similarity between the two transcriptional responses
2. Whether the correlation is driven by a subset of genes or is broadly distributed
3. What the correlation strength implies about shared versus distinct biology

## Shared Significant Genes

The following genes are differentially expressed in both contrasts:

{shared_genes}

## Contrast-Specific Genes

The following genes are differentially expressed in only one of the contrasts:

{contrast_specific}

## Output Instructions

Always begin your summary with a heading in this exact format (maintain strict adherence):
**{contrast1}** vs **{contrast2}** | Comparison summary

Then synthesize the above information to explain:
- The overall similarity or divergence between the two contrasts based on fold-change correlation
- Key biological processes or pathways represented by the shared significant genes
- What the contrast-specific genes reveal about condition-unique biology
- How the shared and unique gene signatures collectively inform the biological relationship between the two conditions
- Use quantitative metrics to support your interpretation (cite specific correlation values, logFC, and gene counts when available)
