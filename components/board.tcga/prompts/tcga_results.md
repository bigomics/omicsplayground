## Task

Interpret the biological and clinical significance of the TCGA survival analysis results for the contrast **{{contrast}}** based on the data provided below.

## Analysis Overview

**Contrast:** {{contrast}}

**Phenotype:** {{phenotype}}

**Signature type:** {{signature_type}}

The TCGA analysis correlates the gene expression signature from this contrast with survival outcomes across 32 TCGA cancer types. Each cancer cohort is split into patients positively and negatively correlated with the signature, and survival is compared using the Kaplan-Meier method.

## Top Signature Genes

The following genes define the expression signature used for the TCGA survival analysis:

{{top_genes}}

Focus on identifying:
1. Known cancer-related genes and their potential roles in tumor biology
2. Functional coherence among the top signature genes
3. Potential mechanistic links between the gene signature and patient survival

## Survival Results Across Cancer Types

{{survival_results}}

## Output Instructions

Always begin your summary with a heading in this exact format (maintain strict adherence):
**{{contrast}}** | TCGA survival summary

Then synthesize the above information to explain:
- Which cancer types show significant survival associations with this gene signature
- The biological plausibility of the observed survival effects in specific cancer contexts
- How the top signature genes may contribute to the survival associations
- Whether the survival effects are consistent across cancer types (pan-cancer) or tissue-specific
- Use quantitative metrics to support your interpretation (cite specific p-values, gene fold changes, and cancer types when available)
