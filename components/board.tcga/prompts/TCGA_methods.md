## EXPERIMENT CONTEXT

This analysis uses The Cancer Genome Atlas (TCGA) to evaluate survival associations and gene expression patterns across 32 cancer types in over 10,000 cancer patients.

Experiment: {experiment}

## TCGA-SPECIFIC GUIDANCE

- TCGA is a landmark cancer genomics program that molecularly characterized over 20,000 primary cancer samples spanning 33 cancer types
- Survival analysis is performed using the Kaplan-Meier method, comparing patients whose expression profiles are positively vs. negatively correlated with the user's gene signature
- Each TCGA cohort is dichotomized based on correlation with the input signature (contrast fold-changes or a custom gene list)
- The log-rank test is used to assess statistical significance of survival differences between the two patient groups
- The number of top correlated genes (ntop) determines how many genes from the signature are used for the correlation

## QUANTITATIVE INTERPRETATION GUIDELINES

Survival metrics:
- p-value < 0.01: strong evidence of survival difference; 0.01-0.05 standard significance; 0.05-0.10 marginal; > 0.10 not significant
- Hazard ratio > 1: higher risk in the positively correlated group; < 1 lower risk; ~1 no difference
- Median survival difference: clinically meaningful if > 6 months in most cancer contexts

Gene expression metrics:
- logFC > 1: strong differential expression; 0.5-1 moderate; 0.2-0.5 mild; < 0.2 minimal
- Genes with high absolute fold change in the input contrast are likely drivers of the survival association
- Concordance of gene-level changes across multiple cancer types strengthens the biological interpretation

Integration rules:
- Prioritize cancer types with both significant p-values and biologically plausible survival associations
- Look for tissue-specific patterns: signatures may have opposite survival effects in different cancer types
- Consider sample size: cohorts with fewer than 50 patients per group may yield unreliable survival estimates
- Cross-reference top genes with known cancer driver genes and therapeutic targets
- Pan-cancer consistency of survival effects suggests a fundamental biological mechanism rather than tissue-specific artifact
- Treat marginal p-values or small cohorts as tentative findings requiring validation
