## EXPERIMENT CONTEXT

This analysis compares fold-change signatures between two experimental contrasts to identify shared and divergent transcriptional responses.

Experiment: {{experiment}}

## COMPARISON-SPECIFIC GUIDANCE

- Fold-change correlation analysis compares the log fold-change (logFC) profiles of two contrasts across all shared genes
- High Pearson correlation (r > 0.5) indicates similar transcriptional programs between the two contrasts
- Low or negative correlation suggests distinct or opposing biological responses
- Genes that are significant in both contrasts (shared genes) are the most reliable markers of common biology
- Genes significant in only one contrast (contrast-specific) highlight condition-specific responses
- The comparison uses meta fold-change values aggregated across statistical methods (DESeq2, edgeR, limma)

## QUANTITATIVE INTERPRETATION GUIDELINES

Correlation metrics:
- Pearson r > 0.7: strong positive correlation, highly similar responses
- Pearson r 0.5-0.7: moderate correlation, partially overlapping responses
- Pearson r 0.3-0.5: weak correlation, limited shared biology
- Pearson r < 0.3: negligible correlation, largely independent responses
- Negative r: opposing transcriptional programs between contrasts

Gene-level metrics:
- logFC > 1: strong differential expression; 0.5-1 moderate; 0.2-0.5 mild; < 0.2 minimal
- logFC < -1: strong downregulation; -0.5 to -1 moderate; -0.2 to -0.5 mild; > -0.2 minimal
- Genes with concordant direction (same sign logFC in both contrasts) support shared biological mechanisms
- Genes with discordant direction (opposite sign) highlight divergent responses

Integration rules:
- Prioritize genes with large absolute logFC in both contrasts as key shared drivers
- Look for functional coherence among shared significant genes (e.g., pathway membership)
- Contrast-specific genes may reveal condition-unique biology worth investigating
- Consider the magnitude of correlation alongside the number of shared significant genes
- Treat genes with marginal fold changes cautiously and avoid over-interpretation
