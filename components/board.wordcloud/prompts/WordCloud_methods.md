## EXPERIMENT CONTEXT

This analysis uses keyword enrichment (WordCloud analysis) to identify biological themes from gene set enrichment results.

Experiment: {{experiment}}

## WORDCLOUD-SPECIFIC GUIDANCE

- Keywords are extracted from the names and descriptions of enriched gene sets
- Keyword enrichment is computed by running GSEA on the enrichment score profile for each contrast
- For each keyword, the "keyword set" is the collection of gene sets that contain that keyword in the title/description
- The Normalized Enrichment Score (NES) reflects how strongly a keyword is associated with a contrast
- Keywords with high NES and low adjusted p-value represent dominant biological themes
- Semantic clustering (t-SNE/UMAP) groups co-occurring keywords, revealing related biological processes

## QUANTITATIVE INTERPRETATION GUIDELINES

Keyword enrichment metrics:
- NES > 2.0: very strong enrichment; 1.5-2.0 strong; 1.0-1.5 moderate; < 1.0 weak
- Negative NES indicates enrichment in the opposite direction of the contrast
- padj < 0.01: high confidence; 0.01-0.05 standard; 0.05-0.10 marginal; > 0.10 not significant
- Size: number of gene sets containing the keyword; larger size increases statistical power but may indicate generic terms

Interpretation rules:
- Prioritize keywords with both high absolute NES and low padj
- Group related keywords to identify overarching biological themes
- Consider keyword size: very common keywords (large size) may be less specific
- Keywords that cluster together in t-SNE/UMAP space often represent the same biological process
- Compare enrichment direction (positive vs negative NES) to understand which group shows the effect
- Treat marginal p-values or low NES as tentative and avoid over-interpretation
