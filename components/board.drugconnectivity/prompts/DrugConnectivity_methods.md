## EXPERIMENT CONTEXT

This analysis uses drug connectivity mapping and drug set enrichment analysis (DSEA) to identify drugs whose perturbation signatures match or oppose the experimental gene expression signature.

Experiment: {experiment}

## DRUG CONNECTIVITY-SPECIFIC GUIDANCE

- Drug connectivity analysis correlates experimental gene expression signatures with known drug perturbation profiles from the L1000 database (Connectivity Map / CMap)
- The L1000 database contains gene expression profiles measured after treating cell lines with thousands of small-molecule compounds
- A positive connectivity score (NES) indicates the drug induces a similar transcriptional response to the experiment; a negative score indicates an opposing response
- Drugs with strong negative connectivity scores are candidates for therapeutic repurposing, as they may reverse the disease signature
- Drug Set Enrichment Analysis (DSEA) aggregates multiple perturbation experiments for the same drug using GSEA, producing a single enrichment score per drug
- Mechanism of action (MOA) analysis groups drugs by their known pharmacological class or molecular target to identify enriched drug mechanisms
- Target gene analysis identifies molecular targets shared among the top-ranked drugs

## QUANTITATIVE INTERPRETATION GUIDELINES

Drug connectivity metrics:
- NES > 1.5: strong positive connectivity (drug mimics the experimental signature); 1.0-1.5 moderate; 0.5-1.0 mild; < 0.5 minimal
- NES < -1.5: strong negative connectivity (drug opposes the experimental signature); -1.0 to -1.5 moderate; -0.5 to -1.0 mild; > -0.5 minimal
- p-value < 0.01: high confidence; 0.01-0.05 standard significance; 0.05-0.10 marginal; > 0.10 not significant
- q-value (adjusted p-value): accounts for multiple testing across all drugs tested

MOA and target analysis:
- MOA NES reflects whether drugs sharing a mechanism collectively match or oppose the signature
- Target gene NES indicates whether drugs acting on the same molecular target show concordant connectivity
- Higher absolute NES with lower q-values indicates stronger and more reliable MOA/target enrichment
- MOA classes or targets with few drugs should be interpreted cautiously

Integration rules:
- Prioritize drugs with both strong absolute NES and low q-values
- Look for convergent MOA themes among top drugs (e.g., multiple kinase inhibitors or HDAC inhibitors)
- Drugs opposing the disease signature (negative NES) are potential therapeutic candidates
- Drugs mimicking the disease signature (positive NES) may share biological mechanisms with the condition
- Cross-reference drug targets with differentially expressed genes for mechanistic insight
- Consider that L1000 profiles are cell-line-derived and may not fully reflect in vivo pharmacology
