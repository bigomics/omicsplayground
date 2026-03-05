## EXPERIMENT CONTEXT

This analysis uses feature-level dimensionality reduction (UMAP) to explore co-expression structure among genes and gene sets.

Experiment: {{experiment}}

## FEATURE MAP-SPECIFIC GUIDANCE

- Feature maps cluster genes (or gene sets) based on pairwise co-expression, using covariance as the distance metric
- UMAP (Uniform Manifold Approximation and Projection) is computed from either the normalized log-expression matrix (logCPM) or the log-foldchange matrix (logFC)
- Features that are close together on the UMAP have high pairwise covariance, meaning they tend to change together across conditions
- This is feature-level clustering (grouping genes by behavior) as opposed to sample-level clustering (grouping samples by profile)
- Gene-level maps show individual gene positions; geneset-level maps show pathway/gene-set positions
- Coloring by fold-change signature allows visual comparison of the global effect between different contrasts or phenotype conditions
- The rms(FC) metric is the root-mean-square fold change across all contrasts, summarizing overall variability of a feature

## QUANTITATIVE INTERPRETATION GUIDELINES

Feature map metrics:
- rms.FC > 1.0: highly variable feature across contrasts; 0.5-1.0 moderately variable; 0.2-0.5 mildly variable; < 0.2 stable
- Fold change (logFC) per contrast: positive values indicate upregulation, negative values indicate downregulation
- Standard deviation (sd.X): higher values indicate more variable features across samples

Spatial interpretation:
- Features clustered tightly together are co-regulated (high covariance) and likely belong to the same biological module or pathway
- Features at the periphery with high rms.FC are potential drivers of variation between conditions
- Dense clusters suggest coordinated gene programs; isolated features may have unique regulatory patterns
- When colored by a specific contrast, red/blue gradients reveal which regions of the feature space respond to that condition

Integration rules:
- Identify coherent gene clusters that respond similarly across contrasts
- Look for known gene families or pathway members clustering together as validation
- Features with high fold change at the map periphery are candidate driver genes
- Compare multiple contrast colorings to find features with condition-specific versus pan-condition responses
- Consider the number of nearby features: isolated high-FC features may be less biologically meaningful than clustered ones
