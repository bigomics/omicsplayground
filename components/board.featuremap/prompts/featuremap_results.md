## Task

Interpret the biological significance of the feature map analysis for the contrast **{{contrast}}** based on the feature clustering and expression data provided below.

## Analysis Overview

**Contrast:** {{contrast}}

**Feature level:** {{feature_level}}

The feature map visualizes genes (or gene sets) in a 2D UMAP embedding based on co-expression structure. Features that are nearby on the map have high pairwise covariance and tend to change together across conditions.

## Top Variable Features

The following features show the strongest variation across contrasts (ranked by rms.FC or sd):

{{top_features}}

Focus on identifying:
1. The dominant biological themes among the most variable features
2. Whether the top features form coherent functional groups
3. Which features are most strongly affected by this specific contrast

## Nearby Features on the Map

{{nearby_features}}

## Output Instructions

Always begin your summary with a heading in this exact format (maintain strict adherence):
**{{contrast}}** | Feature map summary

Then synthesize the above information to explain:
- The primary gene modules or co-expression clusters suggested by the feature map
- Which biological processes are most strongly activated or suppressed in this contrast
- How the spatial clustering of features supports a coherent biological interpretation
- Key driver genes or gene sets that appear at the periphery with high fold changes
- Use quantitative metrics to support your interpretation (cite specific logFC values and feature statistics when available)
