## Task

Interpret the biological significance of the biomarker discovery results for the prediction target **{target}** based on the feature importance data provided below.

## Analysis Overview

**Prediction target:** {target}

**Phenotype:** {phenotype}

The biomarker analysis identifies the most discriminative molecular features for classifying or predicting the selected phenotype. Multiple machine learning algorithms (sPLS, elastic nets, random forests, XGBoost) are used, and their importance rankings are combined into a cumulative score.

## Top Biomarker Features

The following features have the highest cumulative importance scores across all machine learning methods:

{top_features}

Focus on identifying:
1. Which features are the strongest and most consistent biomarker candidates
2. Biological functions or pathways represented by the top features
3. Potential mechanistic connections between the top biomarkers and the phenotype

## Model Performance

{model_performance}

## Output Instructions

Always begin your summary with a heading in this exact format (maintain strict adherence):
**{target}** | Biomarker summary

Then synthesize the above information to explain:
- The most promising biomarker candidates and why they are biologically relevant
- How the top features relate to the phenotype being predicted
- Whether the top biomarkers suggest specific biological mechanisms or pathways
- The robustness of the biomarker selection (agreement across methods)
- Use quantitative metrics to support your interpretation (cite specific importance scores and method-level rankings when available)
