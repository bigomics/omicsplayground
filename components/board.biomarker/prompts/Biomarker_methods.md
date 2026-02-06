## EXPERIMENT CONTEXT

This analysis uses machine learning-based feature selection to identify potential biomarkers for classification or prediction of phenotypic outcomes.

Experiment: {experiment}

## BIOMARKER-SPECIFIC GUIDANCE

- Biomarker discovery aims to identify the most discriminative molecular features (genes, proteins, metabolites) that distinguish between phenotypic groups
- Multiple machine learning algorithms are used in parallel: sparse partial least squares (sPLS), elastic nets, random forests, and extreme gradient boosting (XGBoost)
- Results are combined via cumulative ranking across all methods to select robust, consensus biomarkers
- The importance score reflects how consistently a feature is ranked as important across the different algorithms
- A decision tree is fitted using the top features to provide an interpretable classification model

## QUANTITATIVE INTERPRETATION GUIDELINES

Feature importance metrics:
- Cumulative rank: higher values indicate a feature was consistently ranked as important across multiple methods
- Features ranked highly by all methods are the most reliable biomarker candidates
- Features ranked highly by only one method may be algorithm-specific and less generalizable

Machine learning methods used:
- sPLS (sparse Partial Least Squares): identifies features with the strongest covariance with the outcome; good for correlated features
- Elastic net (LASSO + Ridge): performs variable selection via regularization; favors sparse solutions
- Random forest: ensemble of decision trees; captures non-linear relationships and interactions
- XGBoost (Extreme Gradient Boosting): sequential ensemble method; captures complex patterns and feature interactions

Integration rules:
- Prioritize features that appear in the top ranks across multiple methods (high cumulative rank)
- Consider biological plausibility: known disease-associated genes strengthen biomarker candidacy
- Features appearing in the decision tree are particularly informative for classification
- Look for functional coherence among top biomarkers (e.g., members of the same pathway)
- A single dominant biomarker may indicate a strong biological signal, while many weak biomarkers may indicate a distributed or noisy signal
- Consider that biomarker panels (combinations of features) often outperform individual markers
