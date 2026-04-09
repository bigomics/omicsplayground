## **Methods**

WGCNA identified {{n_modules}} co-expression modules (excluding the grey/unassigned bin of {{grey_size}} {{feature_type}}).
Module-trait correlations were computed using Pearson correlation between module eigengenes and each binary trait indicator.
Pathway enrichment was tested against {{n_genesets_tested}} gene sets. Hub genes were ranked by module membership (MM).
Module preservation analysis assessment follows the approach described in Langfelder et al. (2011).

Module strength classification:
strong: ≥10 significant enrichments at q < 0.05 & |trait r| > 0.7, ≥50 genes;
moderate: 1-9 significant enrichments at q < 0.05;
weak: no significant enrichment.

### **Compute Settings**

| Parameter | Value |
|-----------|-------|
| Network type | {{network_type}} |
| Soft-thresholding power | {{power}} |
| Minimum module size | {{min_mod_size}} |
| Merge cut height | {{merge_cut_height}} |
| Minimum kME | {{min_kme}} |
| Features used | {{n_features_wgcna}} {{feature_type}} |
| Samples | {{n_samples}} |
| Gene sets tested | {{n_genesets_tested}} |

## **Bibliography**

- Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics*. 2008;9:559.
- Langfelder P, Luo R, Oldham MC, Horvath S. Is my network module preserved and reproducible? *PLoS Computational Biology*. 2011;7(1):e1001057.
- WGCNAplus R package: https://github.com/bigomics/WGCNAplus

#### **Disclaimer**: AI-generated content. May contain inaccuracies. Independent verification advised.
