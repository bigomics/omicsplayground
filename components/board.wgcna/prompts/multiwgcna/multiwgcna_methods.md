## Methods

Multi-Omics WGCNA was performed on {{n_features_wgcna}} features across {{n_layers}} omics layers and {{n_samples}} samples. WGCNA modules were inferred separately within each layer using {{network_type}} networks with soft-thresholding power {{power}}, minimum module size {{min_mod_size}}, merge cut height {{merge_cut_height}}, and minimum kME {{min_kme}}.

Layer-specific module eigengenes were correlated with phenotype indicators using Pearson correlation. Pathway enrichment was evaluated for each module against {{n_genesets_tested}} gene sets. Hub features were ranked by module membership (MM), defined as the correlation between a feature profile and its module eigengene.

Cross-layer integration was summarized by correlating module eigengenes across layers and by extracting highly correlated cross-layer feature sets. These integrated views were used to identify convergent and opposing multi-omics programs.

### Compute Settings

| Parameter | Value |
|-----------|-------|
| Layers | {{layers}} |
| Soft-thresholding power | {{power}} |
| Minimum module size | {{min_mod_size}} |
| Merge cut height | {{merge_cut_height}} |
| Minimum kME | {{min_kme}} |
| Features used | {{n_features_wgcna}} |
| Samples | {{n_samples}} |
| Gene sets tested | {{n_genesets_tested}} |

### References

- Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics*. 2008;9:559.
- Zhang B, Horvath S. A general framework for weighted gene co-expression network analysis. *Stat Appl Genet Mol Biol*. 2005;4(1).
- Multi-layer integration in OmicsPlayground combines layer-specific WGCNA with cross-layer module correlation and LASAGNA-style graph integration.

_This report was generated with OmicsPlayground (BigOmics, {{date}})._
_Note: AI-generated interpretation may contain inaccuracies and must be independently verified._
