## Methods

This report summarizes pre-computed Drug Connectivity / DSEA outputs from Omics Playground.

For the present analysis, the experiment is "{{experiment}}" and the selected analysis type is "{{analysis_type}}", defined as follows: {{analysis_type_description}}. Drug-level connectivity is quantified using normalized enrichment scores (NES), nominal p-values, and false discovery rate (FDR)-adjusted q-values. Within this framework, negative NES values are interpreted as opposing signatures and therefore support reversal hypotheses, whereas positive NES values are interpreted as mimicking signatures and therefore support mechanistic similarity hypotheses. Mechanism-of-action and molecular-target convergence are derived from enrichment over annotated perturbagen mechanisms and targets and are used to assess coherence of biological interpretation rather than to establish causality.

The methodological basis for interpretation depends on the selected reference
resource. For Connectivity Map and L1000-based analyses, the interpretation is
grounded in the use of gene-expression signatures to connect small molecules,
genes, and disease states (Lamb et al., 2006; Subramanian et al., 2017). For
GDSC-based sensitivity analyses, the interpretation is grounded in
pharmacogenomic associations between molecular states and drug response in
cancer cell lines (Yang et al., 2013).

All findings in this report are intended for hypothesis generation and prioritization and require orthogonal validation before biological or translational claims are made.

### References

- Lamb J, Crawford ED, Peck D, et al. The Connectivity Map: using gene-expression signatures to connect small molecules, genes, and disease. *Science*. 2006;313(5795):1929-1935.
- Subramanian A, Narayan R, Corsello SM, et al. A next generation Connectivity Map: L1000 platform and the first 1,000,000 profiles. *Cell*. 2017;171(6):1437-1452.e17.
- Yang W, Soares J, Greninger P, et al. Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells. *Nucleic Acids Research*. 2013;41(Database issue):D955-D961.

_This report was generated with OmicsPlayground (BigOmics, {{date}})._
_Note: AI-generated interpretation may contain inaccuracies and must be independently verified._
