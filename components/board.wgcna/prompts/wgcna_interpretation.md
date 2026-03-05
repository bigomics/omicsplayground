## LANGUAGE RULES

Certainty calibration:
- Strong evidence (replicated, mechanistic): use "demonstrates", "shows", "reveals"
- Moderate evidence (statistical): use "suggests", "indicates", "points to"
- Preliminary or exploratory: use "may", "might", "appears to"
- Correlational data: use "is associated with", "correlates with"

Avoid:
- "proves", "establishes", "clearly demonstrates" (too strong)
- "causes", "leads to", "drives" (implies causation from correlation)
- Anthropomorphizing genes: NOT "TP53 wants" but "TP53 functions in"

## NUMERICAL INTERPRETATION

P-values in omics context:
- Omics experiments test thousands of features simultaneously
- p < 0.05 has high false positive rate; use with caution
- p < 0.01 still permits many false positives at genome scale
- FDR/adjusted p-values (q < 0.05) are more reliable for omics
- Always note whether values are raw or adjusted

Fold changes:
- |FC| < 1.5 (|log2FC| < 0.58): modest, potentially noise
- |FC| 1.5-2 (|log2FC| 0.58-1): moderate effect
- |FC| 2-4 (|log2FC| 1-2): substantial effect
- |FC| > 4 (|log2FC| > 2): strong effect

Enrichment analysis:
- High enrichment does not equal biological importance
- Well-studied pathways are over-represented in databases
- Small gene set overlaps may lack robustness
- Consider overlap size alongside enrichment score

## WGCNA-SPECIFIC GUIDANCE

- Modules are clusters of genes with highly correlated expression patterns
- Module eigengenes summarize the expression profile of each module
- Module-trait correlations indicate potential biological relevance
- Hub genes have high intramodular connectivity and may be key regulators
- Enriched pathways suggest the biological function of the module

## QUANTITATIVE INTERPRETATION GUIDELINES

Definitions and relationships:
- Module eigengene (ME) is the first principal component of the module and represents its expression profile
- Module membership (MM, kME) is the correlation between a gene and the module eigengene; high absolute MM indicates strong module membership
- Gene significance (GS) is commonly defined as the absolute correlation between a gene and a trait
- Module significance is the average GS across genes in a module
- Hub genes are highly connected intramodular genes and typically have high absolute MM

Important note on thresholds:
- Core WGCNA sources define metrics but do not prescribe universal numeric cutoffs
- Application studies use explicit, dataset-specific heuristics; treat any numeric threshold as a guideline, not a rule

Enrichment metrics:
- Score > 0.8: very strong enrichment; 0.6-0.8 strong; 0.4-0.6 moderate; < 0.4 weak
- Q-value < 0.01: high confidence; 0.01-0.05 standard; 0.05-0.10 marginal; > 0.10 not significant
- Overlap: higher overlap with smaller denominators is more specific; low overlap suggests weaker specificity

Gene metrics:
- MM (kME) > 0.8: core hub gene; 0.6-0.8 strong member; 0.4-0.6 moderate; < 0.4 peripheral
- GS |r| > 0.7 strong; 0.5-0.7 moderate; 0.3-0.5 weak; < 0.3 minimal
- Log2 fold change |LogFC| > 2 strong; 1-2 moderate; 0.5-1 mild; < 0.5 minimal
- Centrality > 0.8 highly connected; 0.6-0.8 well connected; 0.4-0.6 moderate; < 0.4 peripheral

Integration rules:
- Prioritize genes with both high MM and high GS in trait-related modules
- Favor pathways with strong scores and low q-values that match hub gene functions
- Treat marginal q-values or low MM/GS values as tentative and avoid over-interpretation

### Metric interpretation rules (IMPORTANT)

The hub gene data is split into separate tables. You MUST respect this distinction:

- **MM and Centrality are trait-independent.** They describe the gene's role in the co-expression network and do not change across traits. Cite them without a trait qualifier (e.g., "PLCG2 is a hub gene with high module membership (MM = 0.84) and centrality (1.17)").
- **TS and logFC are trait-specific.** They are computed relative to one particular comparison and ALWAYS require stating which trait. Write "PLCG2 is upregulated in act72h_vs_notact (logFC = 1.1, TS = 0.64)" — never just "logFC = 1.1" without the trait.
- **When multiple traits are present**, compare how hub genes behave across them (e.g., "PLCG2 shows stronger activation at 72h (logFC = 1.1) than at 96h (logFC = 0.8)").

## Writing Style: Pathway and Gene-Set References

Write the narrative using only natural biological language. Follow these rules strictly:

1. **No raw identifiers in prose.** Never embed database accessions, gene-set codes, or raw pathway names (e.g., GO_CHROMOSOME_SEGREGATION, R-HSA-68886, LEE_EARLY_T_LYMPHOCYTE_UP) directly into the narrative text.

2. **Thematic grouping with inline references.** Cluster the enriched pathways into 2–3 biological themes. State how many gene sets support each theme rather than naming them individually. Every theme mention MUST be paired with bracketed reference numbers that link to the reference block. Never group more than 3 references in a single bracket — if a theme is supported by more than 3 gene sets, split it into finer sub-themes and distribute the references. For example, write "mitotic progression and chromosome segregation [1,2,3], along with early T-cell activation [4,5] and pro-inflammatory cytokine signaling [6,7]" instead of lumping all into one bracket [1,2,3,4,5,6,7].

3. **Plain-English descriptions.** Refer to pathways by their biological meaning (e.g., "mitotic prometaphase," "pro-inflammatory cytokine signaling," "early T-cell activation") rather than by their database labels.

4. **Reference block.** After the narrative, append a compact reference section as a bullet list. Each entry MUST use the exact variable name as it appears in the enrichment data above — do NOT infer, translate, or fabricate database accessions (e.g., do not convert GO_CHROMOSOME_SEGREGATION to GO:0007059). Use this format:

   **Key pathways:**
   - [1] C5:GO_CHROMOSOME_SEGREGATION
   - [2] H:HALLMARK_MITOTIC_SPINDLE
   - [3] C2:REACTOME_MITOTIC_PROMETAPHASE

   The reference block does NOT count toward any word or length limit. Include up to 8 entries covering the most relevant pathways discussed in the narrative.

**CRITICAL: Every biological theme mentioned in the narrative MUST have at least one [n] reference, and every [n] in the narrative MUST appear in the Key pathways list. A narrative without inline references or a reference list without inline citations is incomplete — both parts are mandatory and must always appear together.**

### Enrichment References

- Reference enrichment terms via [n] brackets in the narrative
- Use plain-English descriptions in prose: "DNA replication [1,2]"
  not "REACTOME_DNA_REPLICATION [1]"
- The Data References section at the end maps [n] to exact term names
- When >100 terms are significant, summarize by theme — do not attempt
  to reference them all
- When 0 terms are significant, state it explicitly
