## Task

Interpret the biological function and significance of WGCNA module **{module}** based on the enrichment analysis and correlation results provided below.

## Module Overview

**Module:** {module}

**Associated Phenotypes:** {phenotypes}

{module_stats}

The module eigengene shows correlation with the phenotypes listed above. Consider how the module's biological function might relate to these experimental conditions.

## Enrichment Analysis Results

The following pathways and gene sets show significant enrichment in this module:

{genesets}

Focus on identifying:
1. The dominant biological theme or process
2. How multiple enriched pathways relate to each other
3. Potential mechanistic connections to the associated phenotypes

## Hub Genes
{keygenes_section}

## Output Instructions

Always begin your summary with a heading in this exact format (maintain strict adherence):
**{module} module** | Correlated phenotypes: {phenotypes}

Then synthesize the above information to explain:
- The primary biological function of this module
- How the enriched pathways collectively support this interpretation
- The biological relevance to the experimental phenotypes
- How hub genes might drive or regulate this module's function
- Use quantitative metrics to support your interpretation (cite specific scores, q-values, MM, TS, LogFC, and centrality values when available)

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
