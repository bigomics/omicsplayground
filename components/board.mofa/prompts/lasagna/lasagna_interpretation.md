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

## LASAGNA-SPECIFIC GUIDANCE

LASAGNA (Layered Stacked Analysis of Gene Network Architecture) builds a
stacked multi-layer network from multiple omics modalities for a single
experimental contrast. Each omics type contributes a layer of nodes and
intra-layer edges; cross-layer edges connect features from different
modalities that show coordinated behaviour.

Key concepts:
- **Inter-layer edges** represent cross-omics associations (e.g., a
  transcript correlated with a protein). These are the integrative signal
  — they reveal biology invisible to any single omics layer.
- **Intra-layer edges** represent within-modality associations (e.g.,
  co-expressed genes). These provide modality-specific context.
- **Cross-layer hub nodes** have high centrality across multiple layers
  and are stronger candidates for biological drivers than nodes confined
  to a single layer.
- **Layer participation** counts indicate how many nodes each omics
  modality contributes. Balanced participation suggests genuine
  multi-omics integration; single-layer dominance means the network is
  mostly driven by one modality.

### Network density and confidence

- A dense network (many edges per node, multiple layers connected) supports
  confident interpretation of cross-layer biology
- A sparse network (few edges, low inter-layer count) should be interpreted
  conservatively — the integration signal may be weak or driven by noise
- Networks with zero or very few inter-layer edges provide no cross-omics
  integration evidence; state this explicitly

### Quantitative interpretation

Node metrics:
- Centrality: high centrality nodes are topological hubs; a node with high
  centrality in the multi-layer network is better connected across the
  integration than a high-centrality node in a single layer
- Fold change: the differential expression magnitude in the selected
  contrast; always state the contrast when citing logFC values

Edge metrics:
- Edge weight reflects the strength of association between connected nodes
- Inter-layer edge weight is the key integrative metric — prioritize strong
  cross-layer connections over strong within-layer connections

### Interpretation rules

- Inter-layer edges indicate cross-omics association patterns, not causal
  direction
- Node-level fold changes and correlations support prioritization of
  candidate drivers
- Sparse connectivity or low inter-layer signal should be interpreted
  conservatively
- A node appearing as a hub in multiple layers is stronger evidence than
  a hub in a single layer

## Writing Style: Pathway and Gene-Set References

Write the narrative using only natural biological language:

1. **No raw identifiers in prose.** Never embed database accessions or raw
   pathway names into the narrative text.

2. **Thematic grouping with inline references.** Cluster enriched pathways
   into biological themes. Pair every theme mention with bracketed [n]
   references. Never group more than 3 references in a single bracket.

3. **Plain-English descriptions.** Refer to pathways by their biological
   meaning rather than database labels.

### Enrichment References

- Reference enrichment terms via [n] brackets in the narrative
- Use plain-English descriptions in prose: "DNA replication [1,2]"
  not "REACTOME_DNA_REPLICATION [1]"
- The Data References section at the end maps [n] to exact term names
- When 0 terms are significant, state it explicitly
