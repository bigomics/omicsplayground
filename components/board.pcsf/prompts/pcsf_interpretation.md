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

## PCSF-SPECIFIC GUIDANCE

- The Prize-Collecting Steiner Forest algorithm identifies high-confidence subnetworks connecting differentially expressed genes through known protein-protein interactions
- Genes are assigned "prizes" based on their differential expression (fold change or correlation), and the algorithm finds an optimal subnetwork that balances including high-prize nodes against the cost of connecting them
- The backbone network is constructed from STRING and GRAPHITE pathway databases, representing known physical and functional protein interactions
- Steiner nodes are genes not in the original differential expression set but recruited by the algorithm as connectors; these often represent hidden regulators or signaling intermediaries
- Terminal nodes are the original differentially expressed genes that anchor the network
- Hub genes are identified by centrality measures (e.g., page-rank) and represent potential driver genes that coordinate the biological response
- The network topology (clusters, bridges, hubs) reveals functional modules and their interconnections

## QUANTITATIVE INTERPRETATION GUIDELINES

Network metrics:
- logFC > 1: strong differential expression; 0.5-1 moderate; 0.2-0.5 mild; < 0.2 minimal
- logFC < -1: strong downregulation; -0.5 to -1 moderate; -0.2 to -0.5 mild; > -0.2 minimal
- High centrality score: gene is topologically important (potential driver or hub)
- Steiner nodes with high centrality: likely hidden regulators not detected by expression alone

Interpretation rules:
- Hub genes with high centrality AND high |logFC| are the strongest candidates for driver genes
- Steiner nodes (not differentially expressed themselves) that appear as hubs suggest hidden regulatory mechanisms
- Clusters of connected genes often represent coherent functional modules or pathways
- The ratio of Steiner to terminal nodes indicates how much the algorithm relied on known interactions to connect the network
- Enriched pathways identified among network nodes provide functional context for the network modules
- Consider both the direction (up/down) and magnitude of fold changes when interpreting network modules

## Language Rules

- Use language like "suggests", "is consistent with", and "prioritizes for follow-up".
- Do not claim proven causality or clinical efficacy.
- Keep uncertainty explicit when evidence is sparse.

## Actionability Rules

- Recommend concrete, feasible follow-ups tied to reported hubs/pathways.
- Include at least one orthogonal validation strategy.
