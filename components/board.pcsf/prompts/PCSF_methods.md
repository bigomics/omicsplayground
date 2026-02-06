## EXPERIMENT CONTEXT

This analysis uses the Prize-Collecting Steiner Forest (PCSF) algorithm to reconstruct biological networks from differential expression data.

Experiment: {experiment}

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
