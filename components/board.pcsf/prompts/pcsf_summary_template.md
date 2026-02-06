## Task

Interpret the biological significance of the Prize-Collecting Steiner Forest (PCSF) network analysis results for the contrast **{contrast}** based on the network data provided below.

## Analysis Overview

**Contrast:** {contrast}

**Phenotype:** {phenotype}

The PCSF analysis reconstructs a subnetwork of interacting genes connecting the most differentially expressed genes for this contrast, using protein-protein interaction databases as the backbone.

## Network Summary

{network_summary}

## Hub Genes (by centrality)

The following genes have the highest centrality scores in the PCSF network, indicating topological importance as potential driver genes:

{hub_genes}

Focus on identifying:
1. Which hub genes are most likely to be functional drivers of the observed phenotype
2. Whether Steiner nodes (recruited connectors) point to hidden regulatory mechanisms
3. How the network topology relates to the biological contrast

## Enriched Pathways in Network

{network_pathways}

## Output Instructions

Always begin your summary with a heading in this exact format (maintain strict adherence):
**{contrast}** | PCSF network summary

Then synthesize the above information to explain:
- The key hub genes and their potential roles as drivers of the biological response
- How the network topology (hubs, connectors, clusters) reveals functional organization
- The biological significance of Steiner nodes recruited into the network
- How enriched pathways among network genes support a coherent biological interpretation
- Use quantitative metrics to support your interpretation (cite specific logFC, centrality scores, and pathway statistics when available)
