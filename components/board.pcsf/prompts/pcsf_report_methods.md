## Methods

The present analysis was conducted in the context of the experiment **{{experiment}}**.
Network reconstruction was performed using the Prize-Collecting Steiner Forest (PCSF)
framework, which identifies parsimonious subnetworks that connect high-prize molecular
entities over an interaction backbone. Differential signal was encoded as node prizes,
while candidate interactions were derived from STRING/GRAPHITE-based prior knowledge
networks. For each contrast, the PCSF optimization was solved under board-standard
settings to balance prize retention against edge-cost minimization and network
compactness. The resulting subnetworks contain terminal nodes, corresponding to entities
directly supported by the observed differential signal, and Steiner nodes, corresponding
to inferred connectors introduced by the optimization to bridge high-prize regions.
For each reconstructed network, we quantified structural properties including node count,
edge count, connected-component structure, and graph density. Candidate driver
prioritization was performed using centrality-based hub ranking to identify topologically
influential genes. To contextualize direct versus inferred biology, we examined
terminal-to-Steiner composition across reconstructed solutions. Pathway-level support was
then assessed by evaluating overlap between network genes and geneset-level fold-change
signals, enabling a mechanistic summary that links topological structure to functional
interpretation. PCSF-derived subnetworks are hypothesis-generating models and should
therefore be interpreted as mechanistic priors rather than definitive causal maps.
Although centrality-informed prioritization is useful for target nomination, high
centrality alone is insufficient to establish causal relevance. Consequently, orthogonal
validation through perturbation assays, replication cohorts, or complementary multi-omics
evidence remains necessary before translational decision-making. This methodological
framework follows the original PCSF formulation [1].

[1] Tuncbag N, Braunstein A, Pagnani A, Huang SSC, Chayes J, Borgs C, Zecchina R,
Fraenkel E. *Simultaneous reconstruction of multiple signaling pathways via the
Prize-Collecting Steiner Forest problem*. Journal of Computational Biology.
2013;20(2):124-136.

_This report was generated with OmicsPlayground (BigOmics, {{date}})._
_Note: AI-generated interpretation may contain inaccuracies and must be independently verified._
