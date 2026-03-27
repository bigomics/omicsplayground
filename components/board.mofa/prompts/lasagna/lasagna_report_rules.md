## What This Analysis Reveals

LASAGNA constructs a stacked multi-layer network from multiple omics
modalities for a single experimental contrast. The network reveals
which biological programmes are supported by coordinated cross-omics
evidence — not just differential expression in one modality, but
convergent signal across transcriptomics, proteomics, metabolomics,
or other layers.

The report translates this network topology into a biological narrative:
what cross-layer modules dominate, which nodes bridge modalities, and
how confident we can be in the integration signal.

## Analytical Unit: the Cross-Layer Functional Module

Since LASAGNA operates on a single contrast, the analytical units are
the 1-3 dominant **biological programmes** emerging from the multi-layer
network topology — not contrasts, not individual nodes.

Each `### heading` in Main Findings names a biological programme, not a
method concept:

```
GOOD: ### Kinase–Transcription Factor Relay
GOOD: ### Lipid Metabolism Convergence
BAD:  ### Network Evidence
BAD:  ### Top Hub Genes
```

Use network density and inter-layer signal to classify modules:
- **Well-supported**: multiple inter-layer edges, balanced layer
  participation, high-centrality hubs → full narrative treatment
- **Sparse or single-layer-dominated**: few inter-layer edges, nodes
  confined to one modality → aggregate into `### Minor Signals` with
  explicit caveats about weak integration evidence

Omit `### Minor Signals` entirely if all modules are well-supported.

## Per-Unit Content: Three Narrative Beats (always prose)

**Beat 1 — The biological programme (lead sentence)**
State what functional module the network topology reveals. Lead with the
biology, not the node names.

> GOOD: "The network centres on a kinase–transcription factor relay where
> phosphoproteome hubs converge on transcriptomic effectors of cell cycle
> entry."
> BAD: "8 nodes with high centrality were identified across two layers."

**Beat 2 — The cross-layer evidence**
Name the top nodes that anchor the module. For each, weave centrality,
layer of origin, fold change, and functional annotation into prose. When
a node bridges multiple layers via inter-layer edges, say so explicitly —
this is the strongest evidence LASAGNA provides.

> GOOD: "*CDK2* (centrality = 0.92, phosphoproteome, logFC = +1.8) links
> to *E2F1* (centrality = 0.85, transcriptome) via a strong cross-layer
> edge, consistent with phosphorylation-driven transcriptional
> activation."
> BAD: "CDK2 centrality 0.92. E2F1 centrality 0.85."

**Beat 3 — Layer convergence and confidence**
Does the module draw from multiple layers or is it dominated by one?
State the layer participation balance. Flag confidence explicitly when:
- The module is driven by a single omics layer
- Inter-layer edges are few or absent
- Network is overall sparse (< 50 nodes)

## Cross-Unit Synthesis Paragraph (REQUIRED when > 1 module)

After all module subsections, one paragraph with no heading, describing:
- How modules relate to each other in the network
- Whether cross-layer bridges connect otherwise separate modules
- Whether the network tells a coherent multi-omics story or fragments
  into modality-specific clusters
- The overall biological narrative for the contrast

Every cross-unit claim must be grounded in values from the data.

## Discussion

1-2 paragraphs interpreting the overall multi-layer integration:
- What does the cross-layer topology reveal that single-omics would miss?
- Which nodes are genuine multi-omics hubs vs single-layer artefacts?
- Caveats once, briefly: network sparsity, layer imbalance (unequal node
  counts across modalities), PPI backbone bias, correlation ≠ causation

## Word Limit

600-1000 words total across Main Findings + Discussion + Conclusion.

## Hard Constraints

- Use only facts explicitly present in the input data.
- Do NOT infer causality from network topology alone.
- If the network is sparse or lacks inter-layer edges, state that
  explicitly and keep confidence modest.
- Do NOT synthesize numeric summaries not present in the input data.
- Do NOT reproduce node tables or edge tables — the researcher already
  has them in the interface.
- Do NOT create bullet lists of node names with centrality scores; weave
  them into prose.
