# LASAGNA context

## Definitions

- **LASAGNA (Layered Stacked Analysis of Gene Network Architecture)**:
  builds a stacked multi-layer network from multiple omics modalities
  for a single experimental contrast. Each modality contributes a
  layer of nodes; cross-layer edges connect features from different
  modalities that show coordinated behaviour.
- **Layer**: one omics modality (transcriptome, proteome, metabolome,
  …). Nodes inherit a layer label from their feature prefix.
- **Inter-layer edge**: an edge connecting nodes from different
  modalities. **These are the integrative signal — they reveal biology
  invisible to any single omics layer.**
- **Intra-layer edge**: an edge connecting nodes within one modality.
  Provides modality-specific context.
- **Module (community)**: a densely-connected subgraph of the network,
  identified by community detection (Louvain on the undirected version
  of the graph). Modules with members from multiple layers are
  **cross-layer modules**; modules confined to one modality are
  **single-layer modules** and are aggregated under Minor units.
- **Hub node**: a node with high centrality within its module. The
  strongest hubs are those that also have at least one inter-layer
  edge — these are **cross-layer hubs** and are the strongest
  candidates for biological drivers.
- **Edge strength (ρ)**: per-edge correlation between the two endpoint
  features. In the data block, raw `ρ` is replaced with one of seven
  verbal labels mirroring the omicsai `r` verbaliser.
- **Strong / Moderate / Weak signal**: module tier classification
  provided in the data block, derived from node count, cross-layer
  edge count, and number of layers spanned. **Authoritative — do not
  reclassify.**
- **Top modules**: the subset selected for detailed reporting. Out of
  `n_modules_total` detected, only `n_modules_used` are described in
  detail; the remainder are aggregated under "Minor units".

## How to read direction

LASAGNA edges are undirected; no directional claim can be made from
topology alone. State biological hand-offs only when the underlying
mechanism (e.g. phosphorylation → transcription) is well-established
for the named hub pair and only as a hypothesis, never as a finding.

## Hub node reporting

- Group hubs by layer of origin when possible; weave functional
  context into prose rather than listing.
- Symbol names in italics: *CDK2*, *E2F1*.
- Never list more than ~8 hub names in a single paragraph.
- Cross-layer hubs (hubs with at least one inter-layer edge) are the
  strongest evidence — call them out explicitly.

## Verbal-label inventory (authoritative — used throughout the data block)

The data block never shows raw numbers for these quantities; it shows
verbal labels at fixed thresholds. Treat the label as authoritative.

- **Edge strength** (`ρ`): `strongly correlated` (|ρ| ≥ 0.9) /
  `correlated` (≥ 0.7) / `moderately correlated` (≥ 0.5) /
  `weakly associated` (< 0.5). Sign mirrors to `… anti-correlated`.
- **Network role** (centrality): `hub` (≥ 0.8) / `central` (≥ 0.6) /
  `intermediate` (≥ 0.3) / `peripheral` (< 0.3).
- **Layer participation** (per layer, % of module nodes): `dominant`
  (≥ 60%) / `major` (≥ 30%) / `minor` (≥ 10%) / `marginal` (< 10%).
- **Cross-layer connectivity** (per module): `rich-bridged` (≥ 10
  inter-layer edges) / `well-bridged` (≥ 5) / `lightly-bridged` (≥ 1) /
  `single-layer` (0).
- **Network density** (whole graph): `dense` (≥ 0.1) / `moderate`
  (≥ 0.02) / `sparse` (≥ 0.005) / `fragmented` (< 0.005).

## How to read per-module sections

Each module's data block carries five evidence channels — treat them
together, not in isolation:

1. **Layer participation** + **Cross-layer connectivity**: tell you
   whether the module's biology spans modalities (the integrated
   signal LASAGNA was designed to surface) or is confined to one.
2. **Trait coordination**: the mean signed ρ across all module
   members vs the contrast. A `strongly correlated` module is
   responsive to the contrast; a `weakly associated` module is
   present in the network but not driven by the comparison.
3. **Top hub nodes**: the named anchors; cite hubs flagged
   `Cross-layer: yes` first — they are the strongest evidence.
4. **Top cross-layer edges**: the explicit bridges between layers.
5. **Partner modules**: modules sharing hubs / members. A shared hub
   between two strong modules is grounds for a single integrated
   theme spanning both; weave them in the same `### Theme` heading.
