# MOFA context

## Definitions

- **MOFA (Multi-Omics Factor Analysis)**: an unsupervised factor model
  that decomposes a stack of omics matrices (sharing the same samples)
  into latent factors. Each factor is a linear combination of features
  across one or more omics views; the factor captures a shared or
  view-specific source of variation.
- **Factor**: a latent variable indexed `Factor1`, `Factor2`, …. Each
  factor has per-feature loadings (weights) and per-sample scores. The
  sign of a factor is arbitrary — direction-of-effect statements must be
  paired with the trait or contrast they refer to.
- **Weight (loading)**: the coefficient mapping a feature to a factor.
  |weight| measures how strongly the feature contributes to the factor;
  the sign indicates the direction relative to the factor axis. The data
  block reports it under the `Contribution` column.
- **View**: one omics modality (transcriptome, proteome, metabolome,
  geneset, …). A factor that has high-|weight| features in more than one
  view is a **cross-view factor** — the strongest evidence of integrated
  multi-omics biology MOFA provides.
- **Variance explained**: per-factor, per-view fraction of variance
  captured. A factor with high variance across multiple views captures
  shared biology; one with high variance in a single view is
  modality-specific.
- **NES (Normalised Enrichment Score)**: GSEA effect size for a pathway
  on a factor's loading vector. The data block reports `Direction &
  strength` (sign + magnitude as a verbal label); the raw NES is not
  shown.
- **padj (adjusted p-value / q-value)**: enrichment significance for a
  pathway. The data block reports this under the `Significance` column
  as a verbal label.
- **Centrality**: per-feature score from the cross-view network — how
  central the feature is to the multi-omics module supporting the
  factor. The data block reports it under the `Network role` column.

## Verbal-label inventory (authoritative — used throughout the data block)

The data block never shows raw numbers for these six quantities; it shows
verbal labels at fixed thresholds. Treat the label as authoritative.

- **Factor-trait correlation** (`r`): `strongly correlated` (|r| ≥ 0.9) /
  `correlated` (≥ 0.7) / `moderately correlated` (≥ 0.5) /
  `weakly associated` (< 0.5). Sign mirrors to
  `… anti-correlated`. Missing → `no trait correlation`.
- **Pathway significance** (`padj`): `highly significant` (padj < 1e-6) /
  `significant` (< 1e-3) / `nominally significant` (< 0.05) /
  `not significant`. Missing → `not tested`.
- **Feature contribution** (`|weight|`): `dominant` (|w| ≥ 1.0) /
  `strong` (≥ 0.5) / `modest` (≥ 0.25) / `minimal` (< 0.25).
- **Pathway direction & strength** (`NES`): `up-strong` / `up-moderate` /
  `up-weak` / `up-nominal` (NES ≥ 0) and the mirrored `down-…` family
  (NES < 0). Magnitude breaks at |NES| = 2.5 / 1.5 / 0.5.
- **Per-view variance explained**: `major` (≥ 30%) / `moderate` (≥ 15%) /
  `minor` (≥ 5%) / `negligible` (< 5%). The raw percentage is shown as
  a parenthetical anchor on the lead view only.
- **Network role** (centrality): `hub` (≥ 0.8) / `central` (≥ 0.6) /
  `intermediate` (≥ 0.3) / `peripheral` (< 0.3).
- **Strong / Moderate / Weak signal**: tier classification provided in
  the data block, derived from significant-pathway count, max |NES|,
  and max |weight|. **Authoritative — do not reclassify.**
- **Top factors**: the subset selected for detailed reporting. Out of
  `n_factors_total` extracted, only `n_factors_used` are described in
  detail; the rest are aggregated under "Minor units".

## How to read direction

A factor's sign is mathematically arbitrary. State direction relative to
a trait or sample group: *"Factor 1 is elevated in treated samples"* —
never just *"Factor 1 is positive"*. The data block carries direction
via trait correlation labels, not via the bare sign of a weight or
score.

## Top-feature reporting

- Group features by view of origin when possible; weave functional
  context into prose rather than listing.
- Symbol names in italics: *IL6*, *TP53*.
- Never list more than ~8 feature names in a single paragraph.
- Prefer named features over feature IDs — if `symbol` is present in
  the data block, use it.

## Pathway-enrichment reporting

- Translate MSigDB / Reactome / GO IDs into plain-English clusters.
  Bad: *"GO_0007049 enriched"*. Good: *"Cell cycle machinery enriched"*.
- A factor without significant enrichment is informative — report it
  as such, do not omit the factor.
- NES sign on a factor's GSEA result inherits the factor's arbitrary
  sign; pair the direction statement with a trait or group, not with
  the bare NES sign.

## Cross-view evidence

- A feature with high |weight| in two or more views is the strongest
  multi-omics signal — call it out explicitly.
- Variance-explained patterns that span views ground "integrated
  biology" claims; single-view variance-dominated factors are
  modality-specific and must be hedged accordingly.

Experiment: {{experiment}}
