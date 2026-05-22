# Drug connectivity context

## What this analysis is

Drug connectivity mapping scores each drug's pre-computed perturbation
signature against the experimental gene expression signature, using
GSEA-style enrichment. Each drug receives a Normalized Enrichment Score
(NES) and significance values (p, q).

Three levels of aggregation are reported per contrast:

- **DSEA (Drug Set Enrichment Analysis)**: per-drug NES — multiple
  perturbation experiments per drug aggregated into a single score.
- **MOA-class enrichment**: NES recomputed over drugs sharing a
  pharmacological class (e.g., "HDAC inhibitor", "EGFR inhibitor").
- **Target enrichment**: NES recomputed over drugs hitting the same
  molecular target (gene-level).

## Databases (analysis backends)

The data block names the active backend in `analysis_type` /
`analysis_type_description`. Two families exist, with **opposite NES
direction semantics**:

- **L1000/activity** and **L1000/gene** — transcriptional connectivity
  (Connectivity Map / LINCS L1000 perturbation profiles). Activity uses
  Transcriptional Activity Scores aggregated across cell lines; gene
  uses per-cell-line 978-landmark log-fold-changes.
- **CTRPv2/sensitivity** and **GDSC/sensitivity** — pharmacological
  sensitivity connectivity. CTRPv2 uses AUC dose-response over ~481
  compounds × ~860 cancer cell lines; GDSC uses IC50 over ~367 clinical
  anti-cancer drugs × ~987 cell lines.

## NES direction semantics

Read the backend before naming therapeutic candidates.

- **L1000** (transcriptional): negative NES → drug *opposes* the
  experimental signature → reversal candidate. Positive NES → drug
  *mimics* it.
- **CTRPv2 / GDSC** (sensitivity): positive NES → experimental state
  *predicts vulnerability* (sensitivity candidate). Negative NES →
  predicted resistance.

## Definitions

- **NES**: signed enrichment score; semantics depend on backend (above).
- **Annotation coverage**: fraction of tested drugs that carry MOA /
  target metadata. Low coverage weakens MOA / target interpretation
  and is reported per contrast.
- **Tier (strong / moderate / weak / data-limited)**: per-contrast
  classification provided in the data block, derived from significant
  drug count, max |NES|, MOA / target hits and annotation coverage.
  **Authoritative — do not reclassify.**

## Quantitative thresholds

- |NES| ≥ 1.5 strong; 1.0–1.5 moderate; 0.5–1.0 mild; < 0.5 minimal.
- q < 0.05 significant (headline claims).
- p < 0.05 with q ≥ 0.05 nominal (mention with hedging).
- Otherwise unsupported (suppress unless explicitly required).
- Numbers live in the data block, not the prose.

## Evidence hierarchy

When building the narrative for each contrast, walk this hierarchy:

1. **Supported MOA terms** (from the "Interpretation evidence summary"
   lines in the data block) — primary frame per contrast.
2. **Corroborating targets** — validate the mechanism at gene level.
3. **Exemplar drugs** — anchor the mechanism with one or two named
   compounds matching a supported MOA term.
4. **Raw top drugs** — exploratory only; use when no MOA signal exists.

## Caveats to keep visible

- Connectivity reflects cell-line perturbation signatures, not direct
  clinical efficacy.
- Dose, exposure time and cell context of reference profiles may diverge
  from the experimental biology.
- Annotation coverage gates MOA / target claims — when coverage is low,
  flag it before recommending mechanisms.
