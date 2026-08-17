# Methylome Profiler

A Dashboard board for Illumina methylation array data, registered as
`MODULE.methylome` and shown for methylomics datasets. It answers per-sample
questions that need no experimental groups — epigenetic age, cell composition,
methylome character — and runs a covariate-adjusted EWAS when groups exist.

Run it locally with `make methylome` (see the bottom of this file).

## Capabilities

| Screen | Capability | Notes |
|---|---|---|
| **Sample ledger** | Per-sample bimodality, missingness, imprinted-DMR drift | Flagging is cohort-relative (3 MAD), not an absolute cut-off |
| | DNAm age per sample | |
| | Predicted sex, including aneuploidies | Refuses to predict when X/Y probes are absent |
| | Beta density per sample | Bimodality check against the expected 0.2 / 0.8 peaks |
| **Epigenetic age** | 10 clocks: Horvath, skinHorvath, Hannum, BNN, BLUP, EN, PedBE, Wu, Levine (PhenoAge), TL | via `methylclock`; falls back to the 5 in `wateRmelon` if absent |
| | Select by family or individually; coverage floor; recompute on demand | A clock below the floor is **withheld**, not estimated from partial probes |
| | DNAm age vs chronological age per clock, with *r* | |
| | Pairwise clock agreement | |
| | Age acceleration by phenotype, clock and phenotype both selectable | Residual of DNAm age on chronological age; t-test for two levels, Kruskal-Wallis above |
| | **Intrinsic** acceleration, adjusted for cell composition | Optional. Without it, a cohort whose groups differ in blood composition reports that as an accelerated methylome |
| | Per-clock coverage table | Names each clock and why it was withheld |
| **Methylome character** | Mean beta by relation to CpG island (island / shore / shelf / open sea) | |
| | Mean beta by position in gene (TSS1500 → 3'UTR) | |
| | Global mean beta, hypermethylated fraction | |
| **Cell composition** | Houseman reference projection, 9 panels: adult blood, 4 × cord blood, saliva, DLPFC | Estimator ships inside `methylclock` — no extra dependency |
| | Per-sample stacked bars, ordered by phenotype | |
| | Each cell type compared across groups, with a test | The confounding check |
| | Proportions table, CSV export | Feeds the EWAS model as covariates |
| **EWAS — Manhattan & hits** | **limma fitted in-app**, with an outcome selector | Replaces the unadjusted result stored at upload |
| | **Two-group contrast or continuous exposure** — age, BMI, pack-years, dose | One picker, two optgroups. For a continuous outcome Δβ is the slope: methylation change per unit |
| | Covariates from any sample column; cell-composition adjustment | Rank-deficiency and residual-df checks refuse a bad design |
| | Probe masking: common-SNP and cross-reactive | 52,116 SNP probes on 450K (90,084 on EPIC); 53,498 cross-reactive |
| | Fitted formula, n per group and masked count printed | Visible on every EWAS sub-tab |
| | Manhattan with adjustable threshold (FDR or nominal p), optional \|Δβ\| filter | Threshold line drawn at the largest passing p, omitted when nothing passes |
| | Hits table: gene, island context, gene region, **Δβ**, direction, p, FDR | Δβ on the beta scale, not the M-value logFC |
| **EWAS — Regions & pathways** | DMR calling | `dmrff` — runs on the summary statistics already produced |
| | **Region detail plot**: methylation vs genomic position, for the region clicked in the table | Per-sample lines, group (or tertile) means, island track, flanking context. Shows whether a region is coherent or one CpG dragging its neighbours |
| | Gene-set testing corrected for probes-per-gene | `missMethyl::gometh` |
| **EWAS — QQ & context** | QQ plot with λ, plus bacon bias and inflation | bacon is **diagnostic only** and does not adjust p-values |
| | Island-context enrichment of hits vs the tested background | Background is the probes tested in this contrast, not the whole array |
| | Per-CpG beta stripcharts of the top hits | Fixed [0,1] axis so the absolute methylation level stays visible |

Every panel is a platform `PlotModule` / `TableModule`, so the info and methods
popovers, maximize, and PNG / PDF / SVG / CSV download come from the platform
rather than being reimplemented.

## Array support

450K and EPIC are both supported, and **the array is inferred from the probe
IDs**, not from metadata: `pgx$meth_type` is accepted by `pgx.createPGX` but is
not stored on the pgx, so anything keyed on it silently falls through to 450K.
EPIC keeps roughly 453k of the 450K probes and adds ~413k of its own, so the
share of probes that exist only on EPIC identifies the array — and unlike a
probe count, that still works on a subset.

The array selection drives the annotation manifest (and therefore SNP masking)
and the `array.type` passed to `gometh`. `test/test-arrays.Rscript` checks the
whole chain on both: the real 450K dataset, and an EPIC dataset synthesised
from the EPIC manifest with planted signal.

**EPIC v2 is detected and refused, not analysed.** Its probe IDs carry
replicate suffixes (`cg00000029_TC21`) that match neither manifest, so the
EPIC-only fraction comes out at zero and the array would otherwise be reported
as 450K — annotation, SNP masking and the gene-set background all degrading
with nothing said. The suffix is detected instead and everything needing
annotation refuses with an explicit message. Supporting it properly needs
`IlluminaHumanMethylationEPICv2anno.20a1.hg38`, which is not installed and is
not in the configured repositories, plus logic to collapse replicate probes.

## Design decisions worth knowing

- **The model is fitted in the app, not read from the pgx.** The stored result
  is an unadjusted two-group comparison; in bulk tissue that is usually a test
  of cell composition rather than of the phenotype (Jaffe & Irizarry 2014,
  *Genome Biol* 15:R31).
- **Δβ is computed, never taken from `meta.fx`.** The stored effect size is a
  logFC on M-values, which is not interpretable as a change in methylation — a
  logFC of 1 means something different at β = 0.05 than at β = 0.5.
- **`dmrff` rather than `DMRcate`.** dmrff runs on summary statistics and its
  whole dependency list is `parallel` + `ggplot2`; DMRcate would pull bsseq,
  Gviz, biomaRt, missMethyl and an ExperimentHub package that downloads at
  runtime.
- **The |Δβ| filter defaults to off.** True differences in a population EWAS
  are typically under 5%, so a hard effect-size cut deletes the real signal
  (Mansell et al. 2019, *BMC Genomics* 20:366).
- **bacon reports, it does not correct.** Applying it to an unadjusted model
  would launder confounding into well-calibrated-looking p-values.
- **Expensive steps are explicit.** Clocks, deconvolution, the model, DMRs and
  gene sets each run on a button, not on every settings change.
- **Nothing picks its own variable.** The acceleration panel used to take the
  first sample column with 2–6 levels, which on a sheet whose first such column
  is the slide or plate silently produced a t-test against batch under an
  "age acceleration by phenotype" heading. Clock and phenotype are now both
  chosen; a panel that cannot name its variable refuses instead of guessing.
- **Covariates that do not vary are dropped by name.** A model fitted on the
  subset carrying a value for the outcome can find a covariate constant within
  that subset — the drop is reported, not silently absorbed into the intercept.

## Out of scope

These need raw IDAT intensities, which a beta matrix does not carry — β is a
ratio and the total intensity divides out. They are not backlog items:

- copy-number inference, and everything that depends on it
- detection p-values, failed-probe masking, per-sample call rate
- background correction and dye-bias equalisation (Noob, ssNoob)
- read-level heterogeneity (PDR, epipolymorphism) — needs bisulfite sequencing

Also unavailable: GrimAge, whose coefficients are not public.

## Files

| File | Contents |
|---|---|
| `R/methylome_ui.R`, `R/methylome_server.R` | Shell and wiring |
| `R/methylome_utils.R` | Shared helpers and the datatype guard |
| `R/methylome_ledger.R` | Sample ledger |
| `R/methylome_age.R` | Clocks |
| `R/methylome_character.R` | Genomic context profiles |
| `R/methylome_deconv.R` | Cell composition |
| `R/methylome_model.R` | Array inference, probe masking, the limma fit |
| `R/methylome_ewas.R`, `R/methylome_ewas_inputs.R` | EWAS panels and settings |
| `R/methylome_regions.R` | DMRs and gene-set testing |
| `inst/masking/cross_reactive_probes.csv` | Chen 2013, McCartney 2016, Zhou non-SNP mask |
| `test/app.R` | Standalone harness for the screens |
| `test/build-test-pgx.Rscript` | Builds the demo dataset from GSE43976 |
| `test/test-arrays.Rscript` | 450K and EPIC verification |
| `test/test-model.Rscript` | Contrast and continuous fits, EPICv2 refusal, region helper |

## Running it

```bash
make methylome                                  # stages the demo dataset if it can find it
make methylome pgx=/path/to/dataset.pgx         # stage a specific dataset
make methylome port=3939
```

Then: **Library** → load a methylomics dataset. Loading one lands directly on
this board's Sample ledger; the six screens sit under **Methylome** in the
Dashboard menu. Which menu groups appear for methylation data is set by
`MODULES_METHYLOMICS` in `etc/OPTIONS`.

Optional dependencies, all degrading gracefully when absent: `methylclock`
(10 clocks and deconvolution — falls back to 5 clocks and no composition),
`dmrff` (regions), `missMethyl` (gene sets), `bacon` (inflation diagnostics).
None are currently declared in the Dockerfile.
