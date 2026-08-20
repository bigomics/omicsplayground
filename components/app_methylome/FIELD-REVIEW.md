# Methylome Profiler — field review (v4)

What published methylation papers actually do and plot, what this app can currently
reproduce, and what is missing.

## Where we sit today (2026-08-19)

Methylomics is now its **own app**, not a board: loading a methylation dataset
bypasses the Dashboard entirely and lands on the six Methylome screens, and
Dashboard/Studio/Copilot are disabled for the datatype. Three of the seven
items PLAN.md listed as "missing, real work" were built in the same round —
**SVA, the multi-level F-test EWAS and the EWAS Catalog lookup** — so the
scientific gap list is materially shorter than in v3. Epigenetic clocks are
fitted once at dataset creation rather than per session, the methylation code
now lives in a separate package (`playbase.epigenetics`), and a bug that had
silently produced **no clock ages at all on every platform-uploaded dataset**
is fixed.

Against that: **`dmrff`, `missMethyl` and `bacon` are still absent from the
Docker build**, so DMR calling, gene-set testing and the inflation diagnostic
are all dead in production — the "Regions & pathways" sub-tab has nothing to
show at all. That is the single largest gap between this document and what a
deployed user actually sees, and it is a one-line-per-package fix. Everything
below is written against the code as of 2026-08-19; where a capability exists
in the repo but not in the container, it says so.

**Version history.** v1 was compiled 2026-08-17 from four parallel surveys:
population/cohort EWAS, cancer & clinical, clocks & cell composition, and a
code-level audit of the app as it stood then. **v2, same day**, updated the
code-audit and gaps/strengths sections after a round of changes: the app became
a platform board, a volcano plot and a top-hits heatmap were added, seven
defects were fixed, and drugs/connectivity were switched off by default for
methylation data. **v3, also 2026-08-17**, checked in after a further round:
four of v2's seven cheap, high-visibility gaps were built, one long-open
correctness risk was fixed, and a root-cause bug that had silently disabled sex
prediction on every playbase-built dataset was fixed. **v4, 2026-08-19**,
follows the standalone-app conversion, the `playbase.epigenetics` split, three
newly built gap items, two correctness fixes on the upload path, and the
verification that a script and the platform now produce the same numbers. The
literature survey and the conventions section are unchanged from v1.

**On the evidence.** Frequency labels ("near-universal", "common") are the
surveys' judgement from the papers they inventoried, not systematic counts —
each survey said so explicitly. The cancer survey could not access three key
full texts and flagged them. Treat the direction as solid and the precise
percentages as indicative.

---

## 1. Headline

The app is **credible for population EWAS and epigenetic clocks, and essentially
absent for cancer/clinical methylation**. Cancer is where most methylation-array
papers sit, and its figure set barely overlaps with ours. That verdict hasn't
moved since v2 — but the EWAS half got materially stronger this round rather
than merely deeper. Three things the field treats as standard equipment landed
together (commit `75da3d588`):

- **SVA / latent-factor adjustment.** The confounding control available when no
  deconvolution reference exists for the tissue — tumour, buccal, adipose,
  placenta — with the outcome protected in `mod` and the surrogate variables
  named in the printed formula line (`methylome_model.R:386-423`, checkbox at
  `methylome_ewas_inputs.R:52`).
- **Multi-level (moderated F-test) EWAS.** A three-arm categorical phenotype is
  tested properly instead of being redirected to the continuous path, and the
  panels that cannot express an unsigned result refuse rather than approximate:
  the volcano and the DMR caller both say why, bacon is skipped, the hits table
  turns Δβ into an unsigned beta range and adds an `F` column.
- **EWAS Catalog lookup.** Top hits marked known-vs-novel against
  ewascatalog.org, degrading to a message when the container has no egress.

Closest to competent: **clocks**, and now cheaper — the ten-clock fit that used
to cost ~40 s behind a button on every session is computed once at dataset
creation and stored on `pgx$meth$clocks`. Furthest: **cancer
subtyping/prognostic**, unchanged, and in one respect *worse* than v3 claimed:
v3 filed consensus clustering and t-SNE/UMAP as "routing, not backlog" because
`board.clustering` could be pointed at a methylation pgx. Methylomics now loads
no Dashboard board at all, and the compute no longer produces the embeddings or
the `gsetX` matrix that board requires. That routing answer is dead; see §5.

**Since v3**, this is a standalone app (`components/app_methylome`), not a
Dashboard board. The trade is deliberate and it cuts both ways: the "one click
from an unadjusted diffexpr view" hazard v3 flagged is gone, and so is every
free reuse of a platform board.

---

## 2. Coverage by paper archetype

| Archetype | We could produce | We would be missing |
|---|---|---|
| **Single-cohort EWAS** (smoking, BMI, exposure) | Genome-wide overview, QQ + λ + bacon*, top-hits table with catalog annotation, DMRs*, bias-corrected gene sets*, island-context enrichment, per-CpG stripcharts, region detail, volcano, top-hits heatmap, covariate / cell / **SVA** adjustment, **multi-level F-test** | PC-vs-covariate association matrix, replication scatter, forest plots |
| **Consortium meta-analysis** | The per-cohort half only | Forest plots, I²/Cochran's Q, leave-one-cohort-out — i.e. everything that defines the genre. Needs multiple cohorts |
| **Clocks / biological ageing** | DNAm vs chronological age with *r* and MAE, clock-vs-clock correlation, acceleration by group with test, coverage table (probe % and weight % missing), composition bars and by-group | DunedinPACE, PC-clocks, methylation exposure scores, survival, longitudinal trajectories, ICC/replicate reliability |
| **Cancer subtyping / prognostic** | Almost nothing | Annotated heatmap (the centrepiece), consensus clustering, t-SNE/UMAP, Kaplan-Meier, ROC, risk plots, tumour purity, copy number |

\* **In the repo, not in the container.** bacon, DMRs (`dmrff`) and gene-set
testing (`missMethyl::gometh`) are written, guarded and tested, and all three
are absent from `docker/Dockerfile.update` — so on a deployed instance the
EWAS row loses its region calling, its gene-set testing and its inflation
diagnostic. See §4, "Blocked on the deployment image".

Consensus clustering and t-SNE/UMAP in the cancer row were filed as a routing
question in v3. They are back to being build gaps — the standalone app loads no
Dashboard board, and the methylomics compute no longer produces embeddings
(§5). The sex-concordance check and the optional sex-chromosome mask are
QC/design controls that cut across every archetype rather than belonging to one
row's figure list; they show up in §3 and §4, not here.

---

## 3. Strengths

Things this app does that are correct and not universal in the field:

- **The model is refitted in-app** with user-chosen covariates, cell-composition
  adjustment and — since 2026-08-19 — surrogate variables. Most tools display the
  unadjusted result stored at upload. In bulk
  tissue that is usually a test of cell composition rather than of the phenotype
  (Jaffe & Irizarry 2014, *Genome Biol* 15:R31).
- **Latent-factor adjustment is available and honest about itself.** SVA
  estimates surrogate variables from the M-value matrix with the outcome
  protected in `mod`, caps at ten, degrades to the unadjusted fit with a note
  when `num.sv` finds none, and appends `SV1 + SV2 …` to the printed formula so
  which model you are looking at is never a guess (`methylome_model.R:386-423`).
  The info text says what the literature says — report both models — and warns
  that Δβ stays the raw group difference and is not adjusted by anything in the
  design.
- **A three-or-more-level phenotype gets an F-test, and the panels that cannot
  express it refuse.** The moderated F is the right test for "does any group
  differ"; it has no direction, so the volcano refuses with a pointer to the
  stripcharts, region calling refuses because dmrff needs a signed estimate and
  its standard error, bacon is skipped because it models a signed z, the
  Manhattan drops its hypo/hyper colouring, and the Δβ column becomes an
  explicitly unsigned beta range beside an `F` statistic. Nothing is
  approximated into a number that would read as a direction
  (`methylome_ewas.R:115-117,238-259,415-416`, `methylome_regions.R:31-34`,
  `methylome_server.R:135-137`).
- **Hits can be marked known-vs-novel, and the wording refuses to overclaim.**
  `playbase.epigenetics::ewas_catalog_lookup()` queries ewascatalog.org for the
  top 50 hits and returns `NULL` the moment any request fails rather than a
  partial answer — half an answer would read as "these CpGs are unreported",
  the one wrong conclusion the lookup exists to prevent. The column header says
  absence from the catalog is absence from the catalog, not novel biology.
- **Δβ is computed on the beta scale**, never taken from the stored M-value logFC.
  Test on M, report on β is the settled convention; carrying an M logFC into a
  results table is a common error.
- **Enrichment background is the probes actually tested**, not the whole array.
- **bacon reports, it does not correct.** Applying it to an unadjusted model would
  launder confounding into well-calibrated-looking p-values.
- **Clocks below a coverage floor are withheld**, not estimated from partial probes.
- **MAE sits beside *r* on the age-correlation panel, not just *r*.** Median
  absolute error, not mean, so one failed sample doesn't set the headline
  number. On the test cohort (GSE43976), Horvath/Hannum/Levine all correlate
  at *r* ≈ 0.85–0.86, but MAE is 4.2 / 2.2 / 7.1 years — Levine's real offset
  is invisible in *r* alone (`methylome_age.R`, `methylome_plot_agecor_server`).
- **Clock coverage is reported by weight, not just by probe count.**
  `clock_weights()` pulls \|coefficient\| per CpG straight from
  `methylclockData`'s own getters at call time, so nothing here is a constant
  transcribed by hand; the intercept is dropped from the denominator, and BNN
  is excluded outright because a neural net has no weighted sum to divide up.
  The two numbers diverge: on this dataset, dropping the single heaviest
  Horvath CpG costs 0.3% of its probes and over 3% of its weight. Moved to
  `playbase.epigenetics` with the fitting, and reached through thin wrappers
  (`methylome_age.R:168`).
- **Clocks are fitted once at dataset creation, and the app checks the stored
  fit before trusting it.** `pgx.computePGX` calls
  `playbase.epigenetics::compute_clocks()` immediately after the beta
  conversion and stores ages, per-clock probe and weight coverage, the
  producing package and its version. `mp_stored_clocks()` rejects the slot when
  the ages no longer line up with the loaded samples, or when the installed
  package version has moved — a methylclock upgrade can move a coefficient set,
  and a silently stale age is worse than a slow one — and falls back to the
  same live fit as before (`methylome_age.R:122-159`). The one mismatch *not*
  treated as stale is the package having gone away entirely: a stored
  methylclock fit beats refitting on the wateRmelon fallback.
- **The fit is stored with no coverage floor and no clock selection.** Both are
  display choices, applied to the coverage table rather than baked into a
  stored result; verified bit-identical at `min.perc = 0` and `0.8` for every
  clock defined at both.
- **Full published mask set** — Chen 2013, McCartney 2016, Zhou 2017 (53,498 unique
  cross-reactive probes) plus SNP probes (52,116 on 450K, 90,084 on EPIC).
- **A sex-chromosome mask is offered, and deliberately off by default.** SNP
  and cross-reactive probes are unconditionally wrong wherever they sit; X/Y
  probes are only conditionally wrong — pooling them across sexes creates
  spurious hits in a mixed-sex cohort, but masking them out deletes the answer
  in a single-sex cohort or a study of sex itself. 11,648 probes on 450K,
  19,627 on EPIC (`mp_sexchr_probes`, `methylome_model.R:182-201`; checkbox in
  `methylome_ewas_inputs.R`).
- **Array type inferred from probe IDs**, not from metadata that a pgx does not carry.
- **EPIC v2 detected and refused** rather than silently mislabelled as 450K.
- **Predicted sex is checked against the sample sheet, not just displayed.**
  A MISMATCH folds into the ledger's Verdict column — the cheapest sample-swap
  detector an array gives you. Numeric 1/2 sex coding is refused rather than
  guessed, because which digit means male is not standardised across cohorts.
  The logic moved to `playbase.epigenetics` (`sex_check`, `recorded_sex`,
  `has_sexchr`) so the upload QC screen and the ledger cannot drift apart; the
  app keeps only pgx-unpacking wrappers (`methylome_utils.R:44-61`,
  `methylome_ledger.R:23-41`).
- **The per-sample verdict declines to be cohort-relative without a cohort.**
  Flagging is median ± 3 MAD on imprint drift and bimodality, because absolute
  cut-offs do not transfer between platforms or tissues — but below a minimum
  cohort size, or when ties drive the MAD to exactly zero, only the sex
  mismatch (the one absolute check) can raise a CHECK. Previously a zero MAD
  flagged every sample that merely differed from the mode
  (`playbase.epigenetics::sample_qc`).
- **`gometh` rather than a naive hypergeometric** — corrects the probes-per-gene bias
  that otherwise returns the same probe-dense developmental sets whatever the biology.
- **Methylation gets its own app rather than a scoped Dashboard menu.**
  v3's strength here was `MODULES_METHYLOMICS`, which dropped Compare and
  SystemsBio but left Expression and GeneSets open — boards that read the
  stored unadjusted model. That list no longer exists (`etc/OPTIONS:51-53`).
  `MODULES_ACTIVE` is all-`FALSE` for methylomics, so no Dashboard board is
  ever inserted, Dashboard/Studio/Copilot are disabled rather than left
  clickable on shells that were never built for the datatype, and the load
  navigates straight to the Methylome app (`opg_server.R:130-146,251-258,
  386-395`). The routing decision lives in the one observer that already
  navigates after a load, not in a second observer racing it.
- **Script and platform produce the same numbers, and this is now checked.**
  The app's own compute path *is* the script path (`bin/pgxcreate_op.R` →
  `pgx.createPGX` + `pgx.computePGX`, driven by a saved `params.RData`). The
  last divergence was the methylation QC sample removal, which the app used to
  apply itself by pre-filtering counts/samples/contrasts; `drop_samples` now
  travels inside the `preprocess` settings list so `pgx.preprocess()` drops the
  columns *before* normalizing — which matters, because quantile normalization
  is not indifferent to the order (`upload_module_normalization.R:1645`,
  playbase `pgx-preprocess.R:83-86`). Verified against a real UI upload driven
  end to end: same three CSVs, same settings typed by hand into a script,
  identical `X` (max abs diff 0), identical effect sizes, p-values, q-values
  and clock ages (commit `49fff2178`).
- **The compute is scoped to what the app actually shows.** Gene sets,
  `meta.go`, WGCNA and the embeddings were being computed and stored for
  Dashboard screens that are never rendered for methylomics. On 9 samples ×
  17.8k probes the compute drops from ~4.5 min to ~96 s and the pgx from
  13.2 MB to 4.2 MB (commit `4bb9956b0`). `compute_testGenes` is deliberately
  kept: `pgx.checkObject()` requires `gx.meta` and a pgx failing it never
  reaches the Library, and it costs 0.15 s. The cost of this is real and is
  recorded in §5 — no `gsetX` and no `cluster$pos` means no reuse of
  `board.clustering` even in principle.
- **Drugs and connectivity are opt-out, not opt-in, for methylation data, and the
  code says why.** Both correlate the dataset's fold-change against L1000
  gene-*expression* reference signatures after collapsing CpGs to gene
  symbols. A methylation change is not an expression change — promoter
  methylation lowers expression, gene-body methylation often raises it — so
  the correlation compares two different quantities that merely share a symbol
  index. Off by default at upload for methylomics (`upload_module_computepgx.R`);
  still selectable for anyone who wants it.

---

## 4. Gaps

### Blocked on the deployment image — the highest-value open item

**`dmrff`, `missMethyl` and `bacon` are not installed by
`docker/Dockerfile.update`, and are not pulled in by playbase either** —
neither package appears anywhere in playbase's `DESCRIPTION` or in the update
image (checked 2026-08-19). All three are optional-with-a-guard in the app
(`mp_can_dmrff`, `mp_can_gometh`, `mp_bacon`), so nothing crashes; the panels
simply have nothing to show. The damage, per screen:

| Missing | What dies | Screen |
|---|---|---|
| `dmrff` | DMR calling and the region detail plot that reads its table | EWAS → Regions & pathways |
| `missMethyl` | `gometh` gene-set testing — the probes-per-gene bias correction §3 lists as a strength | EWAS → Regions & pathways |
| `bacon` | The bias/inflation diagnostic beside λ | EWAS → QQ & context |

So one of the app's six sidebar screens loses one of its three sub-tabs
entirely, and a second sub-tab loses its inflation diagnostic while keeping the
QQ plot and λ. (An earlier framing of "three of six screens" overstates it —
the three packages land on two sub-tabs of one screen.) All three are ordinary
Bioconductor/CRAN/GitHub installs, and the same round already added the
Illumina manifests and `methylclock`/`methylclockData` to the image
(commit `ea50aebf5`), so the pattern to copy is one line up.

For contrast, `sva` — which the new latent-factor adjustment needs — *is*
covered: it sits in playbase's `Imports`, so `library(playbase)` would fail
without it and it must be present in the base image. That is luck rather than
design; the three above have no such backstop.

### Closed since v3

Three of PLAN.md's seven "missing, real work" items were built on 2026-08-19,
all inside commit `75da3d588` (whose message advertises only the standalone-app
conversion — the SVA, F-test and catalog work is in the diff but not in the
body, which is worth knowing if you go looking): **SVA**, **multi-level F-test
EWAS**, and the **EWAS Catalog lookup**. All three are described in §3 and
tested headlessly in `test/test-model.Rscript`. One further item, **BMIQ**,
turned out never to have been missing — see below.

v2's seven cheap items had already closed four (MAE, the sex-concordance check,
the sex-chromosome mask, clock weight coverage); the other three are in the
dependency table immediately below.

### Blocked on a missing coefficient set

| Gap | Why it matters | Blocked on |
|---|---|---|
| **LUMP tumour purity** | Mean β over 44 immune-unmethylated probes ÷ 0.85 (Aran 2015). Purity confounds cancer clustering | No installed package carries the probe set — `EpiDISH` and `InfiniumPurify` are both absent |
| **Smoking / exposure scores** | Elliott 2014 coefficients are public; used as a *covariate* in modern EWAS, not just an outcome | `EpiSmokEr` and `meffonym` are both absent |
| **DunedinPACE + PC-clocks** | Public coefficients. Lussier 2024 recommends PC versions for cross-array robustness | `DunedinPACE` and `dnaMethyAge` are both absent; `methylclockData` — already a dependency — carries only the nine clocks already in the app, not these |

Checked directly, not assumed: none of `EpiDISH`, `InfiniumPurify`,
`EpiSmokEr`, `meffonym`, `DunedinPACE`, `dnaMethyAge` is installed. All three
were skipped rather than hardcoded from memory — an invented probe list
produces a confident, wrong number, and nothing downstream would flag it as
invented.

### Missing, real work, feasible from a beta matrix

- **Per-probe missingness filter** — a probe with three non-missing samples
  still enters the FDR count beside a fully observed one. Checked directly:
  `mp_fit_ewas()` drops incomplete *samples* (`complete.cases` on the model
  frame, `methylome_model.R:340-341`) but applies no completeness floor to
  probes. This is the smallest item on the list and the only one on it that is
  pure correctness; see PLAN.md §1.
- **Consensus clustering, t-SNE/UMAP** — the cancer centrepieces. v3 filed
  these as "routing, not backlog" via `board.clustering`. **That is no longer
  true and was already only two-thirds true.** Consensus clustering never
  existed anywhere in the platform (`board.clustering` offers a k-cut of one
  hclust/kmeans, no resampling, no consensus matrix), and now the embeddings
  are gone too: methylomics loads no Dashboard board and its compute skips
  `do.cluster` and the gene-set matrix, so `cluster$pos` and `gsetX` — both
  `req()`ed by that board — are simply absent from a methylomics pgx. Building
  either here is a build, not a route. See §5.
- **Survival (KM / Cox)** where the metadata carries follow-up — and that
  clause is the whole feature, not fine print: no dataset the app can currently
  load carries time-to-event columns, so the refusal branch is the only branch
  any user would see. Ruled a separate module, deliberately not built.
- **BMIQ probe-type correction — v3 was wrong to list this as missing.**
  `playbase.epigenetics::normalizeMethylation()` has `method = "BMIQ"` as its
  *default* and is called from `pgx.preprocess()` on the methylomics branch, so
  a dataset uploaded through the platform is already BMIQ-corrected at ingest —
  which is the upstream option PLAN.md §4 called "the correct home". What is
  genuinely missing is a per-analysis toggle, and it is not obvious anyone
  wants one. See the caveat in §6: until 2026-08-19 that call could silently
  return unnormalized data.
- **EWAS Catalog / Atlas lookup** — **built 2026-08-19**, see §3. The Atlas
  half stays unbuilt: NGDC has no stable public API worth depending on.

### Impossible from a beta matrix — not backlog items

- **Copy number.** β = M/(M+U+α) is a ratio; total hybridisation intensity — the
  copy-number signal — divides out and cannot be recovered. `conumee` explicitly
  requires both channels. This is impossible, not degraded.
- Detection p-values, bead counts, pOOBAH masking, control-probe QC.
- funnorm, Noob/ssNoob, dasen, dye-bias correction (all need channel intensities).
- Read-level analyses: PDR, epipolymorphism, methylation haplotypes, UXM
  deconvolution — these need several CpGs on one physical molecule.
- Meta-analysis, forest plots, I² — need multiple cohorts.
- Survival — needs time-to-event follow-up in the metadata.
- mQTL, allele-specific methylation, Mendelian randomisation — need genotypes.
- eQTM, starburst plots — need matched expression.
- GrimAge — coefficients are patent-encumbered; PCGrimAge is the safer public route.

---

## 5. App vs. dashboard: division of labour — settled by removal

**As of 2026-08-19 this question is closed, and closed in the strong
direction.** Methylomics loads no Dashboard board at all: `MODULES_ACTIVE` is
an all-`FALSE` array for the datatype, the removal loop tears down whatever a
previously loaded dataset left behind, and Dashboard/Studio/Copilot are
disabled (`opg_server.R:130-160,251-258`). `MODULES_METHYLOMICS` was deleted
rather than tightened — `etc/OPTIONS:51-53` now carries only a comment saying
there is no module list left to configure.

That resolves the problem v3 flagged: a methylomics user is no longer one click
from the unadjusted diffexpr view or the density-inflated enrichment this app's
EWAS tab exists to avoid. It also closes the escape hatch. The table below is
kept as the record of *why* each board was not worth routing to, and one row
has changed verdict:

**`board.clustering` is no longer reusable, on two independent counts.** It is
not loaded for methylomics, and the compute no longer produces what it
`req()`s: `do.cluster` is `FALSE` and `add.gmt` is `FALSE` for the datatype
(`upload_module_computepgx.R:1143-1148`), so a methylomics pgx now has neither
`cluster$pos` nor `gsetX`. Anyone who wants PCA/t-SNE/UMAP or a consensus
partition on methylation data has to build it inside this app. That is the
price of the standalone conversion, and it should be paid knowingly.

The original audit, for the record:

| Board | Verdict for methylation | Why |
|---|---|---|
| `board.clustering` (PCA/t-SNE/UMAP/heatmap) | **Reuse as-is** | Model-independent — it clusters whatever matrix it is handed. Unsupervised subtyping views for the cancer archetype are a routing question, not something to build here. |
| `board.expression`, `board.enrichment`, `board.pathway`, `board.signature`, `board.intersection`, `board.connectivity`, `board.drugconnectivity` | **Not appropriate as-is** | All read the stored, unadjusted M-value model computed at upload — exactly the problem this board's in-board refit exists to avoid (§3). `enrich`/`pathway` additionally collapse CpGs to gene symbols with no dedup, so hits inflate with probe density per gene; that is the precise bias `gometh` corrects, and these boards don't apply it. |
| `board.tcga` | **Not appropriate** | Would error outright on CpG IDs. |
| `board.correlation` | **Reuse, but it carries a defect** | Partial correlation is computed on M-values; plain correlation is computed on `2**beta` then `log2(1+x)`. The two are shown side by side with no note that they sit on different scales. |

v3 used this table to argue that `MODULES_METHYLOMICS` kept two unsafe menu
groups open. That argument is spent — the menu is gone. The verdicts still hold
for anyone tempted to re-route methylation through a Dashboard board later.

---

## 6. Code audit: fixed and open

The seven defects v1 found are fixed. Checked again at their current location:

| Fixed | Where |
|---|---|
| **bacon bias sign.** Now fed `fit$t`, the already-signed moderated t, instead of rebuilding a z-score as `qnorm(p/2)*sign(logFC)` — `qnorm(p/2)` is always negative, so that construction flipped every sign | `methylome_server.R` (`r_ewas`, `bacon = mp_bacon(fit$t)`) |
| **Chromosome selector for the beta boxplot.** `r.chromosome()` is now read; it was wired into the module signature but never called | `epigenomics_plot_boxplot_beta.R:165-172` |
| **X/Y chromosomes silently undrawable.** `mp_sort_chroms()` replaces `sort(as.numeric(...))`, which turned `"X"`/`"Y"` into `NA` | `methylome_server.R` |
| **Manhattan threshold line disagreeing with its own points.** `mp_ewas_hline()` now derives the line from `mp_ewas_sig()`, the same rule that colours the points and counts the hits, instead of the q cut-off alone | `methylome_ewas.R:47-67` |
| **Ledger's false coverage-check claim.** No longer says "first clock with complete probe coverage"; now says the ledger's DNAm age is "a single quick estimate with no coverage check" and points to the Epigenetic age screen for the real one | `methylome_ui.R` |
| **Composition bars auto-picking their grouping variable.** The by-group plot now takes an `r.pheno` reactive wired to the `comp_pheno` selectInput; the "first column with 2-6 levels" heuristic is gone | `methylome_deconv.R:133-145`, `methylome_server.R` |
| **`test/app.R` dead file.** Removed, along with `test/run.sh` | (deleted) |

One more speed fix in the same vein, not in v1's list: the ledger used to pay
for all five wateRmelon clocks to display one number, so `mp_clocks()` gained a
`method` argument and the ledger asked for Horvath alone. Overtaken since —
with the clocks precomputed, the ledger reads Horvath off the stored fit and
only falls back to a live single-clock fit when there is no stored one
(`methylome_ledger.R:23-41`).

Two more, fixed since v2:

| Fixed | Where |
|---|---|
| **The EWAS cell-adjustment checkbox could fail open.** Flagged as "still open" in v2: if `methylclock` is not installed, `mp_can_deconv()` is `FALSE`, `r_cells()` returns `NULL`, and ticking "Adjust for cell composition" silently fitted the *unadjusted* model with nothing in the model summary to say so. It now validates and refuses, with a message pointing at the Cell composition screen | `methylome_server.R:96-101` (`r_ewas` eventReactive), `methylome_deconv.R:19-26` (`mp_can_deconv`) |
| **Sex prediction was dead board-wide.** `mp_has_xy()` matched `^(chr)?[XY]$` against `pgx$genes$chr`, but that column is a *cytoband* (`"Xp22.2"`) in every playbase-built pgx, never a bare chromosome name. The match returned `FALSE` universally, so sex was never predicted and the ledger's "Sex predicted"/"Sex check" columns always read "not available" — found and fixed while implementing the sex check above. Now matches the band's leading chromosome and requires ≥100 X/Y probes before it will predict at all, so a handful of stray probes still correctly declines rather than guessing confidently | now `playbase.epigenetics::has_sexchr()`, wrapped at `methylome_utils.R:49-51` |

That cytoband trap is not a one-off: `mp_ewas_annotate()`, the function behind
the Manhattan plot and the results table's chromosome column, already has to
strip the same cytoband down to a bare chromosome for the same reason
(`methylome_ewas.R:29-31`, "`pgx$genes$chr` is a cytoband... the arm and band
are stripped as `board.epigenomics` does"). `mp_has_xy()` was a new function
that re-derived the same assumption independently and got it wrong the first
time. Worth treating as a class of bug this codebase keeps producing — any
new code that reads `pgx$genes$chr` expecting `"chrX"` rather than
`"Xp22.2"` — not as an isolated mistake.

### Fixed since v3 — two silent-wrong-answer bugs on the upload path

Both of these sat outside this app's tree, and both are the kind that produce a
confident number rather than an error.

| Fixed | Where |
|---|---|
| **"Append symbol to feature ID" was rewriting probe IDs, and it silently disabled the epigenetic clocks on every platform-uploaded dataset.** `convert.hugo` turns `cg00000029` into `cg00000029_RBL2`. Every methylation lookup — clock coefficient sets, the Illumina manifests, the masking lists — keys on the bare ID, so a dataset uploaded through the Shiny form came back with zero clock coverage and no ages at all. Now forced off for methylomics in both places, and corrected in `settings` so the record says what happened rather than what was asked for | app: `upload_module_computepgx.R:1155-1159`; playbase: `pgx-compute.R:397-400` |
| **`meth_type = NULL` — the functions' own default — was fatal, and in one caller the crash was swallowed.** `!meth_type %in% c(...)` is `logical(0)` when `meth_type` is `NULL`, so `if (c1 \| c2)` throws "argument is of length zero". The Shiny upload always passes a value, so only script callers hit it — and `pgx.preprocess()` wraps `normalizeMethylation()` in `try(silent = TRUE)`, so **a script that omitted `meth_type` got unnormalized data back with no error and no message.** Fixed in playbase `pgx-annot.R` and in both copies inside `playbase.epigenetics::normalizeMethylation()` | playbase `933a299d`, playbase-epigenetics `59a34cb` |

The first of these is the reason "the clocks work" and "the clocks work on a
dataset you uploaded yourself" were different claims until 2026-08-19. Reported
from the verification run (not recorded in a commit body, so take the exact
figures as indicative): 0 → 81 of 90 clock-age cells populated on the 9-sample
test dataset, Horvath coverage 0 → 0.992.

### Open

Small, and none of it currently produces a wrong answer a user would see — but
none of it is fixed either.

- **`pgx$genes$human_ortholog` carries some wrong values** from the gprofiler
  annotation backend — `TKT → DDR2` is the example that was actually observed.
  Not caused by, and not fixed by, the human short-circuit below.
- **`meth_type` is accepted by `pgx.createPGX` but never stored on the pgx.**
  It is passed to `getProbeAnnotation()` (`pgx-compute.R:463`) and then
  forgotten, so anything downstream keyed on it silently falls through to 450K.
  The app defends against this by re-inferring the array from probe IDs
  (`mp_array_type`, `methylome_model.R:117-140`), which is the right defence,
  but the underlying hole is still there for any other consumer.
- **A dead guard in `pgx-compute.R`'s `remove.xy.probes` branch.**
  `intersect(c("chr","map"), colnames(pgx$genes))[1]` is length 1 even when
  nothing matches — it is `NA` — so `if (length(kk) > 0)` always fires, and if
  both columns are absent the next line indexes by `NA` and errors instead of
  degrading. Practically unreachable today because playbase-built pgx objects
  always carry `chr` (`pgx-compute.R:534-545`).
- **The next round of the epigenetics split is in this app, not in playbase.**
  `methylome_model.R` (504 lines) holds `MP_ANNO_PKG`, the manifest
  attach-and-read, array inference and probe masking — knowledge
  `playbase.epigenetics` already owns via `annotate_methylomics()`. Two
  independent derivations of "which manifest for which array" is exactly the
  shape of the cytoband bug above.
- **EPIC v2 is detected and refused, not analysed.** Supporting it needs an
  annotation package that is neither installed nor in the configured
  repositories, plus logic to collapse replicate probes.
- **`getHumanOrtholog` short-circuits for human — a speed fix, not a
  correctness one.** Returning the uppercased symbol for hsapiens → hsapiens
  takes the largest single stage of `createPGX` from ~23 s to 0.006 s on 23k
  features. It first looked like a correctness fix too; measured end to end it
  changes 179 of 17,793 orthologs, 123 of them casing only. Worth having, worth
  not overselling (playbase `1cf67e27`).

### Unverified

- **Volcano and top-hits heatmap have still never been visually confirmed to
  render.** Both are wired end to end — settings inputs
  (`methylome_ewas_inputs.R`) to UI (`methylome_ui.R:293,305`) to server
  (`methylome_server.R:352,354`) to the plotting functions
  (`methylome_ewas.R:407-538`). The app *has* since been driven in a real
  browser, but for the upload-parity check (commit `49fff2178`), which ended at
  the numbers; no record says anyone opened the EWAS tab and looked at these
  two panels.
- **The three new EWAS controls are headless-tested, not eye-tested.** The SVA
  checkbox, the categorical outcome optgroup and the EWAS Catalog button all
  have assertions in `test/test-model.Rscript` (SVA on both the two-group and
  the anova branch, the F-test path against a fabricated tertile factor, the
  catalog parser against a canned JSON fixture), but the settings-panel
  placement, the volcano's refusal message on an F-test, the `-` in the
  Direction column and the `Reported` column render are all browser-only and
  unchecked.
- **Nobody has confirmed the deployed container can reach ewascatalog.org.**
  The lookup degrades to a message on any failure and cannot hang
  (`curl` timeout 5 s, stop at the first failure), so the risk is a dead
  button, not a wedged app — but "behind ShinyProxy the app may or may not have
  outbound network" is still an open question, and no `requireNamespace` can
  answer it.
- **Whether the ideogram errors on a pgx without `$dma` — still not confirmed.**
  The code reads `dma = pgx$dma` and its only use is guarded by
  `if (!is.null(dma))` (`epigenomics_plot_methylIdeogram.R:109,141`), which
  reads as safe for a plain list or environment where `$` on a missing name
  returns `NULL` rather than erroring. But that is inference from a read, not
  a run: the test PGX builder always sets `dma`
  (`test/build-test-pgx.Rscript:49`), so no test exercises the missing case,
  and this audit did not execute the app against one.
- **The `DT::formatStyle("Sex check", ...)` styling is not browser-verified.**
  The table renders a MISMATCH row in a distinct colour and leaves "ok" and
  "no record" alike (`methylome_ledger.R:61-63`); that this actually shows up
  correctly in a live DataTable, in both the compact and modal render paths,
  has not been checked by eye.
- **The MISMATCH path has never fired on real data.** The one dataset this
  audit can run against, `GSE43976-methyl-mini.pgx`, carries exactly 8
  sex-chromosome probes — below the 100-probe floor — so `mp_has_xy()`
  correctly declines and the sex check never gets far enough to disagree with
  anything. The MISMATCH logic is exercised by unit tests with a fabricated
  prediction (`test-model.Rscript`, "sex check" section), not by a live board
  on a dataset with real X/Y coverage.

---

## 7. Conventions worth knowing

**Effect-size thresholds differ by field, and getting this backwards is visible.**
In population EWAS a Δβ of 1–5% is a real, publishable finding; a hard effect-size
cut deletes the signal (Mansell 2019). In cancer, effects are enormous — a silenced
promoter goes β 0.05 → 0.80 — and \|Δβ\| ≥ 0.2 or 0.3 is used as a noise filter.
A cancer paper headlining Δβ = 0.03 will be assumed to have a purity or batch
problem. Our \|Δβ\| filter defaults to off, which is right for EWAS; it is unit-blind
for continuous outcomes, where Δβ is a per-unit slope.

**Significance.** 1×10⁻⁷ is the de facto 450K convention (Saffari 2018 puts the
permutation estimate at 2.4×10⁻⁷; Mansell 2019 at 9×10⁻⁸ for EPIC). FDR is used as
primary when the study is underpowered or the signal diffuse. Cancer papers tighten
BH to 0.01 rather than borrowing GWAS thresholds.

**Cell adjustment is not automatically correct.** Braun 2019 found imputed Houseman
counts *reduced* replicated CpGs by 35–61% for age while barely mattering for BMI
and smoking; Crimmins 2021 found algorithmic estimates moved 15–35 coefficients
where measured flow cytometry did not. Reviewers now want adjusted and unadjusted
side by side — an argument for making that a first-class comparison rather than a
checkbox.

**Acceleration definitions are routinely botched.** EAA is the residual of DNAm age
on chronological age. IEAA adds a specific seven-cell-type set and is defined only
for Horvath. EEAA is *not* "Hannum adjusted for cells" — it deliberately up-weights
age-related cell counts before residualising, and is defined only for Hannum. Note
that `methylclock`'s own IEAA/EEAA column labels do not match these published
definitions; this app does not surface those columns, and computes its own residual.
Our "intrinsic" option residualises on chronological age plus the selected
deconvolution reference for the selected clock — defensible, but not literally IEAA,
and the label should say what it does rather than borrow the term.

---

## 8. Sources

Population EWAS: Joehanes 2016 *Circ Cardiovasc Genet*; Joubert 2016 *AJHG*;
Küpers 2019 *Nat Commun*; Merid 2020 *Genome Med*; Wahl 2017 *Nature*;
Hillary 2023 *PLOS Med*; Jovanova 2018 *JAMA Psychiatry*; Gruzieva 2017 *EHP*;
Campagna 2021 *Clin Epigenetics*; Braun 2019; Saffari 2018; Mansell 2019;
van Iterson 2017 (bacon); Suderman (dmrff); Phipson 2016 (gometh).

Cancer/clinical: Capper 2018 *Nature*; Koelsche 2021 *Nat Commun*;
Noushmehr 2010 *Cancer Cell*; Hinoue 2012 *Genome Res*; Ceccarelli 2016 *Cell*;
Janssels 2023 *Clin Epigenetics*; Pinto 2025 *Mol Oncol*; Wang 2021
*Clin Epigenetics*; Zheng 2017 (InfiniumPurify); Aran 2015 (LUMP);
Bady 2012 (MGMT-STP27); Moss 2018 *Nat Commun*; Loyfer 2023 *Nature*;
Zhou 2017 *NAR* (masks).

Clocks/composition: Horvath 2013; Hannum 2013; Levine 2018 (PhenoAge);
Lu 2019/2022 (GrimAge); Belsky 2022 (DunedinPACE); Higgins-Chen 2022 (PC-clocks);
Crimmins 2021; Lussier 2024 *Clin Epigenetics*; Bell 2019 *Genome Biol*;
Houseman 2012; Salas 2018/2022; Koestler 2016 (IDOL); Zheng 2018 (CellDMC);
Jaffe & Irizarry 2014; Elliott 2014; Bollepalli 2019 (EpiSmokEr).

§5 (app vs. dashboard) draws on a separate, code-level audit of
`board.clustering`, `board.expression`, `board.enrichment`, `board.pathway`,
`board.signature`, `board.intersection`, `board.connectivity`,
`board.drugconnectivity`, `board.tcga` and `board.correlation`, not on the
literature surveys above.

§4's "blocked on the deployment image" and §6's open list were re-checked
against the code on 2026-08-19: `docker/Dockerfile.update`, playbase's
`DESCRIPTION`, and `git log master..feat/methylome-app` plus the two
`playbase` / `playbase-epigenetics` branches. Where a figure comes from a
commit body it is attributed; where it comes from a verification run that left
no commit record, it says so.
