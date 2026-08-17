# Methylome Profiler — field review (v3)

What published methylation papers actually do and plot, what this board can currently
reproduce, and what is missing. v1 was compiled 2026-08-17 from four parallel surveys:
population/cohort EWAS, cancer & clinical, clocks & cell composition, and a
code-level audit of the app as it stood then. **v2, same day**, updated the
code-audit and gaps/strengths sections after a round of changes: the app became
a platform board, a volcano plot and a top-hits heatmap were added, seven
defects were fixed, and drugs/connectivity were switched off by default for
methylation data. **v3, also 2026-08-17**, checks in after a further round: four
of v2's seven cheap, high-visibility gaps are now built, one long-open
correctness risk is fixed, and a root-cause bug that had silently disabled sex
prediction on every playbase-built dataset is fixed. The literature survey and
the conventions section are unchanged from v1.

**On the evidence.** Frequency labels ("near-universal", "common") are the
surveys' judgement from the papers they inventoried, not systematic counts —
each survey said so explicitly. The cancer survey could not access three key
full texts and flagged them. Treat the direction as solid and the precise
percentages as indicative.

---

## 1. Headline

The board is **credible for population EWAS and epigenetic clocks, and essentially
absent for cancer/clinical methylation**. Cancer is where most methylation-array
papers sit, and its figure set barely overlaps with ours. That verdict hasn't
moved since v2 — this round added depth, not breadth. The clocks tab now
reports what the ageing literature actually reports (MAE beside *r*, weight-based
coverage beside probe-count coverage), and the ledger's identity check —
predicting sex from the array and comparing it to the sample sheet — went from
silently dead to working, because of a bug that had nothing to do with the sex
check's own logic and everything to do with a chromosome column being a
cytoband string ("Xp22.2"), not a chromosome name.

One correctness risk from v2's "still open" list is also closed: the EWAS
model can no longer silently drop a requested cell-composition adjustment. It
now refuses, with a message, rather than fitting an unadjusted model and
saying nothing about it.

Closest to competent: **clocks**, more so than in v2 — all three near-universal
panels are present (DNAm-vs-chronological scatter, clock-vs-clock correlation,
acceleration by group), and two of the archetype's three standing gaps (MAE,
weight-based coverage) are gone. Furthest: **cancer subtyping/prognostic**,
unchanged — this round's work went to EWAS and clocks, and the one
cancer-adjacent item on last round's cheap list, LUMP purity, stayed on the
shelf because its coefficient set does not exist in any installed package.

**Since v1**, this is a platform board (`components/board.methylome`), not a
standalone app — no more DEVMODE gate, no more launcher tile. That changes
which of the remaining gaps are actually missing versus a routing decision in
a board that already exists elsewhere in the platform; see §5.

---

## 2. Coverage by paper archetype

| Archetype | We could produce | We would be missing |
|---|---|---|
| **Single-cohort EWAS** (smoking, BMI, exposure) | Genome-wide overview, QQ + λ + bacon, top-hits table, DMRs, bias-corrected gene sets, island-context enrichment, per-CpG stripcharts, region detail, volcano, top-hits heatmap | PC-vs-covariate association matrix, replication scatter, forest plots |
| **Consortium meta-analysis** | The per-cohort half only | Forest plots, I²/Cochran's Q, leave-one-cohort-out — i.e. everything that defines the genre. Needs multiple cohorts |
| **Clocks / biological ageing** | DNAm vs chronological age with *r* and MAE, clock-vs-clock correlation, acceleration by group with test, coverage table (probe % and weight % missing), composition bars and by-group | DunedinPACE, PC-clocks, methylation exposure scores, survival, longitudinal trajectories, ICC/replicate reliability |
| **Cancer subtyping / prognostic** | Almost nothing | Annotated heatmap (the centrepiece), consensus clustering, t-SNE/UMAP, Kaplan-Meier, ROC, risk plots, tumour purity, copy number |

Consensus clustering and t-SNE/UMAP in the cancer row are a routing question
now, not a build gap — see §5. The sex-concordance check and the optional
sex-chromosome mask are QC/design controls that cut across every archetype
rather than belonging to one row's figure list; they show up in §3 and §4,
not here.

---

## 3. Strengths

Things this board does that are correct and not universal in the field:

- **The model is refitted in-board** with user-chosen covariates and cell-composition
  adjustment. Most tools display the unadjusted result stored at upload. In bulk
  tissue that is usually a test of cell composition rather than of the phenotype
  (Jaffe & Irizarry 2014, *Genome Biol* 15:R31).
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
  `mp_clock_weights()` pulls \|coefficient\| per CpG straight from
  `methylclockData`'s own getters at call time, so nothing here is a constant
  transcribed by hand; the intercept is dropped from the denominator, and BNN
  is excluded outright because a neural net has no weighted sum to divide up.
  The two numbers diverge: on this dataset, dropping the single heaviest
  Horvath CpG costs 0.3% of its probes and over 3% of its weight
  (`methylome_age.R:168-200`).
- **Full published mask set** — Chen 2013, McCartney 2016, Zhou 2017 (53,498 unique
  cross-reactive probes) plus SNP probes (52,116 on 450K, 90,084 on EPIC).
- **A sex-chromosome mask is offered, and deliberately off by default.** SNP
  and cross-reactive probes are unconditionally wrong wherever they sit; X/Y
  probes are only conditionally wrong — pooling them across sexes creates
  spurious hits in a mixed-sex cohort, but masking them out deletes the answer
  in a single-sex cohort or a study of sex itself. 11,648 probes on 450K,
  19,627 on EPIC (`mp_sexchr_probes`, `methylome_model.R:177-197`; checkbox in
  `methylome_ewas_inputs.R`).
- **Array type inferred from probe IDs**, not from metadata that a pgx does not carry.
- **EPIC v2 detected and refused** rather than silently mislabelled as 450K.
- **Predicted sex is checked against the sample sheet, not just displayed.**
  `mp_sex_check()` folds a MISMATCH into the ledger's Verdict column — the
  cheapest sample-swap detector an array gives you. Numeric 1/2 sex coding is
  refused rather than guessed, because which digit means male is not
  standardised across cohorts (`methylome_utils.R:118-146`,
  `methylome_ledger.R:32-51`).
- **`gometh` rather than a naive hypergeometric** — corrects the probes-per-gene bias
  that otherwise returns the same probe-dense developmental sets whatever the biology.
- **The menu is scoped for methylation, not just gated behind DEVMODE.**
  `MODULES_METHYLOMICS` (`etc/OPTIONS`) drops Compare and SystemsBio from the
  sidebar for methylomics datasets — cross-dataset comparison and drug/cell/PPI
  tooling built around gene-level signatures. Whether the two modules it keeps
  open, Expression and GeneSets, are actually safe is a separate question — §5.
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

### Closed since v2

v2 listed seven items here as cheap and high-visibility. Four are built (§3):
MAE, the sex-concordance check, the sex-chromosome mask, and clock weight
coverage. The other three turned out not to be cheap — each needs a published
CpG coefficient set this environment does not have — so they moved to the
table below rather than sit on a list that implied they were merely unwritten.

### Blocked on a missing dependency

| Gap | Why it matters | Blocked on |
|---|---|---|
| **LUMP tumour purity** | Mean β over 44 immune-unmethylated probes ÷ 0.85 (Aran 2015). Purity confounds cancer clustering | No installed package carries the probe set — `EpiDISH` and `InfiniumPurify` are both absent |
| **Smoking / exposure scores** | Elliott 2014 coefficients are public; used as a *covariate* in modern EWAS, not just an outcome | `EpiSmokEr` and `meffonym` are both absent |
| **DunedinPACE + PC-clocks** | Public coefficients. Lussier 2024 recommends PC versions for cross-array robustness | `DunedinPACE` and `dnaMethyAge` are both absent; `methylclockData` — already a dependency — carries only the nine clocks already in the board, not these |

Checked directly, not assumed: none of `EpiDISH`, `InfiniumPurify`,
`EpiSmokEr`, `meffonym`, `DunedinPACE`, `dnaMethyAge` is installed. All three
were skipped rather than hardcoded from memory — an invented probe list
produces a confident, wrong number, and nothing downstream would flag it as
invented.

### Missing, real work, feasible from a beta matrix

- **SVA / latent-factor adjustment** — the only route to confounding control in
  tissues with no deconvolution reference (tumour, buccal, adipose, placenta).
  Flagged independently by multiple surveys.
- **Multi-level / F-test EWAS** — a true three-arm categorical contrast is still
  refused; the model special-cases two-group contrasts or a continuous
  variable only, and the UI's own error message redirects a multi-level
  phenotype toward the continuous path, which is not the same test.
- **Consensus clustering, t-SNE/UMAP** — the cancer centrepieces, and still
  nothing this board builds. Not a build gap, though: `board.clustering` is
  model-independent and already usable on a methylation pgx (PCA/t-SNE/UMAP/
  heatmap) — routing, not backlog. See §5.
- **Survival (KM / Cox)** where the metadata carries follow-up.
- **BMIQ probe-type correction** — needs no IDATs, `wateRmelon` is already a
  dependency. The README's "Out of scope" section is still scoped to
  IDAT-requiring things and does not acknowledge this one.
- **Per-probe missingness filter** — a probe with three non-missing samples
  currently enters the FDR count beside a fully observed one.
- **EWAS Catalog / Atlas lookup** to mark hits novel vs previously reported.

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

## 5. Board vs. dashboard: division of labour

v1 raised "reuse vs build" as an open question for the cancer-archetype gaps.
Since then a separate audit of the platform's other boards answered it for
some of them, and surfaced a problem with how the methylomics menu is scoped
today.

| Board | Verdict for methylation | Why |
|---|---|---|
| `board.clustering` (PCA/t-SNE/UMAP/heatmap) | **Reuse as-is** | Model-independent — it clusters whatever matrix it is handed. Unsupervised subtyping views for the cancer archetype are a routing question, not something to build here. |
| `board.expression`, `board.enrichment`, `board.pathway`, `board.signature`, `board.intersection`, `board.connectivity`, `board.drugconnectivity` | **Not appropriate as-is** | All read the stored, unadjusted M-value model computed at upload — exactly the problem this board's in-board refit exists to avoid (§3). `enrich`/`pathway` additionally collapse CpGs to gene symbols with no dedup, so hits inflate with probe density per gene; that is the precise bias `gometh` corrects, and these boards don't apply it. |
| `board.tcga` | **Not appropriate** | Would error outright on CpG IDs. |
| `board.correlation` | **Reuse, but it carries a defect** | Partial correlation is computed on M-values; plain correlation is computed on `2**beta` then `log2(1+x)`. The two are shown side by side with no note that they sit on different scales. |

**This cuts against §3's "scoped for methylation" strength.**
`MODULES_METHYLOMICS` keeps Expression and GeneSets open for methylomics
datasets, and both menu groups route straight into boards from the "not
appropriate" row above: Expression bundles Differential expression, Correlation
analysis and Find biomarkers (`board.expression`, `board.correlation`);
GeneSets bundles Geneset Enrichment, Test geneset and Pathway analysis
(`board.enrichment`, `board.signature`, `board.pathway`) alongside Word cloud
(`components/app_opg/R/opg_ui.R`, `full_menu_tree`). A methylomics user is one
click away from the unadjusted diffexpr view and the density-inflated
enrichment this board's own EWAS tab exists to avoid. Keeping the menu open
was a deliberate call; the boards behind it are not yet methylation-safe, and
that's a gap the menu doesn't communicate.

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
for all five wateRmelon clocks to display one number. `mp_clocks()` gained a
`method` argument (default `"all"`, still what the Epigenetic age screen
needs); the ledger now asks for `"horvath"` alone.
(`methylome_utils.R:59-61`, `methylome_ledger.R:31`)

Two more, fixed since v2:

| Fixed | Where |
|---|---|
| **The EWAS cell-adjustment checkbox could fail open.** Flagged as "still open" in v2: if `methylclock` is not installed, `mp_can_deconv()` is `FALSE`, `r_cells()` returns `NULL`, and ticking "Adjust for cell composition" silently fitted the *unadjusted* model with nothing in the model summary to say so. It now validates and refuses, with a message pointing at the Cell composition screen | `methylome_server.R:96-101` (`r_ewas` eventReactive), `methylome_deconv.R:19-26` (`mp_can_deconv`) |
| **Sex prediction was dead board-wide.** `mp_has_xy()` matched `^(chr)?[XY]$` against `pgx$genes$chr`, but that column is a *cytoband* (`"Xp22.2"`) in every playbase-built pgx, never a bare chromosome name. The match returned `FALSE` universally, so sex was never predicted and the ledger's "Sex predicted"/"Sex check" columns always read "not available" — found and fixed while implementing the sex check above. Now matches the band's leading chromosome and requires ≥100 X/Y probes before it will predict at all, so a handful of stray probes still correctly declines rather than guessing confidently | `methylome_utils.R:100-105` |

That cytoband trap is not a one-off: `mp_ewas_annotate()`, the function behind
the Manhattan plot and the results table's chromosome column, already has to
strip the same cytoband down to a bare chromosome for the same reason
(`methylome_ewas.R:29-31`, "`pgx$genes$chr` is a cytoband... the arm and band
are stripped as `board.epigenomics` does"). `mp_has_xy()` was a new function
that re-derived the same assumption independently and got it wrong the first
time. Worth treating as a class of bug this codebase keeps producing — any
new code that reads `pgx$genes$chr` expecting `"chrX"` rather than
`"Xp22.2"` — not as an isolated mistake.

No load-bearing code-audit item is currently open. What's below is what
hasn't been checked, not what is known wrong.

### Unverified

- **Volcano and top-hits heatmap have never been visually confirmed to render.**
  Both are wired end to end — settings inputs (`methylome_ewas_inputs.R`) to UI
  (`methylome_ui.R:289,308`) to server (`methylome_server.R:316,318`) to the
  plotting functions (`methylome_ewas.R:373-499`) — but this audit is a code
  read, not a running session, and nobody has opened the board and looked at
  the pixels.
- **Whether the ideogram errors on a pgx without `$dma` — still not confirmed.**
  The code reads `dma = pgx$dma` and its only use is guarded by
  `if (!is.null(dma))` (`epigenomics_plot_methylIdeogram.R:109,141`), which
  reads as safe for a plain list or environment where `$` on a missing name
  returns `NULL` rather than erroring. But that is inference from a read, not
  a run: the test PGX builder always sets `dma`
  (`test/build-test-pgx.Rscript:49`), so no test exercises the missing case,
  and this audit did not execute the board against one.
- **The `DT::formatStyle("Sex check", ...)` styling is not browser-verified.**
  The table renders a MISMATCH row in a distinct colour and leaves "ok" and
  "no record" alike (`methylome_ledger.R:75-76`); that this actually shows up
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
definitions; this board does not surface those columns, and computes its own residual.
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

§5 (board vs. dashboard) draws on a separate, code-level audit of
`board.clustering`, `board.expression`, `board.enrichment`, `board.pathway`,
`board.signature`, `board.intersection`, `board.connectivity`,
`board.drugconnectivity`, `board.tcga` and `board.correlation`, not on the
literature surveys above.
