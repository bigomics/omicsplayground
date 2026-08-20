# Implementation plan — the seven "missing, real work" items

Plan for the seven items in FIELD-REVIEW.md §4 "Missing, real work, feasible from a
beta matrix". Compiled 2026-08-17 against the code as it stood on
`feat/methylome-app`. Every dependency claim below was checked with
`requireNamespace()` in this environment, not assumed; every consumer claim was
traced to a line, not remembered. Where §4 turned out to be wrong or optimistic,
this document says so at the item rather than hiding it in a caveats section.

## Status, 2026-08-19

Three items are built, three were deliberately not built, one is open. Each
item below opens with its own status line; this is the summary.

| # | Item | Status | Settled by |
|---|---|---|---|
| 1 | Per-probe missingness filter | **Open** — deliberately not implemented so far | Verified absent: `mp_fit_ewas()` filters incomplete samples, never probes |
| 2 | Multi-level (F-test) EWAS | **Done**, essentially as planned | `75da3d588`, `methylome_model.R:326-484` + the refusals listed at the item |
| 3 | SVA / latent-factor adjustment | **Done**, as planned | `75da3d588`, `methylome_model.R:380-423` |
| 4 | BMIQ probe-type correction | **Overtaken — it was already there** | `playbase.epigenetics::normalizeMethylation()` defaults to BMIQ and runs at ingest |
| 5 | Consensus clustering + t-SNE/UMAP | **Consensus deliberately skipped; the t-SNE/UMAP routing answer is now dead** | Standalone-app conversion + `do.cluster = FALSE` for methylomics |
| 6 | Survival (KM / Cox) | **Deliberately skipped** — ruled a separate module | The deferral argument at the item held |
| 7 | EWAS Catalog lookup | **Done** | `75da3d588` + `playbase.epigenetics::ewas_catalog_lookup()` |

Two things changed underneath the whole plan and are worth reading before any
item: the Methylome board became a **standalone app** that loads no Dashboard
board at all, and the methylation primitives moved out of playbase into
**`playbase.epigenetics`**. Where an item's instructions have been overtaken by
either, the item says so instead of leaving directions that would mislead.

### The dependency table has a hole the original run could not see

Verified package status, one run, this machine, 2026-08-17 — and re-checked
2026-08-19, still accurate **for this machine**:

| Present | Absent |
|---|---|
| `sva`, `missMethyl`, `dmrff`, `wateRmelon` (with `BMIQ`), `RPMM`, `survival`, `Rtsne`, `uwot`, `umap`, `bacon`, `limma`, `minfi`, `methylclock(Data)`, `httr`, `jsonlite`, `curl`, `impute` | `ConsensusClusterPlus`, `M3C`, `ENmix`, `survminer` |

**`requireNamespace()` on a dev machine was the wrong question.** `dmrff`,
`missMethyl` and `bacon` are installed here and are installed nowhere in the
build: not in `docker/Dockerfile.update`, not in playbase's `DESCRIPTION`. So
three features this plan treated as working already are dead in production, and
that outranks every item below. `sva` is fine only by accident — it sits in
playbase's `Imports`, so the base image must carry it. See FIELD-REVIEW.md §4,
"Blocked on the deployment image".

---

## 1. Per-probe missingness filter

> **Status 2026-08-19: OPEN.** Deliberately not implemented. Verified absent —
> `mp_fit_ewas()` drops incomplete *samples* via `complete.cases()` on the model
> frame (`methylome_model.R:340-341`) and applies no completeness floor to
> probes at all. Nothing since 2026-08-17 has touched this, and nothing has
> weakened the argument for it: it is still the smallest diff on the list and
> still the only item that is pure correctness rather than a new capability.
> The plan below is current except for line numbers — the insertion point moved
> from line 330 to `methylome_model.R:376-378`, where `B` is subset and handed
> to `betaToM()`, and the masking block just above it is now
> `methylome_model.R:370-373`.

**What it is.** Not a package — a row filter inside the model fit. limma tolerates
NAs per probe, so today a probe observed in 3 of 95 samples is tested on 3 points
and its p-value enters the same BH correction as a fully observed probe's. The
standard remedy is a completeness floor: drop probes observed in fewer than a
fixed fraction of the fitted samples (90–95% is the field's usual line), and for a
group contrast additionally require a minimum per group, because 90% overall can
still mean one group contributes two samples.

**Dependency.** None. Base R (`rowSums(is.na(B))`).

**Where.** `mp_fit_ewas()` only — `methylome_model.R`, immediately after `B` is
subset to the fitted samples (line 330, before `betaToM`). The masking block just
above (lines 323–326) is the template: compute the dropped set, count it, report
the count in the returned list beside `masked`, and let the model summary
(`output$ewas_model`, `methylome_server.R:115-130`) print it. No new UI control:
a fixed floor (probe kept when ≥90% of fitted samples are observed and, on the
contrast path, ≥3 per group) reported in the summary line beats a numeric input
nobody will know how to set. The dbeta computation (lines 339, 356) and the DMR
caller inherit the filter for free because both read the fitted table.

**What it breaks.** Nothing. Every consumer reads `fit$table` or `fit$t`, which
just get shorter. On the test dataset the filter is a no-op — `GSE43976-methyl-mini`
carries zero NAs — which is exactly why it must be tested with injected NAs, not
with the dataset as shipped.

**Effort and risk.** Half a day including tests. The one design question is
whether the floor applies before or after masking; after, so the counts in the
summary don't double-report. Nothing makes this not worth doing.

**Testing.** In `test-model.Rscript`: copy the pgx, set 20% of one probe's values
and 95% of another's to NA, refit; assert the mostly-missing probe is absent from
`fit$table`, the 20% probe is present, the reported drop count is right, and a
fully-NA-free refit returns the same rows as today. Fully headless; nothing needs
a browser beyond the summary line's wording.

---

## 2. Multi-level (F-test) EWAS

> **Status 2026-08-19: DONE**, and built essentially as specified below —
> commit `75da3d588`, whose message advertises only the standalone-app
> conversion. Read the diff, not the body, if you go looking.
>
> What landed, checked against the plan: the `nlevels == 2` validate is gone
> and `anova_fit <- !continuous && nlevels(dat$.group) >= 3` branches the fit
> (`methylome_model.R:353-355`); `mp_is_categorical()` and a categorical
> optgroup make a 3-level phenotype reachable at all; `topTable` over the group
> coefficients yields the `F` column, surfaced as `tab$Fstat`
> (`methylome_model.R:454-483`); `dbeta` becomes `max group mean − min group
> mean` on the beta scale (`:461`); bacon is skipped on the anova branch
> (`methylome_server.R:135-137`); the Manhattan drops direction colouring for a
> single accent (`methylome_ewas.R:115-117`); the volcano refuses with the
> "no single direction" message and points at the stripcharts (`:415-416`); the
> hits table prints `-` for Direction and relabels the column "Beta range"
> (`:238-259`); the DMR caller refuses (`methylome_regions.R:31-34`).
> `test/test-model.Rscript` covers the fit, the 2-level regression, the
> collapse-to-2-levels regression and SVA on the anova branch.
>
> The one thing worth re-reading below is the consumer table: it is still the
> map of which panel means what under an F-test, and it was accurate.

**What it is.** The moderated F-test limma already contains: fit `~ 0 + .group`
with k ≥ 3 levels plus covariates, then `topTable(fit2, coef = <the k group
columns>, ...)` — with multiple coefficients topTable returns an `F` column and
its p-value, testing "any difference among the groups" per probe. No new package;
this is the third branch `mp_fit_ewas()` refuses to grow (the validate at
`methylome_model.R:268-269`, whose message redirects a categorical phenotype to
the continuous path — the wrong test, as §4 says).

**Dependency.** `limma`, present.

**Where.** Three places, all extending existing code rather than sitting beside it:

- `mp_fit_ewas()` (`methylome_model.R:245-386`): relax the `nlevels == 2`
  validate to `>= 2` and branch. The k = 2 path is untouched. The k ≥ 3 path
  fits `~ 0 + .group + covars`, moderated F via `topTable` over the k group
  coefficients, and returns `anova = TRUE` in the list.
- The outcome picker. Today it offers only pgx contrast columns and continuous
  variables (`methylome_server.R:27-53`); playbase-built contrasts are two-group
  by construction, so a 3-level phenotype is currently *unreachable*, not just
  refused. Add a third optgroup "Categorical variables" from a
  `mp_categorical_vars()` sibling of `mp_continuous_vars()` (3–8 non-NA levels,
  ≥3 samples each), and teach the routing (`mp_is_continuous` gets a categorical
  counterpart) that such a selection takes the group path via
  `factor(pgx$samples[[sel]])` instead of `mp_contrast_groups()`.
- `r_ewas` (`methylome_server.R:80-112`): skip `mp_bacon()` when `meta$anova` —
  bacon models signed z/t statistics and an F is neither; feeding it one would
  produce a confident, wrong bias number. The QQ legend already handles a NULL
  bacon (`methylome_ewas.R:159-164`).

**What it breaks — the traced consumer list.** This is most of the work, and the
news is better than the framing suggests. Because the continuous path already
returns 3-level tertile `strata`, every per-sample panel was generalised to k
levels months ago and needs *nothing*:

| Consumer | Reads | k ≥ 3 behaviour | Action |
|---|---|---|---|
| Stripcharts (`methylome_ewas.R:291-364`) | `strata` | loops over `levels(grp)`, palette ramps to k | none |
| Hits heatmap (`:444-498`) | `strata` | group bar and legend already k-level | none |
| Region detail (`methylome_regions.R:151-234`) | `strata` | per-level mean lines, k-safe | none |
| QQ / lambda (`:148-171`) | `p`, `bacon` | fine on F-test p-values | guard NULL bacon (above) |
| Enrichment context, gometh | sig list | direction-free already | none |
| Manhattan (`:96-135`) | `dbeta` sign for hit colour (line 116) | **breaks**: no sign exists | one accent colour when `meta$anova`; `r.ewas()` is in scope |
| Volcano (`:382-426`) | signed `dbeta` x-axis | **breaks**: axis is a two-group concept | refuse with a message pointing at the stripcharts ("a k-group difference has no single direction"), or plot the beta *range* — refusal first, range plot only if someone asks |
| Hits table (`:229-274`) | `Delta beta`, `Direction` | **breaks** | `dbeta` := max group mean − min group mean on the beta scale ("beta range" — honest, unsigned, keeps the `min_dbeta` filter meaningful); `Direction` prints `-`, which `styleEqual` already leaves unstyled |
| `mp_ewas_sig` (`:47-54`) | `abs(dbeta)` | works on the range unchanged | none |
| DMR caller (`methylome_regions.R:25-63`) | `t` for `se = |est/t|` | **breaks**: no single estimate/se exists for an F-test, dmrff cannot run | `validate(need(!anova, ...))` — refuse, don't approximate. Return `t = NULL` from the anova branch so nothing downstream mistakes an F for a t |
| bacon (`methylome_ewas.R:76-83`) | `t` | **breaks** (wrong model) | skipped in `r_ewas`, see above |

So `dbeta` keeps its slot but changes meaning (range, not signed difference), the
volcano and DMR caller refuse with messages, the Manhattan drops direction
colouring, and everything else runs as-is. The table gains an `F` column in
place of `logFC_M` (which has no single value across k coefficients).

**Effort and risk.** 2–3 days. The risk is semantic, not technical: a user who
doesn't read the axis change may interpret the beta range as a signed effect.
Mitigation is the model summary line saying "F-test across k groups" and the
tooltip on `min_dbeta`. Not building it keeps the current state where the UI's
own error message sends users to a wrong test — this is the strongest candidate
on the list.

**Testing.** No natural 3-level factor exists in the test pgx (`cell_type` and
`smoking` are two-level), so the test fabricates one — `cut(age, 3)` written into
a copy's `$samples`. Assert: the fit succeeds with `anova = TRUE`; `groups` has 3
entries; `strata` has 3 levels named by sample; `table$p` is a valid p (the
fabricated age tertiles must rediscover ELOVL2-adjacent signal, same trick as the
continuous test at `test-model.Rscript:46-49`); `dbeta` ≥ 0 everywhere; `t` is
NULL; a covariate still enters the formula; and collapsing the factor to 2 levels
reproduces (up to contrast parameterisation) today's two-group p-values — the
regression that proves the old path is untouched. Browser-only: the Manhattan
recolour, the volcano refusal message, and the `-` Direction column.

---

## 3. SVA / latent-factor adjustment

> **Status 2026-08-19: DONE**, as specified — same commit `75da3d588`.
> `mp_fit_ewas()` gained `n_sv` (0 off, −1 estimate, positive to force), the
> insertion point is where the plan put it, `mod0` is the design minus the
> outcome columns with an intercept restored when the `~ 0 + .group`
> parameterisation leaves none, the SVs cbind onto the design and their names
> onto `used` so the printed formula reads `+ SV1 + SV2 …`, `num.sv` returning
> zero degrades to the unadjusted fit with a note rather than erroring, and the
> whole thing is wrapped so a degenerate design refuses instead of crashing the
> reactive (`methylome_model.R:380-423`, returned as `n_sv` at `:496`). The
> checkbox is at `methylome_ewas_inputs.R:52`, passed through `r_ewas` at
> `methylome_server.R:129`, and its info text carries the two caveats the plan
> asked for — report both models, and Δβ stays raw.
>
> One deviation, deliberate: the cap is applied inside the branch and there is
> no numeric input for `n_sv`. The checkbox is on/off, i.e. `n_sv = -1L` or `0L`.

**What it is.** `sva::sva(dat = M, mod, mod0, n.sv)` with `n.sv` from
`sva::num.sv(M, mod, method = "be")` — surrogate variables estimated from the
M-value matrix with the outcome protected in `mod`, then appended to the design
as covariates. The only confounding control available when no deconvolution
reference exists for the tissue, which is §4's whole point.

**Dependency.** `sva` — present, verified.

**Where.** Extends `mp_fit_ewas()`. A `n_sv` argument (default 0 = off; -1 =
estimate via num.sv, capped at ~10). The insertion point is after `M` exists and
after the design is built (`methylome_model.R:330-332`): `mod` is the current
design, `mod0` the same design minus the outcome columns (`.x`, or the `.group`
columns), SVs cbind onto `design` and their names onto `used` so the formula line
prints `+ SV1 + SV2 ...` — the summary honesty the board already trades on. On
the contrast path the same augmented design flows into `contrasts.fit`
(`makeContrasts` takes `levels = design`, so the new columns ride along); on the
continuous path the dbeta refit on `B` (line 339) uses the augmented design too,
so the reported slope is SV-adjusted, consistent with the p-value beside it.
UI: one checkbox "Adjust for latent factors (SVA)" in the hits settings group
(`methylome_ewas_inputs.R:33-66`), passed through `r_ewas`
(`methylome_server.R:102`) like `mask` is. It runs inside the existing Run EWAS
button — no new gating needed, but wrap the fit in `shiny::withProgress`.

**What it breaks.** Nothing structural — the return shape is unchanged and every
consumer is indifferent to what sat in the design. Two real risks instead:

- *Contrast-path dbeta stays raw.* `dbeta` for a two-group contrast is the plain
  difference of group means (line 356-357), unadjusted by covariates — already
  true today for cell fractions and covariates, so SVA doesn't make it worse,
  but the gap between an SV-adjusted p and a raw dbeta widens when the SVs soak
  up real structure. Worth one sentence in the info text, not a code change.
- *Runtime and stability.* On the 23k-probe test set sva is seconds; on a full
  450k × n matrix `num.sv` + `sva` is a couple of minutes of SVD and IRW
  iterations. Acceptable behind the existing button with progress, but the plan
  should not promise interactivity. sva also errors on some degenerate designs;
  wrap in tryCatch and refuse with the message rather than crashing the reactive.

**Effort and risk.** 1–2 days. Scientific risk is the standard one — with the
outcome protected, sva can still absorb signal correlated with the outcome
(e.g. cell composition when the phenotype drives composition). That is inherent
to the method, the literature's answer is "report both models", and the formula
line already makes which-model-this-is visible. Not a reason to skip it.

**Testing.** Assert: `n_sv = 2` on the smoking contrast fits, formula contains
"SV1", `ncol` bookkeeping holds (residual-df validate still passes), p-values
differ from the unadjusted fit, and `n_sv = 0` is bit-identical to today's
output (the regression). `num.sv` returning 0 must degrade to the unadjusted
fit with the summary saying so, not error. All headless; only the checkbox
placement needs eyes.

---

## 4. BMIQ probe-type correction

> **Status 2026-08-19: NOT BUILT, and the premise was wrong.** Deliberately not
> implemented in-app — and while checking why, it turned out **option 1 has
> been live all along**. `playbase.epigenetics::normalizeMethylation()` takes
> `method = "BMIQ"` as its *default*, runs `wateRmelon::BMIQ()` per sample with
> `design.v` from the manifest's `Type` column, and is called from
> `pgx.preprocess()` on the methylomics branch (playbase
> `pgx-preprocess.R:183`). So a dataset uploaded through the platform is
> already BMIQ-corrected at ingest, which is exactly the home this plan called
> correct. The in-app overlay of option 2 would have been throwaway plumbing —
> the risk the plan itself named — and it is right that nobody built it.
>
> Two live caveats, neither of which the plan could have known:
>
> - Until 2026-08-19 that call could **silently return unnormalized data**.
>   `normalizeMethylation()`'s own default is `meth_type = NULL`, which made
>   `!meth_type %in% c(...)` a `logical(0)` and the guarding `if` throw — and
>   `pgx.preprocess()` wraps the call in `try(silent = TRUE)`. A script that
>   omitted `meth_type` got its input matrix back with no error and no message.
>   Fixed in both copies (playbase-epigenetics `59a34cb`, playbase `933a299d`).
> - The runtime concern in this item is real and now lives at ingest rather
>   than behind a button, where a user has no progress bar and no way to skip it.
>
> What genuinely remains unbuilt is a *per-analysis* BMIQ toggle, and it is not
> obvious anyone wants one. FIELD-REVIEW §4 has been corrected accordingly.

**What it is.** Teschendorff 2013's beta-mixture quantile normalisation of type II
probes onto the type I distribution — `wateRmelon::BMIQ(beta.v, design.v)`, and
the generic *does* have a plain-vector method (`beta.v="ANY"`; verified with
`showMethods`), so no IDATs and no MethyLumiSet are needed. `design.v` (1 = type
I, 2 = type II) comes from the manifest's `Type` column via `mp_manifest()`. Its
EM engine is `RPMM::blc` — also installed.

**Dependency.** `wateRmelon` + `RPMM`, both present. §4 is right that this is
unblocked.

**Where — and here §4's "cheap" framing is wrong.** BMIQ is a *normalisation*:
it sits upstream of everything, not beside the model. Two honest options:

1. **Upstream, in playbase at pgx build time.** The correct home — every board
   then sees corrected betas — but it is a playbase change, out of this
   component's tree, and re-normalising can't be toggled per analysis.
2. **In-board, as an explicit corrected-beta overlay.** A button ("Apply BMIQ")
   computes the corrected matrix once per dataset into a `reactiveVal`; a wrapper
   reactive `PGX_adj` serves the pgx with `$X` swapped, and *every* module wired
   in `methylome_server.R` takes `PGX_adj` instead of `PGX`. The landscape
   modules that receive the raw `pgx` argument (lines 276-285) need the same
   swap or an explicit note that they show uncorrected betas.

Option 2 is what this plan can actually deliver, but the plumbing touches every
server wiring line, and the cache must invalidate on `r_dataset()` like the
clocks do.

**What it breaks.** By design, everything downstream changes numerically: beta
densities, clock estimates (Horvath's own advice is BMIQ first, so this is a
*fix* for clocks), deconvolution, dbeta magnitudes, and — mildly — the EWAS
itself. The per-probe test compares one probe across samples, all the same type,
so BMIQ's within-sample redistribution moves it only second-order; the real value
is cross-probe-type comparability (densities, heatmap absolute scale, clocks).
The break to worry about is *cost*: BMIQ fits a three-state mixture per sample.
Minutes for the 95-sample test set at 23k probes; on a full 450k cohort, an hour
is possible. `nfit` (probes subsampled for the fit) mitigates; still, this must
be a button with progress and a stale-notice, the pattern `r_clockset` already
implements (`methylome_server.R:159-183`).

**Effort and risk.** 2–3 days for option 2, mostly wiring and cache correctness,
plus real runtime risk. What could make it not worth doing in-board: if the
preprocess-single-source refactor lands BMIQ in playbase's methylation ingest,
option 1 gets it for free and the in-board toggle becomes dead plumbing. Worth a
decision before building, not after.

**Testing.** Headless: corrected matrix same dim/rownames, all values in [0,1],
type I probes nearly unchanged, type II beta distribution's modes move toward the
type I modes (compare mode positions before/after), idempotence not required but
determinism per seed is. Clocks on corrected betas still correlate with age at
r ≥ the uncorrected value. Browser-only: the density overlay and the stale
notice.

---

## 5. Consensus clustering + t-SNE/UMAP — the scepticism item

> **Status 2026-08-19: consensus clustering deliberately SKIPPED. The
> t-SNE/UMAP half of this item is no longer a routing question — the route was
> removed.** This is the item the architecture change hit hardest, so read the
> status note rather than the plan.
>
> The scepticism below was right and is now more so. Everything this item said
> routing genuinely covered has stopped being true:
>
> - `MODULES_METHYLOMICS` no longer exists. `MODULES_ACTIVE` is an all-`FALSE`
>   array for methylomics (`opg_server.R:144`), so `board.clustering` is never
>   loaded for a methylation dataset.
> - The pgx no longer carries what that board `req()`s. The methylomics compute
>   sets `do.cluster <- FALSE` and `add.gmt <- FALSE`
>   (`upload_module_computepgx.R:1143-1148`), so there is no `cluster$pos` with
>   precomputed pca/tsne/umap and no `gsetX`. The `etc/OPTIONS:63` reference
>   below is stale; the file now carries a comment where the list used to be.
>
> So the split this item recommended — "t-SNE/UMAP: close as routing" — cannot
> be taken. Sample embeddings on methylation data are now a build, in this app,
> or they do not exist. The build itself is small (PCA/t-SNE/UMAP of samples
> never touches probe IDs, which is why the routing argument worked in the
> first place); what changed is who has to do it.
>
> Consensus clustering is unchanged and still deliberately skipped:
> `ConsensusClusterPlus` is still not installed, the platform still has no
> consensus clustering anywhere, and the deciding argument below still holds —
> no cancer datasets in the library. Build on demand, not on archetype.

**§4's claim checked.** "Routing, not backlog" is *two-thirds* true, and the
false third is the headline word.

What routing genuinely covers: `board.clustering` is in `MODULES_METHYLOMICS`
(`etc/OPTIONS:63`) already, and the test pgx carries everything it `req()`s —
`gsetX` (5126 × 95), `families`, and precomputed `cluster$pos` with
pca/tsne/umap in 2d and 3d, built at upload. PCA/t-SNE/UMAP of samples renders
fine on a methylation pgx; sample-level plots never touch probe IDs. The
X/Y-filter checkbox even survives the cytoband trap — its regex (`^X|^Y`,
`clustering_server.R:323`) happens to match `"Xp22.2"`.

What routing does not cover:

- **Consensus clustering does not exist anywhere in the platform.**
  `board.clustering` offers a k-cut of a single hclust/kmeans
  (`hm_clustk`, `clustering_server.R:543-549`) — no resampling, no consensus
  matrix, no delta-CDF, no stability evidence. The Capper-style cancer
  centrepiece is *not* a routing question; §4 is wrong to fold it into one.
  And `ConsensusClusterPlus` is **not installed** (checked; `M3C` also absent).
  It is a small, pure-R Bioconductor package — an ordinary install, not a
  LUMP-class blocker — but it is a dependency decision someone must make.
- **The heatmap labels probes as gene symbols** —
  `probe2symbol(..., fill_na = TRUE)` (`clustering_plot_splitmap.R:245`) — so a
  methylation user sees duplicate symbol rows (several probes per gene) and
  ~24% unmapped fillers on this dataset. Usable, mildly confusing, not wrong.
- **The cluster-annotation panel correlates against `gsetX`**
  (`clustering_server.R:635`), the symbol-collapsed, probe-density-inflated
  geneset matrix that §5 of the review condemns for the enrichment boards. The
  annotation half of the clustering board carries the exact bias this board's
  gometh path exists to avoid.

**Plan.** Split the item. t-SNE/UMAP: close as routing — at most a pointer in
this board's info text; do not rebuild embeddings here. Consensus clustering: a
real build if and when a cancer-archetype user exists. Shape: install
`ConsensusClusterPlus`; new panel (own nav_panel beside the EWAS tab or under
Methylome character), `ConsensusClusterPlus(d = beta of top-N most-variable
probes (SD, N ≈ 5000), maxK = 6, reps = 500, pItem = 0.8)` behind an explicit
button — minutes-scale, the `regions_val` pattern (`methylome_server.R:332-344`)
is the template. Outputs: consensus heatmap for chosen k, delta-CDF curve, and
the per-sample class written where the phenotype pickers can see it. Consumes
nothing from `mp_fit_ewas`; sits fully beside it. Effort 3–4 days.
Testing: headless assertions that k=2..4 returns a partition of all samples,
the known monocytes/whole_blood split is recovered at k=2 (ARI ≈ 1 on this
dataset), and results are deterministic under a set seed; the heatmap itself is
browser-only.

**What could make it not worth doing:** no cancer datasets in the library. The
review itself says cancer is where the papers are, but this board's users so far
are EWAS/clocks-shaped. Build on demand, not on archetype.

---

## 6. Survival (KM / Cox) — the contradiction resolved

> **Status 2026-08-19: deliberately SKIPPED, and reclassified.** The deferral
> argument at the end of this item is the one that won: every dataset the app
> can currently load would see only the refusal message. It was additionally
> ruled a **separate module** rather than a Methylome screen — survival on the
> loaded dataset's own follow-up is not methylation-specific, and building it
> behind this app's datatype guard would put it out of reach of every other
> datatype that has the same need.
>
> The design below stays on file and is still the design to build from,
> including the "never guess the time/event columns from names" rule, which is
> the part most likely to be got wrong by whoever picks it up.

**Resolving §4 vs the impossible list.** Both entries are correct and the
condition is the whole feature: survival is *feasible from a beta matrix* if and
only if the sample sheet carries time-to-event columns, and *impossible* when it
doesn't — which is the common case; the board's own test dataset
(`cell_type, age, sex, ms_status, smoking, slide, geo_title`) has no follow-up
at all. §4 should have carried its own condition inline; as written the two
lists contradict each other only because §4's entry dropped the "where the
metadata carries follow-up" clause into fine print.

**Is it this board's job?** If built anywhere, yes here: `board.tcga` errors on
CpG IDs (§5) and answers a different question — the user's *signature* against
TCGA's cohorts, not the user's *own* follow-up. Nothing else in the platform
does KM on the loaded dataset's metadata.

**What it would be.** `survival::survfit(Surv(time, event) ~ group)` +
`survival::survdiff` (log-rank) for the KM panel, `survival::coxph` with
age/sex/chosen covariates for the adjusted hazard ratio. `survival` is
installed; `survminer` is not and is not needed — `plot(survfit)` in base
graphics is the house style anyway. The grouping variable is the interesting
design choice, and the honest one is to reuse what the board already computes:
the EWAS `strata` (r_ewas), a median-split of a picked CpG's beta, the
acceleration residual from the clocks tab, or a consensus class from item 5.
Column detection: never guess from names alone — two selectInputs (time: any
positive numeric column; event: any 0/1 or two-level column) populated only when
candidates exist, and a `validate(need(...))` refusal ("this dataset carries no
follow-up columns") otherwise, which is what every current user would see.

**Where.** A new panel + settings group, own file (`methylome_survival.R`),
sitting fully beside `mp_fit_ewas` — it consumes `r_ewas()$meta$strata` at most.
Cheap to compute; no button needed beyond the pickers.

**What it breaks.** Nothing — pure addition. The risk is statistical: a
median-split KM on a CpG chosen *because* it was a top hit is circular and will
look impressive on noise. The info text must say the p-value is descriptive, not
a validated prognostic claim.

**Effort.** 2 days. **But argue for deferral** — see the sequence. A feature
whose refusal branch is the only branch any current dataset can exercise is
speculation; the day a dataset with follow-up lands in the library is the day
this earns its build.

**Testing.** Entirely via fabricated columns: add exponential times and a 0/1
event to a pgx copy with a planted group effect; assert `survfit` strata sizes,
log-rank p < 0.05 on the planted effect and ~uniform on a permuted one, Cox HR
sign matches the plant, and the no-follow-up refusal fires on the unmodified
pgx. The KM rendering is browser-only.

---

## 7. EWAS Catalog lookup

> **Status 2026-08-19: DONE** — same commit `75da3d588`, with the pure half
> living in `playbase.epigenetics` rather than in this app, which is a better
> home than the plan proposed: `ewas_catalog_lookup()` and its parser
> `ewas_catalog_traits()` are pure functions with no Shiny in them, so the
> tests are headless and the upload path could reuse them.
>
> Deviations from the plan, all small: `curl` directly rather than `httr`
> (`curl` is already a `playbase.epigenetics` Import and `httr` would have been
> a new one); the function returns `NULL` on the *first* failure rather than a
> partial frame, for the reason its own docs give — half an answer reads as
> "these CpGs are unreported", the one wrong conclusion the lookup exists to
> prevent; and traits are summarised as `"Smoking (12); Age (3); +2 more"`
> counting distinct `study_id`s. The button, the `reactiveVal`, the
> invalidate-on-new-`r_ewas` and the optional `Reported` column on the hits
> table are all where the plan put them (`methylome_ewas_inputs.R:72`,
> `methylome_server.R:363-387`, `methylome_ewas.R:233,260-264`).
>
> **Still unanswered, and it is the deploy-time question this item flagged:**
> nobody has confirmed the container has outbound network. The failure mode is
> a button that reports "not reachable", not a hang, so this is low-risk — but
> do not put the feature in release notes until someone has clicked it on a
> deployed instance.

**What it is.** Marking hits known-vs-novel against the EWAS Catalog's REST API:
`GET https://ewascatalog.org/api/?cpg=cg05575921` returns JSON (`results` array;
verified reachable from this machine, HTTPS — plain http 301-redirects). One
request per CpG; no documented batch endpoint, so cap at the top ~50 hits.
EWAS Atlas (NGDC) has no stable public API worth depending on; the Catalog is
the one to use, and the fallback is its downloadable full-catalog TSV (~100 MB,
a data-shipping decision this plan does not recommend for a first cut).

**Dependency.** `httr` + `jsonlite`, both present. The real dependency is
**network egress from the deployed container**, which no `requireNamespace` can
check: behind ShinyProxy the app may or may not reach the outside world, and the
feature must degrade to a message ("catalog not reachable"), never a hang — set
`httr::timeout(5)` per request and stop at the first failure.

**Where.** Beside the model, not in it. An action button in the hits settings
group ("Look up hits in EWAS Catalog"), the `enrich_val` reactiveVal pattern
(`methylome_server.R:355-368`): on click, query the top-N significant CpGs with
`withProgress`, store `probe → traits (n studies)` in a reactiveVal, invalidate
on a new `r_ewas()`. The hits table (`methylome_table_hits_server`,
`methylome_ewas.R:229-274`) takes an optional `r.catalog` reactive and, when
non-NULL, appends a `Reported` column ("Smoking (12); Age (3)" / "novel").
A tiny pure function `mp_catalog_lookup(cpgs)` holds the HTTP so it tests
headless.

**What it breaks.** Nothing; a NULL reactive leaves the table exactly as today.

**Effort and risk.** 1 day. Risks: the API's uptime and undocumented rate
limits (50 serial GETs is polite enough); the deploy-time egress question;
and interpretive — "not in the catalog" means "not in the catalog", not novel
biology, and the column header tooltip should say so.

**Testing.** Headless: `mp_catalog_lookup` parsing asserted against a canned
JSON fixture (no network in CI); URL construction; empty-result and
timeout paths return the degrade value. One optional live smoke test guarded by
`curl::has_internet()`, skipped otherwise. The button flow and column render are
browser-only.

---

## What is actually left to do (2026-08-19)

The sequence below is the 2026-08-17 recommendation, kept for the record. Most
of it has happened, and not in this order — items 2, 3 and 7 landed together on
2026-08-19. What remains, in the order it is worth doing:

1. **Add `dmrff`, `missMethyl` and `bacon` to `docker/Dockerfile.update`.** Not
   in this plan at all, because the plan asked `requireNamespace()` on a dev
   machine. It outranks everything below: three finished features are dead in
   production without it, and the fix is three lines beside the ones that
   already install the Illumina manifests and `methylclock`.
2. **Per-probe missingness filter** (item 1) — the only unbuilt item here whose
   case is unchanged and unarguable. Half a day.
3. **Sample embeddings**, if anyone wants unsupervised views on methylation
   data. Item 5's routing answer is gone; PCA/t-SNE/UMAP of samples would now
   have to be built in this app. Small, but it is a build now, not a route.
4. **Consensus clustering** (item 5) and **survival** (item 6) stay parked on
   the same conditions as before: a cancer dataset in the library, and a
   follow-up-bearing dataset respectively. Survival is additionally scoped as a
   separate module, not a Methylome screen.

Also worth folding into whatever comes next, from FIELD-REVIEW §6's open list:
`methylome_model.R` still holds manifest reading, array inference and probe
masking that `playbase.epigenetics` already owns — the next round of the split
is in this app, not in playbase.

### The 2026-08-17 sequence, for the record

1. **Per-probe missingness filter** — smallest diff, pure correctness, zero UI.
   Do first; it also touches the exact lines the next two items edit, so land it
   before the model file gets busy.
2. **Multi-level F-test** — the only item where the board currently *misdirects*
   users to a wrong test rather than merely lacking a feature, and the traced
   consumer list shows the per-sample panels are already k-level-safe: the work
   is the model branch, the picker, and three honest refusals.
3. **SVA** — installed, extends the same function, high scientific value, and
   its return shape breaks nothing. After the F-test lands, the SV columns must
   also append on the anova branch — another reason for this ordering.
4. **EWAS Catalog** — one day, self-contained, visible in the flagship table.
   Do it fourth because it depends on nothing above and blocks nothing below;
   confirm container egress before promising it in release notes.
5. **BMIQ** — installed and real, but decide the playbase-vs-in-board question
   first: if the preprocess refactor will normalise at ingest, the in-board
   overlay is throwaway plumbing. If in-board wins, budget for the runtime and
   reuse the clock-cache pattern.
6. **Consensus clustering** — build only when a cancer dataset actually enters
   the library, and install `ConsensusClusterPlus` then. Do **not** rebuild
   t-SNE/UMAP here in any case: routing to `board.clustering` genuinely covers
   the embeddings, and the parts it covers badly (symbol row labels, the
   gsetX-based annotation panel) are that board's defects to fix, not this
   board's features to duplicate.
   **— No longer true. Do not follow this line.** Methylomics loads no
   Dashboard board and its pgx has no `cluster$pos` and no `gsetX`; see the
   status note on item 5.
7. **Survival — argue against building now.** Every dataset the board can
   currently load would see only the refusal message; a feature that can only
   refuse is a maintenance liability wearing a menu entry. Keep the design above
   on file and build it the week a follow-up-bearing dataset exists. If product
   pressure forces it earlier, build the KM panel only and skip Cox — half the
   surface, same refusal logic.

Corrections to §4 found while planning, for the next review revision:
"consensus clustering + t-SNE/UMAP is routing, not backlog" is wrong for the
consensus half — the platform has no consensus clustering and the package is not
installed; the survival entry needs its metadata condition stated where the
claim is, not two sections later; and BMIQ, while genuinely unblocked, is a
normalisation with app-wide plumbing and minutes-to-hours of runtime, not a
peer of the cheap items it was listed beside.

**All three were carried into FIELD-REVIEW v4 (2026-08-19)**, with the BMIQ one
sharpened further: it is not merely unblocked, it is already running at ingest.
Two corrections went the other way, from the code back into this plan — the
routing answer for embeddings no longer exists, and `requireNamespace()` on a
dev machine was never evidence about the deployment image.
