# App bugs found while building the compute-comparer parity harness

Two pre-existing bugs in the `master` upload flow, found by driving the real UI with
Playwright and diffing the resulting `params.RData` against a script-built one. Neither
is caused by the preprocess refactor; both are independent of it and affect users today.

Harness and full parity results: `dev/compute-comparer/`.

---

## 1. The chosen batch-correction method is silently reset to ComBat

**Severity: user-visible wrong results.** A user who picks limma, RUV, SVA or NPM can get
a dataset computed with **ComBat** instead, with no warning.

`components/board.upload/R/upload_module_normalization.R:253-257`, inside the
batch-correction comparison reactive:

```r
methods <- c("ComBat", "limma", "RUV", "SVA", "NPM")
if (ncol(X0) > 100) methods <- methods[methods != "NPM"]
shiny::updateSelectInput(
  session,
  "bec_method",
  choices = methods          # <-- no `selected = input$bec_method`
)
```

`updateSelectInput()` with `choices` and no `selected` makes Shiny select the **first**
choice. The first choice is `"ComBat"`. Every time this comparison re-runs (it is
reactive on the data and on `bec_param`), the user's selection is discarded.

Downstream, `bc_method` (line 1357) reads `input$bec_method`, which now says ComBat, and
that is what reaches `params$batch.correct.method` and the compute.

### When it fires

The enclosing reactive (`results_correction_methods`) depends on `cleanX()`, `imputedX()`,
`r_samples()`, `r_contrasts()`, `input$bec_param` and `input$bec_full_features` — but
**not** on `input$bec_method`. So selecting a method does not by itself retrigger the
reset. It fires whenever the comparison re-runs, i.e. when the user:

- changes the **batch parameter** (`bec_param`) — the most likely case, since picking the
  parameter right after picking the method is the natural order;
- toggles **"Use all features for BC preview"**;
- changes any upstream normalization / imputation / outlier setting.

**Observed:** driving the UI to select `limma` or `NPM` produced
`params$batch.correct.method = "ComBat"` on **every** attempt, across separate browser
sessions, under all five of:

1. `Shiny.setInputValue`, method then batch parameter;
2. `Shiny.setInputValue`, batch parameter then method;
3. re-asserting the method immediately before triggering compute;
4. driving `selectize.setValue()` (the path a real click takes) instead of `setInputValue`;
5. (4) plus verifying the widget held the value for 15 s with no drift, and excluding the
   BC inputs from the pre-compute re-assert so nothing re-invalidated the comparison.

In case 5 the **client widget stably reported `limma`** while the app still computed and
recorded `ComBat`. So the client and the server disagree: the selectize widget holds the
user's choice, `input$bec_method` server-side does not.

### What is and is not established

Established: `updateSelectInput(session, "bec_method", choices = methods)` with no
`selected=` is incorrect and will reset the selection to the first choice; and the app
records `ComBat` regardless of what the widget displays, reproducibly, from five
different client strategies.

**Not** established: the exact propagation path that leaves client and server
out of sync. It could not be pinned from a scripted client. **The cheapest next step is a
human clicking the dropdown manually**: select `limma`, pick a batch parameter, hit
Compute, and read `params$batch.correct.method` from the resulting `raw_*/params.RData`.
That distinguishes "scripted-client artefact" from "every user is affected" in about a
minute. Until then, treat the severity as unconfirmed but the code defect as real.

Note a `ComBat` scenario "passes" trivially: ComBat is the value the reset lands on, so a
test that only ever checks ComBat will never catch this. The harness's `bc-combat`
scenario passed for exactly this reason and should not be read as coverage.

### Unrelated: ComBat errors on the `human-symbol` fixture

With `cluster` as the batch parameter, the run log shows:

```
Error in sva::ComBat(cX, batch = b, mod = mod) :
  At least one covariate is confounded with batch!
Coefficients not estimable: batch1 batch2
```

`cluster` is confounded with `group` in that fixture, so it is a poor batch variable.
Any batch-correction scenario needs a batch column independent of the contrast design;
this is a fixture-choice caveat, not an app defect.

Note this also means a `ComBat` scenario "passes" trivially: ComBat is the value the
reset lands on, so a test that only ever checks ComBat will never catch this.

**Fix:** preserve the selection.

```r
shiny::updateSelectInput(
  session, "bec_method",
  choices  = methods,
  selected = if (input$bec_method %in% methods) input$bec_method else methods[1]
)
```

Note the `ncol(X0) > 100` branch drops NPM from the list, so a user who selected NPM on a
small dataset and then grows it needs the fallback — hence the `%in%` guard.

---

## 2. Every wizard lock/unlock branch is dead code

**Severity: no data corruption, but the wizard never locks.** Users can page past a step
whose upload failed validation, and the "please finish this step" alerts never appear.

`components/board.upload/R/upload_server.R` compares the wizard input against **step ids**:

```r
req(input$upload_wizard == "step_counts")        # line ~918
req(input$upload_wizard == "step_samples")       # line ~932
if (input$upload_wizard == "step_comparisons")   # line ~894
```

But `input$upload_wizard` carries the step **title**, not the id. The installed
`wizardR::wizard_step()` signature is:

```r
wizard_step <- function(..., step_title = NULL, session = getDefaultReactiveDomain())
```

There is no `step_id` (nor `server`) parameter. `upload_server.R:587-640` passes both
anyway, so they fall into `...` and are rendered as stray HTML attributes. The Shiny
input binding's `getValue()` returns `$(el).attr("data-title")`
(`wizardR/inst/main.js`), i.e. `"Step 1: Upload counts"`, `"Step 4: QC/BC"`, `"Compute!"`.

So `"Step 1: Upload counts" == "step_counts"` is always FALSE. Affected:

- the four alert branches at `upload_server.R:882-905`
- the lock/unlock observers at `upload_server.R:916-960`

Combined with `lock_start = FALSE`, the wizard is **never locked at any point**.

**Fix:** compare against the titles, or give the steps a real id. The cleanest is to add a
`step_id` parameter to `wizardR::wizard_step()` that emits `data-step-id`, and have the
binding return that — otherwise the app's ids and titles must be kept in sync by hand.

A minimal in-app fix without touching wizardR:

```r
STEP_TITLES <- c(
  step_counts      = "Step 1: Upload counts",
  step_samples     = "Step 2: Upload samples",
  step_comparisons = "Step 3: Create comparisons",
  step_qc          = "Step 4: QC/BC",
  step_compute     = "Compute!"
)
req(identical(input$upload_wizard, STEP_TITLES[["step_counts"]]))
```

---

## Not a bug, but worth knowing

**`playbase` main and edgy annotate genes differently.** For the same 7439 features in
the same order, `main` keeps multi-mapped symbols joined (`A2M;IGHA2`,
`CRYBG1;AURKB;SLC45A2`) while `edgy` collapses to the first (`A2M`, `CRYBG1`). 465
`symbol`, 467 `gene_title` and 25 `human_ortholog` values differ on the `human-symbol`
fixture. This is the dominant source of `pgx$genes` differences between the two branches
and is unrelated to the preprocess work — worth confirming it is intentional.

**QRILC imputation is stochastic.** Two `pgx.preprocess()` calls on identical input in
one session differ by mean relative difference ~0.022. Any test asserting bit-identical
output under `impute_method = "QRILC"` (or `MinProb`) will flake. `SVD2` is deterministic.
