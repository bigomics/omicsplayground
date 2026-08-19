# Driving the Omics Playground upload flow with Playwright

Knowledge file for `drive_app.mjs`. Everything here was read off the `master`
(= `compute-comparer`) tree; line numbers are from that branch.

Goal of the automation: get the app to write `params.RData` for a chosen dataset +
preprocessing scenario, and copy that file out before the app deletes it. We do **not**
wait for the app's own compute — see "Why we stop at params.RData".

---

## 1. Launching the app

```bash
R_LIBS_USER=dev/compute-comparer/libs/pb-main Rscript dev/run_app_headless.R
```

- `dev/run_app_headless.R` is `run_app.R` with `launch.browser = FALSE` — the right one
  for automation.
- Serves on `0.0.0.0:3838`.
- It `setwd()`s into `components/app/R` and sources `global.R`, `ui.R`, `server.R`.
- `etc/OPTIONS` has `AUTHENTICATION = none`, so there is **no login step**. If that ever
  changes, the driver needs a login leg.
- `R_LIBS_USER` precedes site-library in `.libPaths()`, so it decides which playbase the
  app runs against. This is how Path A is pinned to playbase `main`.

Readiness check: poll `http://localhost:3838` until it returns 200. R takes ~30-60s to
source everything, so give it a generous timeout.

## 2. Page structure

The upload board is mounted as `UploadUI("upload")` (`components/app/R/ui.R:577`), so
**every input id in this flow is prefixed `upload-`**, and module children add another
segment:

| Module | Prefix |
|---|---|
| upload board itself | `upload-` |
| counts preview | `upload-counts_preview-` |
| samples preview | `upload-samples_preview-` |
| contrasts preview | `upload-contrasts_preview-` |
| normalization / QC panel | `upload-checkqc-` |
| compute panel | `upload-compute-` |

CSS-selecting these needs escaping (`#upload-checkqc-normalize` is fine, but ids with
dots would not be). Prefer `page.locator('#id')` with plain ids as above.

## 3. The landing tab

`components/board.upload/R/upload_ui.R`. A `tabsetPanel(id = "upload-upload_tabs")` with
two tabs; we want the default one, `value = "upload"`.

| Input | id | Notes |
|---|---|---|
| Data type | `upload-selected_datatype` | plain `selectInput`. Values: `RNA-seq`, `mRNA microarray`, `proteomics`, `scRNA-seq`, `methylomics`, `metabolomics`, `multi-omics` |
| Organism | `upload-selected_organism` | **server-side selectize** (`upload_server.R:85`, `server = TRUE`) |
| Start | `upload-start_upload` | actionButton, opens the wizard modal |

The organism widget is populated by AJAX as you type, so there is no static
`<option>` list to click. Set it directly:

```js
await page.evaluate(() =>
  Shiny.setInputValue('upload-selected_organism', 'Human', {priority: 'event'}));
```

`observeEvent(input$selected_organism)` (`upload_server.R:760`) just copies it into
`upload_organism()`, so a direct set is fully equivalent to typing.

## 4. The wizard

Built with `wizardR` (`upload_server.R:587-654`), rendered into `uiOutput("upload-upload_wizard")`.

Five steps, in order:

| # | `step_id` | `step_title` (= `data-title`) | Content |
|---|---|---|---|
| 1 | `step_counts` | `Step 1: Upload counts` | counts.csv |
| 2 | `step_samples` | `Step 2: Upload samples` | samples.csv |
| 3 | `step_comparisons` | `Step 3: Create comparisons` | contrasts.csv |
| 4 | `step_qc` | `Step 4: QC/BC` | normalization panel (this is what we vary) |
| 5 | `step_compute` | `Compute!` | methods + dataset name |

It is `modal = TRUE`, so the whole thing lives inside a Bootstrap modal with id
`wizard-modal-upload-upload_wizard`.

### DOM (from `wizardR/inst/Wizard-JS/wizard.min.js`)

```
.wizard
  .wizard-nav
  .wizard-content
    .wizard-step[data-step][data-title]   <- .active marks the current one
  .wizard-buttons
    .wizard-btn.btn.next
    .wizard-btn.btn.prev
    .wizard-btn.btn.finish
```

Advance with `.wizard-btn.btn.next`, go back with `.prev`.

### GOTCHA: `input$upload_wizard` is the TITLE, not the step_id

`wizard_step()` in the installed wizardR takes only `(..., step_title, session)` — it has
no `step_id` or `server` parameter. `upload_server.R` passes both anyway, so they fall
into `...` and land as stray HTML attributes on the div. The Shiny input binding's
`getValue()` returns `$(el).attr("data-title")`, i.e. the **title string**.

Consequence: every `req(input$upload_wizard == "step_counts")` /
`== "step_samples"` / `== "step_comparisons"` / `== "step_compute"` branch in
`upload_server.R` (the alerts at 882-905 and the lock/unlock observers at 916-960) is
comparing a title against an id and **never fires**. Combined with `lock_start = FALSE`,
the wizard is effectively **never locked**.

For the driver that means:

- Do **not** wait for an unlock signal — it will never come.
- The next button is always clickable, including before the upload has been validated.
  So the driver must gate on its own readiness signal (see §6) or it will happily walk
  past a step whose data has not loaded, and `params.RData` will be built from nothing.

This is a pre-existing app bug, out of scope for the harness. It is recorded here
because it changes how the driver has to wait.

### The finish button

`wizardR/inst/utils.js`:

```js
$(document).on('click', '.wizard-btn.btn.finish', function () {
    ... modal.hide();
    Shiny.setInputValue(modalId, "wizard_finished", {priority: "event"});
});
```

and `upload_module_computepgx.R:1005` reacts to exactly that:

```r
shiny::observeEvent(upload_wizard(), {
  if (!is.null(upload_wizard()) && upload_wizard() != "wizard_finished") return(NULL)
  ...
})
```

So compute can be triggered either by clicking `.wizard-btn.btn.finish`, or directly:

```js
await page.evaluate(() =>
  Shiny.setInputValue('upload-upload_wizard', 'wizard_finished', {priority: 'event'}));
```

The direct form is preferred for the harness: it skips Wizard-JS's own DOM bookkeeping,
and the observer re-validates everything itself anyway (`iv$is_valid()` on
`selected_name`/`selected_description`, plus a `probetype()` readiness gate at
`upload_module_computepgx.R:1020-1030`). Click the real button only when testing the UI
path itself.

## 5. File uploads

`fileInputArea` widgets, standard Shiny file inputs underneath:

| Step | id |
|---|---|
| counts | `upload-counts_preview-counts_csv` |
| samples | `upload-samples_preview-samples_csv` |
| contrasts | `upload-contrasts_preview-contrasts_csv` |

Use `setInputFiles` on the underlying `input[type=file]`:

```js
await page.locator('#upload-counts_preview-counts_csv').setInputFiles(countsPath);
```

Wait for the input with `{ state: 'attached' }` and a long timeout before calling
`setInputFiles` — a step's UI is rendered server-side and can easily take longer than
`setInputFiles`' 30s default to appear (step 3 builds its comparison UI from the sample
table first).

### GOTCHA: step 3 has no file input until you switch it to upload mode

`show_comparison_builder` is a `reactiveVal` initialised to **TRUE**
(`upload_server.R:32`), and the contrasts module renders the `fileInputArea` only in the
`if (!show_comparison_builder())` branch (`upload_module_preview_contrasts.R:244`). So
step 3 opens on the *interactive comparison builder* and
`#upload-contrasts_preview-contrasts_csv` **does not exist in the DOM at all** — any wait
for it runs to timeout.

Click "Upload comparisons file" first — `input$goUploadComparison`
(`upload_module_preview_contrasts.R:109`, observer at 306) sets the flag FALSE and
renders the upload area:

```js
await page.click('#upload-contrasts_preview-goUploadComparison');
await settle(page);
// only now does the contrasts file input exist
```

Steps 1 and 2 have no such gate; their upload areas render immediately.

Shiny uploads asynchronously, then runs `pgx.checkINPUT` server-side and writes a
`CHECKS_OUTPUT` line into `raw_dir` (`upload_utils.R:301-315`). Wait for the preview
table to render before advancing.

Note `upload_module_preview_counts.R:784` treats a file literally named `params.RData`
in the counts upload as a special re-import path — irrelevant here, but do not name
fixture files that.

## 6. Waiting reliably

Since the wizard never locks, use these signals instead:

1. **Shiny idle**: `document.querySelectorAll('.shiny-busy').length === 0`, plus
   `document.documentElement.classList.contains('shiny-busy')`. Necessary but not
   sufficient (there are gaps between reactive flushes) — require it to stay quiet for
   a second or two.
2. **Preview table populated**: each step renders a table once its CSV is accepted.
3. **Data-driven check**: the most reliable gate is a round-trip through the server —
   e.g. wait until `#upload-checkqc-normalization_method` is *attached* on step 4, which
   only happens once counts/samples/contrasts have all been accepted.

Avoid fixed sleeps; the app's reactive graph is slow and variable under load.

### GOTCHA: a step's controls do not exist until you are ON that step

Shiny's `suspendWhenHidden` defaults to TRUE, so an output inside a `display: none`
wizard step is **suspended and never renders**. The QC panel's `uiOutput` is in step 4,
which is hidden while you are on step 3 — so

```js
await nextStep(page, '#upload-checkqc-normalization_method', 'step 3');  // DEADLOCK
```

waits for a control that cannot appear until the click that the wait itself is blocking.
Always click next **first**, then wait for the next step's controls:

```js
await nextStep(page, null, 'step 3 (contrasts)');
await page.waitForSelector('#upload-checkqc-normalization_method', { state: 'attached' });
```

The same applies to `#upload-compute-selected_name` on step 5.

### GOTCHA: never wait on `.recalculating`

Measured on the welcome tab: `.shiny-busy` is 0 but **`.recalculating` sits permanently
at 33**. Outputs belonging to tabs that were never shown never finish recalculating, so
any idle-check that includes `.recalculating` blocks forever. Use `.shiny-busy` only.

### GOTCHA: `selectInput` is a selectize widget, so the `<select>` is never "visible"

Shiny's `selectInput` defaults to `selectize = TRUE`: it renders a `.selectize-control`
and sets `display: none` on the original `<select>`. Playwright's `waitForSelector`
defaults to `state: 'visible'`, so

```js
await page.waitForSelector('#upload-selected_datatype');           // hangs forever
await page.waitForSelector('#upload-checkqc-normalization_method'); // hangs forever
```

Use `{ state: 'attached' }` for any select id, or gate on a genuinely visible control
such as the `#upload-start_upload` action button. `setInputFiles` needs no visibility, so
the file inputs are fine either way.

## 6b. Overlays you must clear first

Two things cover the page and silently eat clicks. Both present as a `page.click`
timeout on an element that `waitForSelector` just reported as visible.

**1. The sign-in splash.** Even with `AUTHENTICATION = none`, `AuthenticationModule.R:29-40`
shows `splashLoginModal()` with a single click-through button labelled **"Sure I am!"**
(id `ns("login_emailSubmit")`). Nothing is reachable until it is dismissed. Match it by
label rather than id so the auth namespace can move:

```js
await page.getByRole('button', { name: 'Sure I am!' }).first().click();
```

**2. The chirp promo modal.** `ENABLE_CHIRP = TRUE` in `etc/OPTIONS` (with
`etc/chirp_data.csv`) pops a Bootstrap modal — a video carousel, e.g. "Drug Discovery
Data Analysis" — over the upload board. It appears *after* the board renders, so a
`waitForSelector` on `#upload-start_upload` succeeds and the click that follows still
times out.

Dismiss both generically before each real click, and **spare `wizard-modal-*`** — that
one is the wizard itself:

```js
document.querySelectorAll('.modal.show').forEach((m) => {
  if (m.id?.startsWith('wizard-modal-')) return;
  bootstrap.Modal.getInstance(m)?.hide()
    ?? m.querySelector('.btn-close, [data-bs-dismiss=modal]')?.click();
});
document.querySelectorAll('.swal-overlay--show-modal .swal-button').forEach((b) => b.click());
```

Also strip any orphaned `.modal-backdrop` once no `.modal.show` remains — a leftover
backdrop still intercepts pointer events.

Turning chirp off in `etc/OPTIONS` would also work, but that is a tracked repo file;
dismissing client-side keeps the harness free of app-config edits.

## 6c. Navigating to the upload board

The board is a `bigdash` `bigTabItem("upload-tab", ...)`. Its nav trigger is

```html
<a data-target="upload-tab" class="dropdown-item cursor-pointer tab-trigger">Upload new</a>
```

which lives inside a **closed navbar dropdown**, so it is not visible and `page.click()`
fails with "Element is not visible". jQuery fires the bound handler regardless
(`bigdash/srcjs/sidebar.js:207` → `toggleTabs(target)`), which is exactly what bigdash's
own `bigdash-select-tab` message handler does:

```js
await page.evaluate(() => window.$('.tab-trigger[data-target="upload-tab"]').trigger('click'));
```

Verified working: the third `#big-tabs .big-tab` flips from `d-none` to `display: block`.

## 7. Step 4 — the normalization panel (what the scenarios vary)

`upload_module_normalization.R:1055-1180`. All inside a `bslib::accordion` with
`multiple = FALSE`, so **only one panel is open at a time** — opening "Normalization"
collapses "Missing values". Inputs inside a collapsed accordion panel still exist in the
DOM and their values are still bound, so setting them via JS works regardless of which
panel is visible. Clicking them does not, unless the panel is expanded first.

| Input | id (`upload-checkqc-` +) | Type | Choices / range | Default |
|---|---|---|---|---|
| Treat zero as NA | `zero_as_na` | checkbox | | `FALSE` |
| Remove NA rows | `filtermissing` | checkbox | | `FALSE` |
| NA threshold | `filterthreshold` | select | `0.1`, `0.2`, `0.5`, `3`, `-0.5` | `0.2` |
| Impute NA | `impute` | checkbox | | `DEFAULTS$qc$impute` |
| Impute method | `impute_method` | select | `SVD2`, `QRILC`, `MinProb`, `Perseus` | `SVD2` |
| Normalize data | `normalize` | checkbox | | `TRUE` (`FALSE` for Olink/NULISA) |
| Normalization method | `normalization_method` | select | datatype-dependent, see below | first choice |
| Reference gene | `ref_gene` | **server-side selectize** | rownames of counts | — |
| Remove outliers | `remove_outliers` | checkbox | | `FALSE` |
| Outlier threshold | `outlier_threshold` | slider | 1-12 step 1 | `6` |
| Batch correct | `batchcorrect` | checkbox | | `FALSE` |
| BC method | `bec_method` | select | `ComBat`, `limma`, `NPM`, `RUV`, `SVA` | `SVA` |
| BC parameter | `bec_param` | selectize, multiple | **reactive** — sample table columns | — |

`normalization_method` choices depend on `upload_datatype()`:

- RNA-seq / microarray: `CPM`, `CPM+quantile`, `TMM`, `quantile`, `maxMedian`, `maxSum`, `reference`
- proteomics / metabolomics: `maxMedian`, `maxSum`, `quantile`, `reference`
- methylomics: `BMIQ`, `quantile`
- multi-omics: `median`

Several of these sit behind `conditionalPanel` (`filterthreshold` under
`input.filtermissing == true`, `impute_method` under `input.impute == true`,
`normalization_method` under `input.normalize == true`, `ref_gene` under
`normalization_method == 'reference'`, `outlier_threshold` under
`input.remove_outliers == true`, `bec_method`/`bec_param` under
`input.batchcorrect == true`). `conditionalPanel` is client-side only — the server reads
`input$...` whether or not the panel is visible.

### Interaction strategy

- **Checkboxes**: real click is fine and exercises the real path, but the accordion panel
  must be open. Setting via JS is equivalent and order-independent — prefer JS.
- **selectInput**: `selectElement.value = v` + `dispatchEvent(new Event('change'))`, or
  just `Shiny.setInputValue`.
- **server-side selectize** (`ref_gene`, `bec_param`) and **slider**
  (`outlier_threshold`): no practical DOM path — use `Shiny.setInputValue`.

Uniform helper:

```js
const setInput = (page, id, value) => page.evaluate(
  ([id, value]) => Shiny.setInputValue(id, value, {priority: 'event'}),
  [id, value]);
```

Caveat: `Shiny.setInputValue` sets the server-side value but does not update the widget's
own display, and a later `updateSelectInput`/`updateSliderInput` from the server can
overwrite it. Set inputs **after** the step has settled, and re-assert them right before
triggering compute. The harness verifies what actually landed by diffing the captured
`params.RData` against the scenario (gate 1) — so a silently-reverted input shows up as a
test failure, not a silent wrong result.

### `filterthreshold` is a character string

`selectInput` returns `"0.2"`, `"3"`, `"-0.5"` etc. as **character**, and both
`upload_module_normalization.R:94` and `playbase::pgx.preprocess` compare it with
`f >= 1` / `f < 0` without coercing. That is string comparison in R. It happens to give
the intended answer for all five offered values (`"3" >= "1"` TRUE, `"0.2" >= "1"` FALSE,
`"-0.5" < "0"` TRUE), and both sides do the same thing, so it is not a parity risk — but
`run_script.R` must pass the **string**, not a number, to stay on the identical code path.

## 8. Step 5 — the compute panel

`upload_module_computepgx.R`, ids prefixed `upload-compute-`:

| Input | id | Notes |
|---|---|---|
| Dataset name | `selected_name` | required (`shinyvalidate::sv_required`) |
| Description | `selected_description` | required |
| Gene test methods | `gene_methods` | defaults from `GENETEST.SELECTED()` |
| Geneset methods | `gset_methods` | defaults from `GENESET.SELECTED()` |
| Extra methods | `extra_methods` | forced to include `meta.go`, `infer` |
| Run extra | `do_extra` | if FALSE, `extra.methods <- c()` |
| Filter methods | `filter_methods` | checkbox group: `append.symbol`, `proteingenes`, `remove.unknown`, `average.duplicated`, `remove.xy.probes`, `excl.immuno`, `excl.xy`, `remove.notexpressed` |
| Exclude genes | `exclude_void` / `exclude_genes` | |
| Time series | `dotimeseries` | |

Both name and description are `sv_required` — leave either blank and the compute observer
bails with a "Missing required fields" alert instead of writing `params.RData`.

Note `GENETEST.METHODS()`/`GENETEST.SELECTED()` (`upload_module_computepgx.R:47-93`)
branch on `sum(is.na(countsX()))`. So **the normalization settings change which DE methods
are selected by default**. For a clean parity comparison, pin `gene_methods` and
`gset_methods` explicitly in every scenario rather than letting them float — otherwise a
scenario that introduces NAs silently also changes the method list, and gate 3 would be
comparing two different analyses.

## 9. Capturing `params.RData`

`upload_module_computepgx.R:1239-1252`:

```r
path_to_params <- file.path(raw_dir(), "params.RData")
saveRDS(params, file = path_to_params)
script_path <- normalizePath(file.path(get_opg_root(), "bin", "pgxcreate_op.R"))
tmpdir <- normalizePath(raw_dir())
new.job <- list(process = processx::process$new(
  "Rscript", args = c(script_path, tmpdir, OPG), supervise = TRUE, ...))
```

`raw_dir` is created by `create_raw_dir()` (`upload_server.R:252-257`):

```r
raw_dir <- tempfile(pattern = prefix, tmpdir = file.path(PGX.DIR, "USER_INPUT"))
```

with `PGX.DIR <- file.path(OPG, "data")` (`components/app/R/global.R:64`). So the file
lands at:

```
data/USER_INPUT/raw_<random>/params.RData
```

### The cleanup race

`on_process_completed()` deletes the whole directory on success
(`upload_module_computepgx.R:1462-1469`):

```r
if (grepl("raw_", raw_dir)) {
  if (length(list.files(raw_dir, pattern = "ERROR_")) == 0) {
    unlink(raw_dir, recursive = TRUE)
  }
}
```

The compute takes 30-60 min, so there is a huge window — but the driver must still poll
for `data/USER_INPUT/raw_*/params.RData` and copy it out, rather than assume it persists.
Poll every ~2s from the moment compute is triggered.

Also copy the sibling artifacts while there: `CHECKS_OUTPUT`, `processx-error.log`,
`processx-output.log` — they are the first place to look when a scenario misbehaves.

### Why we stop at params.RData

The app does not compute the pgx itself; it shells out to `bin/pgxcreate_op.R` with
`params.RData`. Re-running that same script on the same file reproduces the app's
`.pgx` exactly. So once `params.RData` is captured we tear the app down, and run both
pipelines ourselves under controlled seeds and library paths. This buys:

- no 30-60 min in-browser wait per scenario,
- no cleanup race on the output,
- a seeded, comparable Path A instead of an uncontrolled background job,
- `params.RData` itself as the sharpest diff artifact (gate 1 names the divergent knob
  before any compute runs).

## 10. Teardown

Kill the app process after capture. The `processx` compute child is supervised
(`supervise = TRUE`), so it dies with the parent — which is what we want, since we
re-run it ourselves.

Between scenarios, clear `data/USER_INPUT/` so the poller cannot pick up a stale
`raw_*` directory from a previous run.
