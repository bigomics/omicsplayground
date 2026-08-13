# Reproducing the compute pipeline from a script (edgy)

Knowledge file for `run_script.R`. Describes the `params.RData` contract, the
`preprocess` option mapping, and every known place where the `master` app pipeline and
`playbase::pgx.preprocess` can diverge.

---

## 1. The entry point

`bin/pgxcreate_op.R` is the whole compute. It is what the Shiny app shells out to, and
what a script should call:

```bash
R_LIBS_USER=dev/compute-comparer/libs/pb-edgy \
  Rscript bin/pgxcreate_op.R <temp_dir> <OPG_root>
```

- `argv[1]` = `temp_dir`, a directory containing `params.RData` (an `saveRDS`'d list —
  despite the `.RData` name it is read with `readRDS`).
- `argv[2]` = `OPG`, the omicsplayground root; used only to `scan()` `VERSION` into
  `pgx$versions$omicsplayground_version`.
- Falls back to `PARAMS.yml` if `params.RData` is absent.
- Output: `save(pgx, file = <params$pgx.save.folder>/<params$name>.pgx)`, or into
  `temp_dir` if that folder does not exist.

### master vs edgy

The two versions of this file are identical **except** for one line in the
`pgx.createPGX()` call:

```r
  X = params$countsX,
  preprocess = params$preprocess,   # <- edgy only
```

That single line is the whole refactor at this level. Everything else — the
`pgx.computePGX()` call, the save logic — is the same.

## 2. `pgx.createPGX(preprocess=)`

playbase edgy, `R/pgx-compute.R:221-228`:

```r
## Opt-in: build X from raw counts via the shared preprocessing pipeline.
if (is.null(X) && !is.null(preprocess) && datatype != "scRNA-seq") {
  message("[pgx.createPGX] building X via pgx.preprocess()")
  pp <- pgx.preprocess(counts, samples = samples, contrasts = contrasts,
                       annot = annot_table, options = preprocess)
  ...
}
```

Three gates, all of which must hold or preprocessing is silently skipped:

1. `X` must be **NULL**. If you pass both `X` and `preprocess`, `X` wins and
   `preprocess` is ignored — no warning. This is the single easiest way to write a
   parity test that silently proves nothing.
2. `preprocess` must be non-NULL.
3. `datatype != "scRNA-seq"` — scRNA has its own path
   (`pgx.createSingleCellPGX`, which has no `X` parameter at all and always rebuilds
   from counts).

## 3. The `preprocess` options list

The authoritative mapping lives in the app on edgy, as the `preprocess` reactive in
`components/board.upload/R/upload_module_normalization.R:1504-1520`. Transcribed:

| Option | Source in the app | Default in `pgx.preprocess` |
|---|---|---|
| `datatype` | `upload_datatype()` | `"RNA-seq"` |
| `is_npx` | `is.olink() \|\| is.nulisa()` | `FALSE` |
| `zero_as_na` | `input$zero_as_na` | `FALSE` |
| `normalize` | `input$normalize` | `TRUE` |
| `norm_method` | `input$normalization_method` | `"CPM"` |
| `ref_gene` | `input$ref_gene` if method is `reference`, else NULL | `NULL` |
| `filter_missing` | `input$filtermissing` | `FALSE` |
| `filter_threshold` | `input$filterthreshold` | `3` |
| `impute` | `input$impute` | `FALSE` |
| `impute_method` | `input$impute_method` | `"SVD2"` |
| `remove_outliers` | `input$remove_outliers` | `FALSE` |
| `outlier_threshold` | `input$outlier_threshold` | `3` |
| `meth_type` | `meth_type()` | `NULL` |
| `outlier_methods` | *(not exposed in the UI)* | `c("z.correlation","z.distance","z.features")` |

Note the defaults inside `pgx.preprocess` are **not** the app's defaults
(`remove_outliers` threshold is 6 in the UI vs 3 here; `filter_threshold` is `"0.2"` in
the UI vs `3` here). Scenarios must set every field explicitly rather than relying on
either set of defaults.

`filter_threshold` must be passed as the **character** string the `selectInput` yields
(`"0.1"`, `"0.2"`, `"0.5"`, `"3"`, `"-0.5"`). Both implementations compare it with
`f >= 1` / `f < 0` without coercion, i.e. string comparison. It gives the intended answer
for all five offered values, and both sides do it identically — but only if the script
also passes a string.

### Batch correction is NOT in `preprocess`

It runs later, inside `pgx.createPGX`, via `batch.correct.method` and `batch.pars`
(top-level params). The `X` produced by both the app's normalization module and
`pgx.preprocess` is **pre-batch-correction**, which is what makes them comparable.

## 4. The `params` contract

Assembled at `components/board.upload/R/upload_module_computepgx.R:1162-1220`. The fields
`bin/pgxcreate_op.R` actually reads:

**createPGX**: `organism`, `counts`, `countsX`, `preprocess` (edgy), `norm_method`,
`samples`, `contrasts`, `azimuth_ref`, `dotimeseries`, `name`, `datatype`,
`datatype_subtype`, `probe_type`, `description`, `metadata`, `creator`,
`batch.correct.method`, `batch.pars`, `covariates`, `dma`, `remove.xy.probes`,
`meth_type`, `prune.samples`, `filter.genes`, `exclude.genes`, `only.known`,
`average.duplicated`, `only.proteincoding`, `only.hugo`, `convert.hugo`,
`custom.geneset`, `max.genesets`, `annot_table`, `settings`, `sc_compute_settings`.

**computePGX**: `max.genes`, `gx.methods`, `gset.methods`, `custom_fc`, `extra.methods`,
`use.design`, `prune.samples`, `do.cluster`, `cluster.contrasts`, `pgx.save.folder`,
`libx.dir`, (edgy also: `ai_features`).

**Not read by the script** but present in the app's params: `date`, `ETC`, `email`,
`sendSuccessMessageToUser` (a closure — it will be in the captured app params and must be
on the comparer's known-diff list).

### Fields a script gets wrong by default

Found by diffing a script-built `params` against one captured from a live browser
session (gate 1). All four are silent — nothing errors, the values just differ.

1. **`contrasts` are a LABEL matrix, not the raw CSV.** The upload flow runs
   `playbase::contrasts.convertToLabelMatrix(contrasts, samples)`, so `params$contrasts`
   holds `"notact"` / `"act12h"` group labels. An old-style `contrasts.csv` on disk holds
   `+1`/`-1`. Read it and pass it through unconverted and you feed a different encoding
   into `createPGX` than the app does. (`pgx.preprocess` converts internally for its
   group-wise NA filter, so this does *not* show up as an `X` difference — only a
   `params` one, and then downstream in `createPGX`.)
2. **`custom.geneset` is `list(gmt = NULL, info = NULL)`, not `NULL`** — the app
   initialises it even when no GMT is uploaded.
3. **`sc_compute_settings` is a list of 8 FALSE flags, not `NULL`** — `sc_compute_settings.PARS`
   is built unconditionally, so a bulk run still carries
   `nfeature_threshold`, `mt_threshold`, `hb_threshold`, `compute_supercells`,
   `regress_mt`, `regress_hb`, `regress_ribo`, `regress_ccs`, all FALSE.
4. **`datatype_subtype` is `"MS"` for plain proteomics**, not `NULL` — the app resolves
   proteomics to `Olink NPX` / `Nulisa NPQ` / `MS` from `is.olink()`/`is.nulisa()`. Every
   other datatype stays `NULL`.
5. **`norm_method` is `"skip_normalization"` when normalization is off**, not the
   still-selected method (`upload_module_normalization.R:1339-1343`):

   ```r
   norm_method <- reactive({
     m <- input$normalization_method
     if (!input$normalize) m <- "skip_normalization"
     m
   })
   ```

   This is the *reported/stored* value only — it goes into `params$norm_method` and
   `params$settings$norm_method`. `preprocess$norm_method` must still carry the **real**
   method name, because `pgx.preprocess` decides whether to normalize from the separate
   `normalize` flag but still reads the method for its `grepl("CPM|TMM")` log2-prior
   branch. Setting `preprocess$norm_method = "skip_normalization"` would silently change
   the prior. The same encoding pattern applies to `remove_outliers`, which reports
   `"no_outlier_removal"` instead of the threshold.

Only `ETC`, `email` and `sendSuccessMessageToUser` are legitimately app-only (the last is
a closure and cannot be reproduced or compared).

`settings` is a nested list the app fills with
`imputation_method`, `bc_method`, `remove_outliers`, `norm_method`, `custom_fc`. It is
carried into the pgx for later re-compute; it does not drive preprocessing.

## 5. Building `counts` for Path B — the trap

On `master`, `params$counts` is `countsRT()`, which resolves to the normalization
module's `cleanX()$counts` — i.e. **already NA-filtered and outlier-removed**. It is not
raw.

So Path B must **not** copy the app's counts. It must start from the raw CSV and let
`pgx.preprocess` re-derive the reduction. Passing reduced counts *plus* `preprocess`
would double-apply filtering, and `detectOutlierSamples` re-run on an
already-outlier-stripped matrix can flag a *different* sample — the step is not
idempotent.

The correct Path B input is the raw counts matrix as it exists just before the
normalization module, i.e. `checked_samples_counts()$COUNTS`, plus
`annot_table = checked_annot()$matrix`.

### Upstream steps to replicate

Between the CSV on disk and that boundary the app applies:

1. `playbase::pgx.checkINPUT(counts, "COUNTS")` — validation/repair.
2. The `e29` log-scale back-conversion (`upload_server.R:371-380`): if `checkINPUT`
   flags `e29`, the data is on a log scale and gets `2^x`'d back to linear.
3. `sum_treps` replicate summation (`upload_module_preview_samples.R:56-66`) —
   interactive, only when the user asks for it.

`opg-exampledata/check.R:128-134` already implements (1) and (2) exactly:

```r
res <- playbase::pgx.checkINPUT(counts, "COUNTS")
if ("e29" %in% names(res$checks)) {
  res$df <- 2**res$df
  if (min(res$df, na.rm = TRUE) >= 1) res$df <- res$df - 1
  counts <- res$df
}
```

Reuse that block verbatim. (3) is skipped — the fixtures have no replicates, and the
scenarios do not enable it.

## 5b. Environment: playdata must be `main`

Both playbase `main` and `edgy` call `playdata::PHOSPHOSITE` (`R/pgx-annot.R`). The
system playdata (0.1.0, branch `feat/reactome-diagrams`) does **not** ship it, so
`pgx.createPGX` dies with:

```
Error: 'PHOSPHOSITE' is not an exported object from 'namespace:playdata'
```

`playdata@main` (0.1.2) has it. The harness installs that into `libs/pb-main` and
`libs/pb-edgy` only, leaving the system playdata alone. Note this affects **both** paths
— it is an environment gap, not a version-pairing issue.

Checking for it needs an actual access, not `getNamespaceExports()`: `PHOSPHOSITE` is a
lazy-loaded dataset, so it never appears in the namespace export list even when present.

```r
tryCatch({playdata::PHOSPHOSITE; "OK"}, error = function(e) "FAIL")
```

## 6. Known implementation deltas (master app vs `pgx.preprocess`)

These were found by reading the two implementations side by side. They are the concrete
hypotheses the harness exists to confirm or refute — gates 1 and 2 will catch each of
them.

### 6a. `detectOutlierSamples` — signature differs, results should not

`master` app (`upload_module_normalization.R:181`) calls:

```r
res <- playbase::detectOutlierSamples(X, plot = FALSE)
```

In playbase **main** the signature is `detectOutlierSamples(X, plot, par)` — there is no
`methods` argument. It always computes `z.correlation`, `z.distance`, `z.features` and
takes their `rowMeans` (`playbase-main/R/pgx-outlier.R:7-44`).

In playbase **edgy** the signature is `detectOutlierSamples(X, methods = ..., plot, ...)`
where the default is `c("z.correlation","z.distance","z.features","z.isoforest")[1:3]`,
and `pgx.preprocess` additionally pins `outlier_methods` to those same three with the
comment that a NULL would switch the isolation forest on behind the caller's back.

So both sides compute the same three z-scores and the pin is deliberately parity-preserving.
**Expected to match** — but it depends on edgy's refactored implementation computing each
score identically to main's inline version, which only a real run confirms.

### 6b. NA-handling in the outlier flag

```r
# master app
is.outlier <- (res$z.outlier > threshold)
# pgx.preprocess
is.outlier <- !is.na(res$z.outlier) & (res$z.outlier > opt$outlier_threshold)
```

With any non-finite z-score, master propagates NA into `if (any(is.outlier) && ...)`,
which errors with "missing value where TRUE/FALSE needed"; edgy treats it as
not-an-outlier. Only reachable when the z-scores are degenerate, so most scenarios will
not hit it. Reachable via a scenario combining `zero_as_na` with outlier removal.

### 6c. counts/X realignment when feature names are duplicated

```r
# master app cleanX (unconditional)
kk <- intersect(rownames(X), rownames(counts))
X <- X[kk, , drop = FALSE]; counts <- counts[kk, , drop = FALSE]

# pgx.preprocess (only when they actually differ)
if (!identical(rownames(X), rownames(counts))) { ... }
```

Two differences:

- master realigns **unconditionally**; `pgx.preprocess` is a no-op when the rownames
  already match.
- `intersect()` returns **unique** names, and indexing a matrix by a duplicated name
  returns only the *first* matching row. So on a counts matrix with duplicated feature
  names, master silently collapses duplicates here, while `pgx.preprocess`'s `%in%`
  branch keeps every duplicate row.

For a dataset with duplicated rownames this is a **row-count divergence**, which gate 2
reports directly. Neither `mini-example` (917 features) nor `human-symbol` (7439
features) has duplicates, so the default scenario matrix does **not** exercise this.
Add a fixture with duplicated gene symbols to test it deliberately.

## 7. Running Path B

```r
params <- list(
  counts     = <raw counts, post checkINPUT/e29>,
  countsX    = NULL,                 # MUST be NULL or preprocess is ignored
  preprocess = <scenario options>,
  samples = ..., contrasts = ..., annot_table = ...,
  ...                                # rest per §4
)
saveRDS(params, file.path(run_dir, "params.RData"))
```

then

```bash
R_LIBS_USER=.../libs/pb-edgy \
  Rscript /home/massagno/bigomics/GitHub/omicsplayground-edgy/bin/pgxcreate_op.R \
          <run_dir> /home/massagno/bigomics/GitHub/omicsplayground-edgy
```

Path A is run the same way but against the **master** worktree and `pb-main`, on the
`params.RData` captured from the browser. Both sides are therefore executed by us, under
the same seed, rather than one of them coming from an uncontrolled background job.

### Seeding

There is no `set.seed` anywhere in `pgx-compute.R`, so `pgx.computePGX` is not
reproducible run-to-run (t-SNE, UMAP, clustering, WGCNA). To seed the top-level process
without editing the script:

```bash
Rscript -e 'set.seed(42); source("bin/pgxcreate_op.R")' --args <run_dir> <OPG>
```

`commandArgs(trailingOnly = TRUE)` picks up everything after `--args`, so the script sees
the same arguments it would normally get.

This does **not** fix draws inside forked or `BiocParallel` workers. Gate 3 must therefore
establish an A-vs-A noise floor first, and judge parity on the deterministic slots (`X`,
`counts`, `gx.meta` fold changes and p-values, `gset.meta`) with the stochastic ones
reported as informational.
