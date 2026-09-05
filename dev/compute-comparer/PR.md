# fix: stop the upload preview silently dropping duplicated features

## The bug

`upload_module_normalization.R`'s `cleanX` realigned `counts` to `X` with:

```r
kk <- intersect(rownames(X), rownames(counts))
X      <- X[kk, , drop = FALSE]
counts <- counts[kk, , drop = FALSE]
```

`intersect()` returns **unique** names, and indexing a matrix by a duplicated name
returns only the **first** matching row. So on any dataset with duplicated features the
preview silently collapsed them. It also ran **unconditionally**, reducing the data even
when nothing needed realigning.

Measured on datasets we already ship in `opg-exampledata`:

| dataset | features | duplicated | after old `cleanX` |
|---|---|---|---|
| `salmon-ENSOKI` | 999 | 30 | **969** |
| `horse` | 31217 | 14462 | ~16755 |
| `lipidomics_underscore` | 180 | 22 | 158 |

Duplicated gene symbols are normal in RNA-seq, so this is not an edge case.

## Why it matters now

On `edgy` the compute path no longer uses the module's `X`: `upload_module_computepgx.R`
sends **raw counts + `preprocess` settings** and `pgx.createPGX()` rebuilds `X` through
`playbase::pgx.preprocess()`. `pgx.preprocess` does **not** have this bug — it only
realigns when the rownames actually differ, and uses a `%in%` subset when they are
duplicated.

So the two disagreed: the **preview** showed collapsed data while the **compute** used
the full matrix. That gap matters beyond cosmetics, because the module's `countsX` still
drives:

- default DE-method selection (`sum(is.na(countsX))`, `upload_module_computepgx.R:53-89`)
- the NA warnings in the compute UI (lines ~621, 657, 674)

## The fix

Make `cleanX` mirror `playbase::pgx.preprocess()`, which is what the compute actually
runs — realign only when needed, and preserve duplicates:

```r
if (!identical(rownames(X), rownames(counts))) {
  if (anyDuplicated(rownames(counts))) {
    counts <- counts[rownames(counts) %in% rownames(X), , drop = FALSE]
  } else {
    counts <- counts[rownames(X), , drop = FALSE]
  }
}
```

Also folded in, from the same parity audit: the outlier flag is now NA-safe,
matching `pgx.preprocess`:

```r
is.outlier <- !is.na(res$z.outlier) & (res$z.outlier > threshold)
```

Previously a non-finite z-score propagated `NA` into `if (any(is.outlier) && ...)`, which
errors with *"missing value where TRUE/FALSE needed"*.

## Verification

`dev/compute-comparer/test_cleanX_parity.R` — 12 assertions, no Shiny or browser needed:

- duplicates preserved when normalization changes nothing
- a whole dropped feature removed, `counts` still row-aligned to `X`
- BMIQ-style reorder followed correctly (unique names)
- no-op when rownames already match
- end-to-end against the real `playbase::pgx.preprocess` on a duplicated fixture
  (999 → 999 with the fix; 999 → 969 without)
- one assertion pins a **documented limit**: if a normalizer ever dropped only *part* of
  a duplicated feature, neither this code nor `pgx.preprocess` would stay aligned. Not
  reachable today (normalization is column-wise, so it keeps every row of a feature or
  none), and `createPGX` compares rownames element-wise so it would error rather than
  corrupt — but the assertion fails loudly if that ever changes.

```
R_LIBS_USER=<playbase-edgy-lib> Rscript dev/compute-comparer/test_cleanX_parity.R
```

`dev/compute-comparer/test_preprocess_sweep.R` — 45 checks, seconds to run, no browser.
Runs the app's full `imputedX -> normalizedX -> cleanX` chain (with this fix) against
`playbase::pgx.preprocess()` over 9 fixtures x 5 option sets. **45/45 pass**, covering the
datatype branches the browser harness cannot reach (the app rejects some fixtures'
organisms):

| fixture | branch exercised | result |
|---|---|---|
| `methylation` | `mToBeta` + `normalizeMethylation` | 50000x95 |
| `multiomics-pxmx` | `is.mox`, per-datatype prior | 7445x13, filter -> 5661 |
| `olink-uniprot` | `is_npx` 2^x back-transform | 90x494, outliers -> 471 |
| `metabolomics-kegg` | metabolomics norm | 429x220, outliers -> 214 |
| `human-dups` / `salmon-ENSOKI` / `lipidomics_underscore` | **duplicated features** | all preserved |

`human-dups` gives 7479 -> 7479 with the fix, where the old `cleanX` gave 7439.

```
R_LIBS_USER=<playbase-edgy-lib> Rscript dev/compute-comparer/test_preprocess_sweep.R
```

Note the sweep transcribes the module's chain rather than calling it (it is Shiny
reactive code). The browser harness is what proves the transcription is faithful; the
sweep is what makes breadth affordable.

`dev/compute-comparer/fixtures/make_human_dups.R` generates the Human/RNA-seq fixture
(`human-symbol` + 40 duplicated symbols). It exists because the app rejects
`salmon-ENSOKI`'s organism (*Oncorhynchus kisutch* is not in its namespace list), so the
shipped duplicate-bearing datasets could not drive the upload UI. The generated CSVs are
**not committed** — the repo `.gitignore` excludes `*.csv`, and a 2 MB fixture does not
belong in git. Run the script once before the sweep:

```
cd dev/compute-comparer && Rscript fixtures/make_human_dups.R
```

## Behaviour change

Datasets with duplicated feature names now keep every feature through the preview, where
before some were silently discarded. That is a **change in what the preview reports** (and
in the DE-method defaults derived from it), and it makes the preview agree with what
compute has been doing on `edgy` all along. Duplicate handling proper still happens
downstream in `createPGX` (`average.duplicated` /
`counts.mergeDuplicateFeatures`), which is where it belongs.

## Context

Found while building `dev/compute-comparer/`, a harness that drives the real upload UI
with Playwright, captures the `params.RData` the app hands to `bin/pgxcreate_op.R`, and
diffs it against a script-built one. Across 20 scenarios — every normalization method,
all three NA-filter encodings, outlier removal, imputation — `pgx.preprocess()` already
reproduces the app's `X` bit-for-bit. Duplicated feature names were the one case no
fixture covered, and the one case that diverged.

See `dev/compute-comparer/README.md` and `docs/FINDINGS.md`.

## Not included

Two further findings are written up in `docs/FINDINGS.md` but deliberately **not** patched
here, to keep this diff to one concern:

1. **The chosen batch-correction method is reset to ComBat.**
   `upload_module_normalization.R:253` calls `updateSelectInput(session, "bec_method",
   choices = methods)` with no `selected=`, so Shiny selects the first choice. The app
   recorded ComBat across five different scripted attempts to select limma/NPM. Needs a
   human to confirm severity by clicking manually before fixing.
2. **Every wizard lock/unlock branch is dead code.** `input$upload_wizard` carries the
   step *title*, not `step_id`, so `req(input$upload_wizard == "step_counts")` and friends
   never fire and the wizard never locks.
