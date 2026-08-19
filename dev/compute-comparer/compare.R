#!/usr/bin/env Rscript
##
## The three parity gates. See docs/script-pipeline.md for what each one proves.
##
##   R_LIBS_USER=libs/pb-edgy Rscript compare.R <scenario-id> 1   # params + preprocess seam
##   Rscript compare.R <scenario-id> 2                            # createPGX matrices
##   Rscript compare.R <scenario-id> 3                            # full .pgx
##
## Gate 1 MUST run under pb-edgy: it calls pgx.preprocess(), which playbase main
## does not have. Gates 2 and 3 only read saved objects, so either lib works.
##
## Appends one row per gate to runs/report.csv.
##

suppressMessages(library(playbase))

HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) < 2) stop("usage: compare.R <scenario-id> <gate 1|2|3> [tolerance]")
scenario_id <- argv[1]
gate <- as.integer(argv[2])
TOL <- if (length(argv) >= 3) as.numeric(argv[3]) else 1e-8

A <- file.path(HERE, "runs", scenario_id, "A")
B <- file.path(HERE, "runs", scenario_id, "B")

## Missing inputs are a harness problem, not a parity result. Say which side is
## missing and exit 2, so a batch run can tell "not captured yet" from "diverged".
need <- function(...) {
  for (f in c(...)) {
    if (!file.exists(f)) {
      cat(sprintf("== gate %d %s: MISSING INPUT %s\n", gate, scenario_id, f))
      cat("   (run `make a SC=", scenario_id, "` / `make b SC=", scenario_id,
        "` first)\n", sep = "")
      quit(status = 2)
    }
  }
}
if (gate == 1) need(file.path(A, "params.RData"), file.path(B, "params.RData"))
if (gate == 2) need(file.path(A, "createpgx.rds"), file.path(B, "createpgx.rds"))

results <- list()
record <- function(check, status, detail = "") {
  results[[length(results) + 1]] <<- data.frame(
    scenario = scenario_id, gate = gate, check = check,
    status = status, detail = substr(detail, 1, 400), stringsAsFactors = FALSE
  )
  cat(sprintf("  [%-5s] %-28s %s\n", status, check, substr(detail, 1, 160)))
}

## Compare two matrices: shape first (a dim/dimname mismatch is the interesting
## failure -- it means filtering or outlier removal diverged), then values.
cmp_matrix <- function(name, x, y) {
  if (is.null(x) && is.null(y)) return(record(name, "SKIP", "both NULL"))
  if (is.null(x) || is.null(y)) return(record(name, "FAIL", "one side is NULL"))
  if (!identical(dim(x), dim(y))) {
    return(record(name, "FAIL", sprintf("dim %s vs %s",
      paste(dim(x), collapse = "x"), paste(dim(y), collapse = "x"))))
  }
  if (!identical(rownames(x), rownames(y))) {
    d <- setdiff(rownames(x), rownames(y))
    return(record(name, "FAIL", sprintf("rownames differ (%d only in A, e.g. %s)",
      length(d), paste(utils::head(d, 3), collapse = ","))))
  }
  if (!identical(colnames(x), colnames(y))) {
    return(record(name, "FAIL", sprintf("colnames differ: %s vs %s",
      paste(utils::head(colnames(x), 4), collapse = ","),
      paste(utils::head(colnames(y), 4), collapse = ","))))
  }
  eq <- all.equal(unname(as.matrix(x)), unname(as.matrix(y)), tolerance = TOL)
  if (isTRUE(eq)) {
    record(name, "PASS", sprintf("%s identical", paste(dim(x), collapse = "x")))
  } else {
    md <- suppressWarnings(max(abs(as.matrix(x) - as.matrix(y)), na.rm = TRUE))
    record(name, "FAIL", sprintf("%s; max|diff|=%.3g", paste(eq, collapse = "; "), md))
  }
}

## Run pgx.preprocess under a fixed seed, so stochastic imputers are reproducible.
withr_seed <- function(seed, pb) {
  set.seed(seed)
  playbase::pgx.preprocess(
    counts = pb$counts, samples = pb$samples, contrasts = pb$contrasts,
    annot = pb$annot_table, options = pb$preprocess
  )$X
}

## Strip per-run metadata that can never match, wherever it is nested:
##   timings -- wall-clock durations
##   myid    -- a UUID regenerated per object (graph/visNetwork environments), e.g.
##              5973916f-e6ab-... vs 8aed1a3a-7bad-...
## Environments are converted to lists so both sides normalise identically (an
## environment otherwise compares by reference and always differs).
## Environments are REPLACED by a marker, not expanded: pgx$meta.go$graph holds
## graph objects that reference themselves, so converting them to lists makes
## all.equal walk a cycle forever ("Component 10: Component 1: Component 10: ...").
## They are object handles rather than data -- the data they wrap is compared
## through the surrounding slots -- so collapsing both sides to the same marker is
## the honest treatment.
META_KEYS <- c("timings", "myid")
scrub <- function(x, depth = 0) {
  if (depth > 12) return("<max-depth>")
  if (is.environment(x)) return("<environment>")
  if (is.list(x)) {
    nm <- names(x)
    if (!is.null(nm)) x <- x[!nm %in% META_KEYS]
    return(lapply(x, scrub, depth = depth + 1))
  }
  x
}

cmp_obj <- function(name, x, y) {
  eq <- all.equal(x, y, tolerance = TOL)
  if (isTRUE(eq)) record(name, "PASS", "") else record(name, "FAIL", paste(eq, collapse = "; "))
}

## ===================================================================== gate 1
if (gate == 1) {
  cat("== gate 1:", scenario_id, "== params contract + preprocess seam\n")
  pa <- readRDS(file.path(A, "params.RData"))
  pb <- readRDS(file.path(B, "params.RData"))

  ## Which app are we comparing against?
  ##   master: sends a precomputed countsX, and counts/annot_table are the REDUCED
  ##           cleanX outputs -- so those legitimately differ from a script's raw ones.
  ##   edgy:   sends RAW counts + `preprocess` and countsX = NULL, i.e. exactly what a
  ##           script sends. counts/annot_table/preprocess then become MUST-MATCH, which
  ##           makes this a direct test of run_script.R's params contract.
  app_is_edgy <- is.null(pa$countsX) && !is.null(pa$preprocess)
  cat("   app mode:", if (app_is_edgy) "edgy (raw counts + preprocess)" else
    "master (precomputed countsX)", "\n")

  ## Session-bound or per-run metadata. Not data, and not reproducible by a script:
  ## sendSuccessMessageToUser and ai_features$reports$credentials are CLOSURES built
  ## from the live Shiny session (get_ai_credentials(session)), so they can never
  ## compare equal. ai_features only drives optional AI report generation, not the
  ## computed matrices -- a script leaves it NULL unless it wants those reports.
  run_meta <- c("name", "description", "creator", "date",
    "pgx.save.folder", "libx.dir", "ETC", "email", "sendSuccessMessageToUser",
    "ai_features")
  expected_diff <- if (app_is_edgy) {
    ## Only run metadata may differ. Everything else, including the preprocess
    ## options and the raw counts, has to be identical.
    run_meta
  } else {
    c("countsX",     # A only: the app's precomputed X. B sends NULL by design.
      "preprocess",  # B only: master has no such field.
      "counts",      # A is cleanX()$counts (reduced); B is raw. See sec.5 of the docs.
      "annot_table", # A is subset alongside the NA filter; B is raw.
      "settings",    # encoded differently; compared field-wise below
      run_meta)
  }

  keys <- union(names(pa), names(pb))
  must_match <- setdiff(keys, expected_diff)
  mismatched <- character(0)
  for (k in must_match) {
    if (!isTRUE(all.equal(pa[[k]], pb[[k]], tolerance = TOL))) mismatched <- c(mismatched, k)
  }
  if (length(mismatched) == 0) {
    record("params.shared_fields", "PASS", sprintf("%d fields match", length(must_match)))
  } else {
    for (k in mismatched) {
      record(paste0("params$", k), "FAIL", sprintf("A=%s | B=%s",
        paste(utils::head(format(pa[[k]]), 4), collapse = ","),
        paste(utils::head(format(pb[[k]]), 4), collapse = ",")))
    }
  }
  a_only <- setdiff(names(pa), names(pb)); b_only <- setdiff(names(pb), names(pa))
  record("params.field_inventory",
    if (all(c(a_only, b_only) %in% expected_diff)) "PASS" else "WARN",
    sprintf("A-only: %s | B-only: %s", paste(a_only, collapse = ","), paste(b_only, collapse = ",")))

  ## --- the sharp test -------------------------------------------------
  ## Does the edgy pgx.preprocess reproduce the master app's normalization
  ## module byte for byte, from the same raw counts and the same settings?
  ## This isolates the seam completely: no createPGX, no compute, no RNG.
  if (app_is_edgy) {
    ## The edgy app hands createPGX the same raw counts + options a script does, so
    ## there is no second X to compare -- one implementation, not two. What is worth
    ## asserting is that the options themselves match field by field.
    pa_pp <- pa$preprocess; pb_pp <- pb$preprocess
    keys <- union(names(pa_pp), names(pb_pp))
    bad <- Filter(function(k) !isTRUE(all.equal(pa_pp[[k]], pb_pp[[k]])), keys)
    if (length(bad) == 0) {
      record("preprocess.options_match", "PASS",
        sprintf("%d options identical (app and script feed createPGX the same call)",
          length(keys)))
    } else {
      for (k in bad) {
        record(paste0("preprocess$", k), "FAIL",
          sprintf("A=%s | B=%s", paste(format(pa_pp[[k]]), collapse = ","),
            paste(format(pb_pp[[k]]), collapse = ",")))
      }
    }
  } else if (!"pgx.preprocess" %in% getNamespaceExports("playbase")) {
    record("preprocess.equivalence", "SKIP",
      "playbase has no pgx.preprocess -- rerun gate 1 under R_LIBS_USER=libs/pb-edgy")
  } else {
    pp <- playbase::pgx.preprocess(
      counts = pb$counts, samples = pb$samples, contrasts = pb$contrasts,
      annot = pb$annot_table, options = pb$preprocess
    )

    ## QRILC / MinProb draw from a truncated normal, so an unseeded call differs run
    ## to run. That is NOT inherent non-determinism though -- measured: set.seed(42)
    ## twice gives identical output, and 42 vs 99 differs. So rather than tolerating
    ## a fuzzy match, assert the real property: under a fixed seed the result is
    ## reproducible. The app's countsX came from an unseeded session and cannot be
    ## reproduced after the fact, so that one comparison stays informational.
    stochastic_impute <- c("QRILC", "MinProb")
    is_stochastic <- isTRUE(pb$preprocess$impute) &&
      pb$preprocess$impute_method %in% stochastic_impute

    if (is_stochastic) {
      pp_a <- withr_seed(42, pb)
      pp_b <- withr_seed(42, pb)
      pp_c <- withr_seed(99, pb)
      record("preprocess.seeded_reproducible",
        if (identical(pp_a, pp_b) && !identical(pp_a, pp_c)) "PASS" else "FAIL",
        sprintf("%s: same seed identical=%s, different seed differs=%s",
          pb$preprocess$impute_method, identical(pp_a, pp_b), !identical(pp_a, pp_c)))
      obs <- mean(abs(pa$countsX - pp$X)) / mean(abs(pp$X))
      record("preprocess.X_vs_app_countsX", "INFO",
        sprintf("%s unseeded in the app session; rel.diff %.4g. Seed it to compare exactly.",
          pb$preprocess$impute_method, obs))
    } else {
      cmp_matrix("preprocess.X_vs_app_countsX", pa$countsX, pp$X)
    }
    cmp_matrix("preprocess.counts_vs_app", pa$counts, pp$counts)
  }
}

## ===================================================================== gate 2
if (gate == 2) {
  cat("== gate 2:", scenario_id, "== pgx.createPGX output\n")
  ga <- readRDS(file.path(A, "createpgx.rds"))
  gb <- readRDS(file.path(B, "createpgx.rds"))
  cmp_matrix("createPGX.counts", ga$counts, gb$counts)
  cmp_matrix("createPGX.X", ga$X, gb$X)
  cmp_obj("createPGX.samples", ga$samples, gb$samples)
  cmp_obj("createPGX.contrasts", ga$contrasts, gb$contrasts)
  cmp_obj("createPGX.genes", ga$genes, gb$genes)
}

## ===================================================================== gate 3
if (gate == 3) {
  cat("== gate 3:", scenario_id, "== full .pgx\n")
  find_pgx <- function(d) {
    f <- list.files(d, pattern = "\\.pgx$", full.names = TRUE)
    if (length(f) == 0) stop("no .pgx in ", d)
    f[1]
  }
  ga <- playbase::pgx.load(find_pgx(A))
  gb <- playbase::pgx.load(find_pgx(B))

  ## Run-to-run noise, not divergence: pgx.computePGX seeds nothing and forked
  ## workers ignore the top-level seed. Reported, never failed.
  stochastic <- c("tsne2d", "tsne3d", "umap2d", "umap3d", "cluster", "wgcna",
    "cluster.genes", "cluster.gsets", "mofa", "deconv", "drugs", "connectivity")
  ## Metadata that legitimately differs between two runs.
  ignored <- c("filename", "date", "version", "versions", "timings", "name",
    "description", "creator", "this.date")

  keys <- intersect(names(ga), names(gb))
  record("pgx.slot_inventory",
    if (setequal(names(ga), names(gb))) "PASS" else "WARN",
    sprintf("A-only: %s | B-only: %s",
      paste(setdiff(names(ga), names(gb)), collapse = ","),
      paste(setdiff(names(gb), names(ga)), collapse = ",")))

  for (k in setdiff(keys, ignored)) {
    eq <- all.equal(scrub(ga[[k]]), scrub(gb[[k]]), tolerance = TOL)
    if (isTRUE(eq)) {
      record(paste0("pgx$", k), "PASS", "")
    } else if (k %in% stochastic) {
      record(paste0("pgx$", k), "INFO", paste(utils::head(eq, 2), collapse = "; "))
    } else {
      record(paste0("pgx$", k), "FAIL", paste(utils::head(eq, 3), collapse = "; "))
    }
  }
}

## ==================================================================== report
df <- do.call(rbind, results)
rep_file <- file.path(HERE, "runs", "report.csv")
write.table(df, rep_file, sep = ",", row.names = FALSE,
  col.names = !file.exists(rep_file), append = file.exists(rep_file), qmethod = "double")

n_fail <- sum(df$status == "FAIL")
cat(sprintf("\n== gate %d %s: %d pass, %d fail, %d info/warn/skip -> %s\n",
  gate, scenario_id, sum(df$status == "PASS"), n_fail,
  sum(!df$status %in% c("PASS", "FAIL")), rep_file))
quit(status = if (n_fail > 0) 1 else 0)
