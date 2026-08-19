#!/usr/bin/env Rscript
##
## Unit check: the normalization module's cleanX realignment must agree with
## playbase::pgx.preprocess(), which is what the compute actually runs.
##
##   R_LIBS_USER=libs/pb-edgy:$HOME/R/x86_64-pc-linux-gnu-library/4.3 \
##     Rscript test_cleanX_parity.R
##
## No Shiny, no browser: the realignment is pure matrix code, so the app's own
## realign_counts_to_X() is loaded straight out of the module file and compared
## against the old intersect() behaviour on inputs with duplicated feature names.
##
## Set OPG_EXAMPLEDATA to point at your opg-exampledata checkout.
##
## Exits non-zero if any assertion fails.
##

HERE <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
## Fixture roots are overridable so this is not tied to one developer's checkout.
EXAMPLEDATA <- Sys.getenv("OPG_EXAMPLEDATA",
  "/home/massagno/bigomics/GitHub/opg-exampledata")

ok <- TRUE
check <- function(label, cond, detail = "") {
  cat(sprintf("  [%-4s] %-46s %s\n", if (isTRUE(cond)) "PASS" else "FAIL", label, detail))
  if (!isTRUE(cond)) ok <<- FALSE
}

## ---------------------------------------------------------------- the two impls

## BEFORE: unconditional intersect(). Returns unique names, so indexing by it
## keeps only the FIRST row of each duplicated feature.
realign_old <- function(X, counts) {
  kk <- intersect(rownames(X), rownames(counts))
  list(X = X[kk, , drop = FALSE], counts = counts[kk, , drop = FALSE])
}

## AFTER: loaded FROM THE MODULE ITSELF, so this is a real regression test -- revert
## realign_counts_to_X() in the app and these assertions fail. (A local copy of the
## implementation would pass no matter what the app does.)
##
## Only that one top-level assignment is evaluated: the rest of the file is Shiny
## module code that needs a session to run.
MODULE <- file.path(HERE, "..", "..", "components", "board.upload", "R",
  "upload_module_normalization.R")
if (!file.exists(MODULE)) stop("cannot find the normalization module at ", MODULE)
mod_env <- new.env(parent = globalenv())
found <- FALSE
for (e in parse(MODULE)) {
  if (is.call(e) && length(e) >= 3 &&
    as.character(e[[1]]) %in% c("<-", "=") &&
    identical(as.character(e[[2]]), "realign_counts_to_X")) {
    eval(e, mod_env)
    found <- TRUE
  }
}
if (!found) {
  stop("realign_counts_to_X() not found in ", MODULE,
    " -- was the realignment inlined back into cleanX? This test must exercise the ",
    "app's own implementation, not a copy.")
}
realign_new <- function(X, counts) {
  list(X = X, counts = mod_env$realign_counts_to_X(X, counts))
}
cat("loaded realign_counts_to_X() from the app module\n")

mk <- function(names, ncol = 3) {
  matrix(seq_len(length(names) * ncol), nrow = length(names),
    dimnames = list(names, paste0("s", seq_len(ncol))))
}

cat("== duplicated feature names, nothing dropped by normalization\n")
n <- c("A", "B", "B", "C", "C", "C")
X <- mk(n); counts <- mk(n)
o <- realign_old(X, counts); w <- realign_new(X, counts)
check("old collapses duplicates", nrow(o$X) == 3, sprintf("%d -> %d rows", nrow(X), nrow(o$X)))
check("new preserves every row", nrow(w$X) == 6 && nrow(w$counts) == 6,
  sprintf("%d -> %d rows", nrow(X), nrow(w$X)))
check("new keeps X untouched", identical(w$X, X))
check("new keeps counts aligned to X", identical(rownames(w$counts), rownames(w$X)))

cat("\n== duplicated names, normalization dropped a whole feature\n")
## Realistic drop: every row of "D" goes, the duplicated ones are untouched.
X2 <- mk(c("A", "B", "B", "C", "C", "C")); counts2 <- mk(c("A", "B", "B", "C", "C", "C", "D"))
w2 <- realign_new(X2, counts2)
check("dropped feature removed", !("D" %in% rownames(w2$counts)),
  paste(rownames(w2$counts), collapse = ","))
check("counts stays row-aligned to X", identical(rownames(w2$counts), rownames(w2$X)))

cat("\n== KNOWN LIMIT: partially dropping a duplicated feature\n")
## If a normalizer kept only SOME rows of a duplicated name, the %in% subset would
## keep all of them and counts would no longer align to X. Neither this code nor
## playbase::pgx.preprocess() handles that -- and it is not reachable today, because
## normalization is column-wise: it either keeps every row or drops a feature
## entirely. Asserted here so the assumption is visible and fails loudly if a future
## normalizer breaks it. createPGX also compares rownames element-wise, so such a
## case would error rather than silently corrupt.
X3 <- mk(c("A", "B", "B", "C")); counts3 <- mk(c("A", "B", "B", "C", "C", "C"))
w3p <- realign_new(X3, counts3)
check("documented limit still holds (counts NOT aligned here)",
  nrow(w3p$counts) != nrow(w3p$X),
  sprintf("X=%d counts=%d -- partial-duplicate drops are unsupported by design",
    nrow(w3p$X), nrow(w3p$counts)))

cat("\n== unique names, normalization reordered (methylation BMIQ case)\n")
X3 <- mk(c("C", "A", "B")); counts3 <- mk(c("A", "B", "C"))
w3 <- realign_new(X3, counts3)
check("counts follows X's order", identical(rownames(w3$counts), rownames(X3)),
  paste(rownames(w3$counts), collapse = ","))

cat("\n== no-op when rownames already match\n")
X4 <- mk(c("A", "B", "C")); counts4 <- mk(c("A", "B", "C"))
w4 <- realign_new(X4, counts4)
check("counts returned unchanged", identical(w4$counts, counts4))

## ------------------------------------------------- against the real pgx.preprocess
if (requireNamespace("playbase", quietly = TRUE) &&
  "pgx.preprocess" %in% getNamespaceExports("playbase")) {
  cat("\n== end-to-end against playbase::pgx.preprocess on a real duplicated fixture\n")
  fx <- file.path(EXAMPLEDATA, "salmon-ENSOKI")
  if (dir.exists(fx)) {
    counts <- playbase::read_counts(file.path(fx, "counts.csv"))
    samples <- playbase::read_samples(file.path(fx, "samples.csv"))
    contrasts <- playbase::read_contrasts(file.path(fx, "contrasts.csv"))
    contrasts <- playbase::contrasts.convertToLabelMatrix(contrasts, samples)
    ndup <- sum(duplicated(rownames(counts)))
    pp <- playbase::pgx.preprocess(counts, samples, contrasts, NULL,
      list(datatype = "RNA-seq", normalize = TRUE, norm_method = "CPM", impute = TRUE))
    check("fixture really has duplicates", ndup > 0, sprintf("%d duplicated names", ndup))
    check("pgx.preprocess keeps all features", nrow(pp$X) == nrow(counts),
      sprintf("%d -> %d", nrow(counts), nrow(pp$X)))
    check("old cleanX would have dropped them",
      nrow(realign_old(pp$X, pp$counts)$X) < nrow(pp$X),
      sprintf("%d -> %d", nrow(pp$X), nrow(realign_old(pp$X, pp$counts)$X)))
    check("new cleanX agrees with pgx.preprocess",
      identical(realign_new(pp$X, pp$counts)$counts, pp$counts))
  } else {
    cat("  [SKIP] fixture not found:", fx, "\n")
  }
} else {
  cat("\n  [SKIP] playbase with pgx.preprocess not on .libPaths()\n")
}

cat(if (ok) "\nALL CHECKS PASSED\n" else "\nFAILURES\n")
quit(status = if (ok) 0 else 1)
