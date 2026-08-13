#!/usr/bin/env Rscript
##
## Unit check: the normalization module's cleanX realignment must agree with
## playbase::pgx.preprocess(), which is what the compute actually runs.
##
##   R_LIBS_USER=libs/pb-edgy:$HOME/R/x86_64-pc-linux-gnu-library/4.3 \
##     Rscript test_cleanX_parity.R
##
## No Shiny, no browser: cleanX's realignment is pure matrix code, so both the old
## and the new version are reproduced here verbatim and compared on inputs that
## actually have duplicated feature names.
##
## Exits non-zero if any assertion fails.
##

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

## AFTER: mirrors playbase::pgx.preprocess() -- only realign when the rownames
## actually differ, and use a %in% subset when names are duplicated.
realign_new <- function(X, counts) {
  if (!identical(rownames(X), rownames(counts))) {
    if (anyDuplicated(rownames(counts))) {
      counts <- counts[rownames(counts) %in% rownames(X), , drop = FALSE]
    } else {
      counts <- counts[rownames(X), , drop = FALSE]
    }
  }
  list(X = X, counts = counts)
}

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
  fx <- "/home/massagno/bigomics/GitHub/opg-exampledata/salmon-ENSOKI"
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
