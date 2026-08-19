#!/usr/bin/env Rscript
##
## Broad R-level sweep: run playbase::pgx.preprocess() against the app's
## imputedX -> normalizedX -> cleanX chain (transcribed verbatim, with the cleanX
## fix applied) across many fixtures x option sets.
##
##   R_LIBS_USER=libs/pb-edgy:$HOME/R/... Rscript test_preprocess_sweep.R
##
## Why this exists alongside the browser harness: driving the real UI costs ~10 min
## per scenario and cannot reach datatypes whose organism the app rejects. This runs
## the same comparison in seconds, so it can cover the branches the browser matrix
## never got to (methylomics, multi-omics, metabolomics, duplicated features).
##
## It is NOT a replacement for the browser runs -- the transcription below could
## drift from the real module. The browser harness is what proves the transcription
## is faithful; this is what gives breadth cheaply.
##

suppressMessages(library(playbase))

ok <- TRUE
rows <- list()
note <- function(fixture, opts, status, detail = "") {
  cat(sprintf("  [%-4s] %-22s %-34s %s\n", status, fixture, opts, detail))
  if (status == "FAIL") ok <<- FALSE
  rows[[length(rows) + 1]] <<- data.frame(fixture, opts, status, detail,
    stringsAsFactors = FALSE)
}

## ------------------------------------------- the app's chain, cleanX fix applied
## Transcribed from components/board.upload/R/upload_module_normalization.R
## (imputedX / normalizedX / cleanX), with input$* replaced by opt$*.
app_chain <- function(counts, samples, contrasts, annot, opt) {
  counts[which(is.nan(counts))] <- NA
  counts[which(is.infinite(counts))] <- NA
  if (isTRUE(opt$is_npx)) counts <- 2**counts
  if (any(counts < 0, na.rm = TRUE)) counts <- pmax(counts, 0)
  if (isTRUE(opt$zero_as_na)) counts[which(counts == 0)] <- NA

  is.mox <- playbase::is.multiomics(rownames(counts))
  if (is.mox) {
    X <- counts
    for (dt in unique(sub(":.*", "", rownames(X)))) {
      ii <- grep(paste0("^", dt, ":"), rownames(counts))
      prior <- if (dt != "gx") playbase::getPrior(counts[ii, ]) else 1
      X[ii, ] <- log2(counts[ii, ] + prior)
    }
  } else if (opt$datatype == "methylomics") {
    X <- playbase::mToBeta(counts); prior <- 0
  } else {
    prior0 <- playbase::getPrior(counts)
    prior <- ifelse(grepl("CPM|TMM", opt$norm_method), 1, prior0)
    X <- log2(counts + prior)
  }

  if (sum(is.na(X)) > 0 && isTRUE(opt$filter_missing)) {
    f <- opt$filter_threshold
    sc <- playbase::contrasts.convertToLabelMatrix(contrasts, samples)
    grp <- apply(sc, 1, paste, collapse = "_")
    if (f >= 1) {
      gs <- tapply(seq_len(ncol(counts)), grp, function(i) rowSums(!is.na(counts[, i, drop = FALSE])))
      sel <- apply(do.call(cbind, gs), 1, max, na.rm = TRUE) >= 3
    } else if (f < 0) {
      ga <- tapply(seq_len(ncol(counts)), grp, function(i) rowMeans(!is.na(counts[, i, drop = FALSE])))
      sel <- apply(do.call(cbind, ga), 1, max, na.rm = TRUE) >= 0.5
    } else {
      sel <- rowMeans(is.na(X)) <= f
    }
    X <- X[which(sel), , drop = FALSE]
    counts <- counts[which(sel), , drop = FALSE]
    if (!is.null(annot)) annot <- annot[which(sel), , drop = FALSE]
  }

  if (any(is.na(X)) && isTRUE(opt$impute)) {
    X <- if (is.mox) playbase::imputeMissing.mox(X, method = opt$impute_method) else
      playbase::imputeMissing(X, method = opt$impute_method)
  }

  if (isTRUE(opt$normalize)) {
    if (opt$datatype == "multi-omics") {
      X <- playbase::normalizeMultiOmics(X)
    } else if (opt$datatype == "methylomics") {
      nX <- try(playbase::normalizeMethylation(X, opt$norm_method, opt$meth_type), silent = TRUE)
      ## matches the module (which now guards try-error too)
      if (!inherits(nX, "try-error") && !is.null(nX)) X <- nX
    } else {
      X <- playbase::normalizeExpression(X, method = opt$norm_method,
        ref = opt$ref_gene, prior = prior)
    }
  }

  ## cleanX -- WITH THE FIX (see PR.md). The old code was:
  ##   kk <- intersect(rownames(X), rownames(counts)); X <- X[kk,]; counts <- counts[kk,]
  if (!identical(rownames(X), rownames(counts))) {
    if (anyDuplicated(rownames(counts))) {
      counts <- counts[rownames(counts) %in% rownames(X), , drop = FALSE]
    } else {
      counts <- counts[rownames(X), , drop = FALSE]
    }
  }
  if (isTRUE(opt$remove_outliers)) {
    if (sum(is.na(X)) > 0) {
      X <- if (is.mox) playbase::imputeMissing.mox(X, method = "SVD2") else
        playbase::imputeMissing(X, method = "SVD2")
    }
    res <- playbase::detectOutlierSamples(X, plot = FALSE)
    is.outlier <- !is.na(res$z.outlier) & (res$z.outlier > opt$outlier_threshold)
    if (any(is.outlier) && !all(is.outlier)) {
      X <- X[, which(!is.outlier), drop = FALSE]
      counts <- counts[, colnames(X), drop = FALSE]
    }
  }
  list(counts = counts, X = X)
}

## ------------------------------------------------------------------- fixtures
E <- Sys.getenv("OPG_EXAMPLEDATA", "/home/massagno/bigomics/GitHub/opg-exampledata")
H <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
FIX <- list(
  list(id = "human-symbol",  dir = file.path(E, "human-symbol"),          dt = "RNA-seq"),
  list(id = "human-dups",    dir = file.path(H, "fixtures/human-dups"),   dt = "RNA-seq"),
  list(id = "salmon-dups",   dir = file.path(E, "salmon-ENSOKI"),         dt = "RNA-seq"),
  list(id = "mini-example",  dir = file.path(E, "mini-example"),          dt = "proteomics"),
  list(id = "lipid-dups",    dir = file.path(E, "lipidomics_underscore"), dt = "metabolomics"),
  list(id = "metab-kegg",    dir = file.path(E, "metabolomics-kegg"),     dt = "metabolomics"),
  list(id = "olink",         dir = file.path(E, "olink-uniprot"),         dt = "proteomics", npx = TRUE),
  list(id = "methylation",   dir = file.path(E, "methylation"),           dt = "methylomics"),
  list(id = "multiomics",    dir = file.path(E, "multiomics-pxmx"),       dt = "multi-omics")
)

base_opts <- function(dt, npx = FALSE) list(
  datatype = dt, is_npx = npx, zero_as_na = FALSE, normalize = TRUE,
  norm_method = if (dt == "RNA-seq") "CPM" else if (dt == "methylomics") "BMIQ" else
    if (dt == "multi-omics") "median" else "maxMedian",
  ref_gene = NULL, filter_missing = FALSE, filter_threshold = "0.2",
  impute = TRUE, impute_method = "SVD2", remove_outliers = FALSE,
  outlier_threshold = 6, meth_type = NULL,
  outlier_methods = c("z.correlation", "z.distance", "z.features")
)

VARIANTS <- list(
  list(lab = "default",            mod = list()),
  list(lab = "normalize=off",      mod = list(normalize = FALSE)),
  list(lab = "zero_as_na",         mod = list(zero_as_na = TRUE, impute = FALSE)),
  list(lab = "filter_missing=0.2", mod = list(zero_as_na = TRUE, impute = FALSE,
                                              filter_missing = TRUE, filter_threshold = "0.2")),
  list(lab = "outliers@3",         mod = list(remove_outliers = TRUE, outlier_threshold = 3))
)

for (fx in FIX) {
  if (!dir.exists(fx$dir)) { note(fx$id, "-", "SKIP", "fixture not present"); next }
  cnt <- tryCatch(playbase::read_counts(file.path(fx$dir, "counts.csv")), error = function(e) NULL)
  smp <- tryCatch(playbase::read_samples(file.path(fx$dir, "samples.csv")), error = function(e) NULL)
  ctr <- tryCatch(playbase::read_contrasts(file.path(fx$dir, "contrasts.csv")), error = function(e) NULL)
  if (is.null(cnt) || is.null(smp) || is.null(ctr)) { note(fx$id, "-", "SKIP", "unreadable"); next }
  ctr <- tryCatch(playbase::contrasts.convertToLabelMatrix(ctr, smp), error = function(e) ctr)
  ndup <- sum(duplicated(rownames(cnt)))

  for (v in VARIANTS) {
    opt <- utils::modifyList(base_opts(fx$dt, isTRUE(fx$npx)), v$mod)
    a <- tryCatch(app_chain(as.matrix(cnt), smp, ctr, NULL, opt), error = function(e) e)
    b <- tryCatch(playbase::pgx.preprocess(cnt, smp, ctr, NULL, opt), error = function(e) e)
    lab <- sprintf("%s%s", v$lab, if (ndup > 0) sprintf(" [%d dup]", ndup) else "")
    if (inherits(a, "error") || inherits(b, "error")) {
      msg <- conditionMessage(if (inherits(a, "error")) a else b)
      note(fx$id, lab, if (inherits(a, "error") && inherits(b, "error")) "SKIP" else "FAIL",
        substr(msg, 1, 70))
      next
    }
    same_dim <- identical(dim(a$X), dim(b$X))
    same_val <- isTRUE(all.equal(unname(a$X), unname(b$X), tolerance = 1e-8))
    note(fx$id, lab, if (same_dim && same_val) "PASS" else "FAIL",
      sprintf("%s vs %s%s", paste(dim(a$X), collapse = "x"), paste(dim(b$X), collapse = "x"),
        if (same_dim && !same_val) "  VALUES DIFFER" else ""))
  }
}

df <- do.call(rbind, rows)
cat(sprintf("\n%d checks: %d PASS, %d FAIL, %d SKIP\n", nrow(df),
  sum(df$status == "PASS"), sum(df$status == "FAIL"), sum(df$status == "SKIP")))
quit(status = if (ok) 0 else 1)
