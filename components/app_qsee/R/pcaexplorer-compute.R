##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Computation for the "PCA explorer" Qsee board (formerly the standalone
## app_pcaexplorer app, now merged in as an extra board -- see
## app_qsee/shiny/qsee_pcaexplorer_ui.R / _server.R and
## app_qsee/R/pcaexplorer-plots.R).
##

#' Drop degenerate phenotype columns (constant, or all-unique like a
#' sample-id column).
#'
#' Qsee's own rY() does not do this filtering by itself. Without it,
#' qsee_pcaexplorer_plot_ftest_vs_pcdim() errors on an all-unique column
#' ("length of 'dimnames' [2] not equal to array extent"), so this
#' filtering is required for correctness.
qsee_pcaexplorer_filter_pheno <- function(Y) {
  ngrp <- apply(Y, 2, function(y) length(table(y[!is.na(y)])))
  sel.valid <- ngrp > 1 & ngrp < colSums(!is.na(Y))
  Y[, sel.valid, drop = FALSE]
}

#' PCA computation for the PCA explorer board. `X` is expected in log2
#' space (matching `playbase::normalizeExpression()`'s contract) and `Y`
#' should already be filtered with [qsee_pcaexplorer_filter_pheno()].
qsee_pcaexplorer_compute <- function(X, Y) {
  normX <- playbase::normalizeExpression(X, method = "CPM+quantile",
    ref = NULL, prior = 0)

  cX <- normX - rowMeans(normX, na.rm = TRUE) ## important
  sel <- which(rowMeans(is.na(cX)) == 0) ## only complete rows
  nv <- min(10, dim(cX) - 1)
  res <- irlba::irlba(cX[sel, ], nv = nv)
  rownames(res$u) <- rownames(cX)[sel]
  rownames(res$v) <- colnames(cX)
  colnames(res$u) <- paste0("PC", 1:ncol(res$u))
  colnames(res$v) <- paste0("PC", 1:ncol(res$v))

  H <- playbase::expandPhenoMatrix(Y, drop.ref = FALSE)
  res$R <- stats::cor(H, res$v)

  res$X <- normX
  res$Y <- Y
  res$H <- H

  res
}

#' Range (min/max/step/value) to stretch the "Min. maxFC" slider to, once
#' PCA is ready. Returns NULL when there isn't enough data to compute it
#' (slider is then left at its static default).
qsee_pcaexplorer_minfc_range <- function(res) {
  normX <- res$X
  H <- res$H
  if (is.null(normX) || is.null(H)) {
    return(NULL)
  }
  genes <- intersect(rownames(res$u), rownames(normX))
  if (length(genes) < 2L) {
    return(NULL)
  }
  max_fc <- qsee_pcaexplorer_feature_max_fc(normX[genes, , drop = FALSE], H)
  max_fc <- max_fc[is.finite(max_fc)]
  if (!length(max_fc)) {
    return(NULL)
  }
  hi <- as.numeric(stats::quantile(max_fc, 0.99, na.rm = TRUE))
  if (!is.finite(hi) || hi <= 0) hi <- max(max_fc, na.rm = TRUE)
  if (!is.finite(hi) || hi <= 0) hi <- 1
  hi <- round(hi, 2)
  step <- max(0.01, round(hi / 40, 2))
  list(min = 0, max = hi, step = step, value = 0)
}

## ---------------------------------------------------------------------------
## Feature-PCA compute helpers (used by both the slider range above and the
## feature-PCA plots in pcaexplorer-plots.R)
## ---------------------------------------------------------------------------

#' Handles NA. X: vars x samples (may contain NAs). Y: groups x samples
#' (binary, no NAs assumed). Returns a vars x groups matrix of
#' group1-minus-group0 means.
qsee_pcaexplorer_fast_diff_in_means <- function(X, Y) {
  valid_mask <- !is.na(X)

  X_zeroed <- X
  X_zeroed[!valid_mask] <- 0

  sum_1 <- X_zeroed %*% t(Y)
  sum_total <- matrix(rowSums(X_zeroed), nrow = nrow(X), ncol = nrow(Y))
  sum_0 <- sum_total - sum_1

  n1 <- valid_mask %*% t(Y)
  n_total <- matrix(rowSums(valid_mask), nrow = nrow(X), ncol = nrow(Y))
  n0 <- n_total - n1

  mean_1 <- sum_1 / n1
  mean_0 <- sum_0 / n0
  mean_1[n1 == 0] <- NA
  mean_0[n0 == 0] <- NA

  mean_1 - mean_0
}

## Per-feature max |fold-change| vs expanded phenotype (same metric used to
## color feature PCA points). Returns a named numeric vector.
qsee_pcaexplorer_feature_max_fc <- function(normX, H) {
  if (is.null(normX) || is.null(H) || !nrow(normX) || !ncol(H)) {
    return(numeric(0))
  }
  rx <- qsee_pcaexplorer_fast_diff_in_means(normX, t(H))
  apply(abs(rx), 1, function(z) {
    z <- z[is.finite(z)]
    if (!length(z)) return(NA_real_)
    max(z, na.rm = TRUE)
  })
}

## Keep features with maxFC >= min_fc (drops the middle / weak-effect cloud).
qsee_pcaexplorer_filter_features_by_maxfc <- function(genes, max_fc, min_fc = 0) {
  min_fc <- suppressWarnings(as.numeric(min_fc))
  if (!is.finite(min_fc) || min_fc <= 0) {
    return(genes)
  }
  mf <- max_fc[genes]
  keep <- genes[is.finite(mf) & mf >= min_fc]
  keep
}

## Condition / phenotype direction vectors on the feature PC plane (2D or 3D).
## colMeans(pos * pmax(cor,0)^2), then unit-normalized (length is applied
## later via var.scale * score_span in the plotting function).
qsee_pcaexplorer_feature_condition_dirs <- function(pos, rx, numarrows = Inf) {
  if (is.null(pos) || is.null(rx) || !nrow(pos) || !ncol(rx)) {
    return(NULL)
  }
  genes <- intersect(rownames(pos), rownames(rx))
  if (length(genes) < 2L) {
    return(NULL)
  }
  pos <- as.matrix(pos[genes, , drop = FALSE])
  rx <- as.matrix(rx[genes, , drop = FALSE])
  nd <- ncol(pos)

  D <- matrix(NA_real_, ncol(rx), nd,
    dimnames = list(colnames(rx), colnames(pos)))
  for (i in seq_len(ncol(rx))) {
    w <- pmax(rx[, i], 0, na.rm = TRUE)**2
    if (!any(w > 0, na.rm = TRUE)) next
    dir <- colMeans(pos * w, na.rm = TRUE)
    nrm <- sqrt(sum(dir**2))
    if (!is.finite(nrm) || nrm <= 0) next
    D[i, ] <- dir / nrm
  }
  D <- D[stats::complete.cases(D), , drop = FALSE]
  if (!nrow(D)) {
    return(NULL)
  }
  nkeep <- suppressWarnings(as.integer(numarrows))
  if (is.finite(nkeep) && nkeep >= 0L && nkeep < nrow(D)) {
    if (nkeep == 0L) {
      return(D[FALSE, , drop = FALSE])
    }
    r2 <- rowSums(D^2)
    D <- D[utils::head(order(-r2), nkeep), , drop = FALSE]
  }
  D
}
