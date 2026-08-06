qsee_imputation_plot_pca <- function(res, Y, ph) {
  pcaX <- Filter(Negate(is.null), res$pcaX)
  cc <- 1 + as.integer(factor(Y[, ph]))
  par(mfrow = c(2, 3), cex = 1, mar = c(5, 5, 2, 1))
  for (name in names(pcaX)) {
    pos <- pcaX[[name]]
    plot(pos, main = name, pch = 19, col = cc, xlab = "PC1", ylab = "PC2", cex = 2)
    text(pos, labels = rownames(pos), cex = 0.8, pos = 3)
  }
}

qsee_imputation_plot_heatmaps <- function(res) {
  impX <- res$impX
  cor_hclust <- function(X) {
    corX <- stats::cor(X, use = "pairwise")
    corX[is.na(corX)] <- 0
    fastcluster::hclust(stats::as.dist(1 - corX), method = "ward.D2")
  }
  X <- impX[[1]]
  sel <- rowMeans(is.na(X)) < 0.50 & rowMeans(is.na(X)) > 0
  X <- X[sel, , drop = FALSE]
  X <- X[order(-matrixStats::rowSds(X, na.rm = TRUE)), , drop = FALSE]
  X <- head(X, 200)
  if (!nrow(X)) return(invisible(NULL))
  genes <- rownames(X)
  ii <- cor_hclust(t(X))$order
  jj <- cor_hclust(X)$order
  par(mfrow = c(2, 3), cex = 1, mar = c(5, 1, 2, 6))
  for (name in names(impX)) {
    X1 <- impX[[name]][genes, , drop = FALSE][ii, jj, drop = FALSE]
    playbase::gx.imagemap(X1, clust = FALSE, main = name, cex.main = 1.2)
  }
}

qsee_imputation_plot_histograms <- function(res) {
  impX <- res$impX; ii <- which(!is.na(impX[[1]])); jj <- which(is.na(impX[[1]]))
  par(mfrow = c(2, 3), cex = 1, mar = c(5, 5, 3, 1))
  for (name in names(impX)) {
    X <- impX[[name]]
    hist(X[ii], breaks = 80, main = name, border = "grey60", xlab = "intensity (log2)")
    if (length(jj) && any(!is.na(X[jj]))) hist(X[jj], breaks = 100, add = TRUE, col = "red", border = NA)
  }
}

qsee_imputation_plot_distributions <- function(res, Y, ph) {
  X <- res$impX[[1]]; cc <- 1 + as.integer(factor(Y[, ph]))
  x.avg <- rowMeans(X, na.rm = TRUE); x.nar <- rowMeans(is.na(X)); x.avg2 <- cut(x.avg, breaks = 20)
  x.nar2 <- tapply(seq_len(nrow(X)), x.avg2, function(i) mean(is.na(X[i, , drop = FALSE])))
  aa <- sort(unique(as.numeric(gsub(".*,|\\]", "", as.character(x.avg2)))))
  par(mfrow = c(2, 3), cex = 1, mar = c(5, 5, 3, 1))
  barplot(rbind(x.nar2, 1 - x.nar2), beside = FALSE, names.arg = aa, las = 1,
    col = c("red", "lightgrey"), border = "grey60", xlab = "average intensity (log2)", ylab = "missing ratio")
  title("missingness vs. intensity", cex.main = 1.2)
  plot(x.avg, x.nar, xlab = "average intensity (log2)", ylab = "missing ratio"); title("missingness vs. intensity", cex.main = 1.2)
  if (any(x.nar > 0)) hist(x.nar[x.nar > 0], breaks = 40, main = "missing ratio histogram", xlab = "missing ratio", cex.main = 1.2)
  qq <- seq(0, 1, 0.05); qsum <- sapply(qq, function(q) sum(x.nar <= q))
  plot(qq, qsum, xlab = "max. missing ratio threshold", ylab = "nr. features"); title("nr. features vs. threshold", cex.main = 1.2)
  par(mar = c(7, 5, 3, 1)); barplot(colMeans(is.na(X)), las = 3, cex.names = 0.85, col = cc, ylab = "missing ratio")
  title("missingness per sample", cex.main = 1.2)
}

qsee_imputation_plot_validation <- function(res) {
  impX <- res$impX; V <- res$val_set
  if (!length(V)) { par(mfrow = c(1, 1)); plot.new(); text(0.5, 0.5, "Please add MAR/MNAR missing values."); return(invisible(NULL)) }
  par(mfrow = c(2, 3), cex = 1.1, mar = c(4, 4, 3, 1), mgp = c(2.5, 0.8, 0))
  for (i in 2:length(impX)) {
    jj <- V[, 1]; if (length(jj) > 1e5) jj <- head(jj, 1e5)
    plot(impX[[i]][jj], V[, 2], pch = 19, cex = 0.4, xlim = c(-5, 15), main = names(impX)[i], xlab = "imputed value", ylab = "actual value")
    legend("bottomright", paste("R =", round(cor(impX[[i]][jj], V[, 2]), 2)), bty = "o", bg = "white", cex = 0.9)
  }
}
