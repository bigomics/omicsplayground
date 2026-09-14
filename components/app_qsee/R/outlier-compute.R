# Outlier diagnostics for the Qsee board.
#
# This file owns PCA and heatmap preparation for outlier visualization.
# Numerical scoring and imputation use canonical matrix families.

qsee_outlier_compute <- function(X, Y, progress = NULL) {
  Xc <- X[stats::complete.cases(X), , drop = FALSE]
  if (!is.null(progress)) progress$set(message = "Detecting outliers...", value = 0.2)
  outliers <- playbase::pp.removeOutliers(
    Xc,
    threshold = Inf,
    methods = c("z.correlation", "z.distance", "z.features", "z.isoforest")
  )$scores
  if (!is.null(progress)) progress$set(message = "Computing PCA...", value = 0.55)
  cX <- X - rowMeans(X, na.rm = TRUE)
  if (any(is.na(cX))) cX <- .opg_impute(cX, method = "SVD2")
  pca <- svd(cX, nu = 0, nv = 2)$v[, 1:2, drop = FALSE]
  rownames(pca) <- colnames(X)
  colnames(pca) <- c("PC1", "PC2")
  if (!is.null(progress)) progress$set(message = "Preparing heatmap...", value = 0.85)
  nmax <- min(400L, nrow(X))
  heatX <- X[head(order(-matrixStats::rowSds(X, na.rm = TRUE)), nmax), , drop = FALSE]
  list(X = X, Xc = Xc, Y = Y, outliers = outliers, pca = pca, heatX = heatX)
}
