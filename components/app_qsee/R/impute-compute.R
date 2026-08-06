qsee_imputation_compute <- function(rawX, marlevel = 0, mnarlevel = 0, progress = NULL) {
  if (!is.null(progress)) progress$set(message = "Normalizing...", value = 0.1)
  normX <- playbase::normalizeExpression(rawX, method = "CPM+quantile", prior = 0)
  val_set <- NULL

  if (mnarlevel > 0) {
    q001 <- quantile(normX, probs = 0.01, na.rm = TRUE)
    sdx <- mean(matrixStats::colSds(scale(t(normX), scale = FALSE), na.rm = TRUE))
    jj <- which(!is.na(normX))
    probx <- pnorm(-abs(normX[jj] - q001), mean = 0, sd = sdx)
    ii <- sample(jj, mnarlevel * length(normX), prob = probx, replace = TRUE)
    val_set <- rbind(val_set, cbind(ii, normX[ii]))
    normX[ii] <- NA
  }
  if (marlevel > 0) {
    jj <- sample(which(!is.na(normX)), marlevel * length(normX))
    val_set <- rbind(val_set, cbind(jj, normX[jj]))
    normX[jj] <- NA
  }

  message("[qsee_imputation_compute] missing values: ", round(100 * mean(is.na(normX)), 2), "%")
  if (!is.null(progress)) progress$set(message = "Computing imputation methods...", value = 0.4)
  impX <- list("no imputation" = normX)
  for (method in c("SVD2", "QRILC", "MinProb", "Perseus")) {
    impX[[method]] <- if (playbase::is.multiomics(rownames(normX))) {
      playbase::imputeMissing.mox(normX, method)
    } else {
      playbase::imputeMissing(normX, method)
    }
  }

  if (!is.null(progress)) progress$set(message = "Computing PCA...", value = 0.8)
  pcaX <- lapply(impX, function(X) {
    cX <- X - rowMeans(X, na.rm = TRUE)
    sel <- which(complete.cases(cX))
    if (length(sel) < 2 || ncol(cX) < 2) return(NULL)
    pos <- irlba::irlba(cX[sel, , drop = FALSE], nv = 2)$v
    dimnames(pos) <- list(colnames(cX), c("PC1", "PC2"))
    pos
  })
  list(normX = normX, impX = impX, pcaX = pcaX, val_set = val_set)
}
