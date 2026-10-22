qsee_filtering_compute <- function(X, Y) {
  X <- X[which(rowSums(is.na(X)) == 0), , drop = FALSE]
  cX <- X - rowMeans(X, na.rm = TRUE)
  nv <- min(10, dim(cX) - 1)
  res <- irlba::irlba(cX, nv = nv)
  rownames(res$u) <- rownames(cX)
  rownames(res$v) <- colnames(cX)
  colnames(res$u) <- paste0("PC", seq_len(ncol(res$u)))
  colnames(res$v) <- paste0("PC", seq_len(ncol(res$v)))
  H <- playbase::expandPhenoMatrix(Y, drop.ref = FALSE)
  res$R <- stats::cor(H, res$v)
  res$X <- X
  res$Y <- Y
  res$H <- H
  res
}
