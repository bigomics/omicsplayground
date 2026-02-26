##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

multiomics_ai_experiment_label <- function(pgx, experiment = NULL) {
  experiment %||% pgx$name %||% pgx$description %||% "omics experiment"
}

multiomics_ai_drop_interaction_contrasts <- function(contrasts) {
  contrasts <- contrasts %||% character(0)
  sort(unique(contrasts[!grepl("^IA", contrasts)]))
}

multiomics_ai_top_idx_by_abs <- function(n, ...) {
  vals <- list(...)
  if (length(vals) == 0) return(integer(0))
  len <- length(vals[[1]])
  if (len == 0) return(integer(0))

  score <- rep(0, len)
  for (v in vals) {
    if (length(v) != len) next
    a <- abs(as.numeric(v))
    a[is.na(a)] <- -Inf
    score <- score + a
  }

  head(order(score, decreasing = TRUE, na.last = TRUE), as.integer(n))
}
