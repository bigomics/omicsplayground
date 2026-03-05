##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# Drug Connectivity AI Report — Extraction Helpers
# =============================================================================

dc_get_method_data <- function(pgx, method) {
  if (is.null(pgx$drugs) || is.null(method) || !nzchar(method)) return(NULL)
  pgx$drugs[[method]]
}

dc_get_contrasts <- function(dr) {
  if (is.null(dr) || is.null(dr$X)) return(character(0))
  ct <- colnames(dr$X)
  sort(ct[!grepl("^IA:", ct)])
}

dc_parse_tokens <- function(x) {
  x <- as.character(x %||% "")
  x[is.na(x)] <- ""
  x <- enc2utf8(x)
  out <- lapply(x, function(s) {
    trimws(strsplit(s, split = "[\\|;,]")[[1]])
  })
  out
}

#' Build DSEA table for one contrast/method from pre-computed PGX matrices
dc_build_dsea_table <- function(pgx, method, contrast, only_annotated = FALSE) {
  dr <- dc_get_method_data(pgx, method)
  if (is.null(dr) || is.null(contrast) || !contrast %in% colnames(dr$X)) return(NULL)

  nes <- round(dr$X[, contrast], 4)
  pv <- round(dr$P[, contrast], 4)
  qv <- round(dr$Q[, contrast], 4)
  drug <- rownames(dr$X)

  annot <- dr$annot
  if (is.null(annot)) {
    annot <- tryCatch(
      read.csv(file.path(FILESX, "cmap/L1000_repurposing_drugs.txt"),
        sep = "\t", comment.char = "#"
      ),
      error = function(e) NULL
    )
    if (!is.null(annot) && "pert_iname" %in% colnames(annot)) {
      rownames(annot) <- annot$pert_iname
    }
  }

  nes[is.na(nes)] <- 0
  qv[is.na(qv)] <- 1
  pv[is.na(pv)] <- 1

  if (is.null(annot)) {
    dt <- data.frame(
      drug = drug,
      NES = nes,
      pval = pv,
      padj = qv,
      moa = "",
      target = "",
      stringsAsFactors = FALSE
    )
  } else {
    jj <- match(toupper(drug), toupper(rownames(annot)))
    moa_col <- if ("moa" %in% colnames(annot)) "moa" else NA_character_
    target_col <- if ("target" %in% colnames(annot)) "target" else NA_character_
    aa <- data.frame(
      moa = if (!is.na(moa_col)) annot[jj, moa_col] else "",
      target = if (!is.na(target_col)) annot[jj, target_col] else "",
      stringsAsFactors = FALSE
    )
    dt <- data.frame(drug = drug, NES = nes, pval = pv, padj = qv, aa, stringsAsFactors = FALSE)
  }

  dt <- dt[order(-abs(dt$NES)), , drop = FALSE]
  rownames(dt) <- dt$drug

  if (isTRUE(only_annotated)) {
    sel <- which((dt$moa %||% "") != "" | (dt$target %||% "") != "")
    dt <- dt[sel, , drop = FALSE]
  }

  dt
}

dc_build_moa_enrichment <- function(dsea_table, field = c("moa", "target"), nperm = 5000) {
  field <- match.arg(field)
  if (is.null(dsea_table) || !is.data.frame(dsea_table) || nrow(dsea_table) == 0) return(NULL)

  terms.list <- dc_parse_tokens(dsea_table[[field]])
  names(terms.list) <- rownames(dsea_table)
  terms <- setdiff(unique(unlist(terms.list)), c(NA, "", " ", "NA", "N/A"))
  if (length(terms) == 0) return(NULL)

  gmt <- lapply(terms, function(g) {
    names(which(sapply(terms.list, function(t) g %in% t)))
  })
  names(gmt) <- terms

  rnk <- dsea_table$NES
  names(rnk) <- rownames(dsea_table)

  res <- suppressWarnings(
    tryCatch(
      fgsea::fgsea(gmt, rnk, nperm = nperm),
      error = function(e) NULL
    )
  )
  if (is.null(res) || !is.data.frame(res) || nrow(res) == 0) return(NULL)

  res[order(-abs(res$NES)), , drop = FALSE]
}

dc_support_bucket <- function(padj, pval = NA_real_) {
  if (!is.na(padj) && padj < 0.05) return("significant")
  if (!is.na(pval) && pval < 0.05) return("nominal")
  "unsupported"
}

dc_enrichment_direction_subset <- function(dt, direction = c("opposing", "mimicking")) {
  direction <- match.arg(direction)
  if (is.null(dt) || !is.data.frame(dt) || nrow(dt) == 0) return(dt[0, , drop = FALSE])
  sel <- if (direction == "opposing") dt$NES < 0 else dt$NES > 0
  out <- dt[sel, , drop = FALSE]
  if (nrow(out) == 0) return(out)
  out$support <- vapply(seq_len(nrow(out)), function(i) {
    dc_support_bucket(out$padj[i], out$pval[i])
  }, character(1))
  out$support_rank <- match(out$support, c("significant", "nominal", "unsupported"))
  out[order(out$support_rank, out$padj, out$pval, -abs(out$NES), out$pathway), , drop = FALSE]
}

dc_summary_terms <- function(dt, direction = c("opposing", "mimicking"), n = 3L) {
  dt <- dc_enrichment_direction_subset(dt, direction = direction)
  if (is.null(dt) || nrow(dt) == 0) {
    empty <- if (is.data.frame(dt)) dt[0, , drop = FALSE] else data.frame()
    return(list(
      overall = empty,
      supported = empty,
      significant = empty,
      nominal = empty,
      unsupported = empty
    ))
  }

  significant <- dt[dt$support == "significant", , drop = FALSE]
  nominal <- dt[dt$support == "nominal", , drop = FALSE]
  unsupported <- dt[dt$support == "unsupported", , drop = FALSE]
  supported <- if (nrow(significant) > 0) {
    significant
  } else if (nrow(nominal) > 0) {
    nominal
  } else {
    dt[0, , drop = FALSE]
  }

  list(
    overall = head(dt, n),
    supported = head(supported, n),
    significant = head(significant, n),
    nominal = head(nominal, n),
    unsupported = head(unsupported, n)
  )
}

dc_corroborating_targets <- function(dsea_table,
                                     moa_target,
                                     moa_terms,
                                     direction = c("opposing", "mimicking"),
                                     n = 3L) {
  direction <- match.arg(direction)
  if (is.null(dsea_table) || !is.data.frame(dsea_table) || nrow(dsea_table) == 0 ||
      is.null(moa_target) || !is.data.frame(moa_target) || nrow(moa_target) == 0 ||
      is.null(moa_terms) || length(moa_terms) == 0) {
    return(data.frame())
  }

  moa_tokens <- dc_parse_tokens(dsea_table$moa)
  target_tokens <- dc_parse_tokens(dsea_table$target)
  keep <- vapply(moa_tokens, function(tokens) any(tokens %in% moa_terms), logical(1))
  if (!any(keep)) return(moa_target[0, , drop = FALSE])

  supported_targets <- unique(unlist(target_tokens[keep]))
  supported_targets <- setdiff(supported_targets, c("", " ", NA, "NA", "N/A"))
  if (length(supported_targets) == 0) return(moa_target[0, , drop = FALSE])

  dt <- dc_enrichment_direction_subset(moa_target, direction = direction)
  if (is.null(dt) || nrow(dt) == 0) return(dt)
  dt <- dt[dt$pathway %in% supported_targets, , drop = FALSE]
  if (nrow(dt) == 0) return(dt)
  head(dt[order(dt$support_rank, dt$padj, dt$pval, -abs(dt$NES), dt$pathway), , drop = FALSE], n)
}

dc_annotation_confidence <- function(frac_annotated) {
  if (is.na(frac_annotated)) return("unknown")
  if (frac_annotated < 0.20) return("very low")
  if (frac_annotated < 0.30) return("low")
  if (frac_annotated < 0.60) return("moderate")
  "high"
}

dc_build_moa_summary <- function(moa_class, n = 3L) {
  list(
    opposing = dc_summary_terms(moa_class, direction = "opposing", n = n),
    mimicking = dc_summary_terms(moa_class, direction = "mimicking", n = n)
  )
}

dc_build_target_summary <- function(dsea_table, moa_target, moa_summary, n = 3L) {
  opposing_terms <- moa_summary$opposing$supported$pathway %||% character(0)
  mimicking_terms <- moa_summary$mimicking$supported$pathway %||% character(0)

  list(
    opposing = list(
      supported = dc_corroborating_targets(
        dsea_table = dsea_table,
        moa_target = moa_target,
        moa_terms = opposing_terms,
        direction = "opposing",
        n = n
      ),
      overall = dc_summary_terms(moa_target, direction = "opposing", n = n)
    ),
    mimicking = list(
      supported = dc_corroborating_targets(
        dsea_table = dsea_table,
        moa_target = moa_target,
        moa_terms = mimicking_terms,
        direction = "mimicking",
        n = n
      ),
      overall = dc_summary_terms(moa_target, direction = "mimicking", n = n)
    )
  )
}

dc_recommended_exemplars <- function(dsea_table,
                                     moa_terms,
                                     direction = c("opposing", "mimicking"),
                                     n = 3L) {
  direction <- match.arg(direction)
  if (is.null(dsea_table) || !is.data.frame(dsea_table) || nrow(dsea_table) == 0 ||
      is.null(moa_terms) || length(moa_terms) == 0) {
    return(data.frame())
  }

  dt <- dsea_table[!is.na(dsea_table$NES), , drop = FALSE]
  dt <- if (direction == "opposing") {
    dt[dt$NES < 0, , drop = FALSE]
  } else {
    dt[dt$NES > 0, , drop = FALSE]
  }
  if (nrow(dt) == 0) return(dt)

  moa_tokens <- dc_parse_tokens(dt$moa)
  keep <- vapply(moa_tokens, function(tokens) any(tokens %in% moa_terms), logical(1))
  dt <- dt[keep, , drop = FALSE]
  if (nrow(dt) == 0) return(dt)

  dt <- dt[!(is.na(dt$moa) | dt$moa == ""), , drop = FALSE]
  if (nrow(dt) == 0) return(dt)

  dt$support_rank <- vapply(seq_len(nrow(dt)), function(i) {
    match(dc_support_bucket(dt$padj[i], dt$pval[i]), c("significant", "nominal", "unsupported"))
  }, integer(1))

  head(dt[order(dt$support_rank, dt$padj, -abs(dt$NES), dt$drug), , drop = FALSE], n)
}

#' Extract deterministic data for one contrast/method context
extract_drugconnectivity_context_data <- function(pgx,
                                                  contrast,
                                                  method,
                                                  dsea_table = NULL,
                                                  moa_class = NULL,
                                                  moa_target = NULL,
                                                  only_annotated = FALSE,
                                                  n_top = 15L) {
  if (is.null(dsea_table)) {
    dsea_table <- dc_build_dsea_table(
      pgx = pgx,
      method = method,
      contrast = contrast,
      only_annotated = only_annotated
    )
  }

  if (is.null(dsea_table) || !is.data.frame(dsea_table) || nrow(dsea_table) == 0) {
    return(NULL)
  }

  if (is.null(moa_class)) moa_class <- dc_build_moa_enrichment(dsea_table, "moa")
  if (is.null(moa_target)) moa_target <- dc_build_moa_enrichment(dsea_table, "target")

  experiment <- pgx$name %||% pgx$description %||% "omics experiment"
  has_annot <- (dsea_table$moa != "" | dsea_table$target != "")
  has_annot[is.na(has_annot)] <- FALSE

  top_opposing <- head(dsea_table[order(dsea_table$NES), , drop = FALSE], n_top)
  top_mimicking <- head(dsea_table[order(-dsea_table$NES), , drop = FALSE], n_top)
  top_abs <- head(dsea_table[order(-abs(dsea_table$NES)), , drop = FALSE], n_top)
  moa_summary <- dc_build_moa_summary(moa_class, n = 3L)
  target_summary <- dc_build_target_summary(
    dsea_table = dsea_table,
    moa_target = moa_target,
    moa_summary = moa_summary,
    n = 3L
  )
  exemplars <- list(
    opposing = dc_recommended_exemplars(
      dsea_table = dsea_table,
      moa_terms = moa_summary$opposing$supported$pathway %||% character(0),
      direction = "opposing",
      n = 3L
    ),
    mimicking = dc_recommended_exemplars(
      dsea_table = dsea_table,
      moa_terms = moa_summary$mimicking$supported$pathway %||% character(0),
      direction = "mimicking",
      n = 3L
    )
  )
  frac_annotated <- mean(has_annot)

  caveats <- c(
    "Connectivity reflects L1000 cell-line perturbation signatures, not direct clinical efficacy.",
    "Dose, cell context, and perturbation time in L1000 may differ from study biology."
  )

  list(
    context_id = paste(contrast, method, sep = "::"),
    contrast = contrast,
    method = method,
    experiment = experiment,
    n_drugs_tested = nrow(dsea_table),
    n_drugs_sig = sum(dsea_table$padj < 0.05, na.rm = TRUE),
    n_pos = sum(dsea_table$NES > 0, na.rm = TRUE),
    n_neg = sum(dsea_table$NES < 0, na.rm = TRUE),
    max_abs_NES = max(abs(dsea_table$NES), na.rm = TRUE),
    top_opposing = top_opposing,
    top_mimicking = top_mimicking,
    top_abs = top_abs,
    moa_class = moa_class,
    moa_target = moa_target,
    moa_summary = moa_summary,
    target_summary = target_summary,
    exemplars = exemplars,
    reliability = list(
      has_annotations = any(has_annot),
      frac_annotated = frac_annotated,
      annotation_confidence = dc_annotation_confidence(frac_annotated),
      n_sig_moa_classes = if (is.data.frame(moa_class)) sum(moa_class$padj < 0.05, na.rm = TRUE) else 0L,
      n_sig_targets = if (is.data.frame(moa_target)) sum(moa_target$padj < 0.05, na.rm = TRUE) else 0L
    ),
    caveats = caveats
  )
}
