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
    reliability = list(
      has_annotations = any(has_annot),
      frac_annotated = mean(has_annot),
      n_sig_moa_classes = if (is.data.frame(moa_class)) sum(moa_class$padj < 0.05, na.rm = TRUE) else 0L,
      n_sig_targets = if (is.data.frame(moa_target)) sum(moa_target$padj < 0.05, na.rm = TRUE) else 0L
    ),
    caveats = caveats
  )
}
