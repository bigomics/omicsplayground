##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

multiwgcna_ai_aux_layers <- function() {
  c("gset", "gs", "pheno", "ph")
}

multiwgcna_ai_unwrap_layers <- function(mwgcna) {
  if (!is.null(mwgcna$layers)) mwgcna$layers else mwgcna
}

multiwgcna_ai_data_layers <- function(mwgcna) {
  if (is.null(mwgcna)) return(character(0))
  mwgcna <- multiwgcna_ai_unwrap_layers(mwgcna)
  layers <- setdiff(names(mwgcna), multiwgcna_ai_aux_layers())
  layers[vapply(layers, function(layer) {
    w <- mwgcna[[layer]]
    !is.null(w$me.genes) && length(w$me.genes) > 0
  }, logical(1))]
}

multiwgcna_module_ref <- function(layer, module) {
  paste(layer, module, sep = "::")
}

multiwgcna_split_module_ref <- function(module_ref) {
  parts <- strsplit(module_ref %||% "", "::", fixed = TRUE)[[1]]
  if (length(parts) < 2) {
    return(list(layer = NA_character_, module = module_ref %||% ""))
  }
  list(
    layer = parts[1],
    module = paste(parts[-1], collapse = "::")
  )
}

multiwgcna_ai_module_choices <- function(mwgcna) {
  layers <- multiwgcna_ai_data_layers(mwgcna)
  refs <- unlist(lapply(layers, function(layer) {
    w <- mwgcna[[layer]]
    mods <- names(sort(lengths(w$me.genes), decreasing = TRUE))
    mods <- mods[!grepl("grey$", tolower(mods))]
    stats::setNames(
      vapply(mods, function(mod) multiwgcna_module_ref(layer, mod), character(1)),
      paste(layer, mods, sep = " | ")
    )
  }), use.names = TRUE)
  refs[!duplicated(unname(refs))]
}

multiwgcna_select_top_enrichment <- function(sig_gse, n) {
  if (!is.data.frame(sig_gse) || nrow(sig_gse) == 0) return(sig_gse)

  term_col <- intersect(c("geneset", "pathway", "term"), colnames(sig_gse))[1]
  if (is.na(term_col)) return(utils::head(sig_gse, n))

  go_id <- ifelse(
    grepl("GO_\\d{5,}", sig_gse[[term_col]]),
    sub(".*?(GO_\\d{5,}).*", "\\1", sig_gse[[term_col]]),
    NA_character_
  )
  dup <- duplicated(go_id) & !is.na(go_id)
  sig_gse <- sig_gse[!dup, , drop = FALSE]
  utils::head(sig_gse, n)
}

multiwgcna_classify_artifact <- function(layer_wgcna, module) {
  if (is.null(layer_wgcna)) return(list(confidence = "none", reason = NULL))
  tryCatch(classify_artifact(layer_wgcna, module), error = function(e) {
    list(confidence = "none", reason = NULL)
  })
}

multiwgcna_module_links <- function(mwgcna, module_ref, min_abs_r = 0.6, top_n = 5L) {
  ref <- multiwgcna_split_module_ref(module_ref)
  layers <- setdiff(multiwgcna_ai_data_layers(mwgcna), ref$layer)
  base <- mwgcna[[ref$layer]]
  if (is.null(base)) return(data.frame())

  base_me <- base$datME %||% base$net$MEs
  if (is.null(base_me) || !ref$module %in% colnames(base_me)) return(data.frame())

  x <- base_me[, ref$module]
  rows <- list()
  idx <- 1L

  for (layer in layers) {
    other <- mwgcna[[layer]]
    other_me <- other$datME %||% other$net$MEs
    if (is.null(other_me)) next

    mods <- colnames(other_me)
    mods <- mods[!grepl("grey$", tolower(mods))]
    for (mod in mods) {
      y <- other_me[, mod]
      r <- tryCatch(stats::cor(x, y, use = "pairwise.complete.obs"), error = function(e) NA_real_)
      if (is.na(r) || abs(r) < min_abs_r) next
      rows[[idx]] <- data.frame(
        partner_ref = multiwgcna_module_ref(layer, mod),
        layer = layer,
        module = mod,
        r = as.numeric(r),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  if (length(rows) == 0) return(data.frame())
  out <- do.call(rbind, rows)
  out <- out[order(-abs(out$r)), , drop = FALSE]
  utils::head(out, as.integer(top_n))
}

multiwgcna_extract_module_data <- function(mwgcna, module_ref, pgx, max_traits = 3L) {
  ref <- multiwgcna_split_module_ref(module_ref)
  layer <- ref$layer
  module <- ref$module

  if (is.na(layer) || !layer %in% names(mwgcna)) return(NULL)
  wgcna <- mwgcna[[layer]]
  annot <- wgcna$annot %||% pgx$genes

  top <- tryCatch(
    playbase::wgcna.getTopGenesAndSets(
      wgcna,
      annot = annot,
      ntop = 40,
      level = "gene",
      rename = "gene_title"
    ),
    error = function(e) list(pheno = list(), genes = list(), sets = list())
  )

  M <- tryCatch(playbase:::wgcna.get_modTraits(wgcna), error = function(e) NULL)

  phenotypes <- "None"
  if (module %in% names(top$pheno)) {
    phenotypes <- paste(top$pheno[[module]], collapse = ", ")
  }

  top_trait <- ""
  top_r <- NA_real_
  if (!is.null(M) && module %in% rownames(M)) {
    cors <- M[module, ]
    cors <- cors[!is.na(cors)]
    if (length(cors) > 0) {
      top_idx <- which.max(abs(cors))
      top_trait <- names(cors)[top_idx]
      top_r <- cors[top_idx]
    }
  }

  traits <- character(0)
  if (!identical(phenotypes, "None") && nzchar(phenotypes)) {
    traits <- trimws(strsplit(phenotypes, ",\\s*")[[1]])
    traits <- traits[nzchar(traits) & traits != "None"]
    traits <- utils::head(traits, max_traits)
  }
  if (length(traits) == 0 && nzchar(top_trait)) traits <- top_trait

  eigengene_profile <- NULL
  peak_condition <- ""
  datME <- wgcna$datME %||% wgcna$net$MEs
  groups <- pgx$samples$group %||% NULL
  if (!is.null(datME) && !is.null(groups) && module %in% colnames(datME)) {
    group_means <- tapply(datME[, module], groups, mean)
    eigengene_profile <- group_means
    peak_condition <- names(which.max(abs(group_means)))
  }

  gse <- wgcna$gsea[[module]] %||% wgcna$gse[[module]]
  n_sig <- 0L
  n_total <- 0L
  if (is.data.frame(gse) && nrow(gse) > 0) {
    qcol <- intersect(c("q.value", "padj"), colnames(gse))[1]
    n_total <- nrow(gse)
    if (!is.na(qcol)) {
      n_sig <- sum(gse[[qcol]] < 0.05, na.rm = TRUE)
    }
  }

  trait_stats <- list()
  for (tr in traits) {
    gs <- tryCatch(
      playbase::wgcna.getGeneStats(wgcna, module = module, trait = tr, plot = FALSE),
      error = function(e) NULL
    )
    if (is.data.frame(gs) && nrow(gs) > 0) {
      trait_stats[[tr]] <- gs
    }
  }

  cross_links <- multiwgcna_module_links(mwgcna, module_ref)
  cross_genes <- tryCatch(
    playbase::wgcna.getModuleCrossGenes(
      mwgcna,
      ref = NULL,
      multi = TRUE,
      ngenes = 100,
      modules = module
    )[[1]],
    error = function(e) NULL
  )

  list(
    layer = layer,
    module = module,
    module_ref = module_ref,
    annot = annot,
    size = length(wgcna$me.genes[[module]] %||% character(0)),
    phenotypes = phenotypes,
    traits = traits,
    fallback_genes = top$genes[[module]] %||% character(0),
    fallback_sets = top$sets[[module]] %||% character(0),
    eigengene_profile = eigengene_profile,
    peak_condition = peak_condition,
    top_trait = top_trait,
    top_r = top_r,
    gse = gse,
    n_sig = n_sig,
    n_total = n_total,
    trait_stats = trait_stats,
    cross_links = cross_links,
    cross_genes = cross_genes,
    artifact = multiwgcna_classify_artifact(wgcna, module)
  )
}
