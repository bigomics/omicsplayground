##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

MULTIWGCNA_PROMPTS_DIR <- file.path(OPG, "components/board.wgcna/prompts/multiwgcna")

multiwgcna_overview_nsig <- function(md) md$n_sig %||% 0L

multiwgcna_build_hub_table <- function(md, ntop_genes = 20L) {
  trait_name <- if (length(md$traits) > 0) md$traits[1] else md$top_trait
  gs <- md$trait_stats[[trait_name]] %||% NULL
  if (!is.data.frame(gs) || nrow(gs) == 0) return(NULL)

  if (!"moduleMembership" %in% colnames(gs)) gs$moduleMembership <- NA_real_
  if (!"centrality" %in% colnames(gs)) gs$centrality <- NA_real_
  if ("moduleMembership" %in% colnames(gs)) {
    gs <- gs[order(-abs(gs$moduleMembership)), , drop = FALSE]
  }

  gs <- utils::head(gs, ntop_genes)
  features <- gs$feature
  symbols <- resolve_symbols(features, md$annot)
  funcs <- resolve_functions(features, md$annot)

  data.frame(
    symbol = symbols,
    MM = gs$moduleMembership,
    TS = if ("traitSignificance" %in% colnames(gs)) gs$traitSignificance else NA_real_,
    logFC = if ("foldChange" %in% colnames(gs)) gs$foldChange else NA_real_,
    centrality = gs$centrality,
    func = funcs,
    stringsAsFactors = FALSE
  )
}

multiwgcna_build_report_tables <- function(mwgcna, pgx,
                                           n_modules = 12L,
                                           ntop_enrichment = 12L,
                                           ntop_genes = 20L,
                                           include_contrasts = TRUE) {
  layers <- multiwgcna_ai_data_layers(mwgcna)
  refs <- unname(multiwgcna_ai_module_choices(mwgcna))
  refs <- refs[!duplicated(refs)]

  extracted <- lapply(refs, function(ref) multiwgcna_extract_module_data(mwgcna, ref, pgx))
  names(extracted) <- refs
  extracted <- extracted[!vapply(extracted, is.null, logical(1))]

  score <- vapply(extracted, function(md) {
    abs(md$top_r %||% 0) * 5 + multiwgcna_overview_nsig(md) + nrow(md$cross_links %||% data.frame())
  }, numeric(1))
  score[is.na(score)] <- 0
  selected_refs <- names(sort(score, decreasing = TRUE))
  selected_refs <- utils::head(selected_refs, as.integer(n_modules))

  overview_rows <- lapply(selected_refs, function(ref) {
    md <- extracted[[ref]]
    data.frame(
      Layer = md$layer,
      Module = md$module,
      Features = as.character(md$size),
      Peak_Condition = md$peak_condition %||% "",
      Top_Trait = md$top_trait %||% "",
      Trait_r = if (is.na(md$top_r)) "NA" else sprintf("%+.2f", md$top_r),
      Sig_Enrichments = as.character(md$n_sig),
      Cross_Layer_Links = as.character(nrow(md$cross_links %||% data.frame())),
      stringsAsFactors = FALSE
    )
  })
  overview_df <- do.call(rbind, overview_rows)
  overview_table <- paste(omicsai::omicsai_format_mdtable(overview_df), collapse = "\n")

  module_blocks <- vapply(selected_refs, function(ref) {
    md <- extracted[[ref]]
    lines <- c(sprintf("### %s | %s (%d features)", md$layer, md$module, md$size))

    if (!is.null(md$eigengene_profile)) {
      profile_str <- paste(
        sprintf(
          "%s=%s",
          names(md$eigengene_profile),
          omicsai::omicsai_format_num(md$eigengene_profile, 2)
        ),
        collapse = ", "
      )
      lines <- c(lines, paste0("Eigengene profile: ", profile_str))
    }

    if (nzchar(md$top_trait) && !is.na(md$top_r)) {
      lines <- c(lines, sprintf("Top trait: %s (r=%+.2f)", md$top_trait, md$top_r))
    }

    if (is.data.frame(md$cross_links) && nrow(md$cross_links) > 0) {
      link_str <- paste(
        sprintf("%s | %s (r=%+.2f)", md$cross_links$layer, md$cross_links$module, md$cross_links$r),
        collapse = "; "
      )
      lines <- c(lines, paste0("Cross-layer partners: ", link_str))
    } else {
      lines <- c(lines, "Cross-layer partners: none above |r| >= 0.60")
    }
    lines <- c(lines, "")

    gse <- md$gse
    qcol <- intersect(c("q.value", "padj"), colnames(gse %||% data.frame()))[1]
    term_col <- intersect(c("geneset", "pathway", "term"), colnames(gse %||% data.frame()))[1]
    score_col <- intersect(c("score", "NES"), colnames(gse %||% data.frame()))[1]
    overlap_col <- intersect(c("overlap", "size"), colnames(gse %||% data.frame()))[1]

    if (is.data.frame(gse) && nrow(gse) > 0 && !is.na(term_col)) {
      if (!is.na(qcol)) gse <- gse[order(gse[[qcol]]), , drop = FALSE]
      sig <- if (!is.na(qcol)) gse[!is.na(gse[[qcol]]) & gse[[qcol]] < 0.05, , drop = FALSE] else gse
      sig <- multiwgcna_select_top_enrichment(sig, ntop_enrichment)
      if (nrow(sig) > 0) {
        enrich_df <- data.frame(
          Term = sig[[term_col]],
          Score = if (!is.na(score_col)) sig[[score_col]] else NA_real_,
          `q-value` = if (!is.na(qcol)) sig[[qcol]] else NA_real_,
          Overlap = if (!is.na(overlap_col)) as.character(sig[[overlap_col]]) else "",
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        lines <- c(
          lines,
          sprintf("Top enrichment (%d significant of %d):", md$n_sig, md$n_total),
          omicsai::omicsai_format_mdtable(enrich_df, formatters = list(
            Score = function(x) omicsai::omicsai_format_num(x, 2),
            `q-value` = omicsai::omicsai_format_pvalue
          ))
        )
      } else {
        lines <- c(lines, "No significant enrichment (all q > 0.05)")
      }
    } else {
      lines <- c(lines, "No enrichment table available")
    }
    lines <- c(lines, "")

    hub_df <- multiwgcna_build_hub_table(md, ntop_genes = ntop_genes)
    if (!is.null(hub_df) && nrow(hub_df) > 0) {
      trait_label <- if (length(md$traits) > 0) md$traits[1] else md$top_trait
      logfc_label <- if (nzchar(trait_label)) paste0("logFC (vs ", trait_label, ")") else "logFC"
      lines <- c(
        lines,
        "Hub features:",
        omicsai::omicsai_format_mdtable(
          hub_df,
          col_labels = c(
            symbol = "Feature",
            centrality = "Centrality",
            logFC = logfc_label,
            func = "Known function"
          ),
          formatters = list(
            MM = function(x) omicsai::omicsai_format_num(x, 2),
            TS = function(x) omicsai::omicsai_format_num(x, 2),
            logFC = function(x) omicsai::omicsai_format_num(x, 1),
            centrality = function(x) omicsai::omicsai_format_num(x, 2)
          )
        )
      )
    } else if (length(md$fallback_genes) > 0) {
      lines <- c(lines, paste0("Hub features: ", paste(utils::head(md$fallback_genes, 10), collapse = ", ")))
    }

    if (is.data.frame(md$cross_genes) && nrow(md$cross_genes) > 0) {
      cg <- utils::head(md$cross_genes, 5)
      gene_col <- intersect(c("gene", "feature"), colnames(cg))[1]
      rho_col <- intersect(c("rho", "cor"), colnames(cg))[1]
      if (!is.na(gene_col)) {
        cg_lines <- sprintf(
          "%s%s",
          cg[[gene_col]],
          if (!is.na(rho_col)) paste0(" (rho=", omicsai::omicsai_format_num(cg[[rho_col]], 2), ")") else ""
        )
        lines <- c(lines, paste0("Cross-layer feature evidence: ", paste(cg_lines, collapse = "; ")))
      }
    }

    paste(lines, collapse = "\n")
  }, character(1))

  module_detail <- paste(module_blocks, collapse = "\n\n")

  contrasts_str <- "(no contrasts available)"
  if (include_contrasts && !is.null(pgx$model.parameters$contr.matrix)) {
    cm <- pgx$model.parameters$contr.matrix
    contrast_names <- colnames(cm)
    group_names <- rownames(cm)
    ct_lines <- character(0)
    for (cn in contrast_names) {
      pos <- group_names[cm[, cn] > 0]
      neg <- group_names[cm[, cn] < 0]
      ct_lines <- c(ct_lines, sprintf("- %s: %s vs %s", cn, paste(pos, collapse = "+"), paste(neg, collapse = "+")))
    }
    if (length(ct_lines) > 0) contrasts_str <- paste(ct_lines, collapse = "\n")
  }

  cor_rows <- list()
  idx <- 1L
  for (ref in selected_refs) {
    md <- extracted[[ref]]
    if (!is.data.frame(md$cross_links) || nrow(md$cross_links) == 0) next
    for (i in seq_len(nrow(md$cross_links))) {
      cor_rows[[idx]] <- sprintf(
        "%s | %s <-> %s | %s: r = %+.2f",
        md$layer,
        md$module,
        md$cross_links$layer[i],
        md$cross_links$module[i],
        md$cross_links$r[i]
      )
      idx <- idx + 1L
    }
  }
  cor_rows <- unique(cor_rows)
  module_cors <- if (length(cor_rows) > 0) paste(cor_rows, collapse = "\n") else "(no strong cross-layer correlations detected)"

  template <- omicsai::omicsai_load_template(file.path(MULTIWGCNA_PROMPTS_DIR, "multiwgcna_report_data.md"))
  text <- omicsai::omicsai_substitute_template(template, list(
    experiment = multiomics_ai_experiment_label(pgx),
    organism = pgx$organism %||% "unknown",
    n_samples = as.character(nrow(pgx$samples %||% data.frame())),
    sample_groups = paste(unique(pgx$samples$group %||% "unknown"), collapse = ", "),
    layers = paste(layers, collapse = ", "),
    n_features = as.character(nrow(pgx$X %||% data.frame())),
    wgcna_params = paste0(
      "power=", paste(unique(vapply(layers, function(layer) as.character(mwgcna[[layer]]$power %||% mwgcna[[layer]]$net$power %||% "NA"), character(1))), collapse = ","),
      ", minModSize=", paste(unique(vapply(layers, function(layer) as.character(mwgcna[[layer]]$minModSize %||% "NA"), character(1))), collapse = ","),
      ", mergeCutHeight=", paste(unique(vapply(layers, function(layer) as.character(mwgcna[[layer]]$mergeCutHeight %||% "NA"), character(1))), collapse = ",")
    ),
    overview_table = overview_table,
    module_detail = module_detail,
    contrasts = contrasts_str,
    module_cors = module_cors
  ))

  list(
    text = text,
    data = extracted[selected_refs]
  )
}

multiwgcna_build_summary_params <- function(mwgcna, module_ref, pgx) {
  md <- multiwgcna_extract_module_data(mwgcna, module_ref, pgx)
  if (is.null(md)) return(NULL)

  gse <- md$gse
  ss <- ""
  qcol <- intersect(c("q.value", "padj"), colnames(gse %||% data.frame()))[1]
  term_col <- intersect(c("geneset", "pathway", "term"), colnames(gse %||% data.frame()))[1]
  score_col <- intersect(c("score", "NES"), colnames(gse %||% data.frame()))[1]
  overlap_col <- intersect(c("overlap", "size"), colnames(gse %||% data.frame()))[1]

  if (is.data.frame(gse) && nrow(gse) > 0 && !is.na(term_col)) {
    if (!is.na(qcol)) gse <- gse[order(gse[[qcol]]), , drop = FALSE]
    sig <- if (!is.na(qcol)) gse[!is.na(gse[[qcol]]) & gse[[qcol]] < 0.05, , drop = FALSE] else utils::head(gse, 20)
    top_gse <- utils::head(sig, 20)
    if (nrow(top_gse) > 0) {
      gset_df <- data.frame(
        Pathway = top_gse[[term_col]],
        Score = if (!is.na(score_col)) top_gse[[score_col]] else NA_real_,
        `Q-value` = if (!is.na(qcol)) top_gse[[qcol]] else NA_real_,
        Overlap = if (!is.na(overlap_col)) as.character(top_gse[[overlap_col]]) else "",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      ss <- omicsai::collapse_lines(
        "**Top Enriched Pathways (q < 0.05):**",
        omicsai::omicsai_format_mdtable(gset_df, formatters = list(
          Score = function(x) omicsai::omicsai_format_num(x, 2),
          `Q-value` = omicsai::omicsai_format_pvalue
        )),
        sep = "\n\n"
      )
    }
  }
  if (ss == "" && length(md$fallback_sets) > 0) {
    ss <- paste0("- ", md$fallback_sets, collapse = "\n")
  }

  hub_df <- multiwgcna_build_hub_table(md, ntop_genes = 15L)
  keygenes_section <- ""
  if (!is.null(hub_df) && nrow(hub_df) > 0) {
    keygenes_section <- omicsai::collapse_lines(
      "### Trait-independent metrics",
      omicsai::omicsai_format_mdtable(
        hub_df[, c("symbol", "MM", "centrality", "func"), drop = FALSE],
        col_labels = c(symbol = "Feature", centrality = "Centrality", func = "Known function"),
        formatters = list(
          MM = function(x) omicsai::omicsai_format_num(x, 2),
          centrality = function(x) omicsai::omicsai_format_num(x, 2)
        )
      ),
      "### Trait-specific metrics",
      omicsai::omicsai_format_mdtable(
        hub_df[, c("symbol", "TS", "logFC"), drop = FALSE],
        col_labels = c(symbol = "Feature"),
        formatters = list(
          TS = function(x) omicsai::omicsai_format_num(x, 2),
          logFC = function(x) omicsai::omicsai_format_num(x, 1)
        )
      ),
      sep = "\n\n"
    )
  } else if (length(md$fallback_genes) > 0) {
    keygenes_section <- paste(utils::head(md$fallback_genes, 15), collapse = ", ")
  }

  cross_layer_section <- "No strong cross-layer partners detected."
  if (is.data.frame(md$cross_links) && nrow(md$cross_links) > 0) {
    cross_layer_section <- paste(
      sprintf("%s | %s (r=%+.2f)", md$cross_links$layer, md$cross_links$module, md$cross_links$r),
      collapse = "\n"
    )
  }

  stat_lines <- c(
    "**Module Statistics:**",
    paste0("- **Layer:** ", md$layer),
    paste0("- **Size:** ", md$size, " features")
  )
  if (nzchar(md$top_trait) && !is.na(md$top_r)) {
    stat_lines <- c(stat_lines, paste0("- **Top trait:** ", md$top_trait, " (r=", omicsai::omicsai_format_num(md$top_r, 2), ")"))
  }
  stat_lines <- c(stat_lines, paste0("- **Cross-layer partners above |r| >= 0.60:** ", nrow(md$cross_links %||% data.frame())))

  list(
    layer = md$layer,
    module = md$module,
    phenotypes = md$phenotypes,
    experiment = multiomics_ai_experiment_label(pgx),
    genesets = ss,
    keygenes_section = keygenes_section,
    cross_layer_section = cross_layer_section,
    module_stats = paste(stat_lines, collapse = "\n")
  )
}

multiwgcna_rank_modules <- function(mwgcna, pgx) {
  refs <- unname(multiwgcna_ai_module_choices(mwgcna))
  classification <- list()
  artifact_flags <- list()

  for (ref in refs) {
    md <- multiwgcna_extract_module_data(mwgcna, ref, pgx)
    if (is.null(md)) next

    max_r <- abs(md$top_r %||% 0)
    n_cross <- nrow(md$cross_links %||% data.frame())
    size <- md$size %||% 0L

    if (md$n_sig >= 8 && max_r >= 0.6 && n_cross >= 1 && size >= 25) {
      tier <- "strong"
    } else if (md$n_sig >= 1 || max_r >= 0.5 || n_cross >= 1) {
      tier <- "moderate"
    } else {
      tier <- "weak"
    }

    classification[[ref]] <- list(
      tier = tier,
      n_sig = md$n_sig,
      max_r = max_r,
      size = size,
      n_total = md$n_total
    )

    if (!is.null(md$artifact) && !identical(md$artifact$confidence, "none")) {
      artifact_flags[[ref]] <- md$artifact
    }
  }

  artifact <- names(Filter(function(x) identical(x$confidence, "probable"), artifact_flags))
  strong <- names(Filter(function(x) identical(x$tier, "strong"), classification))
  moderate <- names(Filter(function(x) identical(x$tier, "moderate"), classification))
  weak <- names(Filter(function(x) identical(x$tier, "weak"), classification))

  list(
    strong = strong,
    moderate = moderate,
    weak = weak,
    artifact = artifact,
    artifact_flags = artifact_flags,
    order = c(strong, moderate, weak),
    grey_size = NA_integer_,
    classification = classification
  )
}

multiwgcna_build_methods <- function(mwgcna, pgx) {
  layers <- multiwgcna_ai_data_layers(mwgcna)
  template <- paste(readLines(file.path(MULTIWGCNA_PROMPTS_DIR, "multiwgcna_methods.md"), warn = FALSE), collapse = "\n")

  power_vals <- unique(vapply(layers, function(layer) {
    as.character(mwgcna[[layer]]$power %||% mwgcna[[layer]]$net$power %||% "NA")
  }, character(1)))
  minmod_vals <- unique(vapply(layers, function(layer) {
    as.character(mwgcna[[layer]]$minModSize %||% "NA")
  }, character(1)))
  merge_vals <- unique(vapply(layers, function(layer) {
    as.character(mwgcna[[layer]]$mergeCutHeight %||% "NA")
  }, character(1)))
  minkme_vals <- unique(vapply(layers, function(layer) {
    as.character(mwgcna[[layer]]$minKME %||% "0.3")
  }, character(1)))
  net_vals <- unique(vapply(layers, function(layer) {
    as.character(mwgcna[[layer]]$net$networkType %||% mwgcna[[layer]]$networkType %||% "signed")
  }, character(1)))

  n_features <- sum(vapply(layers, function(layer) {
    length(unlist(mwgcna[[layer]]$me.genes))
  }, numeric(1)))

  n_genesets <- max(vapply(layers, function(layer) {
    w <- mwgcna[[layer]]
    mods <- setdiff(names(w$me.genes), names(w$me.genes)[grepl("grey$", tolower(names(w$me.genes)))])
    first_mod <- mods[1]
    if (length(first_mod) == 0) return(0)
    nrow(w$gsea[[first_mod]] %||% w$gse[[first_mod]] %||% data.frame())
  }, numeric(1)), na.rm = TRUE)

  params <- list(
    n_features_wgcna = as.character(n_features),
    n_layers = as.character(length(layers)),
    layers = paste(layers, collapse = ", "),
    n_samples = as.character(nrow(pgx$samples %||% data.frame())),
    network_type = paste(net_vals, collapse = ", "),
    power = paste(power_vals, collapse = ", "),
    min_mod_size = paste(minmod_vals, collapse = ", "),
    merge_cut_height = paste(merge_vals, collapse = ", "),
    min_kme = paste(minkme_vals, collapse = ", "),
    n_genesets_tested = as.character(n_genesets),
    date = format(Sys.Date(), "%Y-%m-%d")
  )

  omicsai::omicsai_substitute_template(template, params)
}
