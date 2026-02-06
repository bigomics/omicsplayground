##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Build AI prompt parameters for WGCNA module summary
#'
#' Extracts parameters from WGCNA results to populate the prompt template.
#' Replicates exact extraction logic from playbase::wgcna.describeModules()
#' (see playbase/R/pgx-wgcna.R lines 5340-5361)
#'
#' @param wgcna WGCNA results object
#' @param module Character; module name (e.g., "blue", "turquoise")
#' @param annot Data frame; gene annotation (optional)
#' @param ntop Integer; number of top genes/sets to include (default 40)
#' @param multi Logical; TRUE for multi-omics WGCNA
#'
#' @return Named list with template parameters (board_params):
#'   module, phenotypes, experiment, genesets, keygenes_section
wgcna_build_ai_params <- function(wgcna,
                                  module,
                                  annot = NULL,
                                  ntop = 40,
                                  multi = FALSE) {
  fmt_num <- function(x, digits = 2) {
    ifelse(is.na(x), "NA", sprintf(paste0("%.", digits, "f"), x))
  }

  # Get top genes and sets using playbase function

# (mirrors playbase/R/pgx-wgcna.R lines 5280-5288)
  if (multi) {
    top <- playbase::wgcna.getMultiTopGenesAndSets(
      wgcna,
      annot = annot,
      ntop = ntop,
      level = "gene",
      rename = "gene_title"
    )
  } else {
    top <- playbase::wgcna.getTopGenesAndSets(
      wgcna,
      annot = annot,
      ntop = ntop,
      level = "gene",
      rename = "gene_title"
    )
  }

  # Extract phenotypes - format as comma-separated list
  pp <- "None"
  if (module %in% names(top$pheno)) {
    pp <- paste(top$pheno[[module]], collapse = ", ")
  }

  # Extract genesets with quantitative metrics (fallback to names only)
  ss <- ""
  gse <- NULL
  if (!is.null(wgcna$gse) && !is.null(wgcna$gse[[module]])) {
    gse <- wgcna$gse[[module]]
  }
  if (is.data.frame(gse) && nrow(gse) > 0 &&
    all(c("geneset", "score", "q.value", "overlap") %in% colnames(gse))) {
    gse <- gse[order(gse$q.value, gse$score, na.last = TRUE), , drop = FALSE]
    sig <- gse[!is.na(gse$q.value) & gse$q.value < 0.05, , drop = FALSE]
    top_gse <- head(sig, 8)
    if (nrow(top_gse) > 0) {
      pathways <- top_gse$geneset
      gset_table <- paste0(
        "| Pathway | Score | Q-value | Overlap |\n",
        "|---------|-------|---------|---------|\n",
        paste(
          sprintf(
            "| %s | %s | %s | %s |",
            pathways,
            fmt_num(top_gse$score, 2),
            fmt_num(top_gse$q.value, 3),
            top_gse$overlap
          ),
          collapse = "\n"
        )
      )

      score_vals <- gse$score
      score_vals <- score_vals[!is.na(score_vals)]
      score_range <- if (length(score_vals) > 0) {
        paste0(
          fmt_num(min(score_vals), 2), "-",
          fmt_num(max(score_vals), 2),
          " (median: ", fmt_num(median(score_vals), 2), ")"
        )
      } else {
        "NA"
      }
      n_sig <- sum(gse$q.value < 0.05, na.rm = TRUE)

      ss <- paste0(
        "**Top Enriched Pathways (q < 0.05):**\n\n",
        gset_table, "\n\n",
        "**Score range:** ", score_range, "  \n",
        "**Total significant pathways:** ", n_sig, " of ", nrow(gse), " tested"
      )
    }
  }
  if (ss == "" && !is.null(top$sets[[module]])) {
    pathways <- top$sets[[module]]
    ss <- paste0("- ", pathways, collapse = "\n")
  }

  # Extract key genes with metrics — multi-trait approach
  # WGCNA modules are trait-independent (MM, centrality are network properties).
  # TS and logFC are trait-dependent, so we compute them for ALL correlated traits
  # to give the AI a complete picture across phenotypes.
  keygenes_section <- ""
  traits <- NULL
  if (is.character(pp) && pp != "None" && nzchar(pp)) {
    traits <- trimws(strsplit(pp, ",\\s*")[[1]])
    traits <- traits[nzchar(traits) & traits != "None"]
    traits <- head(traits, 3) # cap to avoid very wide tables
  }

  # Collect gene stats for each correlated trait
  trait_stats_list <- list()
  if (length(traits) > 0) {
    for (tr in traits) {
      gs <- tryCatch(
        playbase::wgcna.getGeneStats(
          wgcna,
          module = module,
          trait = tr,
          plot = FALSE
        ),
        error = function(e) NULL
      )
      if (is.data.frame(gs) && nrow(gs) > 0) {
        trait_stats_list[[tr]] <- gs
      }
    }
  }

  if (length(trait_stats_list) > 0) {
    # Use first trait's stats as base for shared columns (MM, centrality)
    base_stats <- trait_stats_list[[1]]
    if (!"moduleMembership" %in% colnames(base_stats)) {
      base_stats$moduleMembership <- NA_real_
    }
    if (!"centrality" %in% colnames(base_stats)) {
      base_stats$centrality <- NA_real_
    }

    # Sort by module membership (trait-independent)
    if ("moduleMembership" %in% colnames(base_stats)) {
      base_stats <- base_stats[order(-abs(base_stats$moduleMembership)), , drop = FALSE]
    } else if ("score" %in% colnames(base_stats)) {
      base_stats <- base_stats[order(-base_stats$score), , drop = FALSE]
    }

    top_genes <- head(base_stats, 8)
    features <- top_genes$feature

    # Resolve gene symbols
    symbols <- features
    if (!is.null(annot)) {
      symbols <- tryCatch(
        playbase::probe2symbol(features, annot, "symbol"),
        error = function(e) features
      )
      symbols <- ifelse(is.na(symbols) | symbols == "", features, symbols)
    }

    # --- Table 1: trait-independent network metrics ---
    network_table <- paste0(
      "### Trait-independent metrics (network structure)\n\n",
      "| Gene | MM | Centrality |\n",
      "|------|----|------------|\n",
      paste(
        sprintf(
          "| %s | %s | %s |",
          symbols,
          fmt_num(top_genes$moduleMembership, 2),
          fmt_num(top_genes$centrality, 2)
        ),
        collapse = "\n"
      )
    )

    # --- Tables 2..N: per-trait gene-trait associations ---
    trait_names <- names(trait_stats_list)
    M <- tryCatch(playbase::wgcna.get_modTraits(wgcna), error = function(e) NULL)

    trait_tables <- vapply(trait_names, function(tr) {
      gs <- trait_stats_list[[tr]]

      # Module-trait correlation for header
      mod_cor <- NA_real_
      if (is.matrix(M) || is.data.frame(M)) {
        mi <- which(rownames(M) %in% c(module, paste0("ME", module)))
        ti <- which(colnames(M) == tr)
        if (length(mi) > 0 && length(ti) > 0) mod_cor <- M[mi[1], ti[1]]
      }
      cor_label <- if (!is.na(mod_cor)) {
        paste0(" (module r = ", fmt_num(mod_cor, 2), ")")
      } else {
        ""
      }

      # Build rows for each hub gene
      rows <- vapply(seq_along(features), function(i) {
        idx <- match(features[i], gs$feature)
        ts_val <- NA_real_
        fc_val <- NA_real_
        if (!is.na(idx)) {
          if ("traitSignificance" %in% colnames(gs)) ts_val <- gs$traitSignificance[idx]
          if ("foldChange" %in% colnames(gs)) fc_val <- gs$foldChange[idx]
        }
        sprintf("| %s | %s | %s |", symbols[i], fmt_num(ts_val, 2), fmt_num(fc_val, 1))
      }, character(1))

      paste0(
        "### Gene-trait associations: ", tr, cor_label, "\n\n",
        "| Gene | TS | logFC |\n",
        "|------|----|-------|\n",
        paste(rows, collapse = "\n")
      )
    }, character(1))

    # Module size
    module_size <- NA_integer_
    if (!is.null(wgcna$me.genes) && !is.null(wgcna$me.genes[[module]])) {
      module_size <- length(wgcna$me.genes[[module]])
    } else {
      module_size <- nrow(base_stats)
    }

    keygenes_section <- paste0(
      "**Hub Genes (ranked by module membership, top ",
      nrow(top_genes), " of ",
      if (!is.na(module_size)) module_size else nrow(base_stats),
      " module genes):**\n\n",
      network_table, "\n\n",
      paste(trait_tables, collapse = "\n\n"),
      "\n\n",
      "**MM:** Module Membership (correlation with eigengene) — trait-independent  \n",
      "**Centrality:** Intramodular connectivity — trait-independent  \n",
      "**TS:** Trait Significance (gene-trait correlation) — specific to each trait  \n",
      "**logFC:** Log2 fold change — specific to each trait"
    )
  } else if (!is.null(top$genes[[module]]) && length(top$genes[[module]]) > 0) {
    genes <- head(top$genes[[module]], 15)
    keygenes_section <- paste0(
      "The following hub genes show high intramodular connectivity:\n\n",
      paste(genes, collapse = ", ")
    )
  }

  # Module-level summary statistics (multi-trait)
  module_stats <- ""
  module_size <- NA_integer_
  if (!is.null(wgcna$me.genes) && !is.null(wgcna$me.genes[[module]])) {
    module_size <- length(wgcna$me.genes[[module]])
  } else if (length(trait_stats_list) > 0) {
    module_size <- nrow(trait_stats_list[[1]])
  }

  stat_lines <- c("**Module Statistics:**")
  if (!is.na(module_size)) {
    stat_lines <- c(stat_lines, paste0("- **Size:** ", module_size, " genes"))
  }

  # Per-trait: module-trait correlation and mean fold change
  if (length(traits) > 0) {
    M <- tryCatch(playbase::wgcna.get_modTraits(wgcna), error = function(e) NULL)
    for (tr in traits) {
      trait_cor <- NA_real_
      if (is.matrix(M) || is.data.frame(M)) {
        module_idx <- which(rownames(M) %in% c(module, paste0("ME", module)))
        trait_idx <- which(colnames(M) == tr)
        if (length(module_idx) > 0 && length(trait_idx) > 0) {
          trait_cor <- M[module_idx[1], trait_idx[1]]
        }
      }
      if (!is.na(trait_cor)) {
        stat_lines <- c(stat_lines, paste0(
          "- **Trait correlation with ", tr, ":** ", fmt_num(trait_cor, 2)
        ))
      }
      # Mean fold change for this trait
      if (tr %in% names(trait_stats_list)) {
        gs <- trait_stats_list[[tr]]
        if ("foldChange" %in% colnames(gs)) {
          mean_fc <- mean(gs$foldChange, na.rm = TRUE)
          if (!is.na(mean_fc)) {
            fc_dir <- ifelse(mean_fc > 0, "upregulated", "downregulated")
            stat_lines <- c(stat_lines, paste0(
              "- **Mean expression change (", tr, "):** ",
              fmt_num(abs(mean_fc), 1), "-fold ", fc_dir
            ))
          }
        }
      }
    }
  }

  module_stats <- paste(stat_lines, collapse = "\n")

  # Get experiment description
  experiment <- wgcna$experiment %||% ""

  list(
    module = module,
    phenotypes = pp,
    experiment = experiment,
    genesets = ss,
    keygenes_section = keygenes_section,
    module_stats = module_stats
  )
}

# -----------------------------------------------------------------------------
# Shiny Module: WGCNA AI Summary
# -----------------------------------------------------------------------------

#' WGCNA AI Summary UI
#'
#' Creates AI summary UI wrapped in PlotModuleUI for consistent OmicsPlayground styling.
#'
#' @param id Shiny module namespace ID
#' @param title Card title
#' @param label Optional label
#' @param info.text Help/info text
#' @param caption Caption text
#' @param height Card height (can be vector for responsive sizing)
#' @param width Card width (can be vector for responsive sizing)
#'
#' @return Shiny UI element
wgcna_ai_summary_ui <- function(id,
                                title = "AI Summary",
                                label = "",
                                info.text = "",
                                caption = "AI-generated module summary.",
                                height = c("100%", TABLE_HEIGHT_MODAL),
                                width = c("auto", "100%")) {
  # PlotModuleUI wrapper for omicsai card integration
  card_wrapper <- function(id, content, options, title, label, info.text,
                           caption, height, width, download.fmt, ...) {
    PlotModuleUI(
      id,
      outputFunc = shiny::htmlOutput,
      title = title,
      label = label,
      info.text = info.text,
      options = options,
      caption = caption,
      height = height,
      width = width,
      download.fmt = download.fmt
    )
  }

  omics.ai::omicsai_summary_card_ui(
    id = id,
    card_wrapper = card_wrapper,
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width
  )
}

#' WGCNA AI Summary Server
#'
#' Server logic for WGCNA AI-generated module summaries.
#' Uses omicsai_summary_card_server with PlotModuleServer integration.
#'
#' @param id Shiny module namespace ID
#' @param wgcna Reactive returning WGCNA results object
#' @param r_module Reactive returning selected module name
#' @param r_annot Reactive returning gene annotation (optional)
#' @param session Shiny session object (required for getUserOption)
#' @param multi Logical; TRUE for multi-omics WGCNA
#' @param cache OmicsAI cache object (optional)
#' @param watermark Logical; add watermark to output
#'
#' @return Result from omicsai_summary_card_server
wgcna_ai_summary_server <- function(id,
                                    wgcna,
                                    r_module,
                                    r_annot = shiny::reactive(NULL),
                                    session,
                                    multi = FALSE,
                                    cache = NULL,
                                    watermark = FALSE) {

  # Build board_params reactive (data to inject into user prompt template)
  params_reactive <- shiny::reactive({
    w <- wgcna()
    module <- r_module()
    annot <- r_annot()
    shiny::req(w, module)

    # Get annotation from wgcna if not provided
    if (is.null(annot) && !is.null(w$annot)) {
      annot <- w$annot
    }

    wgcna_build_ai_params(
      wgcna = w,
      module = module,
      annot = annot,
      ntop = 40,
      multi = multi
    )
  })

  # Load template from board.wgcna prompts directory
  board_template_path <- file.path(
    OPG,
    "components/board.wgcna/prompts/wgcna_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omics.ai::omicsai_load_template(board_template_path)
  })

  # Load WGCNA methods context template
  context_template_path <- file.path(
    OPG,
    "components/board.wgcna/prompts/WGCNA_methods.md"
  )
  context_template <- omics.ai::omicsai_load_template(context_template_path)

  # Build config reactive from user's selected LLM model
  # Includes board-specific context in system prompt
  config_reactive <- shiny::reactive({
    ai_model <- getUserOption(session, "llm_model")
    w <- wgcna()
    shiny::req(ai_model, ai_model != "", w)

    # Build context params for template substitution
    context_params <- list(
      experiment = w$experiment %||% "omics experiment"
    )

    omics.ai::omicsai_config(
      model = ai_model,
      context = context_template,
      context_params = context_params
    )
  })

  # PlotModuleServer wrapper
  plot_server_wrapper <- function(id, func, func2, watermark) {
    PlotModuleServer(
      id,
      plotlib = "generic",
      plotlib2 = "generic",
      func = func,
      func2 = func2,
      renderFunc = shiny::renderUI,
      renderFunc2 = shiny::renderUI,
      pdf.width = 8, pdf.height = 5,
      res = c(75, 100),
      add.watermark = watermark
    )
  }

  # Initialize omicsai card module
  omics.ai::omicsai_summary_card_server(
    id = id,
    params_reactive = params_reactive,
    template_reactive = template_reactive,
    config_reactive = config_reactive,
    plot_server_wrapper = plot_server_wrapper,
    cache = cache,
    watermark = watermark
  )
}
