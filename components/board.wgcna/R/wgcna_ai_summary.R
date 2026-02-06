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
      pathways <- sub(".*:", "", top_gse$geneset)
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
    pathways <- sub(".*:", "", top$sets[[module]])
    ss <- paste0("- ", pathways, collapse = "\n")
  }

  # Extract key genes with metrics (fallback to names only)
  keygenes_section <- ""
  gene_stats <- NULL
  trait <- NULL
  if (is.character(pp) && pp != "None" && nzchar(pp)) {
    trait <- strsplit(pp, ",\\s*")[[1]][1]
  }
  if (!is.null(trait) && trait != "None") {
    gene_stats <- tryCatch(
      playbase::wgcna.getGeneStats(
        wgcna,
        module = module,
        trait = trait,
        plot = FALSE
      ),
      error = function(e) NULL
    )
  }

  if (is.data.frame(gene_stats) && nrow(gene_stats) > 0) {
    if (!"moduleMembership" %in% colnames(gene_stats)) {
      gene_stats$moduleMembership <- NA_real_
    }
    if (!"traitSignificance" %in% colnames(gene_stats)) {
      gene_stats$traitSignificance <- NA_real_
    }
    if (!"foldChange" %in% colnames(gene_stats)) {
      gene_stats$foldChange <- NA_real_
    }
    if (!"centrality" %in% colnames(gene_stats)) {
      gene_stats$centrality <- NA_real_
    }

    if ("moduleMembership" %in% colnames(gene_stats)) {
      gene_stats <- gene_stats[order(-abs(gene_stats$moduleMembership)), , drop = FALSE]
    } else if ("score" %in% colnames(gene_stats)) {
      gene_stats <- gene_stats[order(-gene_stats$score), , drop = FALSE]
    }

    top_genes <- head(gene_stats, 8)
    symbols <- top_genes$feature
    if (!is.null(annot)) {
      symbols <- tryCatch(
        playbase::probe2symbol(top_genes$feature, annot, "symbol"),
        error = function(e) top_genes$feature
      )
      symbols <- ifelse(is.na(symbols) | symbols == "", top_genes$feature, symbols)
    }

    gene_table <- paste0(
      "| Gene | MM | TS | LogFC | Centrality |\n",
      "|------|----|----|-------|------------|\n",
      paste(
        sprintf(
          "| %s | %s | %s | %s | %s |",
          symbols,
          fmt_num(top_genes$moduleMembership, 2),
          fmt_num(top_genes$traitSignificance, 2),
          fmt_num(top_genes$foldChange, 1),
          fmt_num(top_genes$centrality, 2)
        ),
        collapse = "\n"
      )
    )

    module_size <- NA_integer_
    if (!is.null(wgcna$me.genes) && !is.null(wgcna$me.genes[[module]])) {
      module_size <- length(wgcna$me.genes[[module]])
    } else if (is.data.frame(gene_stats)) {
      module_size <- nrow(gene_stats)
    }

    keygenes_section <- paste0(
      "**Hub Genes (ranked by module membership):**\n\n",
      gene_table, "\n\n",
      "**MM:** Module Membership (correlation with eigengene)  \n",
      "**TS:** Trait Significance (correlation with ", trait, ")  \n",
      "**LogFC:** Log2 fold change  \n",
      "**Top ", nrow(top_genes), " of ",
      if (!is.na(module_size)) module_size else nrow(gene_stats),
      " module genes shown**"
    )
  } else if (!is.null(top$genes[[module]]) && length(top$genes[[module]]) > 0) {
    genes <- head(top$genes[[module]], 15)
    keygenes_section <- paste0(
      "The following hub genes show high intramodular connectivity:\n\n",
      paste(genes, collapse = ", ")
    )
  }

  # Module-level summary statistics
  module_stats <- ""
  module_size <- NA_integer_
  if (!is.null(wgcna$me.genes) && !is.null(wgcna$me.genes[[module]])) {
    module_size <- length(wgcna$me.genes[[module]])
  } else if (is.data.frame(gene_stats) && nrow(gene_stats) > 0) {
    module_size <- nrow(gene_stats)
  }

  lines <- c("**Module Statistics:**")
  if (!is.na(module_size)) {
    lines <- c(lines, paste0("- **Size:** ", module_size, " genes"))
  }

  trait_cor <- NA_real_
  if (!is.null(trait) && trait != "None") {
    M <- tryCatch(playbase::wgcna.get_modTraits(wgcna), error = function(e) NULL)
    if (is.matrix(M) || is.data.frame(M)) {
      module_idx <- which(rownames(M) %in% c(module, paste0("ME", module)))
      trait_idx <- which(colnames(M) == trait)
      if (length(module_idx) > 0 && length(trait_idx) > 0) {
        trait_cor <- M[module_idx[1], trait_idx[1]]
      }
    }
  }
  if (!is.na(trait_cor) && !is.null(trait) && trait != "None") {
    lines <- c(lines, paste0("- **Trait correlation:** ", fmt_num(trait_cor, 2), " with ", trait))
  }

  mean_fc <- NA_real_
  if (is.data.frame(gene_stats) && "foldChange" %in% colnames(gene_stats)) {
    mean_fc <- mean(gene_stats$foldChange, na.rm = TRUE)
  }
  if (!is.na(mean_fc)) {
    fc_direction <- ifelse(mean_fc > 0, "upregulated", "downregulated")
    lines <- c(
      lines,
      paste0("- **Mean expression change:** ", fmt_num(abs(mean_fc), 1), "-fold ", fc_direction)
    )
  }
  module_stats <- paste(lines, collapse = "\n")

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
