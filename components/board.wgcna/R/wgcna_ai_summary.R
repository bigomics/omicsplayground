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
#' @param docstyle Character; "short summary" or "long detailed scientific discussion"
#' @param numpar Integer; maximum paragraphs (default 2)
#' @param annot Data frame; gene annotation (optional)
#' @param ntop Integer; number of top genes/sets to include (default 40)
#' @param multi Logical; TRUE for multi-omics WGCNA
#'
#' @return Named list with template parameters:
#'   docstyle, module, phenotypes, experiment, numpar, genesets, keygenes_section
wgcna_build_ai_params <- function(wgcna,
                                  module,
                                  docstyle = "short summary",
                                  numpar = 2,
                                  annot = NULL,
                                  ntop = 40,
                                  multi = FALSE) {
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

  # Extract genesets - strip prefix, semicolon-joined
# (mirrors playbase/R/pgx-wgcna.R lines 5344-5345)
  ss <- ""
  if (!is.null(top$sets[[module]])) {
    ss <- sub(".*:", "", top$sets[[module]])
    ss <- paste(ss, collapse = ";")
  }

  # Extract phenotypes - quoted, semicolon-joined
# (mirrors playbase/R/pgx-wgcna.R lines 5346-5349)
  pp <- ""
  if (module %in% names(top$pheno)) {
    pp <- paste0("'", top$pheno[[module]], "'")
    pp <- paste(pp, collapse = ";")
  }

  # Extract key genes - semicolon-joined
# (mirrors playbase/R/pgx-wgcna.R lines 5352-5355)
  keygenes_section <- ""
  if (!is.null(top$genes[[module]]) && length(top$genes[[module]]) > 0) {
    gg <- paste(top$genes[[module]], collapse = ";")
    keygenes_section <- paste0(
      "\nAfter that, shortly discuss if any of these key genes/proteins/metabolites ",
      "might be involved in the biological function. No need to mention all, just a few. ",
      "Here is the list of key genes/proteins/metabolites: ", gg, "\n"
    )
  }

  # Get experiment description
  experiment <- if (is.null(wgcna$experiment)) "" else wgcna$experiment

  list(
    docstyle = docstyle,
    module = module,
    phenotypes = pp,
    experiment = experiment,
    numpar = as.character(numpar),
    genesets = ss,
    keygenes_section = keygenes_section
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

  # Build params reactive for omicsai (note: docstyle handled by card server)
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
      docstyle = "short summary",  # Will be overridden by card UI input
      numpar = 2,
      annot = annot,
      ntop = 40,
      multi = multi
    )
  })

  # Load template from board.wgcna templates directory
  template_path <- file.path(
    OPG,
    "components/board.wgcna/templates/wgcna_summary_template.md"
  )

  template_reactive <- shiny::reactive({
    omics.ai::omicsai_load_template(template_path)
  })

  # Build config reactive from user's selected LLM model
  # Uses session passed from caller (no moduleServer wrapper needed)
  config_reactive <- shiny::reactive({
    ai_model <- getUserOption(session, "llm_model")
    shiny::req(ai_model, ai_model != "")
    omics.ai::omicsai_config(ai_model)
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

  # Initialize omicsai card module (no double moduleServer wrapping)
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
