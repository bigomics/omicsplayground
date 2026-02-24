##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' WGCNA diagram visual style
#'
#' Returns the style list for \code{omicsai::omicsai_diagram_render()}.
#' Encodes the WGCNA board aesthetic choices: green/red solid edges,
#' box-shaped process nodes, and genes in tooltip only.
#'
#' @return A named list of style parameters
wgcna_diagram_style <- function() {
  list(
    edge_colors = list(positive = "#2E8B57", negative = "#C0392B", association = "#94A3B8"),
    edge_dashes = list(positive = FALSE, negative = FALSE, association = TRUE),
    module_bg        = "#ADD8E6",
    module_border    = "#222222",
    phenotype_bg     = "#FFFFE0",
    phenotype_border = "#222222",
    process_bg       = "#AFEEEE",
    process_border   = "#222222",
    process_shape    = "box",
    show_genes_in_label = FALSE
  )
}

#' Build layered WGCNA diagram prompt
#'
#' Assembles the full prompt for diagram generation by layering:
#' generic schema instructions, board-specific WGCNA rules, species
#' context, and the AI report text.
#'
#' @param report_text Character string with the AI report text
#' @param organism Character string identifying the organism (e.g. "human", "mouse")
#' @param board_root Character string path to the board.wgcna root directory
#'
#' @return Single character string with all prompt layers joined
wgcna_build_diagram_prompt <- function(report_text, organism, board_root) {
  layers <- list()

  ## Layer 1: generic diagram JSON schema instructions
  layers[[1]] <- omicsai::omicsai_instructions("diagram_network")

  ## Layer 2: board-specific WGCNA rules
  layers[[2]] <- tryCatch(
    omicsai::omicsai_load_template("prompts/diagram_wgcna_rules.md", root = board_root),
    error = function(e) ""
  )

  ## Layer 3: species-aware context
  layers[[3]] <- tryCatch(
    omicsai::omicsai_species_prompt(organism),
    error = function(e) ""
  )

  ## Layer 4: the AI report itself (this is the data for the LLM)
  layers[[4]] <- paste("## AI Report\n\n", report_text)

  ## Drop empty layers and join with separator
  layers <- layers[nzchar(layers)]
  paste(layers, collapse = "\n\n---\n\n")
}
