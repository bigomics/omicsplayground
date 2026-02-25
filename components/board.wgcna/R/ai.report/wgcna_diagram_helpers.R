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

#' Build layered WGCNA image prompt
#'
#' Assembles the full prompt for infographic generation by layering:
#' species visual context, cleaned report text, and optional diagram
#' edgelist. Parallel to \code{wgcna_build_diagram_prompt()} but
#' targets the image model rather than the diagram LLM.
#'
#' @param report_text Character string with the AI report text
#' @param organism Character string identifying the organism (e.g. "human", "mouse")
#' @param diagram_edgelist List with \code{$nodes} and \code{$edges} from diagram
#'   result, or \code{NULL} if no diagram has been generated
#'
#' @return Single character string with all prompt layers joined
wgcna_build_image_prompt <- function(report_text, organism, diagram_edgelist = NULL) {
  layers <- list()

  ## Layer 1: species visual context
  layers[[1]] <- omicsai::omicsai_image_species_visual(organism)

  ## Layer 2: cleaned report — strip noise, humanise module names
  clean <- omicsai::omicsai_strip_report_noise(report_text)
  clean <- gsub("\\bME(\\w+)\\b", "Module \\1", clean)
  layers[[2]] <- paste("<report>", clean, "</report>", sep = "\n")

  ## Layer 3: diagram edgelist (if available)
  if (!is.null(diagram_edgelist)) {
    dot <- omicsai::omicsai_edgelist_to_text(diagram_edgelist)
    if (nzchar(dot)) {
      layers[[3]] <- paste("<diagram>", dot, "</diagram>", sep = "\n")
    }
  }

  ## Drop empty layers and join — no markdown separator (unlike the diagram

  ## prompt) because the image model ignores markdown structure; XML tags
  ## on each layer provide sufficient delimitation.
  layers <- layers[nzchar(layers)]
  paste(layers, collapse = "\n\n")
}
