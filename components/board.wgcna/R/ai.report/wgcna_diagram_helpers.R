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
#' @return A named list with \code{node_styles}, \code{edge_styles},
#'   and \code{show_genes_in_label}
wgcna_diagram_style <- function() {
  list(
    node_styles = list(
      module    = list(bg = "#ADD8E6", border = "#222222", shape = "box"),
      phenotype = list(bg = "#FFFFE0", border = "#222222", shape = "ellipse"),
      process   = list(bg = "#AFEEEE", border = "#222222", shape = "box")
    ),
    edge_styles = list(
      positive    = list(color = "#2E8B57", dashes = FALSE, arrows = "arrow", arrow_on = TRUE,  width = 2),
      negative    = list(color = "#C0392B", dashes = FALSE, arrows = "bar",   arrow_on = TRUE,  width = 2),
      association = list(color = "#94A3B8", dashes = TRUE,  arrows = "arrow", arrow_on = FALSE, width = 1.5)
    ),
    show_genes_in_label = FALSE
  )
}

#' Build structured WGCNA diagram prompt
#'
#' Assembles a structured prompt for diagram generation using
#' \code{omicsai::diagram_prompt()} and \code{omicsai::build_prompt()}.
#' The schema (diagram/network) goes in the system section; board-specific
#' WGCNA rules, species context, and the AI report go in the board section.
#'
#' Node and link type names are derived from the style registry and
#' substituted into the board rules template via \code{frag()}.
#'
#' @param report_text Character string with the AI report text
#' @param organism Character string identifying the organism (e.g. "human", "mouse")
#' @param board_root Character string path to the board.wgcna root directory
#'
#' @return Named list with \code{system} and \code{board} character elements
wgcna_build_diagram_prompt <- function(report_text, organism, board_root) {
  style <- wgcna_diagram_style()
  fmt_names <- function(nms) paste(sprintf("- `%s`", nms), collapse = "\n")
  node_names <- fmt_names(names(style$node_styles))
  link_names <- fmt_names(names(style$edge_styles))

  rules_path <- file.path(board_root, "prompts/wgcna_diagram_rules.md")

  p <- omicsai::diagram_prompt(
    role        = omicsai::frag("system_base"),
    task        = omicsai::frag("diagram/network"),
    species     = omicsai::omicsai_species_prompt(organism),
    board_rules = omicsai::frag(rules_path, list(node_names = node_names, link_names = link_names)),
    report      = paste("## AI Report\n\n", report_text)
  )
  omicsai::build_prompt(p)
}

#' Build structured WGCNA image prompt
#'
#' Assembles a structured prompt for infographic generation using
#' \code{omicsai::image_prompt()} and \code{omicsai::build_prompt()}.
#' Returns a system prompt built from role + task and a board prompt
#' built from species visual context, cleaned report text, and optional
#' diagram edgelist.
#'
#' @param report_text Character string with the AI report text
#' @param organism Character string identifying the organism (e.g. "human", "mouse")
#' @param diagram_edgelist List with \code{$nodes} and \code{$edges} from diagram
#'   result, or \code{NULL} if no diagram has been generated
#'
#' @return Named list with \code{system} and \code{board} character elements
wgcna_build_image_prompt <- function(report_text, organism,
                                     diagram_edgelist = NULL,
                                     style_name = NULL,
                                     n_blocks = 1L) {
  ## Species visual context
  species_img <- omicsai::omicsai_image_species_visual(organism)

  ## Cleaned report — strip noise, humanise module names
  clean <- omicsai::omicsai_strip_report_noise(report_text)
  clean <- gsub("\\bME(\\w+)\\b", "Module \\1", clean)

  ## Diagram edgelist (optional)
  edge_text <- NULL
  if (!is.null(diagram_edgelist)) {
    dot <- omicsai::omicsai_edgelist_to_text(diagram_edgelist)
    if (nzchar(dot)) edge_text <- dot
  }

  ## Style and blocks fragments (resolved by build_prompt via .resolve_slot)
  style_frag  <- if (!is.null(style_name) && nzchar(style_name)) {
    omicsai::frag(paste0("image/styles/", style_name))
  }
  blocks_frag <- omicsai::frag(paste0("image/blocks/blocks_", as.integer(n_blocks)))

  p <- omicsai::image_prompt(
    role             = omicsai::frag("system_base"),
    task             = omicsai::frag("image/infographic", params = list(board_name = "WGCNA")),
    species          = species_img,
    style            = style_frag,
    blocks           = blocks_frag,
    report           = clean,
    diagram_edgelist = edge_text
  )
  omicsai::build_prompt(p)
}
