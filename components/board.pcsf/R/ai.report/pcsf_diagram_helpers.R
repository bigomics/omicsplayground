##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' PCSF diagram visual style
#'
#' Returns the style list for \code{omicsai::omicsai_diagram_render()}.
#' Encodes the PCSF board aesthetic choices: blue associations, green supports,
#' red risk edges, and box-shaped process nodes.
#'
#' @return A named list with \code{node_styles}, \code{edge_styles},
#'   and \code{show_genes_in_label}
pcsf_diagram_style <- function() {
  list(
    node_styles = list(
      phenotype = list(bg = "#FEF3C7", border = "#92400E", shape = "ellipse"),
      process   = list(bg = "#DCFCE7", border = "#166534", shape = "box"),
      module    = list(bg = "#E0F2FE", border = "#1E3A8A", shape = "box")
    ),
    edge_styles = list(
      association = list(color = "#2563EB", dashes = FALSE, arrows = "arrow", arrow_on = TRUE,  width = 2),
      supports    = list(color = "#15803D", dashes = FALSE, arrows = "arrow", arrow_on = TRUE,  width = 2),
      risk        = list(color = "#B91C1C", dashes = TRUE,  arrows = "arrow", arrow_on = TRUE,  width = 1.5)
    ),
    show_genes_in_label = FALSE
  )
}

#' Build structured PCSF diagram prompt
#'
#' Assembles a structured prompt for diagram generation using
#' \code{omicsai::diagram_prompt()} and \code{omicsai::build_prompt()}.
#' The schema (diagram/network) goes in the system section; board-specific
#' PCSF rules, species context, and the AI report go in the board section.
#'
#' Node and link type names are derived from the style registry and
#' substituted into the board rules template via \code{frag()}.
#'
#' @param report_text Character string with the AI report text
#' @param organism Character string identifying the organism (e.g. "human", "mouse")
#' @param board_root Character string path to the board.pcsf root directory
#'
#' @return Named list with \code{system} and \code{board} character elements
pcsf_build_diagram_prompt <- function(report_text, organism, board_root) {
  style <- pcsf_diagram_style()
  fmt_names <- function(nms) paste(sprintf("- `%s`", nms), collapse = "\n")
  node_names <- fmt_names(names(style$node_styles))
  link_names <- fmt_names(names(style$edge_styles))

  rules_path <- file.path(board_root, "prompts/pcsf_diagram_rules.md")

  p <- omicsai::diagram_prompt(
    role        = omicsai::frag("system_base"),
    task        = omicsai::frag("diagram/network"),
    species     = omicsai::omicsai_species_prompt(organism),
    board_rules = omicsai::frag(rules_path, list(node_names = node_names, link_names = link_names)),
    report      = paste("## AI Report\n\n", report_text)
  )
  omicsai::build_prompt(p)
}

#' Build structured PCSF image prompt
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
pcsf_build_image_prompt <- function(report_text, organism, diagram_edgelist = NULL) {
  species_img <- omicsai::omicsai_image_species_visual(organism)

  clean <- omicsai::omicsai_strip_report_noise(report_text)

  edge_text <- NULL
  if (!is.null(diagram_edgelist)) {
    dot <- omicsai::omicsai_edgelist_to_text(diagram_edgelist)
    if (nzchar(dot)) edge_text <- dot
  }

  p <- omicsai::image_prompt(
    role             = omicsai::frag("system_base"),
    task             = omicsai::frag("image/infographic", params = list(board_name = "PCSF")),
    species          = species_img,
    report           = clean,
    diagram_edgelist = edge_text
  )
  omicsai::build_prompt(p)
}
