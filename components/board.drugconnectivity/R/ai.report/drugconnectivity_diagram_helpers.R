##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

drugconnectivity_diagram_style <- function() {
  list(
    node_styles = list(
      contrast = list(bg = "#FFFFE0", border = "#B8860B", shape = "ellipse"),
      drug     = list(bg = "#ADD8E6", border = "#1A5276", shape = "box"),
      moa      = list(bg = "#DDA0DD", border = "#6C3483", shape = "diamond"),
      target   = list(bg = "#90EE90", border = "#1E8449", shape = "box")
    ),
    edge_styles = list(
      opposes     = list(color = "#2E8B57", dashes = FALSE, arrows = "arrow", arrow_on = TRUE,  width = 2),
      mimics      = list(color = "#C0392B", dashes = FALSE, arrows = "arrow", arrow_on = TRUE,  width = 2),
      has_moa     = list(color = "#6C3483", dashes = TRUE,  arrows = "arrow", arrow_on = TRUE,  width = 1.5),
      hits        = list(color = "#1E8449", dashes = TRUE,  arrows = "arrow", arrow_on = TRUE,  width = 1.5),
      association = list(color = "#94A3B8", dashes = TRUE,  arrows = "arrow", arrow_on = FALSE, width = 1.5)
    ),
    show_genes_in_label = FALSE
  )
}

drugconnectivity_build_diagram_prompt <- function(report_text, organism, board_root) {
  style <- drugconnectivity_diagram_style()
  fmt_names <- function(nms) paste(sprintf("- `%s`", nms), collapse = "\n")
  node_names <- fmt_names(names(style$node_styles))
  link_names <- fmt_names(names(style$edge_styles))

  layers <- list()
  layers[[1]] <- omicsai::omicsai_instructions("diagram_network")
  layers[[2]] <- tryCatch({
    tpl <- omicsai::omicsai_load_template("prompts/diagram_drugconnectivity_rules.md", root = board_root)
    omicsai::omicsai_substitute_template(tpl, list(node_names = node_names, link_names = link_names))
  }, error = function(e) "")
  layers[[3]] <- tryCatch(omicsai::omicsai_species_prompt(organism), error = function(e) "")
  layers[[4]] <- paste("## AI Report\n\n", report_text)
  layers <- layers[nzchar(layers)]
  paste(layers, collapse = "\n\n---\n\n")
}

drugconnectivity_build_image_prompt <- function(report_text, organism, diagram_edgelist = NULL) {
  layers <- list()
  layers[[1]] <- omicsai::omicsai_image_species_visual(organism)
  layers[[2]] <- paste(
    "<report>",
    omicsai::omicsai_strip_report_noise(report_text),
    "</report>",
    sep = "\n"
  )

  if (!is.null(diagram_edgelist)) {
    edges <- omicsai::omicsai_edgelist_to_text(diagram_edgelist)
    if (nzchar(edges)) {
      layers[[3]] <- paste("<diagram>", edges, "</diagram>", sep = "\n")
    }
  }

  layers <- layers[nzchar(layers)]
  paste(layers, collapse = "\n\n")
}
