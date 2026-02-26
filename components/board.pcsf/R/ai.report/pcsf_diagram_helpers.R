##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

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

pcsf_build_diagram_prompt <- function(report_text, organism, board_root) {
  style <- pcsf_diagram_style()
  fmt_names <- function(nms) paste(sprintf("- `%s`", nms), collapse = "\n")
  node_names <- fmt_names(names(style$node_styles))
  link_names <- fmt_names(names(style$edge_styles))

  layers <- list()
  layers[[1]] <- omicsai::omicsai_instructions("diagram_network")
  layers[[2]] <- tryCatch({
    tpl <- omicsai::omicsai_load_template("prompts/diagram_pcsf_rules.md", root = board_root)
    omicsai::omicsai_substitute_template(tpl, list(node_names = node_names, link_names = link_names))
  }, error = function(e) "")
  layers[[3]] <- tryCatch(omicsai::omicsai_species_prompt(organism), error = function(e) "")
  layers[[4]] <- paste("## AI Report\n\n", report_text)
  layers <- layers[nzchar(layers)]
  paste(layers, collapse = "\n\n---\n\n")
}

pcsf_build_image_prompt <- function(report_text, organism, diagram_edgelist = NULL) {
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
