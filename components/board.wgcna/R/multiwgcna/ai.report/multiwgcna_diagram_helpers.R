##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

multiwgcna_diagram_style <- function() {
  list(
    node_styles = list(
      layer = list(bg = "#E0F2FE", border = "#1E3A8A", shape = "box"),
      module = list(bg = "#DCFCE7", border = "#166534", shape = "ellipse"),
      phenotype = list(bg = "#FEF3C7", border = "#92400E", shape = "box"),
      process = list(bg = "#FCE7F3", border = "#9D174D", shape = "diamond")
    ),
    edge_styles = list(
      positive = list(color = "#2E8B57", dashes = FALSE, arrows = "arrow", arrow_on = TRUE, width = 2),
      negative = list(color = "#C0392B", dashes = FALSE, arrows = "bar", arrow_on = TRUE, width = 2),
      cross_layer = list(color = "#2563EB", dashes = FALSE, arrows = "arrow", arrow_on = TRUE, width = 2),
      association = list(color = "#94A3B8", dashes = TRUE, arrows = "arrow", arrow_on = FALSE, width = 1.5)
    ),
    show_genes_in_label = FALSE
  )
}

multiwgcna_build_diagram_prompt <- function(report_text, organism, board_root, data_tables = NULL) {
  style <- multiwgcna_diagram_style()
  fmt_names <- function(nms) paste(sprintf("- `%s`", nms), collapse = "\n")
  node_names <- fmt_names(names(style$node_styles))
  link_names <- fmt_names(names(style$edge_styles))

  rules_path <- file.path(board_root, "prompts/multiwgcna/multiwgcna_diagram_rules.md")

  p <- omicsai::diagram_prompt(
    role = omicsai::frag("system_base"),
    task = omicsai::frag("diagram/network"),
    species = omicsai::omicsai_species_prompt(organism),
    board_rules = omicsai::frag(rules_path, list(node_names = node_names, link_names = link_names)),
    report = paste("## AI Report\n\n", report_text),
    data = data_tables
  )
  omicsai::build_prompt(p)
}

multiwgcna_build_image_prompt <- function(report_text, organism,
                                          diagram_edgelist = NULL,
                                          style_name = NULL,
                                          n_blocks = 1L) {
  species_img <- omicsai::omicsai_image_species_visual(organism)
  clean <- omicsai::omicsai_strip_report_noise(report_text)
  clean <- gsub("\\b([a-z]{1,3})::ME(\\w+)\\b", "Layer \\1 Module \\2", clean, perl = TRUE)

  edge_text <- NULL
  if (!is.null(diagram_edgelist)) {
    dot <- omicsai::omicsai_edgelist_to_text(diagram_edgelist)
    if (nzchar(dot)) edge_text <- dot
  }

  style_frag <- if (!is.null(style_name) && nzchar(style_name)) {
    omicsai::frag(paste0("image/styles/", style_name))
  }
  blocks_frag <- omicsai::frag(paste0("image/blocks/blocks_", as.integer(n_blocks)))

  p <- omicsai::image_prompt(
    role = omicsai::frag("system_base"),
    task = omicsai::frag("image/infographic", params = list(board_name = "MultiWGCNA")),
    species = species_img,
    style = style_frag,
    blocks = blocks_frag,
    report = clean,
    diagram_edgelist = edge_text
  )
  omicsai::build_prompt(p)
}
