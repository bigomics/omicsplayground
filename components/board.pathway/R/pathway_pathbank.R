##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

pathbankview <- function(pb, val) {
  require(xml2)
  shiny::req(pb)

  url <- paste0("https://www.pathbank.org/view/", pb, "/download?type=simple_vector_image")
  destfile <- tempfile(fileext = ".svg")
  down <- tryCatch(
    {
      download.file(url, destfile)
    },
    error = function(w) {
      NULL
    }
  ) |> is.null()

  if (down) {
    return(NULL)
  }

  # Read the file line by line
  lines <- readLines(destfile)

  # Use gsub to replace the line
  lines <- gsub(
    'xmlns="http://www.w3.org/2000/svg"',
    'xmlns:svg="http://www.w3.org/2000/svg"',
    lines
  )

  # Write the lines back to the file
  writeLines(lines, destfile)

  # Load the SVG
  doc <- read_xml(destfile)

  label_nodes <- xml_find_all(doc, ".//text")
  labels <- xml2::xml_text(label_nodes)

  # Find nodes
  g_nodes <- xml_find_all(doc, ".//g")
  g_nodes_labels <- xml_attr(g_nodes, "data-element-id")

  # Find the 'rect' childs nodes of the g_nodes
  rect_nodes <- xml_find_first(g_nodes, ".//circle")

  if (!is.null(val)) {
    if (sum(names(val) %in% g_nodes_labels) > 0) {
      found_indexes <- which(g_nodes_labels %in% names(val))
      g_nodes_labels <- g_nodes_labels[found_indexes]
      rect_nodes <- rect_nodes[found_indexes]
      val <- val[g_nodes_labels]
      rr <- as.character(round(66 * pmin(1, abs(val / 2.0))**0.5))
      rr <- stringr::str_pad(rr, width = 2, pad = "0")
      colors <- ifelse(val > 0, paste0("#ff0000", rr), paste0("#0055ff", rr))
      xml_attr(rect_nodes, "fill") <- colors
    }
  }

  write_xml(doc, destfile)

  # Read the file line by line
  lines <- readLines(destfile)

  # Use gsub to replace the line
  lines <- gsub(
    'xmlns:svg="http://www.w3.org/2000/svg"',
    'xmlns="http://www.w3.org/2000/svg" xmlns:svg="http://www.w3.org/2000/svg"',
    lines
  )

  # Write the lines back to the file
  writeLines(lines, destfile)

  return(destfile)
}
