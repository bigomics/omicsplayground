##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

wikipathview <- function(wp, val) {
  require(xml2)
  shiny::req(wp)

  isClassic <- FALSE
  url <- paste0("https://www.wikipathways.org/wikipathways-assets/pathways/", wp, "/", wp, ".svg")
  destfile <- tempfile(fileext = ".svg")
  down <- tryCatch(
    {
      download.file(url, destfile)
    },
    error = function(w) {
      tryCatch(
        {
          isClassic <<- TRUE
          url <- paste0("https://classic.wikipathways.org//wpi/wpi.php?action=downloadFile&type=svg&pwTitle=Pathway:", wp)
          download.file(url, destfile)
        },
        error = function(w) {
          return(NULL)
        }
      )
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

  # Find all 'text' elements
  label_nodes <- xml_find_all(doc, ".//text")
  labels <- xml2::xml_text(label_nodes)

  # Find the 'a' parent nodes of the label nodes
  parent_nodes <- lapply(label_nodes, xml_parent)
  parent_paths <- lapply(parent_nodes, xml_path)
  duplicated_parents <- duplicated(parent_paths) | duplicated(parent_paths, fromLast = TRUE)
  if (any(duplicated_parents)) {
    duplicated_label_nodes_indices <- which(duplicated_parents)
  } else {
    duplicated_label_nodes_indices <- NULL
  }
  # Remove duplicated parents
  if (!is.null(duplicated_label_nodes_indices)) {
    label_nodes <- label_nodes[-duplicated_label_nodes_indices]
    labels <- labels[-duplicated_label_nodes_indices]
  }
  a_nodes <- xml_parent(label_nodes)

  # Find the 'rect' children of the 'a' nodes
  if (isClassic) {
    rect_nodes <- xml_find_first(a_nodes, ".//path")
  } else {
    rect_nodes <- xml_find_first(a_nodes, ".//rect")
  }

  if (all(is.na(rect_nodes))) {
    val <- NULL
  }

  if (!is.null(val)) {
    if (sum(names(val) %in% toupper(labels)) > 0) {
      found_indexes <- which(toupper(labels) %in% names(val))
      labels <- labels[found_indexes]
      rect_nodes <- rect_nodes[found_indexes]
      val <- val[toupper(labels)]
      rr <- as.character(round(66 * pmin(1, abs(val / 2.0))**0.5))
      rr <- stringr::str_pad(rr, width = 2, pad = "0")
      colors <- ifelse(val > 0, paste0("#ff0000", rr), paste0("#0055ff", rr))
      if (isClassic) {
        current_style <- xml_attr(rect_nodes, "style")
        new_style <- paste0(current_style, " fill:", colors, ";")
        lapply(seq_along(new_style), function(x) {
          xml_set_attr(rect_nodes[x], "style", new_style[x])
        })
      } else {
        xml_attr(rect_nodes, "fill") <- colors
      }
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
