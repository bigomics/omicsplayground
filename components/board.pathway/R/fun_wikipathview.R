#

wikipathview <- function(wp, val, dir) {
  require(xml2)
  require(fluctuator)

  url <- paste0("https://www.wikipathways.org/wikipathways-assets/pathways/", wp, "/", wp, ".svg")
  destfile <- tempfile(fileext = ".svg")
  down <- tryCatch(
    {
      download.file(url, destfile)
    },
    error = function(w) {
      return(NULL)
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
  a_nodes <- xml_parent(label_nodes)

  # Find the 'rect' children of the 'a' nodes
  rect_nodes <- xml_find_first(a_nodes, ".//rect")
  if (all(is.na(rect_nodes))) {
    val <- NULL
  }

  if (!is.null(val)) {
    if (sum(names(val) %in% labels) == 0) {
      return(NULL)
    }
    found_indexes <- which(labels %in% names(val))
    labels <- labels[found_indexes]
    rect_nodes <- rect_nodes[found_indexes]
    val <- val[labels]
    rr <- as.character(round(66 * pmin(1, abs(val / 2.0))**0.5))
    rr <- stringr::str_pad(rr, width = 2, pad = "0")
    colors <- ifelse(val > 0, paste0("#ff0000", rr), paste0("#0055ff", rr))
    xml_attr(rect_nodes, "fill") <- colors
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

  svg <- fluctuator::read_svg(destfile)

  return(svg)
}
