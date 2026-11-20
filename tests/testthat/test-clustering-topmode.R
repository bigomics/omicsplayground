## test-clustering-topmode.R
##
## Regression guard: the Cluster Samples (heatmap) "Top mode" selector in
## Advanced options must keep the short option names and move the long
## descriptions to the tooltip. The server matches on the plain values
## ("sd", "pca", "marker").

comps_dir <- if (dir.exists("components/ui")) "components" else "../../components"

for (f in list.files(file.path(comps_dir, "ui"), pattern = "\\.R$", full.names = TRUE)) {
  tryCatch(source(f, encoding = "UTF-8", local = TRUE), error = function(e) NULL)
}
source(file.path(comps_dir, "board.clustering", "R", "clustering_ui.R"),
  encoding = "UTF-8", local = TRUE
)

## Find the "Top mode" control by walking the rendered tag tree: the input
## container that carries the tooltip title and holds the hm_topmode select.
find_topmode_container <- function(node) {
  if (inherits(node, "shiny.tag.list")) {
    for (child in as.list(node)) {
      found <- find_topmode_container(child)
      if (!is.null(found)) {
        return(found)
      }
    }
    return(NULL)
  }
  if (!inherits(node, "shiny.tag")) {
    return(NULL)
  }
  attrs <- node$attribs
  is_container <- grepl("shiny-input-container", attrs[["class"]], fixed = TRUE)
  has_select <- !is.null(htmltools::tagQuery(node)$find("select")$length())
  if (is_container && has_select && check_id(node)) {
    return(node)
  }
  for (child in node$children) {
    found <- find_topmode_container(child)
    if (!is.null(found)) {
      return(found)
    }
  }
  NULL
}

check_id <- function(node) {
  ## the select inside belongs to hm_topmode
  sel <- htmltools::tagQuery(node)$find("select")$selectedTags()
  any(vapply(sel, function(s) grepl("hm_topmode", s$attribs[["id"]], fixed = TRUE), logical(1)))
}

topmode_control <- function() {
  selects <- htmltools::tagQuery(ClusteringInputs("x"))$find("select")$selectedTags()
  is_topmode <- vapply(selects, function(s) {
    grepl("hm_topmode", s$attribs[["id"]], fixed = TRUE)
  }, logical(1))
  if (sum(is_topmode) != 1) {
    return(list(title = NULL, options = NULL))
  }
  sel <- selects[[which(is_topmode)]]

  ## options are rendered as a single HTML chunk inside the select
  opts_html <- paste(unlist(sel$children), collapse = "")
  options <- regmatches(opts_html, gregexpr("<option[^>]*>[^<]*</option>", opts_html))[[1]]
  list(
    title = topmode_tooltip_title(),
    options = data.frame(
      value = sub("<option value=\"([^\"]*)\".*", "\\1", options),
      label = sub(".*>([^<]*)</option>", "\\1", options),
      stringsAsFactors = FALSE
    )
  )
}

## tooltip title of the input container that wraps the Top mode select
topmode_tooltip_title <- function() {
  ctl <- find_topmode_container(ClusteringInputs("x"))
  if (is.null(ctl)) NULL else ctl$attribs[["title"]]
}

test_that("Top mode options stay short and keep their server values", {
  ctl <- topmode_control()
  expect_false(is.null(ctl$options))
  expect_setequal(ctl$options$value, c("sd", "pca", "marker"))

  ## short labels only, e.g. no "sd (high variance)" or "pca (top components)"
  expect_false(any(grepl("(", ctl$options$label, fixed = TRUE)))
  expect_true(all(nchar(ctl$options$label) <= 10))
})

test_that("Top mode descriptions are in the tooltip", {
  title <- topmode_tooltip_title()
  expect_type(title, "character")
  expect_match(title, "Specify the criteria for selecting top features")
  for (criterion in c("standard deviation", "principal components", "overexpressed")) {
    expect_match(title, criterion)
  }
})
