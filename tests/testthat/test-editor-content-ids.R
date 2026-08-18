## test-editor-content-ids.R
##
## Regression guard: getEditorContent() must never emit the same input id
## twice, otherwise the browser console fills with
## "[shiny] Duplicate input IDs were found".

.ui_dir <- if (dir.exists("components/ui")) "components/ui" else "../../components/ui"

source(file.path(.ui_dir, "ui-ColorDefaults.R"), encoding = "UTF-8", local = TRUE)
source(file.path(.ui_dir, "ui-EditorContent.R"), encoding = "UTF-8", local = TRUE)

.plot_types <- c(
  "volcano", "heatmap", "barplot", "expression_barplot", "expression_boxplot",
  "correlation", "scatterplot", "featuremap", "enrichment", "clustering",
  "clustering_prism", "grouped_barplot", "gradient", "significance",
  "scatter_highlight", "rank_plot", "correlation_matrix", "scatter_updown",
  "boxplot_methyl"
)

for (.ty in .plot_types) {
  local({
    ty <- .ty
    test_that(paste0("editor content for '", ty, "' has no duplicate ids"), {
      ui <- getEditorContent(
        ty,
        ns = shiny::NS(paste0("board-plot_", ty)),
        ns_parent = shiny::NS("board"),
        title = ty,
        outputFunc = plotly::plotlyOutput
      )
      html <- as.character(ui)
      ids <- regmatches(html, gregexpr('id="[^"]+"', html))[[1]]
      ids <- sub('^id="', "", sub('"$', "", ids))
      dups <- names(which(table(ids) > 1))
      expect_equal(dups, character(0))
    })
  })
}
