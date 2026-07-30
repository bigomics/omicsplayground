## Loadings board: vertical top-loading barplots + PC-split expression heatmap.

PCAexplorer_loadings_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_columns(
    col_widths = c(6, 6),
    bslib::navset_card_tab(
      full_screen = TRUE,
      title = "Loadings",
      bslib::nav_spacer(),
      bslib::nav_panel(
        "plot",
        shiny::plotOutput(ns("loadings"), height = "800px")
      ),
      bslib::nav_panel(
        "table",
        shiny::dataTableOutput(ns("loadings_table"))
      )
    ),
    bslib::card(
      full_screen = TRUE,
      shiny::plotOutput(ns("loadings_heatmap"), height = "800px")
    )
  )
}

PCAexplorer_loadings_module <- function(id, get_pcaX) {
  shiny::moduleServer(
    id,
    function(input, output, session) {

      output$loadings <- shiny::renderPlot({
        res <- get_pcaX()
        shiny::req(res)
        pcaexplorer.plot_loadings_bars(res, npc = 4, ntop = 40)
      })

      ## Full feature × PC loading matrix (res$u)
      output$loadings_table <- shiny::renderDataTable({
        res <- get_pcaX()
        shiny::req(res, res$u)
        U <- as.matrix(res$u)
        if (ncol(U) >= 1L) {
          U <- U[order(-abs(U[, 1])), , drop = FALSE]
        }
        dt <- round(U, 4)

        DT::datatable(
          data = dt,
          class = "compact hover",
          rownames = TRUE,
          extensions = c("Buttons", "Scroller"),
          plugins = "scrollResize",
          selection = list(
            mode = "single",
            target = "row",
            selected = 1
          ),
          options = list(
            dom = "lfrtipB",
            buttons = c('copy', 'csv', 'excel'),
            scroller = TRUE,
            scrollX = TRUE,
            scrollY = "700px",
            scrollResize = TRUE,
            deferRender = TRUE,
            autoWidth = TRUE
          )
        )

      })

      output$loadings_heatmap <- shiny::renderPlot({
        res <- get_pcaX()
        shiny::req(res)
        pcaexplorer.plot_loadings_heatmap(res, npc = 4, ntop = 15)
      })

    }
  )
}
