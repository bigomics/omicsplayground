## PCA pairs board: scatterpair matrix, PC boxplots, F-test vs PC dimension.

PCAexplorer_pairs_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_columns(
    col_widths = c(7, 5),
    bslib::card(full_screen = TRUE, shiny::plotOutput(ns("scatterpairs"), height = "800px")),
    bslib::layout_columns(
      col_widths = 12,
      row_heights = c(1, 1),
      bslib::card(full_screen = TRUE, shiny::plotOutput(ns("pc_boxplots"))),
      bslib::card(full_screen = TRUE, shiny::plotOutput(ns("ftest_vs_pcdim")))
    )
  )
}

PCAexplorer_pairs_module <- function(id, get_pcaX, colorby, ellipse) {
  shiny::moduleServer(
    id,
    function(input, output, session) {

      output$scatterpairs <- shiny::renderPlot({
        res <- get_pcaX()
        ph <- colorby()
        shiny::req(res, ph)
        pcaexplorer.plot_scatterpairs(
          res, colorby_var = ph,
          add.ellipse = isTRUE(ellipse())
        )
      })

      output$pc_boxplots <- shiny::renderPlot({
        res <- get_pcaX()
        ph <- colorby()
        shiny::req(res, ph)
        pcaexplorer.plot_pc_boxplots(res, colorby_var = ph, npc = 4)
      })

      output$ftest_vs_pcdim <- shiny::renderPlot({
        res <- get_pcaX()
        shiny::req(res)
        graphics::par(mar = c(4, 4, 3, 1))
        pcaexplorer.plot_ftest_vs_pcdim(res)
      })

    }
  )
}
