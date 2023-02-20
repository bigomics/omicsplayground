##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

compare_table_corr_score_ui <- function(id, width, height) {
  ns <- shiny::NS(id)

  info_text <- "In this table, users can check mean expression values of features across the conditions for the selected genes."

  TableModuleUI(
    ns("datasets"),
    info.text = info_text,
    width = width,
    height = height,
    title = "Correlation score",
    label = "b"
  )

}

compare_table_corr_score_server <- function(id,
                                            getOmicsScoreTable,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    score_table.RENDER <- shiny::reactive({
      df <- getOmicsScoreTable()
      if (is.null(df)) {
        return(NULL)
      }
      shiny::req(df)
      numeric.cols <- 2:ncol(df)

      DT::datatable(
        df,
        rownames = TRUE, ## escape = c(-1,-2),
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          scrollX = TRUE,
          scrollY = "15vh",
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    })

    score_table.RENDER_modal <- shiny::reactive({
      dt <- score_table.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    score_table <- TableModuleServer(
      "datasets",
      func = score_table.RENDER,
      func2 = score_table.RENDER_modal,
      selector = "none"
    )
    return(score_table)
  })
}
