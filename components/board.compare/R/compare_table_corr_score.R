##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
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
                                            getScoreTable,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    score_table.RENDER <- shiny::reactive({
      df <- getScoreTable()

      shiny::req(df)
      numeric.cols <- 2:ncol(df)
      if ("title" %in% colnames(df)) df$title <- substring(df$title, 1, 60)

      DT::datatable(
        df,
        rownames = TRUE, #
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          scrollX = TRUE,
          scrollY = "15vh",
          scrollResize = TRUE,
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
