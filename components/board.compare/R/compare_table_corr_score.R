##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

compare_table_corr_score_ui <- function(id, label = "", height = c(600, 800)) {
  ns <- shiny::NS(id)

  tableWidget(ns("table"))
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
          scrollY = "70vh",
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    })

    score_table_info <- "In this table, users can check mean expression values of features across the conditions for the selected genes."

    score_table <- shiny::callModule(
      tableModule,
      id = "table",
      func = score_table.RENDER,
      info.text = score_table_info,
      title = tags$div(
        HTML('<span class="module-label">(b)</span>Correlation score')
      ),
      height = c(235, 750),
      width = c("auto", 1600)
    )
    return(score_table)
  })
}
