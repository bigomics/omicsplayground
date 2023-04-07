##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

loading_table_datasets_shared_ui <- function(
  id,
  title,
  info.text,
  caption,
  height,
  width) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    width = width,
    caption = caption,
    height = height,
    title = title
  )
}

loading_table_datasets_shared_server <- function(id,
                                                 rl) {
  moduleServer(id, function(input, output, session) {

    pgxTable_DT <- reactive({
      df <- rl$pgxTableShared_data

      # need this, otherwise there is an error on user logout
      if (length(df$dataset) == 0) df <- NULL

      req(df)

      target1 <- grep("date", colnames(df))
      target2 <- grep("description", colnames(df))
      target3 <- grep("conditions", colnames(df))
      target4 <- grep("dataset", colnames(df))

      DT::datatable(
        df,
        class = "compact hover",
        rownames = TRUE,
        editable = FALSE,
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        fillContainer = TRUE,
        options = list(
          dom = "ft",
          pageLength = 9999,
          scrollX = FALSE,
          scrollY = "55vh",
          deferRender = TRUE,
          autoWidth = TRUE,
          columnDefs = list(
            list(width = "60px", targets = target1),
            list(width = "30vw", targets = target2)
          )
        ) ## end of options.list
      )
    })

    pgxTable.RENDER <- function() {
      pgxTable_DT() %>%
        DT::formatStyle(0, target = "row", fontSize = "12px", lineHeight = "95%")
    }

    pgxTable_modal.RENDER <- function() {
      pgxTable_DT() %>%
        DT::formatStyle(0, target = "row", fontSize = "16px", lineHeight = "95%")
    }

    TableModuleServer(
      "datasets",
      func = pgxTable.RENDER,
      func2 = pgxTable_modal.RENDER,
      selector = "single"
    )
  })
}
