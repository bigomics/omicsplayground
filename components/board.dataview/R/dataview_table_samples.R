##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


dataview_table_samples_ui <- function(
  id,
  width,
  height,
  title,
  info.text,
  caption
) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("datasets"),
    width = width,
    height = height,
    title = title,
    info.text = info.text,
    caption = caption,
    label = "c"
  )
}

dataview_table_samples_server <- function(id,
                                          pgx,
                                          r.samples = reactive(""),
                                          scrollY) {
  moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      shiny::req(pgx$Y, pgx$samples, r.samples())
      samples <- r.samples()
      dt <- pgx$samples[samples, , drop = FALSE]
      dt
    })

    table.RENDER <- function() {
      dt <- table_data()
      req(dt)
      DT::datatable(dt,
        class = "compact hover",
        rownames = TRUE,
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = 1),
        options = list(
          dom = "lfrtip",
          scroller = TRUE,
          scrollX = TRUE,
          scrollY = scrollY,
          scrollResize = TRUE,
          deferRender = TRUE
        )
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    modal_table.RENDER <- function() {
      dt <- table_data()
      req(dt)
      DT::datatable(dt,
        class = "compact hover",
        rownames = TRUE,
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        options = list(
          dom = "lfrtip",
          scroller = TRUE,
          scrollX = TRUE,
          scrollY = SCROLLY_MODAL,
          deferRender = TRUE
        )
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "20px", lineHeight = "70%")
    }

    TableModuleServer(
      "datasets",
      func = table.RENDER,
      func2 = modal_table.RENDER,
      selector = "none"
    )
  }) ## end of moduleServer
} ## end of server
