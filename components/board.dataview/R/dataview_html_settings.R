##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


dataview_html_settings_ui <- function(
  id,
  label = "",
  title,
  height,
  width,
  caption
) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("mod"),
    title = title,
    label = label,
    outputFunc = shiny::htmlOutput,
    outputFunc2 = shiny::htmlOutput,
    caption = caption,
    caption2 = NULL,
    options = NULL,
    download.fmt = NULL,
    width = width,
    height = height
  )
}

dataview_html_settings_server <- function(id,
                                          pgx,
                                          watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    settings_data <- shiny::reactive({
      fields <- c("name", "description", "datatype", "date", "settings",
        "omicsplayground_version")
      info <- playbase::pgx.info(pgx, fields=NULL, format="html")
      return(info)
    })


    info.RENDER <- function() {      
      info <- settings_data()
      div(shiny::HTML(info), class = "gene-info",
        style="font-size: 1em; line-height: 1.3em;")
    }

    PlotModuleServer(
      "mod",
      plotlib = "generic",
      func = info.RENDER,
      renderFunc = shiny::renderUI
    )
  })
}
