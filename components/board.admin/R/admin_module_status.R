##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


admin_module_status_ui <- function(
  id,
  label = "",
  title,
  height,
  width = c("auto", "100%"),
  caption,
  info.text
) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("mod"),
    title = title,
    label = label,
    outputFunc = htmlOutput,
    outputFunc2 = htmlOutput,
    info.text = info.text,
    caption = caption,
    caption2 = NULL,
    options = NULL,
    download.fmt = NULL,
    width = width,
    height = height
  )
}

admin_module_status_server <- function(id, auth) {
  moduleServer(id, function(input, output, session) {
    status_data <- shiny::reactive({
      shiny::req(isTRUE(auth$ADMIN))

      r_version <- paste0(R.version$major, ".", R.version$minor)
      platform <- R.version$platform

      res <- paste0(
        "<p><b>R Version: </b>", r_version, "</p>",
        "<p><b>Platform: </b>", platform, "</p>",
        "<p><b>Working Directory: </b>", getwd(), "</p>"
      )
      res
    })

    status.RENDER <- function() {
      res <- status_data()
      div(shiny::HTML(res), class = "admin-status")
    }

    modal_status.RENDER <- function() {
      res <- status_data()
      div(shiny::HTML(res), class = "admin-status", style = "font-size:1.3em;")
    }

    PlotModuleServer(
      "mod",
      plotlib = "generic",
      plotlib2 = "generic",
      func = status.RENDER,
      func2 = modal_status.RENDER,
      renderFunc = shiny::renderUI,
      renderFunc2 = shiny::renderUI
    )
  }) ## end of moduleServer
}
