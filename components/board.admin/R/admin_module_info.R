##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


admin_module_info_ui <- function(
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

admin_module_info_server <- function(id, auth) {
  moduleServer(id, function(input, output, session) {
    info_data <- shiny::reactive({
      shiny::req(isTRUE(auth$ADMIN))

      res <- paste0(
        "<p><b>Current Admin: </b>", auth$username, "</p>",
        "<p><b>Email: </b>", auth$email, "</p>",
        "<p><b>Admin Status: </b>",
        "<span class='badge bg-success'>Active</span></p>"
      )
      res
    })

    info.RENDER <- function() {
      res <- info_data()
      div(shiny::HTML(res), class = "admin-info")
    }

    modal_info.RENDER <- function() {
      res <- info_data()
      div(shiny::HTML(res), class = "admin-info", style = "font-size:1.3em;")
    }

    PlotModuleServer(
      "mod",
      plotlib = "generic",
      plotlib2 = "generic",
      func = info.RENDER,
      func2 = modal_info.RENDER,
      renderFunc = shiny::renderUI,
      renderFunc2 = shiny::renderUI
    )
  }) ## end of moduleServer
}
