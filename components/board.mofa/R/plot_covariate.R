##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_covariate_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width  = 400) {

  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(
      inputId = ns("expand_conditions"),
      label = "Expand conditions",
      value = TRUE
    )
  )

  
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

mofa_plot_covariate_server <- function(id,
                                       mofa,
                                       watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    
    plot.RENDER <- function() {
      res <- mofa()
      do_collapse <- !input$expand_conditions
      playbase::mofa.plot_covariate_correlation(
        res, collapse = do_collapse)
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(85, 100),
      add.watermark = watermark
    )
  })
}
