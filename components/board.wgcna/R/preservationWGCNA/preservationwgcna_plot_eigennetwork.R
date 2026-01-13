##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

preservationWGCNA_plot_eigenNetwork_ui <- function(
  id,
  title = "",
  info.text = "",
  caption = "",
  label = "",
  height = 400,
  width = 400
) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(
      inputId = ns("top20"),
      label = "Show top 20",
      value = FALSE
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    # options = options,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

preservationWGCNA_plot_eigenNetwork_server <- function(id,
                                                       rwgcna) {
  moduleServer(id, function(input, output, session) {
    plot.RENDER <- function() {
      res <- rwgcna()
      shiny::req(res)

      par(mar = c(0, 3, 3, 1))
      WGCNA::plotEigengeneNetworks(
        res$net$multiMEs,
        setLabels = names(res$net$multiMEs),
        marHeatmap = c(1, 3, 3, 1)
      )
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 8,
      pdf.height = 10,
      res = c(65, 90),
      add.watermark = FALSE
    )
  })
}
