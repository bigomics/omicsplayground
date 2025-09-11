##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_lasagna3D_ui <- function(
    id,
    title = "",
    info.text = "",
    info.references = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::checkboxInput(ns("drawlines"), "draw connections", FALSE)
  )

  PlotModuleUI(
    ns("plot"),
    plotlib = "plotly",
    title = title,
    options = options,
    label = label,
    info.text = info.text,
    info.references = info.references,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

mofa_plot_lasagna3D_server <- function(id,
                                       data,
                                       pgx,
                                       watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot.RENDER <- function() {
      res <- data()
      shiny::req(res$posf)

      graph <- res$graph
      vars <- igraph::V(graph)$value
      names(vars) <- igraph::V(graph)$name

      posf <- playbase::mofa.prefix(res$posf)
      posf <- lapply(posf, function(x) x[(rownames(x) %in% names(vars)), , drop = FALSE])
      edges <- NULL
      if (input$drawlines) {
        edges <- data.frame(
          igraph::as_edgelist(graph),
          weight = igraph::E(graph)$weight
        )
      }

      plt <- playbase::plotly_lasagna(
        pos = posf,
        vars = vars,
        min.rho = 0.01,
        edges = edges,
        num_edges = 20
      )

      plt %>%
        plotly::layout(
          margin = list(l = 10, r = 10, b = 10, t = 10),
          scene = list(
            camera = list(eye = list(x = 2.2, y = 0.8, z = 1))
          )
        )
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      plotlib = "plotly",
      pdf.width = 8, pdf.height = 8,
      res = c(80, 100),
      add.watermark = watermark
    )
  })
}
