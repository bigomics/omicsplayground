##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_lasagna_ui <- function(
    id,
    title = "",
    info.text = "",
    info.references = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plot"),
    plotlib = "plotly",
    title = title,
    label = label,
    info.text = info.text,
    info.references = info.references,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

mofa_plot_lasagna_server <- function(id,
                                     data,
                                     pgx,
                                     input_contrast,
                                     watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function() {
      res <- data()
      shiny::req(res$posf)
      k <- input_contrast()
      shiny::req(k, pgx$X)
      vars <- playbase::pgx.getMetaMatrix(pgx)$fc[,k]
      vars <- sign(vars) * abs(vars)**0.5
      
      plt <- playbase::plotly_lasagna(
        pos = res$posf,
        vars = vars,
        X = pgx$X,
        min.rho = 0.01,
        num_edges = 40
      )

      plt %>%
        plotly::layout(
          margin = list(l=10, r=10, b=10, t=10),
          scene = list(
            camera = list(eye = list(x=2.2, y=0.8, z=1))
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
