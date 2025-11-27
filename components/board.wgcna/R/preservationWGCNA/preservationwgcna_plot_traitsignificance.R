##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

preservationWGCNA_plot_traitsignificance_ui <- function(
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
    shiny::selectInput(
      inputId = ns("ntop"),
      label = "Number of plots",
      choices = c(1, 4, 6, 9, 12, 16),
      selected = 12
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    options = options,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

preservationWGCNA_plot_traitsignificance_server <- function(id,
                                                            rwgcna,
                                                            rtrait) {
  moduleServer(id, function(input, output, session) {
    plot.RENDER <- function(format = 1) {
      res <- rwgcna()
      trait <- rtrait()

      shiny::req(res)
      shiny::req(trait %in% colnames(res$datTraits))

      playbase::wgcna.plotTopModules_multi(
        res,
        nmax = as.integer(input$ntop),
        trait = trait,
        setpar = format
      )
    }

    plot1 <- function() {
      plot.RENDER(format = 1)
    }
    plot2 <- function() {
      plot.RENDER(format = 2)
    }

    PlotModuleServer(
      "plot",
      func = plot1,
      func2 = plot2,
      pdf.width = 8,
      pdf.height = 12,
      res = c(85, 100),
      add.watermark = FALSE
    )
  })
}
