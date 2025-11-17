##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_module_heatmap_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = "",
    height,
    width) {
  ns <- shiny::NS(id)

  opts <- shiny::tagList(
    shiny::radioButtons(ns("plottype"), "Plottype:",
      choices=c("expression","correlation"),
      inline = TRUE),
    shiny::checkboxInput(ns("showtop"), "Show top genes", FALSE)
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    options = opts,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

wgcna_plot_module_heatmap_server <- function(id,
                                             wgcna,
                                             pgx,
                                             selected_module,
                                             selected_trait,
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    render_plot <- function(ntop = 30, xmar=FALSE) {
      res <- wgcna()
      module <- selected_module()
      shiny::req(!is.null(module) & module != "")
      nmax = ntop
      if(!input$showtop) nmax <- -1
      heatmap.mar = c(6,7)
      if(xmar) heatmap.mar = c(8,12)
      par(mar=c(2,2,2,8))
      playbase::wgcna.plotModuleHeatmap(
        wgcna = res,
        module = module,
        type = input$plottype,
        heatmap.mar = heatmap.mar,
        nmax = nmax,
        main = ""
      )
    }

    RENDER <- function() {
      render_plot(ntop = 20, xmar=FALSE)
    }

    RENDER2 <- function() {
      render_plot(ntop = 30, xmar=TRUE)
    }

    PlotModuleServer(
      "plot",
      func = RENDER,
      func2 = RENDER2,
      pdf.width = 8,
      pdf.height = 5,
      res = c(75, 100),
      add.watermark = watermark
    )
  })
}
