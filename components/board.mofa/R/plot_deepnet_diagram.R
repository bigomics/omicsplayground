##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

plot_deepnet_diagram_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)
  
  PlotModuleUI(
    ns("plot"),
    plotlib = "svgPanZoom",
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

plot_deepnet_diagram_server <- function(id,
                                        net,
                                        pgx,
                                        update,
                                        watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function(n=12) {
      update()
      shiny::req(pgx$X)  ## react on pgx change
      net <- isolate(net()) ## do not react everytime
      shiny::req(net)

      dbg("[deepnet_diagram_server] reacted!")
      dbg("[deepnet_diagram_server] names(net$xx) = ", names(net$xx))
      
      svgfile <- playbase::deep.plotNeuralNet(net, outfile=NULL)
      validate(
        need(!is.null(svgfile), "Could not create model diagram")
      )
      img.svg <- readChar(svgfile, nchars = file.info(svgfile)$size)
      pz <- svgPanZoom::svgPanZoom(
        img.svg,
        controlIconsEnabled = TRUE,
        zoomScaleSensitivity = 0.4,
        minZoom = 1,
        maxZoom = 5,
        viewBox = FALSE
      )
      return(pz)
    }


    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      plotlib = "svgPanZoom",
      pdf.width = 12, pdf.height = 5,
      res = c(90, 130),
      add.watermark = watermark
    )


  })
}
