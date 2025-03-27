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
    download.fmt = c("png", "pdf", "svg")
  )
}

plot_deepnet_diagram_server <- function(id,
                                        net,
                                        update,
                                        watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    svgfile <- NULL    
    
    plot.RENDER <- eventReactive({
      list( update(), net() )
    },{      

      if(update()==TRUE || is.null(svgfile)) {
        net <- net() ## do not react everytime      
        progress <- shiny::Progress$new(session, min=0, max=1)
        on.exit(progress$close())
        progress$set(message = paste("Creating diagram..."), value = 0.33)
        svgfile <<- playbase::deep.plotNeuralNet(net, svgfile=NULL)
        update(FALSE)
      }

      dbg("[plot_deepnet_diagram_server] svgfile =", svgfile)

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
    })

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
