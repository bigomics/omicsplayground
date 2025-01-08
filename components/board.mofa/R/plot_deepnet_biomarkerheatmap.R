##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

plot_deepnet_biomarkerheatmap_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = c("100%", TABLE_HEIGHT_MODAL),
    width = c("auto", "100%")
    ) {
  ns <- shiny::NS(id)
  
  options <- tagList(
    shiny::checkboxInput(ns("onlypositive"), "positive only", TRUE)
  )
  
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    options = options,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

plot_deepnet_biomarkerheatmap_server <- function(id,
                                                 net,
                                                 update,
                                                 datatypes,
                                                 watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function(n=12) {
      update()  ## react on updates
      net <- net()
      datatypes <- datatypes()
      playbase::deep.plotBiomarkerHeatmap(
        net, ntop = 24, datatypes = datatypes,
        cexCol = 0.8, show_colnames = FALSE
      ) 
    }

    plot.RENDER2 <- function(n=12) {
      update()  ## react on updates
      net <- net()
      datatypes <- datatypes()
      nsamples <- ncol(net$X[[1]])
      playbase::deep.plotBiomarkerHeatmap(
        net, ntop = 40, datatypes = datatypes,
        rownames_width = 80, rowlab.maxlen = 60,
        cexCol = 0.8, show_colnames = (nsamples<100)
      ) 
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      func2 = plot.RENDER2,      
      pdf.width = 12, pdf.height = 5,
      res = c(90, 120),
      add.watermark = watermark
    )


  })
}
