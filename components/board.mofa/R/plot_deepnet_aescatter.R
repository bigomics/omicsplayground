##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

plot_deepnet_aescatter_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = c("100%", TABLE_HEIGHT_MODAL),
    width = c("auto", "100%"))
{
  ns <- shiny::NS(id)
  
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

plot_deepnet_aescatter_server <- function(id,
                                          net,
                                          update,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function(n=12) {
      update()  ## react on updates
      net <- net()
      
      par(mfrow=c(1,1), mar=c(4,4,2,1))
      playbase::deep.plotAutoEncoderReconstructions(
        net, dtypes="mixed", par=FALSE) 
    }

    plot.RENDER2 <- function(n=12) {
      update()  ## react on updates
      net <- net()
      ntypes <- length(net$X)
      nc <- ceiling(sqrt(ntypes))
      nr <- ceiling(ntypes / nc)
      par(mfrow=c(nr,nc), mar=c(4,5,3,3))
      playbase::deep.plotAutoEncoderReconstructions(
        net, dtypes=NULL, par=FALSE) 
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      func2 = plot.RENDER2,      
      pdf.width = 10, pdf.height = 10,
      res = c(75, 100),
      add.watermark = watermark
    )
  })
}
