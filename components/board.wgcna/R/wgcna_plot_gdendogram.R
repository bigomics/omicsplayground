##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_gdendogram_ui <- function(id, info.text, caption, height, width) {
  ns <- shiny::NS(id)

##  info_text <- "<b>WGCNA gene dendrogram and gene modules.</b>"

  PlotModuleUI(
    ns("plot"),
    title = "Gene dendrogram and gene modules",
    label = "a",
    caption = caption,
    info.text = info.text,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

wgcna_plot_gdendogram_server <- function(id,
                                         wgcna.compute,
                                         labels2rainbow,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    geneDendro.RENDER <- shiny::reactive({

      out <- wgcna.compute()
      net <- out$net

      ## Convert labels to colors for plotting
      mergedColors <- labels2rainbow(net)
      ## Plot the dendrogram and the module colors underneath
      plotDendroAndColors(
        dendro = net$dendrograms[[1]],
        colors = mergedColors[net$blockGenes[[1]]],
        dendroLabels = FALSE, hang = 0.03,
        addGuide = FALSE, guideHang = 0.05,
        marAll = c(0.2, 5, 0.4, 0.2),
        main = NULL
      )
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = geneDendro.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 90),
      add.watermark = watermark
    )
  })
}
