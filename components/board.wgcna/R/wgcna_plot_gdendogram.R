##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_gdendogram_ui <- function(id, label, info.text, caption, height, width) {
  ns <- shiny::NS(id)

  #

  PlotModuleUI(
    ns("plot"),
    title = "Gene dendrogram and gene modules",
    label = label,
    caption = caption,
    info.text = info.text,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

wgcna_plot_gdendogram_server <- function(id,
                                         wgcna.compute,
                                         power,
                                         networktype,
                                         tomtype,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    geneDendro.RENDER <- function() {
#      shiny::req(power())
#      shiny::req(networktype())
#      shiny::req(tomtype())      

      out <- wgcna.compute()
      net <- out$net

      ## playbase::plotDendroFromResults(
      ##   out,
      ##   power = as.numeric(power()),
      ##   networktype = networktype(),
      ##   tomtype = tomtype(),
      ##   nSelect = 1000
      ## ) 
          
      ## Convert labels to colors for plotting
      mergedColors <- playbase::labels2rainbow(net)
      ## Plot the dendrogram and the module colors underneath
      WGCNA::plotDendroAndColors(
        dendro = net$dendrograms[[1]],
        colors = mergedColors[net$blockGenes[[1]]],
        dendroLabels = FALSE, hang = 0.03,
        addGuide = FALSE, guideHang = 0.05,
        marAll = c(0.2, 5, 0.4, 0.2),
        main = NULL
      )

    }

    PlotModuleServer(
      "plot",
      func = geneDendro.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 90),
      add.watermark = watermark
    )
  })
}
