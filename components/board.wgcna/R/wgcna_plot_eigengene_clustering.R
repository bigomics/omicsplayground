##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_eigengene_clustering_ui <- function(id, height, width) {
  ns <- shiny::NS(id)

  info_text <- "<b>eigenClustering</b> <b>(a)</b> ..."

  PlotModuleUI(
    ns("plot"),
    title = "Eigengene clustering",
    label = "a",
    info.text = info_text,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

wgcna_plot_eigengene_clustering_server <- function(id,
                                                   wgcna.compute,
                                                   watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    eigenClustering.RENDER <- shiny::reactive({
      out <- wgcna.compute()

      MET <- out$net$MEs
      if (NCOL(MET) <= 2) MET <- cbind(MET, MET) ## error if ncol(MET)<=2 !!!!

      ## Plot the relationships among the eigengenes and the trait
      WGCNA::plotEigengeneNetworks(
        MET, "",
        marDendro = c(0, 4, 1, 2),
        marHeatmap = c(3, 4, 1, 2), cex.lab = 0.8,
        xLabelsAngle = 90
      )
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = eigenClustering.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(80, 90),
      add.watermark = watermark
    )
  })
}
