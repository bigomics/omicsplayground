##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

consensusWGCNA_plot_preservationDendro_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("dendro"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

consensusWGCNA_plot_preservationHeatmap_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("heatmap"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

consensusWGCNA_plot_preservation_server <- function(id,
                                               mwgcna
                                               ) {
  moduleServer(id, function(input, output, session) {
    
    renderDendro <- function() {
      cons <- mwgcna()
      shiny::req(cons)
      
      par(mfrow=c(2,1), mar=c(10,3,2,0.5))
      multiMEs <- cons$net$multiMEs
      for(k in names(multiMEs)) {
        ME <- multiMEs[[k]]$data
        playbase::wgcna.plotEigenGeneClusterDendrogram(
          wgcna = NULL,
          ME = ME,
          main = k,
          setMargins = FALSE,
          method = "hclust"
        )
      }      
    }

    renderHeatmap <- function() {
      cons <- mwgcna()
      shiny::req(cons)

      multiMEs <- cons$net$multiMEs
      
      par(cex = 0.9)
      WGCNA::plotEigengeneNetworks(
        multiMEs,
        setLabels = names(multiMEs),
        plotDendrograms = FALSE,
        marDendro = c(0,2,2,1)*1.1,
        marHeatmap = c(3,3,2,1)*1.1,
        #zlimPreservation = c(0.5, 1),
        xLabelsAngle = 90
      )
      
    }
    
    PlotModuleServer(
      "dendro",
      func = renderDendro,
      pdf.width = 12,
      pdf.height = 8,
      res = c(75, 110),
      add.watermark = FALSE
    )

    PlotModuleServer(
      "heatmap",
      func = renderHeatmap,
      pdf.width = 12,
      pdf.height = 8,
      res = c(80, 110),
      add.watermark = FALSE
    )

    
  })
}



