##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_gsetmofa_factorCor_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)


  options <- shiny::tagList(
    shiny::checkboxInput(
      inputId = ns("cluster"),
      label = "Cluster heatmap",
      value = TRUE
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
    download.fmt = c("png", "pdf")
  )
}

mofa_plot_gsetmofa_factorCor_server <- function(id,
                                                mofa,
                                                watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function() {
      res <- mofa()      
      K <- t(cor(res$F, res$gset.mofa$F))

      ## if(input$condition) {
      ##   pheno <- input$pheno
      ##   z1 <- res$gset.mofa$Z[pheno,]
      ##   z2 <- res$Z[pheno,]
      ##   W <- abs(outer(z1, z2))
      ##   K <- K * W
      ## }
      
      if(input$cluster) {
        ii <- hclust(dist(K))$order
        jj <- hclust(dist(t(K)))$order
        K <- K[ii,jj]
      }

      ftext <- round(K, digits=2)
      par(mar=c(6,5,1,1))
      WGCNA::labeledHeatmap(
        Matrix = K,
        xLabels = colnames(K), 
        yLabels = rownames(K), 
        textMatrix = ftext,
        cex.text = 0.7,
        cex.lab = 0.9, 
        #  ySymbols = colnames(res$F),
        colorLabels = TRUE, 
        colors = WGCNA::blueWhiteRed(50), 
        setStdMargins = FALSE, 
        zlim = c(-1,1),
        main = ""
      )
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 12, pdf.height = 8,
      res = c(75, 110),
      add.watermark = watermark
    )

    
  })
}



