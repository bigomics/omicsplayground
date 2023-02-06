##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

wgcna_plot_MTrelationships_ui <- function(id, height, width) {
  ns <- shiny::NS(id)

  info_text <- "<b>WGCNA module and trait relationship.</b>"

  moduleTrait_opts = shiny::tagList(
    shiny::checkboxInput(ns("traits_binarize"),"binarize continuous vars", FALSE)
  )

  PlotModuleUI(
    ns("plot"),
    title = "Module-Trait relationships",
    label = "a",
    info.text = info_text,
    options = moduleTrait_opts,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

wgcna_plot_MTrelationships_server <- function(id,
                                              wgcna.compute,
                                              labels2rainbow,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    moduleTrait.RENDER <- shiny::reactive({
      out <- wgcna.compute()
      net <- out$net
      datExpr <- out$datExpr
      datTraits <- out$datTraits
      moduleColors <- labels2rainbow(out$net)
      MEs <- out$net$MEs

      ## Define numbers of genes and samples
      nGenes = ncol(datExpr);
      nSamples = nrow(datExpr);

      if(0) {
        ## Recalculate MEs with color as labels
        MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
        MEs  = orderMEs(MEs0)
      }

      moduleTraitCor = cor(MEs, out$datTraits, use = "pairwise");
      moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

      textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                         signif(moduleTraitPvalue, 1), ")", sep = "")
      textMatrix = signif(moduleTraitCor, 2)

      dim(textMatrix) = dim(moduleTraitCor)
      dim(moduleTraitCor)

      sel1 <- 1:nrow(moduleTraitCor)
      sel2 <- 1:ncol(moduleTraitCor)

      sel2 <- sort(head(order(-colMeans(abs(moduleTraitCor))),40)) ## conditions
      sel1 <- sort(head(order(-rowMeans(abs(moduleTraitCor[,sel2]))),12)) ## eigenvectors

      sel2 <- sort(head(order(-colMeans(pmax(moduleTraitCor,0))),40)) ## conditions
      sel1 <- sort(head(order(-rowMeans(pmax(moduleTraitCor[,sel2],0))),12)) ## eigenv

      message("[moduleTrait.RENDER] sel1 = ",paste(sel1,collapse=" "))
      message("[moduleTrait.RENDER] sel2 = ",paste(sel2,collapse=" "))

      par(mar = c(3, 12, 1.6, 1.5))
      ## Display the correlation values within a heatmap plot
      labeledHeatmap(Matrix = t(moduleTraitCor[sel1,sel2]),
                     yLabels = colnames(out$datTraits)[sel2],
                     xLabels = colnames(MEs)[sel1],
                     xSymbols = colnames(MEs)[sel1],
                     xLabelsAngle = 90,
                     colorLabels = FALSE,
                     colors = greenWhiteRed(50),
                     textMatrix = t(textMatrix[sel1,sel2]),
                     setStdMargins = FALSE,
                     cex.text = 0.7,
                     cex.lab = 0.9,
                     zlim = c(-1,1),
                     main = NULL
      )
    })

    PlotModuleServer(
      "plot",
      func = moduleTrait.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 90),
      add.watermark = watermark
    )
  })
}
