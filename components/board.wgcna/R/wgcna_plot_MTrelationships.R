##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_MTrelationships_ui <- function(
    id,
    title,
    label,
    info.text,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)

  ## moduleTrait_opts <- shiny::tagList(
  ##  shiny::checkboxInput(ns("traits_binarize"), "binarize continuous vars", FALSE)
  ## )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    ## options = moduleTrait_opts,
    height = height,
    width = width,
    download.fmt = NULL # FIXME png and pdf is not working, to avoid crashing other things, we decided to remove it
  )
}

wgcna_plot_MTrelationships_server <- function(id,
                                              wgcna.compute,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    moduleTrait.RENDER <- function() {
      out <- wgcna.compute()
      net <- out$net
      datExpr <- out$datExpr
      datTraits <- out$datTraits
      moduleColors <- playbase::labels2rainbow(out$net)
      MEs <- out$net$MEs

      ## Define numbers of genes and samples
      nGenes <- ncol(datExpr)
      nSamples <- nrow(datExpr)

      moduleTraitCor <- cor(MEs, out$datTraits, use = "pairwise")
      moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

      textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
        signif(moduleTraitPvalue, 1), ")",
        sep = ""
      )
      textMatrix <- signif(moduleTraitCor, 2)

      dim(textMatrix) <- dim(moduleTraitCor)
      dim(moduleTraitCor)

      sel1 <- 1:nrow(moduleTraitCor)
      sel2 <- 1:ncol(moduleTraitCor)

      sel2 <- sort(head(order(-colMeans(abs(moduleTraitCor), na.rm = TRUE)), 40)) ## conditions
      sel1 <- sort(head(order(-rowMeans(abs(moduleTraitCor[, sel2]), na.rm = TRUE)), 12)) ## eigenvectors

      sel2 <- sort(head(order(-colMeans(pmax(moduleTraitCor, 0), na.rm = TRUE)), 40)) ## conditions
      sel1 <- sort(head(order(-rowMeans(pmax(moduleTraitCor[, sel2], 0), na.rm = TRUE)), 12)) ## eigenv


      par(mar = c(3, 12, 1.6, 1.5))
      ## Display the correlation values within a heatmap plot
      WGCNA::labeledHeatmap(
        Matrix = t(moduleTraitCor[sel1, sel2]),
        yLabels = colnames(out$datTraits)[sel2],
        xLabels = colnames(MEs)[sel1],
        xSymbols = colnames(MEs)[sel1],
        xLabelsAngle = 90,
        colorLabels = FALSE,
        colors = greenWhiteRed(50),
        textMatrix = t(textMatrix[sel1, sel2]),
        setStdMargins = FALSE,
        cex.text = 0.7,
        cex.lab = 0.9,
        zlim = c(-1, 1),
        main = NULL
      )
    }

    PlotModuleServer(
      "plot",
      func = moduleTrait.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 90),
      add.watermark = watermark
    )
  })
}
