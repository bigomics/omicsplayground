##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

wgcna_plot_correlation_network_ui <- function(id, height, width) {
  ns <- shiny::NS(id)

  info_text <- "<b>Correlation network.</b> Partial correlation graph centered on module eigen-gene with top most correlated features. Green edges correspond to positive (partial) correlation, red edges to negative (partial) correlation. Width of the edges is proportional to the correlation strength of the gene pair. The regularized partial correlation matrix is computed using the 'graphical lasso' (Glasso) with BIC model selection."

  PlotModuleUI(
    ns("plot"),
    title = "Correlation network",
    label = "b",
    info.text = info_text,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

wgcna_plot_correlation_network_server <- function(id,
                                              wgcna.compute,
                                              selected_module,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    corGraph.RENDER <- shiny::reactive({
      out <- wgcna.compute()

      k <- selected_module()
      shiny::req(k)
      genes <- out$me.genes[[k]]

      dim(out$datExpr)
      xx <- cbind( out$net$MEs[,k,drop=FALSE], out$datExpr[,genes])
      rho1 <- cor(xx, out$net$MEs[,k] )[,1]
      ntop <- min(nrow(xx)-1,20)
      topgg <- names(sort(rho1,decreasing=TRUE))

      if(0) {
        gs0 <- out$gse[out$gse$module==k,]
        gs0 <- head(gs0[order(gs0$p.value),],10)
        gs.genes <- unique(unlist(strsplit(gs0$genes,split="\\|")))
        gs.genes
        topgg <- intersect(topgg, gs.genes) ## only GSET genes???
      }
      topgg <- head(topgg,ntop)

      rho <- Matrix::nearPD(cor(xx[,topgg]))$mat
      me.color <- out$me.colors[k]
      color1 <- me.color
      color1 <- c("white", me.color)[1 + 1*(colnames(rho)==k)]
      size1  <- c(7,10)[1 + 1*(colnames(rho)==k)]

      qgraph::qgraph(rho, graph="glasso", layout="spring", sampleSize=nrow(xx),
                     labels = rownames(rho), color=color1,
                     tuning = 0,  ## gamma for EBIClasso. 0.5=default, 0=BIC
                     vsize = size1, cut=0, maximum=.45,
                     border.width=1.5)
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = corGraph.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 80),
      add.watermark = watermark
    )
  })
}
