##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

enrichment_plot_scatter_ui <- function(id, height, width) {
  ns <- shiny::NS(id)

  info_text <- "A scatter plot of enrichment scores versus expression values across the samples for the gene set selected from the enrichment analysis Table <code>I</code> and the gene selected from the genes Table <code>II</code>."

  PlotModuleUI(
    ns("plot"),
    title = "Enrichment vs. expression",
    label = "d",
    info.text = info_text,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

enrichment_plot_scatter_server <- function(id,
                                           inputData,
                                           gene_selected,
                                           gs_contrast,
                                           subplot.MAR,
                                           gset_selected,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    getcolors <- function(ngs, comp0) {
      ## get colors (what a mess...)
      contr.matrix <- ngs$model.parameters$contr.matrix

      exp.matrix <- ngs$model.parameters$exp.matrix
      xgroup = as.character(ngs$Y$group)

      grp.name <- strsplit(comp0,split="[._ ]vs[._ ]")[[1]]
      grp.name <- c(grp.name, "other")
      xsign <- sign(exp.matrix[,comp0])
      xgroup = grp.name[1*(xsign>0) + 2*(xsign<0) + 1*(xsign==0)]
      table(xgroup)

      names(xgroup) = rownames(ngs$Y)
      table(xgroup)
      samples = names(which(exp.matrix[,comp0]!=0))

      xgroup1 <- xgroup[samples]
      table(xgroup1)
      ngrp <- length(unique(xgroup1))
      grp.klr = c("grey90",rep(RColorBrewer::brewer.pal(12,"Paired"),99)[1:ngrp])
      names(grp.klr) <- c("other",as.character(sort(unique(xgroup1))))
      grp.klr

      xgroup2 <- as.character(xgroup)
      xgroup2[which(!(xgroup %in% xgroup1))] <- "other"
      sample.klr = grp.klr[xgroup2]
      names(sample.klr) <- rownames(ngs$samples)
      table(sample.klr)
      list(samples=sample.klr, group=grp.klr)
    }

    subplot_scatter.RENDER <- shiny::reactive({
      par(mfrow=c(1,1), mgp=c(1.8,0.8,0), oma=c(0,0,0,0.4) )
      par(mar=subplot.MAR)

      ngs <- inputData()
      shiny::req(ngs)


      gene = rownames(ngs$X)[1]
      sel  = gene_selected()
      gset = gset_selected()
      if(is.null(sel)) return(NULL)
      if(is.null(gset)) return(NULL)

      if(is.null(sel) || length(sel)==0) {
        frame()
      } else {
        gene = sel$gene
        gset <- gset[1]
        gx = ngs$X[sel$probe,]
        sx = ngs$gsetX[gset,]
        if( length(gx)==0 || length(sx)==0 ||
            length(gx)!=length(sx) ) {
          frame()
          return(NULL)
        }
        ## get colors
        comp0 = "Th17_mut_2h_VS_Th17_wt_2h_BLA"
        comp0 = gs_contrast()
        klrs = getcolors(ngs, comp0)
        klr = klrs$samples[names(sx)]
        klr = paste0(gplots::col2hex(klr),"99")

        cex1 = c(1.4,0.8,0.3)[cut(length(gx),c(0,100,500,99999))]
        gset1 = breakstring(substring(gset,1,80),32)
        #tt = paste( breakstring(gset,40,80), "\nvs.", gene,"expression")
        tt = paste( breakstring(gset,40,80), " vs. ", gene)
        base::plot( gx, sx, col=klr, main=tt,
                    ylab = "gene set enrichment",
                    xlab = paste(gene,"expression"),
                    cex.lab=1, pch=19, cex=1.0*cex1, cex.main=0.85)
        abline( lm(sx ~ gx), lty=2, lwd=0.7, col="black" )
      }
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = subplot_scatter.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(72,100),
      add.watermark = watermark
    )
  })
}
