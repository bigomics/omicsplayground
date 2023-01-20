##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

compare_plot_gene_corr_ui <- function(id, label='', height=c(600,800)) {

  ns <- shiny::NS(id)

  info_text <- "<b>FC scatter plots.</b> Scatter plots of gene expression scatter values between two contrasts. Scatters that are similar show high correlation, i.e. are close to the diagonal."

  genecorr.opts <- shiny::tagList(
    withTooltip( shiny::selectInput(ns('colorby'),'Color by:',choices=NULL, multiple=FALSE),
                 "Color samples by phenotype.",
                 placement="right", options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = "Gene correlation",
    label = "c",
    info.text = info_text,
    options = genecorr.opts,
    height = c(740,750),
    width=c('auto',900),
    download.fmt = c("png","pdf")
  )
}

compare_plot_gene_corr_server <- function(id,
                                          inputData,
                                          dataset2,
                                          input.contrast1,
                                          input.contrast2,
                                          hilightgenes,
                                          getOmicsScoreTable,
                                          score_table,
                                          contrast1,
                                          watermark=FALSE){
  moduleServer( id, function(input, output, session) {

    shiny::observeEvent(contrast1(), {
      ct <- contrast1()
      shiny::req(ct)
      shiny::updateSelectInput(session, "colorby", choices=ct, selected=ct[1])
    })

    plot_data <- shiny::reactive({

    })

    genecorr.RENDER <- shiny::reactive({

      ngs1 <- inputData()
      ngs2 <- dataset2()

      ct1 <- head(names(ngs1$gx.meta$meta),2)
      ct2 <- head(names(ngs2$gx.meta$meta),2)
      ct1 <- input.contrast1()
      ct2 <- input.contrast2()
      shiny::req(ct1)
      shiny::req(ct2)
      if(!all(ct1 %in% names(ngs1$gx.meta$meta))) return(NULL)
      if(!all(ct2 %in% names(ngs2$gx.meta$meta))) return(NULL)

      gg <- intersect(rownames(ngs1$X),rownames(ngs2$X))
      kk <- intersect(colnames(ngs1$X),colnames(ngs2$X))

      if(length(kk)==0) {
        par(mfrow=c(1,1))
        frame()
        text(0.5,0.6, "Warning: no common samples", col='black')
        text(0.5,0.5, "To compute gene correlation both datasets\nneed to have common samples",
             col='black')
        return()
      }
      if(length(kk) < 10) {
        par(mfrow=c(1,1))
        frame()
        text(0.5,0.6, "Error: too few samples", col='red3')
        text(0.5,0.5, "For gene correlation we need at least 10 common samples", col='red3')
        return()
      }

      ## conform matrices
      X1 <- ngs1$X[gg,kk]
      X2 <- ngs2$X[gg,kk]
      Y1 <- ngs1$samples[kk,]
      Y2 <- ngs2$samples[kk,]

      dset1 <- paste0("[dataset1]  expression")
      dset2 <- paste0("[dataset2]  expression")
      dset1 <- paste0("1: expression")
      dset2 <- paste0("2: expression")

      if(0) {
        F <- pgx.getMetaMatrix(ngs1)$fc
        higenes <- names(head(sort(-rowMeans(F**2)),16))
        higenes <- hilightgenes()
        higenes <- intersect(higenes,rownames(X1))
        higenes <- head(higenes, 16)
      }

      df <- getOmicsScoreTable()
      if(is.null(df)) return(NULL)

      sel <- score_table$rows_all()      ## from module
      shiny::req(sel)
      if(is.null(sel)) return(NULL)

      higenes <- head(rownames(df)[sel],16)
      if(length(higenes)==0) return(NULL)

      ## Set color for points
      klrpal <- rep(1:7,99)
      klrpal <- rep(RColorBrewer::brewer.pal(12,"Paired"),99)

      # colorby="ER_STATUS"
      # colorby = ct1[1]
      colorby <- input$colorby

      if(0) {
        grp <- factor(Y1[,colorby])
        klr1 <- klrpal[as.integer(grp)]
      } else {
        grp <- pgx.getContrastGroups(ngs1, colorby, as.factor=TRUE)
        grp <- grp[colnames(X1)]
        klr1 <- klrpal[as.integer(grp)]
      }

      nc <- ceiling(sqrt(length(higenes)))
      nr <- (length(higenes)-1) %/% nc
      par(mfrow=c(nc,nc), mar=c(2.6,2.3,1.5,0.5)*1, oma=c(3,3,0,0), mgp=c(2.4,0.7,0) )
      i=1
      for(i in 1:length(higenes)) {
        j <- match(higenes[i],rownames(X1))
        base::plot( X1[j,], X2[j,],
                    xlab="", ylab="",
                    pch=20, col=klr1, cex=1.2)
        title(higenes[i], line=0.4, cex.main=1.1)
        if( i%%nc == 1) {
          mtext(dset2, 2, line=2, cex=0.8)
        }
        if((i-1)%/%nc==nr) {
          mtext(dset1, 1, line=2, cex=0.8)
        }

        if(i%%nc == 1) {
          tt <- c("   ",levels(grp))
          legend("topleft", legend=tt,
                 fill = c(NA,klrpal),
                 border = c(NA,"black","black"), bty='n',
                 cex=0.92, box.lwd=0, pt.lwd=0,
                 x.intersp=0.5, y.intersp=0.8)
          legend("topleft", colorby, x.intersp=-0.2,
                 cex=0.92, y.intersp=0.45, bty='n')
        }
      }
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = genecorr.RENDER,
      # csvFunc = plot_data,
      pdf.width = 5, pdf.height = 5,
      res = c(80,90),
      add.watermark = watermark
    )
  }
  )
}