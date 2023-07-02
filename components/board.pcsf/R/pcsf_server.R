##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

PcsfBoard <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    fullH <- 800
    tabH <- "70vh"

    pcsf_info <- div(
      "This PCSF analysis module..."
    )

    observeEvent(input$pcsf_info, {
      showModal(
        modalDialog(
          title = tags$strong("PCSF Network Analysis"),
          pcsf_info,
          easyClose = TRUE,
          size = "xl"
        )
      )
    })

    observe({
      if (is.null(pgx)) {
        return(NULL)
      }
      comparisons <- colnames(pgx$model.parameters$contr.matrix)
      comparisons <- sort(comparisons)
      updateSelectInput(session, "contrast", choices = comparisons,
                        selected = head(comparisons,1))
    })

    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    ## PCSF  analysis
    pcsf_compute <- eventReactive(
    {
      pgx$X
    },
    {

      dbg("[pcsf_server.R:pcsf_compute] reacted!")
      
      NTOP = 4000
      data(STRING, package="PCSF")
      meta <- playbase::pgx.getMetaMatrix(pgx)$fc
      mx <- rowMeans(meta**2)**0.5
      genes <- head(names(mx)[order(-abs(mx))], NTOP)
      
      gsmeta <- playbase::pgx.getMetaMatrix(pgx, level="geneset")$fc
      rx <- rowMeans(gsmeta**2)**0.5
      top.gs <- head(names(sort(-abs(rx))),1000)
      ngmt <- Matrix::rowSums(pgx$GMT[,top.gs]!=0)
      ngmt <- log(1+ngmt)
      
      ## ------------ find gene clusters --------------      
      clust <- playbase::pgx.FindClusters(
        X = scale(t(pgx$X[genes,])),
        #
        method=c("kmeans"), 
        top.sd=1000, npca=40 )

      idx <- clust[['kmeans']][,3]
      names(idx) <- genes

      ## balance number of gene per group??
      table(idx)
      genes <- unlist(tapply(genes, idx[genes], function(gg) head(gg,800)))
      table(idx[genes])
      
      ## ------------ create edge table --------------
      sel <- (STRING$from %in% genes & STRING$to %in% genes)
      table(sel)
      ee <- STRING[sel,]
      genes <- genes[which(genes %in% c(STRING$from, STRING$to))]
      rho <- cor(t(pgx$X[genes,]))
      #
      rho.ee <- pmax(rho[cbind(ee$from, ee$to)],0.001)
      ee$cost <- ee$cost**1 / rho.ee  ## balance PPI with experiment
      #

      ## ------------ create igraph network --------------      
      ppi <- PCSF::construct_interactome(ee)
      gset.wt <- (ngmt[genes]/max(ngmt[genes]))  ## gene set weighting
      terminals <- abs(mx[genes])**1 * (0.1 + gset.wt)**1
      
      list(
        genes = genes,
        meta = meta[genes,],
        idx = idx[genes],
        clust = clust,
        ppi = ppi,
        terminals = terminals
      )
      
    })

    
    pcsf_plot_network_server(
      "pcsf_network",
      pgx,
      pcsf_compute = pcsf_compute,
      pcsf_beta = shiny::reactive(input$pcsf_beta),
      colorby = shiny::reactive(input$colorby),
      contrast = shiny::reactive(input$contrast),
      highlightby = shiny::reactive(input$highlightby),
      watermark = WATERMARK
    )

    pcsf_plot_heatmap_server(
      "pcsf_heatmap",
      pgx,
      pcsf_compute = pcsf_compute,
      watermark = WATERMARK
    )
    

  })
}
