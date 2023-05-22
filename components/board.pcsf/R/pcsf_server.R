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

      if(0) {
        library(playbase)
        library(PCSF)
        source("./util_pcsf.R")
        pgx <- playbase::pgx.load("~/Playground/pgx/example-data.pgx")
        names(pgx$gx.meta$meta)
      }

      dbg("[pcsf_server.R:pcsf_compute] reacted!")
      
      NTOP = 4000
      data(STRING, package="PCSF")
      meta <- playbase::pgx.getMetaMatrix(pgx)$fc
      mx <- rowMeans(meta**2)**0.5
      genes <- head(names(mx)[order(-abs(mx))], NTOP)
      genes <- genes[which(genes %in% c(STRING$from, STRING$to))]
      genes <- genes[order(-abs(mx[genes]))]
      
      gsmeta <- playbase::pgx.getMetaMatrix(pgx, level="geneset")$fc
      rx <- rowMeans(gsmeta**2)**0.5
      top.gs <- head(names(sort(-abs(rx))),1000)
      ngmt <- Matrix::rowSums(pgx$GMT[,top.gs]!=0)
      ngmt <- log(1+ngmt)
      
      ## ------------ find gene clusters --------------      
      clust <- playbase::pgx.FindClusters(
        X = scale(t(pgx$X[genes,])),
        ##  method=c("kmeans","hclust","louvain","meta"),
        method=c("kmeans"), 
        top.sd=1000, npca=40 )

      idx <- clust[['kmeans']][,3]
      names(idx) <- genes

      ## balance number of gene per group??
      table(idx)
      genes <- unlist(tapply(genes, idx[genes], function(gg) head(gg,400)))
      table(idx[genes])
      
      ## ------------ create igraph network --------------
      sel <- (STRING$from %in% genes & STRING$to %in% genes)
      table(sel)
      ee <- STRING[sel,]
      rho <- cor(t(pgx$X[genes,]))
      ##rho.ee <- abs(rho[cbind(ee$from, ee$to)])
      rho.ee <- pmax(rho[cbind(ee$from, ee$to)],0.001)
      ee$cost <- ee$cost**1 / rho.ee  ## balance PPI with experiment
      ##ee$cost <- (1 / rho.ee)**1
      ppi <- PCSF::construct_interactome(ee)
      gset.wt <- (ngmt[genes]/max(ngmt[genes]))  ## gene set weighting
      terminals <- abs(mx[genes])**1 * (0.1 + gset.wt)**1
      net <- PCSF::PCSF(ppi, terminals, w=2, b=1)

      ## remove small clusters...
      cmp <- igraph::components(net)
      sel.kk <- which(cmp$csize > 0.10 * max(cmp$csize))
      net <- igraph::subgraph(net, cmp$membership %in% sel.kk)
      class(net) <- c("PCSF","igraph")
      
      igraph::V(net)$group <- idx[igraph::V(net)$name]
      igraph::V(net)$type  <- idx[igraph::V(net)$name]

      dbg("[pcsf_server.R:pcsf_compute] done!")
      
      list(
        genes = genes,
        meta = meta[genes,],
        clust = clust,
        net = net
      )
      
    })

    
    pcsf_plot_network_server(
      "pcsf_network",
      pgx,
      pcsf_compute = pcsf_compute,
      colorby = shiny::reactive(input$colorby),
      contrast = shiny::reactive(input$contrast),
      watermark = WATERMARK
    )

    pcsf_plot_heatmap_server(
      "pcsf_heatmap",
      pgx,
      pcsf_compute = pcsf_compute,
      colorby = shiny::reactive(input$colorby),
      contrast = shiny::reactive(input$contrast),
      watermark = WATERMARK
    )
    

  })
}
