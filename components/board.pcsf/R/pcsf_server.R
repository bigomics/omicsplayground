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
      updateSelectInput(session, "contrast",
        choices = comparisons,
        selected = head(comparisons, 1)
      )
    })

    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    ## PCSF  analysis
    pcsf_compute.SAVE <- shiny::eventReactive(
      {
        pgx$X
      },
      {
        dbg("[pcsf_server.R:pcsf_compute] reacted!")

        NTOP <- 4000
        data(STRING, package = "PCSF")

        ## get foldchanges
        meta <- playbase::pgx.getMetaMatrix(pgx)$fc

        ## Genes in STRING are based in human genome, thus use
        ## ortholog. For now we collapse to human gene. maybe in
        ## future we may keep original features and use map features
        ## to symbol when creating STRING edges. For now this is
        ## simpler...
        meta <- playbase::rename_by(meta, pgx$genes, "human_ortholog")
        mx <- rowMeans(meta**2, na.rm = TRUE)**0.5
        genes <- head(names(mx)[order(-abs(mx))], NTOP)

        ## compute gene weight for number of gene sets that they are
        ## in. this is used to prioritize genes later.
        gsmeta <- playbase::pgx.getMetaMatrix(pgx, level = "geneset")$fc
        rx <- rowMeans(gsmeta**2, na.rm = TRUE)**0.5
        top.gs <- head(names(sort(-abs(rx))), 1000)

        ## GMT has organism specific symbol, but we need to map to
        ## human ortholog
        G <- pgx$GMT[, top.gs]
        G <- playbase::rename_by2(G, pgx$genes,
          ## from_id = "symbol",
          new_id = "human_ortholog"
        )
        ngmt <- Matrix::rowSums(G != 0, na.rm = TRUE)
        ngmt <- log(1 + ngmt)

        ## ------------ find gene clusters --------------

        ## also X needs to be in human ortholog
        X <- playbase::rename_by(pgx$X, pgx$genes, "human_ortholog")
        X <- X[genes, ]
        clust <- playbase::pgx.FindClusters(
          X = scale(t(X)),
          method = c("kmeans"),
          top.sd = 1000, npca = 40
        )

        ## take the minimum level that has at least 2 clusters
        nclust <- apply(clust[["kmeans"]], 2, function(x) length(table(x)))
        sel <- min(which(nclust > 1))
        idx <- clust[["kmeans"]][, sel]
        names(idx) <- genes
        idx <- idx[!duplicated(names(idx))]

        ## balance number of gene per group??
        genes <- unlist(tapply(names(idx), idx, function(gg) head(gg, 800)))
        genes <- as.vector(genes)
        idx <- idx[genes]

        ## ------------ create edge table --------------
        sel <- (STRING$from %in% genes & STRING$to %in% genes)
        ee <- STRING[sel, ]
        genes <- genes[which(genes %in% c(STRING$from, STRING$to))]

        # Create cost matrix from correlation
        rho <- cor(t(X[genes, ]))
        rho.ee <- pmax(rho[cbind(ee$from, ee$to)], 1e-4)
        ee$cost <- ee$cost**1 / rho.ee ## balance PPI with experiment

        ## ------------ create igraph network --------------
        ppi <- PCSF::construct_interactome(ee)
        ## gene set weighting. give more weight to genes that are in
        ## many gene sets.
        gset.wt <- (ngmt[genes] / max(ngmt[genes]))
        terminals <- abs(mx[genes])**1 * (0.1 + gset.wt)**1

        out <- list(
          genes = genes,
          meta = meta,
          X = X,
          idx = idx,
          clust = clust,
          ppi = ppi,
          terminals = terminals
        )
        dbg("[pcsf_server.R:pcsf_compute] done!")

        return(out)
      }
    )

    pcsf_compute <- shiny::eventReactive(
      {
        list(pgx$X, input$contrast, input$pcsf_beta, input$pcsf_ntop)
      },
      {
        contrast <- input$contrast
        beta <- as.numeric(input$pcsf_beta)
        ntop <- as.integer(input$pcsf_ntop)

        pcsf <- playbase::pgx.computePCSF(
          pgx,
          contrast,
          level = "gene",
          ntop = ntop,
          ncomp = 2,
          beta = 10^beta,
          use.corweight = TRUE,
          dir = "both",
          rm.negedge = TRUE,
          as.name = c("mx")
        )
        if (is.null(pcsf)) {
          validate()
          shiny::validate(
            !is.null(pcsf),
            "No PCSF solution found. Beta value is probably too small. Please adjust beta or increase network size."
          )
          return(NULL)
        }

        pcsf
      }
    )

    pcsf_plot_network_server(
      "pcsf_network",
      pgx,
      pcsf_compute = pcsf_compute,
      r_layout = shiny::reactive(input$layout),
      watermark = WATERMARK
    )

    pcsf_table_centrality_server(
      "centrality_table",
      pgx,
      r_contrast = shiny::reactive(input$contrast),
      r_pcsf = pcsf_compute
    )

    ## pcsf_plot_heatmap_server(
    ##   "pcsf_heatmap",
    ##   pgx,
    ##   pcsf_compute = pcsf_compute,
    ##   watermark = WATERMARK
    ## )
  })
}
