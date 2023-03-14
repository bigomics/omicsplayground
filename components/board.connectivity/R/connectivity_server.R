##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


ConnectivityBoard <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    fullH <- 750 # row height of panel
    tabH <- "70vh"

    cmap_infotext <- strwrap(
      "The <strong>Experiment connectivity</strong> module enables users to
      compare their data to other datasets. For the selected contrast, this
      module provides pairwise correlation plots and/or enrichment plots with
      signatures from other data sets. The <strong>Connectivity map</strong>
      shows the similarity of the contrasts profiles as a t-SNE plot.<br><br>
      <br><br><center><iframe width='500' height='333'
      src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=5'
      frameborder='0' allow='accelerometer; autoplay; encrypted-media;
      gyroscope; picture-in-picture' allowfullscreen></iframe></center>"
    )


    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    shiny::observeEvent(input$cmap_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Connectivity Analysis Board</strong>"),
        shiny::HTML(cmap_infotext),
        easyClose = TRUE, size = "l"
      ))
    })

    ## update choices upon change of data set
    shiny::observe({
      if (is.null(pgx)) {
        return(NULL)
      }
      comparisons <- colnames(pgx$model.parameters$contr.matrix)
      comparisons <- sort(comparisons)
      shiny::updateSelectInput(session, "cmap_contrast",
        choices = comparisons,
        selected = head(comparisons, 1)
      )

      sigdb <- c("A", "B", "C")
      sigdb0 <- dir(SIGDB.DIR, pattern = "sigdb-.*h5")
      sigdb <- names(pgx$connectivity) ## only precomputed inside PGX object??
      sigdb <- sort(intersect(sigdb, sigdb0))
      sel <- sigdb[1]
      shiny::updateSelectInput(session, "cmap_sigdb", choices = sigdb, selected = sel)
    })

    shiny::observe({
      ## reset CMap threshold zero/max
      res <- getConnectivityScores()
      shiny::req(res)
      max <- round(0.999 * max(abs(res$score), na.rm = TRUE), digits = 1)
      max <- round(0.999 * tail(sort(abs(res$score)), 10)[1], digits = 1)
      shiny::updateSliderInput(session, "cmap_scorethreshold", value = 0, max = max)
    })

    ## update choices upon change of chosen contrast
    shiny::observeEvent(input$cmap_contrast, {
      shiny::req(pgx)

      ## reset CMap threshold zero/max
      res <- getConnectivityScores() ## result gets cached
      shiny::req(res)
      max <- round(0.999 * max(abs(res$score), na.rm = TRUE), digits = 1)
      shiny::updateSliderInput(session, "cmap_scorethreshold", value = 0, max = max)
    })

    getCurrentContrast <- shiny::reactive({
      shiny::req(pgx, input$cmap_contrast)
      ct <- input$cmap_contrast
      fc <- pgx$gx.meta$meta[[ct]]$meta.fx
      names(fc) <- rownames(pgx$gx.meta$meta[[ct]])
      gs <- pgx$gset.meta$meta[[ct]]$meta.fx
      names(gs) <- rownames(pgx$gset.meta$meta[[ct]])
      names(fc) <- toupper(names(fc)) ## de-MOUSE
      list(name = ct, fc = fc, gs = gs)
    })

    ## ================================================================================
    ## ========================= REACTIVE FUNCTIONS ===================================
    ## ================================================================================
    cumEnrichmentTable <- shiny::reactive({
      cmap_sigdb <- input$cmap_sigdb
      shiny::req(cmap_sigdb)
      if (!grepl(".h5$", cmap_sigdb)) {
        return(NULL)
      }

      df <- getConnectivityScores()
      if (is.null(df)) {
        return(NULL)
      }
      ii <- connectivityScoreTable$rows_all()
      shiny::req(ii)

      sel <- head(df$pathway[ii], 10)
      sigdb <- input$cmap_sigdb
      F <- getEnrichmentMatrix(sigdb, select = sel)
      if (is.null(F)) {
        return(NULL)
      }

      ## multiply with sign of enrichment
      rho1 <- df$rho[match(colnames(F), df$pathway)]
      F <- t(t(F) * sign(rho1))

      F <- F[order(-rowMeans(F**2)), , drop = FALSE]

      ## add current contrast
      ct <- getCurrentContrast()
      gx <- ct$gs[match(rownames(F), names(ct$gs))]
      names(gx) <- rownames(F)
      gx[is.na(gx)] <- 0
      F <- cbind(gx, F)
      colnames(F)[1] <- ct$name

      F
    })

    getConnectivityFullPath <- function(sigdb) {
      db.exists <- sapply(SIGDB.DIR, function(d) file.exists(file.path(d, sigdb)))
      db.dir <- names(which(db.exists))[1]
      file.path(db.dir, sigdb)
    }

    getConnectivityContrasts <- function(sigdb) {
      if (length(sigdb) == 0 || is.null(sigdb) || sigdb == "") {
        return(NULL)
      }

      db <- getConnectivityFullPath(sigdb)
      cn <- NULL
      if (file.exists(db)) {
        cn <- rhdf5::h5read(db, "data/colnames")
      }
      cn
    }

    getConnectivityMatrix <- function(sigdb, select = NULL, genes = NULL) {
      if (sigdb == "" || is.null(sigdb)) {
        warning("[getConnectivityMatrix] ***WARNING*** sigdb=", sigdb)
        return(NULL)
      }

      db.exists <- sapply(SIGDB.DIR, function(d) file.exists(file.path(d, sigdb)))
      X <- NULL
      if (any(db.exists)) {
        db.dir <- names(which(db.exists))[1]
        if (grepl("csv$", sigdb)) {
          X <- read.csv(file.path(db.dir, sigdb), row.names = 1, check.names = FALSE)
          X <- as.matrix(X)
          X <- X[, colMeans(is.na(X)) < 0.99, drop = FALSE] ## omit empty columns
          if (!is.null(genes)) X <- X[intersect(genes, rownames(X)), , drop = FALSE]
          if (!is.null(select)) X <- X[, intersect(select, colnames(X))]
        }
        if (grepl("h5$", sigdb)) {
          h5.file <- file.path(db.dir, sigdb)
          cn <- rhdf5::h5read(h5.file, "data/colnames")
          rn <- rhdf5::h5read(h5.file, "data/rownames")
          rowidx <- 1:length(rn)
          colidx <- 1:length(cn)
          if (!is.null(genes)) rowidx <- match(intersect(genes, rn), rn)
          if (!is.null(select)) colidx <- match(intersect(select, cn), cn)

          nr <- length(rowidx)
          nc <- length(colidx)

          X <- rhdf5::h5read(h5.file, "data/matrix", index = list(rowidx, colidx))
          rownames(X) <- rn[rowidx]
          colnames(X) <- cn[colidx]
        }
      }
      return(X)
    }

    getEnrichmentMatrix <- function(sigdb, select = NULL, nc = -1) {
      if (sigdb == "" || is.null(sigdb)) {
        warning("[getEnrichmentMatrix] ***WARNING*** sigdb=", sigdb)
        return(NULL)
      }
      if (!grepl("h5$", sigdb)) {
        stop("getEnrichmentMatrix:: only for H5 database files")
        return(NULL)
      }

      h5exists <- function(h5.file, obj) {
        xobjs <- apply(rhdf5::h5ls(h5.file)[, 1:2], 1, paste, collapse = "/")
        obj %in% gsub("^/|^//", "", xobjs)
      }

      db.exists <- sapply(SIGDB.DIR, function(d) file.exists(file.path(d, sigdb)))
      Y <- NULL
      if (any(db.exists)) {
        db.dir <- names(which(db.exists))[1]
        db.dir
        h5.file <- file.path(db.dir, sigdb)
        cn <- rhdf5::h5read(h5.file, "data/colnames")

        has.gs <- h5exists(h5.file, "enrichment/genesets")
        has.gsea <- h5exists(h5.file, "enrichment/GSEA")
        if (!has.gs && has.gsea) {
          return(NULL)
        }

        rn <- rhdf5::h5read(h5.file, "enrichment/genesets")
        rowidx <- 1:length(rn)
        colidx <- 1:length(cn)
        if (!is.null(select)) colidx <- match(intersect(select, cn), cn)
        Y <- rhdf5::h5read(h5.file, "enrichment/GSEA", index = list(rowidx, colidx))
        rownames(Y) <- rn[rowidx]
        colnames(Y) <- cn[colidx]
        sdy <- apply(Y, 1, sd)
        Y <- Y[order(-sdy), ]
      }

      ## cluster genesets into larger groups
      if (nc > 0) {
        hc <- hclust(dist(Y[, ]))
        idx <- paste0("h", cutree(hc, nc))
        Y2 <- tapply(1:nrow(Y), idx, function(i) colMeans(Y[i, , drop = FALSE]))
        Y2 <- do.call(rbind, Y2)
        idx.names <- tapply(rownames(Y), idx, paste, collapse = ",")
        idx.names <- gsub("H:HALLMARK_", "", idx.names)
        idx.names <- gsub("C2:KEGG_", "", idx.names)
        rownames(Y2) <- as.character(idx.names[rownames(Y2)])
        Y <- Y2
      }

      if (nrow(Y) == 0) {
        return(NULL)
      }

      return(Y)
    }

    getSignatureMatrix <- function(sigdb) {
      if (sigdb == "" || is.null(sigdb)) {
        warning("[getSignatureMatrix] ***WARNING*** sigdb=", sigdb)
        return(NULL)
      }

      if (!grepl("h5$", sigdb)) {
        stop("getEnrichmentMatrix:: only for H5 database files")
      }

      db.exists <- sapply(SIGDB.DIR, function(d) file.exists(file.path(d, sigdb)))
      up <- dn <- NULL
      if (any(db.exists)) {
        db.dir <- names(which(db.exists))[1]
        h5.file <- file.path(db.dir, sigdb)
        rhdf5::h5ls(h5.file)
        cn <- rhdf5::h5read(h5.file, "data/colnames")
        dn <- rhdf5::h5read(h5.file, "signature/sig100.dn")
        up <- rhdf5::h5read(h5.file, "signature/sig100.up")
        colnames(dn) <- cn
        colnames(up) <- cn
      }
      list(up = up, dn = dn)
    }

    getConnectivityScores <- shiny::reactive({
      # browser()
      shiny::req(pgx, input$cmap_contrast)
      shiny::validate(shiny::need("connectivity" %in% names(pgx), "no 'connectivity' in object."))

      ntop <- 1000
      sigdb <- input$cmap_sigdb
      shiny::req(sigdb)

      all.scores <- NULL
      if (sigdb %in% names(pgx$connectivity)) {
        all.scores <- pgx$connectivity[[sigdb]]
      } else {
        warning("[getConnectivityScores] ERROR : could not get scores")
        return(NULL)
      }

      ct <- input$cmap_contrast
      if (!ct %in% names(all.scores)) {
        warning("[getConnectivityScores] ERROR : contrast not in connectivity scores")
        return(NULL)
      }

      scores <- as.data.frame(all.scores[[ct]])
      if (input$cmap_abs_score == FALSE) {
        ## put sign back!!!
        scores$score <- scores$score * sign(scores$rho)
      }
      scores <- scores[order(-abs(scores$score)), ]
      scores <- scores[!duplicated(scores$pathway), ]
      rownames(scores) <- scores$pathway

      if (nrow(scores) == 0 || ncol(scores) == 0) {
        warning("[getConnectivityScores] ERROR : scores has zero dimensions")
        return(NULL)
      }

      if (input$cmap_hideclustcontrasts) {
        sel <- grep("cluster[:]", scores$pathway, invert = TRUE)
        scores <- scores[sel, , drop = FALSE]
      }

      ## only those in existing database
      cts <- getConnectivityContrasts(sigdb)
      scores <- scores[which(rownames(scores) %in% cts), , drop = FALSE]

      ## filter on significance
      qsig <- input$connectivityScoreTable_qsig
      scores <- scores[which(scores$padj <= qsig), , drop = FALSE]
      scores <- scores[order(-scores$score), , drop = FALSE]

      no.le <- !("leadingEdge" %in% colnames(scores))
      abs_score <- input$cmap_abs_score
      ntop <- 100

      if (no.le && abs_score == TRUE) {
        ## recreate "leadingEdge" list
        sig <- getSignatureMatrix(sigdb)
        fc <- getCurrentContrast()$fc
        fc <- fc[order(-abs(fc))]

        fc.up <- head(names(fc[fc > 0]), ntop)
        fc.dn <- head(names(fc[fc < 0]), ntop)
        ff <- c(fc.up, fc.dn)
        e1 <- apply(sig$up, 2, function(g) intersect(ff, g))
        e2 <- apply(sig$dn, 2, function(g) intersect(ff, g))
        ee <- mapply(c, e1, e2)
        ee <- ee[match(scores$pathway, names(ee))]
        scores$leadingEdge <- ee
      }
      if (no.le && abs_score == FALSE) {
        ## recreate "leadingEdge" list
        sig <- getSignatureMatrix(sigdb)
        fc <- getCurrentContrast()$fc
        fc <- fc[order(-abs(fc))]

        fc.up <- head(names(fc[fc > 0]), ntop)
        fc.dn <- head(names(fc[fc < 0]), ntop)

        p1 <- apply(sig$up, 2, function(g) intersect(fc.up, g))
        p2 <- apply(sig$dn, 2, function(g) intersect(fc.dn, g))
        pp <- mapply(c, p1, p2)

        n1 <- apply(sig$up, 2, function(g) intersect(fc.dn, g))
        n2 <- apply(sig$dn, 2, function(g) intersect(fc.up, g))
        nn <- mapply(c, n1, n2)

        ee <- vector("list", nrow(scores))
        pos.rho <- which(scores$rho >= 0)
        neg.rho <- which(scores$rho < 0)
        ee[pos.rho] <- pp[match(scores$pathway[pos.rho], names(pp))]
        ee[neg.rho] <- nn[match(scores$pathway[neg.rho], names(nn))]
        scores$leadingEdge <- ee
      }

      ## bail out
      if (nrow(scores) == 0) {
        return(NULL)
      }

      return(scores)
    })

    ## ================================================================================
    ## Correlation score table
    ## ================================================================================

    PERTINFO <- NULL
    pert_info.file <- file.path(FILESX, "GSE92742_Broad_LINCS_pert_info.txt")
    if (file.exists(pert_info.file)) {
      PERTINFO <- read.csv(pert_info.file, sep = "\t", row.names = 1)
    }

    getTopProfiles <- shiny::reactive({
      ## Get profiles of top-enriched contrasts (not all genes...)
      ##
      ##
      df <- getConnectivityScores()

      ii <- 1:100
      sigdb <- "sigdb-archs4.h5"

      ii <- connectivityScoreTable$rows_all()
      shiny::req(ii, input$cmap_sigdb)
      ii <- head(ii, 50) ## 50??
      pw <- df$pathway[ii]

      sigdb <- input$cmap_sigdb
      shiny::req(sigdb)

      fc <- getCurrentContrast()$fc
      ngenes <- 1000
      ngenes <- 500
      var.genes <- head(names(sort(-abs(fc))), ngenes)
      var.genes <- unique(c(var.genes, sample(names(fc), ngenes))) ## add some random
      F <- getConnectivityMatrix(sigdb, select = pw, genes = var.genes)
      pw <- intersect(pw, colnames(F))
      F <- F[, pw, drop = FALSE]
      return(F)
    })

    ## ============================================================================
    ## FC correlation/scatter plots
    ## ============================================================================

    connectivity_plot_cmap_FCFCplots_server(
      "cmap_FCFCplots",
      pgx,
      reactive(input$cmap_contrast),
      getCurrentContrast,
      getTopProfiles,
      getConnectivityScores,
      watermark = WATERMARK
    )

    connectivityScoreTable <- connectivity_table_similarity_scores_server(
      "connectivityScoreTable",
      getConnectivityScores = getConnectivityScores,
      cmap_sigdb = shiny::reactive(input$cmap_sigdb)
    )


    ## ================================================================================
    ## Cumulative FC barplot
    ## ================================================================================

    connectivity_plot_cumFCplot_server(
      "cumFCplot",
      getTopProfiles,
      getConnectivityScores,
      getCurrentContrast
    )

    ## ================================================================================
    ## Cumulative enrichment barplot
    ## ================================================================================

    connectivity_plot_cumEnrichmentPlot_server(
      "cumEnrichmentPlot",
      pgx,
      reactive(input$cmap_sigdb),
      getConnectivityScores,
      connectivityScoreTable,
      getEnrichmentMatrix,
      getCurrentContrast,
      watermark = WATERMARK
    )


    ## =============================================================================
    ## CONNECTIVITY MAP
    ## =============================================================================
    connectivity_plot_connectivityMap_server(
      "connectivityMap",
      pgx,
      reactive(input$cmap_sigdb),
      getConnectivityScores,
      getEnrichmentMatrix
    )

    connectivityScoreTable2 <- connectivity_table_similarity_scores2_server(
      "connectivityScoreTable2",
      getConnectivityScores = getConnectivityScores
    )

    ## -------------------------------------------------------------------------------
    ## Leading-edge graph
    ## -------------------------------------------------------------------------------

    getLeadingEdgeGraph <- connectivity_plot_leadingEdgeGraph_server(
      "leadingEdgeGraph",
      getConnectivityScores,
      connectivityScoreTable,
      getCurrentContrast,
      getTopProfiles
    )

    ## -------------------------------------------------------------------------------
    ## Enrichment graph
    ## -------------------------------------------------------------------------------

    connectivity_plot_enrichmentGraph_server(
      "enrichmentGraph",
      getLeadingEdgeGraph,
      getConnectivityScores,
      connectivityScoreTable,
      getGSETS,
      cumEnrichmentTable
    )

    ## ======================================================================
    ## Pairs
    ## ======================================================================

    ## ----------------------------------------------------------------------
    ## Scatterplot matrix in plotly
    ##
    ## From: https://plot.ly/r/splom/
    ## ----------------------------------------------------------------------

    connectivity_plot_cmapPairsPlot_server(
      "cmapPairsPlot",
      pgx,
      reactive(input$cmap_contrast),
      reactive(input$cmap_sigdb),
      getConnectivityContrasts,
      getCurrentContrast,
      connectivityScoreTable,
      getConnectivityScores,
      getConnectivityMatrix,
      watermark = WATERMARK
    )

    ## =============================================================================
    ## CONNECTIVITY HEATMAP
    ## =============================================================================
    connectivity_plot_connectivityHeatmap_server(
      "connectivityHeatmap",
      getTopProfiles,
      getConnectivityScores,
      getCurrentContrast
    )
  })
} ## end-of-Board
