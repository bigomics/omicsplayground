##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ConnectivityBoard <- function(
    id,
    auth = NoAuthenticationModule(id = "auth", show_modal = FALSE),
    pgx,
    reload_pgxdir = reactive(auth$user_dir)
    ) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    fullH <- 750 # row height of panel
    tabH <- "70vh"

    infotext <- strwrap(
      "The <strong>Experiment connectivity</strong> module enables users to
      compare their data to other datasets. For the selected contrast, this
      module provides pairwise correlation plots and/or enrichment plots with
      signatures from other data sets. The <strong>Connectivity map</strong>
      shows the similarity of the contrasts profiles as a t-SNE plot.<br><br>
      <br><br><center><iframe width='560' height='315' src='https://www.youtube.com/embed/4-2SkBNcTZk?si=m4qEXCuQJo6o-A9o&amp;start=38' title='YouTube video player' frameborder='0' allow='accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share' referrerpolicy='strict-origin-when-cross-origin' allowfullscreen></iframe></center>"
    )

    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Connectivity Analysis Board</strong>"),
        shiny::HTML(infotext),
        easyClose = TRUE, size = "l"
      ))
    })

    ## update choices upon change of data set
    shiny::observeEvent(pgx$model.parameters$contr.matrix, {
      shiny::req(pgx$model.parameters$contr.matrix)
      ## update contrasts
      comparisons <- playbase::pgx.getContrasts(pgx)
      comparisons <- sort(comparisons[!grepl("^IA:", comparisons)])
      shiny::updateSelectInput(
        session,
        "contrast",
        choices = comparisons,
        selected = head(comparisons, 1)
      )
    })

    shiny::observe({
      shiny::req(pgx$X, pgx$connectivity)
      ## update sigdb choices
      my_sigdb <- "datasets-sigdb.h5"
      computed_sigdb <- NULL ## only precomputed inside PGX object??
      if (dir.exists(SIGDB.DIR) && length(pgx$connectivity) > 0) {
        ## only show if we have the libx h5 files available
        libx_sigdb <- dir(SIGDB.DIR, pattern = ".h5$")
        computed_sigdb <- intersect(names(pgx$connectivity), libx_sigdb)
      }
      available_sigdb <- c(my_sigdb, computed_sigdb)
      shiny::updateSelectInput(session, "sigdb", choices = available_sigdb, selected = my_sigdb)
    })

    shiny::observeEvent(pgx$X, {
      shiny::updateTextAreaInput(
        session,
        inputId = "genelist",
        placeholder = tspan("Paste your gene list", js = FALSE)
      )
    })


    ## ================================================================================
    ## =============================  FUNCTIONS =======================================
    ## ================================================================================

    getConnectivityFilename <- function(sigdb) {
      db1 <- file.path(SIGDB.DIR, sigdb)
      db2 <- file.path(auth$user_dir, sigdb)
      if (file.exists(db1)) {
        return(db1)
      }
      if (file.exists(db2)) {
        return(db2)
      }
      return(NULL)
    }

    #' Get the path/folder to the signature database file.
    #'
    #' @param sigdb signature h5 file
    getConnectivityPath <- function(sigdb) {
      db1 <- file.path(SIGDB.DIR, sigdb)
      db2 <- file.path(auth$user_dir, sigdb)
      if (file.exists(db1)) {
        return(SIGDB.DIR)
      }
      if (file.exists(db2)) {
        return(auth$user_dir)
      }
      return(NULL)
    }

    getConnectivityContrasts <- function(sigdb) {
      if (length(sigdb) == 0 || is.null(sigdb) || sigdb == "") {
        return(NULL)
      }
      cpath <- getConnectivityPath(sigdb)
      playbase::sigdb.getConnectivityContrasts(sigdb, path = cpath)
    }

    getConnectivityMatrix <- function(sigdb, select = NULL, genes = NULL) {
      cpath <- getConnectivityPath(sigdb)
      playbase::sigdb.getConnectivityMatrix(sigdb, select = select, genes = genes, path = cpath)
    }

    getEnrichmentMatrix <- function(sigdb, select = NULL, nc = -1) {
      cpath <- getConnectivityPath(sigdb)
      playbase::sigdb.getEnrichmentMatrix(sigdb,
        select = select, path = cpath,
        which = c("gsea", "rankcor")
      )
    }

    getSignatureMatrix <- function(sigdb) {
      cpath <- getConnectivityPath(sigdb)
      playbase::sigdb.getSignatureMatrix(sigdb, path = cpath)
    }


    ## ================================================================================
    ## ========================= REACTIVE FUNCTIONS ===================================
    ## ================================================================================

    getCurrentContrast <- shiny::reactive({
      shiny::req(pgx$gx.meta, pgx$gset.meta, input$contrast)
      ct <- input$contrast
      meta1 <- pgx$gx.meta$meta
      meta2 <- pgx$gset.meta$meta
      has.contrast <- ct %in% names(meta1) && ct %in% names(meta2)
      shiny::req(has.contrast)

      ## convert to human symbols so we can match different organism
      fc <- meta1[[ct]]$meta.fx
      names(fc) <- rownames(meta1[[ct]])
      fc <- playbase::collapse_by_humansymbol(fc, pgx$genes)

      gs <- meta2[[ct]]$meta.fx
      names(gs) <- rownames(meta2[[ct]])

      list(name = ct, fc = fc, gs = gs)
    })

    observe({
      contr <- getCurrentContrast()
      shiny::req(contr)
      ntop <- as.integer(input$genelist_ntop)
      top50 <- head(names(sort(abs(contr$fc), decreasing = TRUE)), ntop)
      top50 <- paste(top50, collapse = " ")
      updateTextAreaInput(session, "genelist", value = top50)
    })

    cumEnrichmentTable <- shiny::reactive({
      sigdb <- input$sigdb

      shiny::req(sigdb, pgx$connectivity)
      if (!grepl(".h5$", sigdb)) {
        return(NULL)
      }

      df <- getConnectivityScores()
      if (is.null(df)) {
        return(NULL)
      }

      ii <- connectivityScoreTable$rows_all()
      shiny::req(ii)

      sel <- head(df$pathway[ii], 10)
      sigdb <- input$sigdb
      F <- getEnrichmentMatrix(sigdb, select = sel)
      if (is.null(F)) {
        return(NULL)
      }

      ## multiply with sign of enrichment
      rho1 <- df$rho[match(colnames(F), df$pathway)]
      F <- t(t(F) * sign(rho1))
      F <- F[order(-rowMeans(F**2, na.rm = TRUE)), , drop = FALSE]

      ## add current contrast
      contr <- getCurrentContrast()
      shiny::req(contr)
      gx <- contr$gs[match(rownames(F), names(contr$gs))]
      names(gx) <- rownames(F)
      gx[is.na(gx)] <- 0
      F <- cbind(gx, F)
      colnames(F)[1] <- contr$name

      return(F)
    })

    compute_connectivity <- shiny::eventReactive(
      {
        auth$user_dir
        pgx$counts
        pgx$name
        reload_pgxdir()
      },
      {
        ## shiny::req(pgx$connectivity)  ##??

        ## COMPUTE HERE??? or in pgxCompute() !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        pgxdir <- auth$user_dir
        sigdb.file <- file.path(pgxdir, "datasets-sigdb.h5")

        need_update <- playbase::pgxinfo.needUpdate(
          pgxdir,
          check.sigdb = TRUE,
          verbose = FALSE
        )

        if (need_update || !file.exists(sigdb.file)) {
          pgx.showSmallModal("Updating your signature database<br>Please wait...")
          info("[compute_connectivity] calling updateDatasetFolder")
          shiny::withProgress(message = "Updating signature database...", value = 0.33, {
            playbase::pgxinfo.updateDatasetFolder(pgxdir, update.sigdb = TRUE)
          })
          shiny::removeModal(session)
        }

        has.user_sigdb <- "datasets-sigdb.h5" %in% names(pgx$connectivity)
        if (need_update || !has.user_sigdb) {
          user.scores <- NULL
          if (file.exists(sigdb.file)) {
            info("[compute_connectivity] re-computing connectivity scores...")
            pgx.showSmallModal("Computing connectivity scores<br>Please wait...")
            shiny::withProgress(message = "Computing connectivity scores...", value = 0.33, {
              user.scores <- playbase::pgx.computeConnectivityScores(
                pgx, sigdb.file,
                ntop = 50,
                contrasts = NULL,
                remove.le = TRUE
              )
            })
            shiny::removeModal(session)
          }
          pgx$connectivity[["datasets-sigdb.h5"]] <- user.scores
          ## save results back?? but what is the real filename?????
          if (!is.null(pgx$filename)) {
            pgx.filepath <- file.path(pgxdir, basename(pgx$filename))
            try(playbase::pgx.save(shiny::reactiveValuesToList(pgx), file = pgx.filepath)) # on board snap test this fails, wrap in try
          }
        }

        ## return connectivity results object
        pgx$connectivity
      },
      ignoreNULL = TRUE,
      ignoreInit = FALSE
    )

    getConnectivityScores <- shiny::reactive({
      pgx.connectivity <- compute_connectivity()

      sigdb <- input$sigdb
      shiny::req(sigdb)
      all.scores <- pgx.connectivity[[sigdb]]

      ct <- "contrast"
      ct <- input$contrast
      if (!ct %in% names(all.scores)) {
        warning("[getConnectivityScores] ERROR : contrast ", ct, " not in connectivity scores")
        return(NULL)
      }

      scores <- as.data.frame(all.scores[[ct]])

      if (input$abs_score == FALSE) {
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

      if (input$hideclustcontrasts) {
        sel <- grep("cluster:|PC[0-9]+:", scores$pathway, invert = TRUE)
        scores <- scores[sel, , drop = FALSE]
      }

      ## only those in existing database
      sigpath <- getConnectivityPath(sigdb)
      cts <- playbase::sigdb.getConnectivityContrasts(sigdb, path = sigpath)
      scores <- scores[which(rownames(scores) %in% cts), , drop = FALSE]
      scores <- scores[order(-scores$score), , drop = FALSE]

      ## compute leading edge
      no.le <- !("leadingEdge" %in% colnames(scores))
      abs_score <- input$abs_score
      ntop <- 100

      contr <- getCurrentContrast()
      shiny::req(contr)
      fc <- contr$fc
      fc <- fc[order(-abs(fc))]
      fc.up <- head(names(fc[fc > 0]), ntop)
      fc.dn <- head(names(fc[fc < 0]), ntop)

      if (no.le && abs_score == TRUE) {
        ## recreate "leadingEdge" list
        sig <- playbase::sigdb.getSignatureMatrix(sigdb, path = sigpath)
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

    getSelectedGenes <- reactive({
      genes <- input$genelist
      genes <- strsplit(genes, split = " ")[[1]]
      genes
    })

    getTopProfiles <- shiny::reactive({
      ## Get profiles of top-enriched contrasts (not all genes...)
      sigdb <- input$sigdb
      shiny::req(sigdb)

      ii <- connectivityScoreTable$rows_all()
      shiny::req(ii, input$sigdb)

      df <- getConnectivityScores()
      pw <- head(df$pathway[ii], 100)

      contr <- getCurrentContrast()
      fc <- contr$fc
      ngenes <- min(500, length(fc))
      top.genes <- head(names(sort(-abs(fc))), ngenes)
      top.genes <- unique(c(top.genes, sample(names(fc), ngenes))) ## add some random

      F <- getConnectivityMatrix(sigdb, select = pw, genes = top.genes)
      pw <- intersect(pw, colnames(F))
      F <- F[, pw, drop = FALSE]
      return(F)
    })

    getSelectedProfiles <- shiny::reactive({
      ## Get profiles of top-enriched contrasts (not all genes...)
      sigdb <- input$sigdb
      shiny::req(sigdb)

      ii <- connectivityScoreTable$rows_all()
      shiny::req(ii, input$sigdb)

      df <- getConnectivityScores()
      pw <- head(df$pathway[ii], 100)

      selected_genes <- getSelectedGenes()
      F <- getConnectivityMatrix(sigdb, select = pw, genes = selected_genes)
      pw <- intersect(pw, colnames(F))
      F <- F[, pw, drop = FALSE]

      return(F)
    })


    ## ============================================================================
    ## FC correlation/scatter plots
    ## ============================================================================
    connectivityScoreTable <- NULL
    connectivityFoldchangeTable <- NULL
    getLeadingEdgeGraph <- NULL

    connectivity_plot_FCFCplots_server(
      "FCFCplots",
      pgx = pgx,
      r_contrast = reactive(input$contrast),
      getCurrentContrast = getCurrentContrast,
      getTopProfiles = getTopProfiles,
      getConnectivityScores = getConnectivityScores,
      watermark = WATERMARK
    )

    connectivityScoreTable <- connectivity_table_similarity_scores_server(
      "connectivityScoreTable",
      getConnectivityScores = getConnectivityScores,
      columns = c("pathway", "score", "rho", "NES", "odd.ratio", "tau"),
      height = "200px" ## scrollY height
    )

    ## ================================================================================
    ## Cumulative FC barplot
    ## ================================================================================

    connectivity_plot_cumFCplot_server(
      "cumFCplot",
      ## getTopProfiles,
      getProfiles = getSelectedProfiles,
      getConnectivityScores = getConnectivityScores,
      getCurrentContrast = getCurrentContrast
    )

    ## ================================================================================
    ## Cumulative enrichment barplot
    ## ================================================================================

    connectivity_plot_cumEnrichmentPlot_server(
      "cumEnrichmentPlot",
      pgx = pgx,
      sigdb = reactive(input$sigdb),
      getConnectivityScores = getConnectivityScores,
      connectivityScoreTable = connectivityScoreTable,
      getEnrichmentMatrix = getEnrichmentMatrix,
      getCurrentContrast = getCurrentContrast,
      watermark = WATERMARK
    )

    ## =============================================================================
    ## CONNECTIVITY MAP
    ## =============================================================================
    connectivity_plot_connectivityMap_server(
      "connectivityMap",
      pgx,
      reactive(getConnectivityFilename(input$sigdb)),
      getConnectivityScores,
      getEnrichmentMatrix
    )

    connectivityFoldchangeTable <- connectivity_table_foldchange_server(
      id = "connectivityFoldchangeTable",
      pgx = pgx,
      getConnectivityScores = getConnectivityScores,
      columns = c("pathway", "score", "rho", "NES", "padj"),
      getProfiles = getSelectedProfiles,
      getConnectivityMatrix = getConnectivityMatrix,
      sigdb = reactive(input$sigdb),
      height = "550px"
    )

    ## -------------------------------------------------------------------------------
    ## Leading-edge graph
    ## -------------------------------------------------------------------------------

    getLeadingEdgeGraph <- connectivity_plot_leadingEdgeGraph_server(
      "leadingEdgeGraph",
      getConnectivityScores = getConnectivityScores,
      connectivityScoreTable = connectivityScoreTable,
      getCurrentContrast = getCurrentContrast,
      getProfiles = getTopProfiles
    )

    ## -------------------------------------------------------------------------------
    ## Enrichment graph
    ## -------------------------------------------------------------------------------

    connectivity_plot_enrichmentGraph_server(
      "enrichmentGraph",
      getLeadingEdgeGraph = getLeadingEdgeGraph,
      getConnectivityScores = getConnectivityScores,
      connectivityScoreTable = connectivityScoreTable,
      cumEnrichmentTable = cumEnrichmentTable
    )

    ## ======================================================================
    ## Scatter SPLOM
    ## ======================================================================

    ## ----------------------------------------------------------------------
    ## Scatterplot matrix in plotly
    ##
    ## From: https://plot.ly/r/splom/
    ## ----------------------------------------------------------------------

    connectivity_plot_scatterPlot_server(
      "scatterPlot",
      pgx = pgx,
      r_sigdb = reactive(input$sigdb),
      getConnectivityContrasts = getConnectivityContrasts,
      getCurrentContrast = getCurrentContrast,
      connectivityScoreTable = connectivityScoreTable,
      getConnectivityScores = getConnectivityScores,
      getConnectivityMatrix = getConnectivityMatrix,
      watermark = WATERMARK
    )

    ## =============================================================================
    ## CONNECTIVITY HEATMAP
    ## =============================================================================
    connectivity_plot_connectivityHeatmap_server(
      "connectivityHeatmap",
      getProfiles = getSelectedProfiles,
      getConnectivityScores = getConnectivityScores,
      getCurrentContrast = getCurrentContrast
    )
  }) ## end of moduleserver
} ## end-of-Board
