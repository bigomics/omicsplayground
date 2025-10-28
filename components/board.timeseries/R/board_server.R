##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


##' TimeSeries board server module
##'
##' .. content for \details{} ..
##' @title
##' @param id
##' @param pgx
##' @return
##' @author kwee
TimeSeriesBoard <- function(id,
                            pgx,
                            labeltype = shiny::reactive("feature")) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    fullH <- 850 ## full height of page

    clust_infotext <-
      '<center><iframe width="1120" height="630" src="https://www.youtube.com/embed/hyDEk_MCaTk"
       title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write;
       encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></center>'

    ## ============================================================================
    ## ======================== OBSERVERS =========================================
    ## ============================================================================

    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Clustering" = list(
        enable = NULL,
        disable = c("contrast")
      ),
      "Statistics" = list(
        enable = NULL,
        disable = c("timefactor", "module", "knn")
      )
    )

    shiny::observeEvent(input$tabs1, {
      bigdash::update_tab_elements(input$tabs1, tab_elements)
    })

    # Board info
    shiny::observeEvent(input$board_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>TimeSeries Board</strong>"),
        shiny::HTML(clust_infotext),
        easyClose = TRUE, size = "xl"
      ))
    })

    shiny::observeEvent(pgx$samples, {
      ## set time variable
      vars <- sort(colnames(pgx$samples))
      timevars <- unique(c(grep(playbase::get_timevars(), vars, value = TRUE, ignore.case = TRUE), vars))
      ## filter timevars to only include those with more than 1 unique value
      valid_timevars <- timevars[sapply(timevars, function(var) {
        var %in% colnames(pgx$samples) && length(unique(pgx$samples[, var])) > 1
      })]
      timevars <- valid_timevars
      shiny::updateSelectInput(session, "timevar", choices = timevars, selected = timevars[1])

      ## set available contrasts
      contrasts <- playbase::pgx.getContrasts(pgx)
      contrasts <- contrasts[!grepl("^IA:", contrasts)]
      shiny::updateSelectInput(session, "contrast", choices = contrasts, selected = contrasts[1])
    })

    ## ===========================================================================
    ## ======================== REACTIVES ========================================
    ## ===========================================================================

    timeseries_full <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(input$timevar)
      ## shiny::req(input$knn)

      X <- pgx$X
      if (any(is.na(X))) X <- playbase::imputeMissing(X, method = "SVD2")
      sdx <- matrixStats::rowSds(X, na.rm = TRUE)
      if (any(sdx == 0)) X <- X + runif(length(X), 0, 1e-5)

      ## pre-filter by SD.
      X <- playbase::mofa.topSD(X, 4000)

      ## collapse by time variable
      timevar <- "time"
      timevar <- input$timevar
      time <- pgx$samples[, timevar]
      timeX <- t(playbase::rowmean(t(X), time))
      cX <- t(scale(t(X)))
      timeZ <- t(playbase::rowmean(t(cX), time))

      ## compute gene clusters. Actually this should be precomputed in
      ## pgx$cluster.genes???
      clust <- playbase::pgx.FindClusters(
        t(timeZ),
        km.sizes = c(4, 6, 9, 12),
        method = "kmeans"
      )[[1]]
      rownames(clust) <- rownames(timeZ)

      knn <- as.integer(input$knn)
      modules <- clust[, paste0("kmeans.", knn)]
      modules <- paste0("M", modules)
      names(modules) <- rownames(clust)

      ## compute geneset enrichment
      gset.rho <- NULL
      if (!is.null(pgx$gsetX) && nrow(pgx$gsetX)) {
        mX <- playbase::rowmean(playbase::rowscale(X), modules)
        gset.rho <- cor(t(pgx$gsetX), t(mX))
      }

      res <- list(X = timeX, Z = timeZ, modules = modules, gset.rho = gset.rho)
      return(res)
    })

    timeseries_filtered <- shiny::reactive({
      res <- timeseries_full()
      shiny::req(res)

      if (input$filtermodules) {
        filtered.modules <- playbase::wgcna.filterColors(
          res$Z,
          res$modules,
          minKME = 0.3,
          mergeCutHeight = 0.05,
          minmodsize = 10,
          ntop = 1000
        )
      } else {
        filtered.modules <- res$modules
      }
      zz <- res$Z[, ]
      xx <- res$X[, ]
      time <- colnames(res$Z)

      ## remove grey group?
      if (1) {
        jj <- which(!filtered.modules %in% c(NA, 0, "---", "grey", "T0"))
        if (length(jj) > 0) {
          zz <- zz[jj, ]
          xx <- xx[jj, ]
          filtered.modules <- filtered.modules[jj]
        }
      }

      ## update selectinput
      modulenames <- sort(unique(filtered.modules))
      shiny::updateSelectInput(
        session,
        "module",
        choices = modulenames,
        selected = modulenames[1]
      )

      res <- list(
        X = xx, Z = zz, time = time,
        modules = filtered.modules,
        gset.rho = res$gset.rho
      )
      return(res)
    })

    ## ============================================================================
    ## =============================== MODULES ====================================
    ## ============================================================================

    TimeSeriesBoard.parcoord_server(
      id = "parcoord",
      pgx = pgx,
      data = timeseries_filtered,
      select_module = reactive(input$module),
      watermark = WATERMARK
    )

    TimeSeriesBoard.clustering_server(
      id = "clustering",
      pgx = pgx,
      data = timeseries_filtered,
      timefactor = reactive(input$timefactor),
      watermark = WATERMARK
    )

    TimeSeriesBoard.enrichment_server(
      id = "enrichment",
      pgx = pgx,
      data = timeseries_filtered,
      select_module = reactive(input$module),
      watermark = WATERMARK
    )

    TimeSeriesBoard.features_server(
      id = "features",
      pgx = pgx,
      data = timeseries_full,
      contrast = shiny::reactive(input$contrast),
      timevar = shiny::reactive(input$timevar),
      watermark = WATERMARK
    )
  }) ## end of moduleServer
} ## end of Board
