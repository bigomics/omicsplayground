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
                            labeltype = shiny::reactive("feature"),
                            board_observers = NULL
                            ) {

  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    fullH <- 850 ## full height of page

    clust_infotext <-
      '<center><iframe width="1120" height="630" src="https://www.youtube.com/embed/hyDEk_MCaTk"
       title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write;
       encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></center>'

    ## ===================================================================================
    ## ======================== OBSERVERS ================================================
    ## ===================================================================================

    my_observers <- list()

    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Clustering" = list(
        enable = NULL,
        disable = c("contrast")
      ),
      "Statistics" = list(
        enable = NULL,
        disable = c("timefactor","module","knn","maxfeatures")
      )
    )

    my_observers[[1]] <- shiny::observeEvent(input$tabs1, {
      bigdash::update_tab_elements(input$tabs1, tab_elements)
    })

    # Board info
    my_observers[[2]] <- shiny::observeEvent(input$board_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>TimeSeries Board</strong>"),
        shiny::HTML(clust_infotext),
        easyClose = TRUE, size = "xl"
      ))
    })

    my_observers[[3]] <- shiny::observeEvent( pgx$samples, {

      ## set time variable
      vars <- sort(colnames(pgx$samples))
      timevars <- unique( c(grep("time|second|minute|day|week|month|year", vars, value=TRUE, ignore.case=TRUE), vars))
      shiny::updateSelectInput(session, "timevar", choices=timevars, selected=timevars[1])

      ## set available contrasts
      contrasts <- playbase::pgx.getContrasts(pgx)
      shiny::updateSelectInput(session, "contrast", choices=contrasts, selected=contrasts[1])
    })

    ## assign to global list of observers. suspend by default.
    # lapply( my_observers, function(b) b$suspend() )
    board_observers[[id]] <- my_observers

    ## ===========================================================================
    ## ======================== REACTIVES ========================================
    ## ===========================================================================
    
    timeseries_full <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(input$timevar)
      shiny::req(input$knn)
      knn=7
      knn <- input$knn

      ## collapse by time variable
      timevar="time"
      timevar <- input$timevar
      time  <- pgx$samples[,timevar]
      X <- pgx$X
      sd <- matrixStats::rowSds(X, na.rm = TRUE)
      if (any(sd == 0)) X <- X + runif(length(X), 0, 1e-5)
      timeX <- t(playbase::rowmean(t(X), time))
      
      cX <- t(scale(t(X)))
      timeZ <- t( playbase::rowmean(t(cX), time))
      clust <- playbase::pgx.FindClusters(t(timeZ), method="kmeans")[[1]]
      rownames(clust) <- rownames(timeZ)
      modules <- clust[,paste0("kmeans.",knn)]
      modules <- paste0("M",modules)
      
      ## compute geneset enrichment
      mX <- playbase::rowmean( playbase::rowscale(X), modules)
      gset.rho <- cor(t(pgx$gsetX), t(mX))
      
      ## update selectinput
      modulenames <- sort(unique(modules))
      shiny::updateSelectInput(
        session,
        "module",
        choices = modulenames,
        selected = modulenames[1]
      )
            
      res <- list(X = timeX, Z = timeZ, modules = modules, gset.rho = gset.rho)
      return(res)

    })
    
    timeseries_filtered <- shiny::reactive({

      res <- timeseries_full()
      shiny::req(res)
      
      ##minKME=0.8;mergeCutHeight=0.15;minmodsize=20;ntop=10      
      filtered.modules <- playbase::wgcna.filterColors(
        res$Z,
        res$modules,
        minKME=0.8,
        mergeCutHeight=0.05,
        minmodsize = 10,
        ntop=100
      )

      ## set remove garbage group?
      if(1) {
        jj <- which(!filtered.modules %in% c(NA,0,"---","grey"))
        zz <- res$Z[jj,]
        xx <- res$X[jj,]        
        filtered.modules <- filtered.modules[jj]
      }

      time <- colnames(zz)
      res <- list(X=xx, Z = zz, time = time, modules = filtered.modules,
        gset.rho = res$gset.rho )
      return(res)

    })
    
    ## ===============================================================================
    ## =============================== MODULES =======================================
    ## ===============================================================================
    
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
#      data = timeseries_full,
      contrast = shiny::reactive(input$contrast),      
      timevar = shiny::reactive(input$timevar),
      watermark = WATERMARK
    )

  }) ## end of moduleServer
} ## end of Board
