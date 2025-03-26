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
      "Time clustering" = list(
        enable = NULL,
        disable = c("contrast","groupvar")
      ),
      "Features" = list(
        enable = NULL,
        disable = c("module","knn","maxfeatures")
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
      vars <- sort(colnames(pgx$samples))
      timevars <- unique( c(grep("time|second|minute|day|week|month|year", vars, value=TRUE, ignore.case=TRUE), vars))
      shiny::updateSelectInput(session, "timevar", choices=timevars, selected=timevars[1])
      groups <- c("<none>",vars)
      shiny::updateSelectInput(session, "groupvar", choices=groups, selected=groups[1])
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
      timevar="week"
      timevar <- input$timevar
      time  <- pgx$samples[,timevar]
      X <- pgx$X
      sd <- matrixStats::rowSds(X, na.rm = TRUE)
      if (any(sd == 0)) X <- X + runif(length(X), 0, 1e-5)
      cX <- t(scale(t(X)))
      timeX <- t( playbase::rowmean(t(cX), time))
      clust <- playbase::pgx.FindClusters(t(timeX), method="kmeans")[[1]]
      rownames(clust) <- rownames(timeX)
      colors <- clust[,paste0("kmeans.",knn)]
      ##colors <- WGCNA::standardColors(435)[colors]
      colors <- paste0("M",colors)
      
      ## update selectinput
      modulenames <- sort(unique(colors))
      shiny::updateSelectInput(
        session,
        "module",
        choices = modulenames,
        selected = head(modulenames, 3)
      )

      ## Statistical methods
      ts.gx.methods <- c("trend.limma", "deseq2.lrt", "edger.lrt", "edger.qlf")
      gx.methods <- intersect(ts.gx.methods, colnames(pgx$gx.meta$meta[[1]]$fc))
      if (length(gx.methods) == 0)
        gx.methods <-  colnames(pgx$gx.meta$meta[[1]]$fc)
      shiny::updateCheckboxGroupInput(
        session,
        "gx_statmethod",
        choices = gx.methods,
        selected = gx.methods[1:3]
      )
      
      res <- list(X = timeX, colors = colors)
      return(res)

    })
    
    timeseries_filtered <- shiny::reactive({

      res <- timeseries_full()
      shiny::req(res)
      
      ##minKME=0.8;mergeCutHeight=0.15;minmodsize=20;ntop=10      
      filtered.colors <- playbase::wgcna.filterColors(
        res$X,
        res$colors,
        minKME=0.8,
        mergeCutHeight=0.05,
        minmodsize = 10,
        ntop=100
      )

      jj <- which(!filtered.colors %in% c(NA,0,"---","grey"))
      xx <- res$X[jj,]
      filtered.colors <- filtered.colors[jj]

      time <- colnames(xx)
      res <- list(X = xx, time = time, colors = filtered.colors)
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

    TimeSeriesBoard.features_server(
      id = "features",
      pgx = pgx,
      data = timeseries_full,
      contrast = shiny::reactive(input$contrast),      
      timevar = shiny::reactive(input$timevar),
      groupvar = shiny::reactive(input$groupvar),
      gx_statmethod = shiny::reactive(input$gx_statmethod),
      watermark = WATERMARK
    )

  }) ## end of moduleServer
} ## end of Board
