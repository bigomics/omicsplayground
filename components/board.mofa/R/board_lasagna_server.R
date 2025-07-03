##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

LasagnaBoard <- function(id, pgx, board_observers = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    fullH <- 700 ## full height of page
    rowH1 <- 250 ## row 1 height
    rowH2 <- 440 ## row 2 height

    infotext <- tspan("<b>Weighted gene co-expression network analysis (WGCNA)</b> is a systems biology method for describing the correlation patterns among genes across microarray samples. Weighted correlation network analysis can be used for finding clusters (modules) of highly correlated genes, for summarizing such clusters using the module eigengene or an intramodular hub gene, for relating modules to one another and to external sample traits (using eigengene network methodology), and for calculating module membership measures. Correlation networks facilitate network based gene screening methods that can be used to identify candidate biomarkers or therapeutic targets.

<p>References:<br>
<ol>
<li>Langfelder, P. and Horvath, S., 2008. WGCNA: an R package for weighted correlation network analysis. BMC bioinformatics, 9(1), p.559.
<li>Zhang, B. and Horvath, S., 2005. A general framework for weighted gene co-expression network analysis. Statistical applications in genetics and molecular biology, 4(1).
</ol>
", js = FALSE)


    ## ============================================================================
    ## ============================ OBSERVERS =====================================
    ## ============================================================================
    my_observers <- list()

    infotext <-
      '<center><iframe width="1120" height="630" src="https://www.youtube.com/embed/rRIRMW_RRS4"
        title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write;
        encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></center>'

    my_observers[[1]] <- shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>LASAGNA Analysis Board</strong>"),
        shiny::HTML(infotext),
        size = "xl",
        easyClose = TRUE
      ))
    })

    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Multi-layer model" = list(disable = c("mpartite_options")),
      "Multi-partite graph" = list(disable = c())
    )

    my_observers[[2]] <- shiny::observeEvent( input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })

    
    ## add to list global of observers. suspend by default.
    my_observers <- my_observers[!sapply(my_observers,is.null)]
    # lapply( my_observers, function(b) b$suspend() )
    if(!is.null(board_observers)) board_observers[[id]] <- my_observers
    
    ## ============================================================================
    ## ============================ REACTIVES =====================================
    ## ============================================================================

    shiny::observeEvent( pgx$mofa, {
      
      shiny::validate( shiny::need( !is.null(pgx$mofa), "missing MOFA slot"))
      ##shiny::validate( shiny::need( !is.null(pgx$mofa$lasagna), "missing LASAGNA slot"))
      
      ## update factors in selectInput
      ct1 <- colnames(pgx$mofa$contrasts)      
      ct2 <- colnames(playbase::pgx.getMetaMatrix(pgx)$fc)
      contrasts <- intersect(ct1, ct2)
      contrasts <- contrasts[!grepl("^IA",contrasts)]  ## no interaction contrasts
      updateSelectInput(session, "contrast", choices = contrasts,
        selected = contrasts[1])
      
      datatypes <- names(pgx$mofa$xx)
      updateSelectInput(session, "datatypes", choices = datatypes,
        selected = datatypes)
      
    }, ignoreNULL=FALSE)

    lasagna_model <- shiny::eventReactive( {
      list( input$updateplots, pgx$X ) 
    }, {
      shiny::validate( shiny::need( !is.null(pgx$mofa), "missing MOFA slot"))
      ##shiny::validate( shiny::need( !is.null(pgx$mofa$lasagna), "missing LASAGNA slot"))
      shiny::validate( shiny::need( pgx$datatype=="multi-omics", "not multi-omics data"))

      shiny::req(pgx$X)
      shiny::req(pgx$mofa)
      
      ##pgx <- pgx.load("~/Playground/omicsplayground/data/mox-brca.pgx")
      ##pgx <- pgx.load("~/Playground/omicsplayground/data/mox-maartenAML.pgx")
      xdata <- playbase::mofa.split_data(pgx$X)
      dt <- c("mir","gx","px")
      dt <- input$datatypes
      if(!is.null(dt) && length(dt)>0) {
        dt <- intersect(dt, names(xdata))
        if(length(dt)) xdata <- xdata[dt]
      }

      moxdata <- list(
        X = xdata,
        samples = pgx$samples,
        contrasts = pgx$contrasts
      )
      dbg("[LasagnaBoard:lasagna_model] creating model...")
      res <- playbase::lasagna.create_model(moxdata, pheno="contrasts",
        ntop=1000, nc=20, add.sink=TRUE, intra=TRUE)
      names(res)
      res$layers
      
      layers <- setdiff(res$layers, c("SOURCE","SINK"))
      sel <- intersect(layers, names(pgx$mofa$posx))
      res$posx <- pgx$mofa$posx[sel]
      res$posf <- pgx$mofa$posf[sel]

      if(0) {
        res.svd <- svd(scale(res$Y))
        pos1 <- res.svd$v[,1:2]
        rownames(pos1) <- colnames(res$Y)
        res$posf[['PHENO']] <- pos1
        pos2 <- res.svd$u[,1:2]
        rownames(pos2) <- rownames(res$Y)
        res$posx[['PHENO']] <- pos2
      }
            
      return(res)
    }, ignoreNULL=FALSE)

    data <- reactive({
      shiny::req(input$contrast)

      res <- lasagna_model()      
      pheno <- input$contrast
      value <- input$node_value      
      
      graph <- playbase::lasagna.set_weights(
        res, pheno,
        value=value, 
        min_rho=0.1,
        max_edges=1000,
        fc.weight=TRUE,
#        sp.weight=FALSE,
        sp.weight = input$sp_weight,        
        prune = FALSE
      )
      res$graph <- graph
      res
    })
    

    ## ==========================================================================
    ## ========================== BOARD FUNCTIONS ===============================
    ## ==========================================================================

    
    ## ==========================================================================
    ## =========================== MODULES ======================================
    ## ==========================================================================
    
    mofa_plot_lasagna3D_server(
      "lasagna",
      data = data,
      pgx = pgx,
      watermark = WATERMARK
    )

    mofa_plot_lasagna_partite_server(
      "lasagnaPartite",
      data = data,
      pgx = pgx,
      input_datatypes = reactive(input$datatypes),
      input_minrho = reactive(input$minrho),
      input_edgesign = reactive(input$edgesign),      
      #input_labeltype = reactive(input$labeltype),
      input_nodevalue = reactive(input$node_value),            
      watermark = WATERMARK
    )
    
    mofa_plot_clustering_server(
      "clusters",
      data = data,
      pgx = pgx,
      type = "features",
      input_contrast = reactive(input$contrast),
      watermark = WATERMARK
    )

    return(NULL)
  })
} ## end of Board
