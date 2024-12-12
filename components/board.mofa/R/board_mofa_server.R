##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


MofaBoard <- function(id, pgx) {
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

    
    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Overview" = list(disable = c("selected_factor","selected_pheno","show_types")),
      "Response" = list(disable = c("show_types","selected_pheno")),
      "Factor" = list(disable = c())
    )
    shiny::observeEvent( input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })

    ## =======================================================================
    ## ======================= PRECOMPUTE FUNCTION ===========================
    ## =======================================================================
    
    mofa <- shiny::eventReactive(input$compute, {

      shiny::req(pgx$X, input$kernel, input$numfactors)
            
      kernel <- input$kernel
      
      ##source("~/Playground/playbase/dev/include.R",chdir=TRUE)
      pgx.showSmallModal(paste("Calculating",kernel,"...<br>please wait"))
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = paste("Calculating",kernel,"..."), value = 0)

      message(paste("[mofa] >>> Calculating ",kernel,"..."))

      ## Select dataset: current PGX or one of the example datasets
      dataset <- input$dataset
      message("[mofa] dataset = ", dataset)
      shiny::req(dataset)

      numfactors <- as.integer(input$numfactors)
      dbg("[MofaBoard] numfactors = ", numfactors)      

      add_gsets <- input$add_gsets
      mofa <- NULL
      has.mofa <- ("mofa" %in% names(pgx))

      dbg("[MofaBoard] names(pgx$mofa) = ", names(pgx$mofa))
      
      if( dataset == "<this pgx>" && has.mofa) {

        dbg("[MofaBoard] using internal MOFA for PGX")
        mofa <- pgx$mofa
        
      } else if( dataset == "<this pgx>" && !has.mofa) {

        dbg("[MofaBoard] NO MOFA!!! calculating MOFA for PGX")
        numfactors <- min( numfactors, min(dim(pgx$X)) )
        
        pgx$mofa <- playbase::pgx.compute_mofa(
          pgx,
          kernel = kernel,
          numfactors = numfactors,
          add_gsets = add_gsets)
        
        
      } else {

        ##--------------------dev only ----------------------
        ##--------------------dev only ----------------------
        ##--------------------dev only ----------------------
        
        dbg("[MofaBoard] using EXAMPLE data")
        
        data <- playbase::mofa.exampledata(dataset)
        numfactors <- min( numfactors, min(dim(data$X[[1]])) )

        ## add geneset layers if asked
        if(add_gsets) {
          gsetX <- playbase::mofa.add_genesets(data$X, GMT=NULL) 
          names(gsetX)
          if(!is.null(gsetX$gset.px)) {
            gsetX$gset.px <- gsetX$gset.px[grep("^GO",rownames(gsetX$gset.px)),]
          }
          if(!is.null(gsetX$gset.gx)) {          
            gsetX$gset.gx <- gsetX$gset.gx[grep("^GO",rownames(gsetX$gset.gx)),]
          }
          if(!is.null(gsetX$gset.mx)) {                    
            gsetX$gset.mx <- gsetX$gset.mx[grep("SMP[0-9]*",rownames(gsetX$gset.mx)),]
          }
          lapply(gsetX,dim)
          data$X <- c( data$X, gsetX)
        }
      
        mofa <- playbase::mofa.compute(
          data$X,
          data$samples,
          contrasts = data$contrasts,
          pheno = NULL,
          kernel = kernel,
          scale_views = TRUE, 
          ntop = 1000,
          max_iter = 200,
          num_factors = numfactors)

        ## LASAGNA computation
        message("computing LASAGNA model...")            
        las.data <- list(
          X = data$X,
          samples = data$samples,
          contrasts = data$contrasts
        )
        lasagna <- playbase::lasagna.create_model(
          las.data, pheno="contrasts", ntop=1000, nc=20,
          use.graphite = FALSE)
        
        ## pre-compute cluster positions
        xx <- playbase::mofa.split_data(lasagna$X)
        xx <- xx[names(data$X)]
        mofa$posx <- playbase::mofa.compute_clusters(xx, along="samples")
        mofa$posf <- playbase::mofa.compute_clusters(xx, along="features")       
        mofa$lasagna <- lasagna
        pgx$mofa <- mofa
      } 
      ##--------------------dev only ----------------------
      ##--------------------dev only ----------------------
      ##--------------------dev only ----------------------

      
      mofa <- pgx$mofa
      names(mofa)
      
      ## update factors in selectInput
      factors <- colnames(mofa$F)
      dtypes <- names(mofa$ww)
      contrasts <- colnames(mofa$contrasts)
      phenotypes <- colnames(mofa$samples)
      updateSelectInput(session, "selected_factor", choices = factors,
                        selected = factors[1])
      updateSelectInput(session, "show_types", choices = dtypes,
                        selected = dtypes)
      shiny::removeModal()
      
      message("[mofa] compute MOFA done!")
      
      mofa
    }, ignoreNULL = TRUE)


    ## ========================================================================
    ## ======================= OBSERVE FUNCTIONS ==============================
    ## ========================================================================

    infotext <-
      '<center><iframe width="1120" height="630" src="https://www.youtube.com/embed/rRIRMW_RRS4"
        title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write;
        encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></center>'

    shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>WGCNA Analysis Board</strong>"),
        shiny::HTML(infotext),
        size = "xl",
        easyClose = TRUE
      ))
    })

    
    ## ========================================================================
    ## ======================= PLOTTING FUNCTIONS =============================
    ## ========================================================================

    ## ----------------------------------------
    ## --------- enrichment table -------------
    ## ----------------------------------------


    ## ========================================================================
    ## =========================== MODULES ====================================
    ## ========================================================================

    # Gene dendrogram and gene modules
    mofa_plot_variance_server(
      "factorxview",
      type = "factorxview",
      mofa = mofa,
      watermark = WATERMARK
    )

    mofa_plot_variance_server(
      "variance_view",
      type = "view",
      mofa = mofa,      
      watermark = WATERMARK
    )
    
    mofa_plot_variance_server(
      "variance_factor",
      type = "factor",
      mofa = mofa,      
      watermark = WATERMARK
    )

    mofa_plot_weights_server(
      "weights",
      mofa = mofa,
      input_factor = reactive(input$selected_factor),
      show_types = reactive(input$show_types),
      watermark = WATERMARK
    )
    
    mofa_plot_enrichment_server(
      "enrichment",
      gsea = reactive({ mofa()$gsea }),
      input_k = reactive(input$selected_factor),
      ntop = 15,
      watermark = WATERMARK
    )

    mofa_plot_loadingheatmap_server(
      "loading_heatmap",
      mofa = mofa,
      input_factor = reactive(input$selected_factor),
      watermark = WATERMARK
    )

    mofa_plot_factortrait_server(
      "factortrait",
      mofa = mofa,
      watermark = WATERMARK
    )

    mofa_plot_factorheatmap_server(
      "factor_heatmap",
      mofa = mofa,
      watermark = WATERMARK
    )

    mofa_plot_dendrogram_server(
      "dendrogram",
      mofa = mofa,      
      watermark = WATERMARK
    )

    mofa_plot_boxplots_server(
      "boxplots",
      mofa = mofa,
      input_factor = reactive(input$selected_factor),
      watermark = WATERMARK
    )

    mofa_plot_factorgraph_server(
      "factorgraph",
      mofa = mofa,
      watermark = WATERMARK
    )

    mofa_plot_modulegraph_server(
      "modulegraph",
      mofa = mofa,
      input_k = reactive(input$selected_factor),
      input_pheno = reactive(input$selected_pheno),
      filter_types = reactive(input$show_types),
      watermark = WATERMARK
    )

    mofa_plot_moduleheatmap_server(
      "integrated_heatmap",
      mofa = mofa,
      input_factor = reactive(NULL),      
      watermark = WATERMARK
    )

    mofa_plot_moduleheatmap_server(
      "module_heatmap",
      mofa = mofa,
      input_factor = reactive(input$selected_factor),      
      show_types = reactive(input$show_types),
      watermark = WATERMARK
    )
    
    mofa_plot_centrality_server(
      "centrality",
      mofa = mofa,
      input_factor = reactive(input$selected_factor),
      ## input_pheno = reactive(input$selected_pheno),
      show_types = reactive(input$show_types),
      watermark = WATERMARK
    )
    
    ## # Table Modules
    ## mofa_table_mgsea_server(
    ##   "gsea_table1",
    ##   mofa = mofa,
    ##   datatype = 1,
    ##   input_factor = reactive(input$selected_factor)            
    ## )


    return(NULL)
  })
} ## end of Board
