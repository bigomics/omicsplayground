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
    
    mofa <- shiny::eventReactive(
      list(
        pgx$X,
        pgx$mofa
      ), {
        shiny::req(pgx$X)
        dbg("[MofaBoard] *** eventReactive triggered ***")

        mofa <- NULL
        has.mofa <- ("mofa" %in% names(pgx)) && !is.null(pgx$mofa)     
        shiny::validate( shiny::need( has.mofa, "No MOFA slot in object. Please recompute MOFA"))
        
#        if(!has.mofa) {
#          shinyalert::shinyalert("Warning","No MOFA slot in object. Please recompute MOFA.")
#          return(NULL)
#        }

        dbg("[MofaBoard] *** using precomputed MOFA results ***")
        dbg("[MofaBoard] names(input) = ", names(input))
        
        mofa <- pgx$mofa
        
        ## update factors in selectInput
        factors <- colnames(mofa$F)
        dtypes <- names(mofa$ww)
        sel.dtypes <- grep("^gset",dtypes,value=TRUE,invert=TRUE)
        contrasts <- colnames(mofa$contrasts)
        phenotypes <- colnames(mofa$samples)
        updateSelectInput(session, "selected_factor", choices = factors,
          selected = factors[1])
        updateSelectInput(session, "show_types", choices = dtypes,
          selected = sel.dtypes)

        dbg("[MofaBoard]  names(mofa) = ", names(mofa))
        dbg("[MofaBoard] eventReactive done!")
        return( mofa )
      }
    )

    shiny::observeEvent(
      list(
        input$compute
      ), {

      shiny::req(pgx$X, input$kernel, input$numfactors)

      dbg("[mofa] input$compute =  ",input$compute)
      if(input$compute==0) return(NULL)
      
      kernel <- input$kernel
      dbg("[mofa] kernel =  ",kernel,"...")      

      pgx.showSmallModal(paste("Calculating",kernel,"...<br>please wait"))
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = paste("Calculating",kernel,"..."), value = 0)
      
      numfactors <- as.integer(input$numfactors)
      numfactors <- min( numfactors, min(dim(pgx$X)) )        
      dbg("[MofaBoard] numfactors = ", numfactors)      

      dbg("[MofaBoard] *** recalculating MOFA ***")
      mofa <- playbase::pgx.compute_mofa(
          pgx,
          kernel = kernel,
          numfactors = numfactors,
          add_gsets = input$add_gsets)
            
      ## ## update factors in selectInput
      ## factors <- colnames(mofa$F)
      ## dtypes <- names(mofa$ww)
      ## sel.dtypes <- grep("^gset",dtypes,value=TRUE,invert=TRUE)
      ## contrasts <- colnames(mofa$contrasts)
      ## phenotypes <- colnames(mofa$samples)
      ## updateSelectInput(session, "selected_factor", choices = factors,
      ##                   selected = factors[1])
      ## updateSelectInput(session, "show_types", choices = dtypes,
      ##                   selected = sel.dtypes)
      
      shiny::removeModal()
      message("[mofa] compute MOFA done!")
      
      pgx$mofa <- mofa  ## should trigger new mofa
      
    }, ignoreNULL = FALSE)


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
      pgx = pgx,
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
