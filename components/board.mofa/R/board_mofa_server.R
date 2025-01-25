##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


MofaBoard <- function(id, pgx, board_observers = NULL) {
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
      "Overview" = list(disable = c("selected_factor","selected_module","show_types")),
      "Response" = list(disable = c("show_types","selected_module")),
      "Weights" = list(disable = c("selected_module")),
      "Enrichment" = list(disable = c("selected_module")),
      "gsetMOFA" = list(disable = c("show_types","selected_factor"))      
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

        mofa <- NULL
        has.mofa <- ("mofa" %in% names(pgx)) && !is.null(pgx$mofa)     
        shiny::validate( shiny::need( has.mofa, "No MOFA slot in object. Please recompute MOFA"))        
        mofa <- pgx$mofa

        if(!"gset.mofa" %in% names(mofa)) {
          progress <- shiny::Progress$new(session, min=0, max=1)
          on.exit(progress$close())
          progress$set(message = paste("Calculating geneset MOFA..."), value = 0.33)
          mofa$gset.mofa <- playbase::mofa.compute_geneset_mofa(
            mofa,
            kernel = input$kernel,
            factorname = "Module",
            #GMT = pgx$GMT,
            numfactors = 20
          )
        }
        
        ## update factors in selectInput
        factors <- colnames(mofa$W)
        dtypes  <- names(mofa$ww)
        sel.dtypes <- grep("^gset",dtypes,value=TRUE,invert=TRUE)
        contrasts <- colnames(mofa$contrasts)
        phenotypes <- colnames(mofa$samples)
        modules <- colnames(mofa$gset.mofa$W)
        updateSelectInput(session, "selected_factor", choices = factors,
          selected = factors[1])
        updateSelectInput(session, "selected_module", choices = modules,
          selected = modules[1])
        updateSelectInput(session, "show_types", choices = dtypes,
          selected = sel.dtypes)

        return( mofa )
      }
    )

    shiny::observeEvent(
      list(
        input$compute
      ), {

      shiny::req(pgx$X, input$kernel, input$numfactors)

      if(input$compute==0) return(NULL)

      if(!is.null(pgx$datatype) && pgx$datatype != "multi-omics") {
        shinyalert::shinyalert("Error", "This is not a multi-omics dataset.")
        return(NULL)
      }
      
      kernel <- input$kernel
      pgx.showSmallModal(paste("Calculating",kernel,"...<br>please wait"))
      progress <- shiny::Progress$new(session, min=0, max=1)
      on.exit(progress$close())


      numfactors <- as.integer(input$numfactors)
      numfactors <- min( numfactors, min(dim(pgx$X)) )        

      dbg("[MofaBoard] *** recalculating MOFA ***")
      progress$set(message = paste("Calculating",kernel,"..."), value = 0.33)
      mofa <- playbase::pgx.compute_mofa(
          pgx,
          kernel = kernel,
          numfactors = numfactors,
          add_gsets = input$add_gsets)

      if(!"gset.mofa" %in% names(mofa)) {
        progress$set(message = paste("Calculating geneset MOFA..."), value = 0.66)
        mofa$gset.mofa <- playbase::mofa.compute_geneset_mofa(
          mofa,
          kernel = kernel,
          factorname = "Module",
          #GMT = pgx$GMT,
          numfactors = 20
        )
      }
      
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
    enrichmentTable_selected <- reactive({

      tbl <- enrichmentTable
      shiny::req(tbl$data())
      
      has.selection <- length(tbl$rownames_selected())>0 
      search_key <- tbl$search()
      has.search <- length(search_key)>0 && search_key[1]!=""
      
      if(has.search && !has.selection) {
        sel <- tbl$rownames_all()
      } else if(has.selection) {
        sel <- tbl$rownames_selected()
      } else {
        sel <- head( tbl$rownames_all(), 20)
      }
      sel
    })


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
      "enrichmentplot",
      mofa = mofa,
      pgx = pgx,
      input_k = reactive(input$selected_factor),
      ntop = 15,
      select = enrichmentTable_selected,
      req.selection = TRUE,      
      watermark = WATERMARK
    )

    mofa_plot_pathbank_server(
      "pathway",
      pgx = pgx,
      sel_pathway = enrichmentTable_selected,
      sel_contrast = reactive(NULL),
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

    mofa_plot_factorcorheatmap_server(
      "factorcorheatmap",
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
    
    mofa_plot_pathwayheatmap_server(
      "pathwayheatmap",
      mofa = mofa,
      input_factor = reactive(NULL),            
      selected = enrichmentTable_selected,
      watermark = WATERMARK
    )

    mofa_plot_centrality_server(
      "centrality",
      mofa = mofa,
      input_factor = reactive(input$selected_factor),
      show_types = reactive(input$show_types),
      watermark = WATERMARK
    )

    mofa_plot_gsetmofa_traitCor_server(
      "gset_traitcor",
      mofa = mofa,
      watermark = WATERMARK
    )

    mofa_plot_gsetmofa_factorCor_server(
      "gset_factorcor",
      mofa = mofa,
      watermark = WATERMARK
    )
    

    ## ------------- Table Modules --------------------------

    mofa_table_genetable_server(
      "mofa_genetable",
      mofa = mofa,
      selected_factor = reactive(input$selected_factor),
      ## datatypes = reactive(input$show_types),
      annot = reactive(pgx$genes)
    )
    
    enrichmentTable <- mofa_table_enrichment_server(
      "mofa_enrichmenttable",
      gsea = reactive(mofa()$gsea),
      selected_factor = reactive(input$selected_factor)            
    )

    mofa_table_gsetmofa_server(
      "gsetmofa_table",
      mofa = mofa,
      selected_module = reactive(input$selected_module)
    )

    return(NULL)
  })
} ## end of Board
