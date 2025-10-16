##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ConsensusWGCNA_Board <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    fullH <- 700 ## full height of page
    rowH1 <- 250 ## row 1 height
    rowH2 <- 440 ## row 2 height

    infotext <- tspan("<b>Multi-Omics WGCNA</b> ...
<p>References:<br>
<ol>
<li>Langfelder, P. and Horvath, S., 2008. WGCNA: an R package for weighted correlation network analysis. BMC bioinformatics, 9(1), p.559.
<li>Zhang, B. and Horvath, S., 2005. A general framework for weighted gene co-expression network analysis. Statistical applications in genetics and molecular biology, 4(1).
</ol>
", js = FALSE)


    ## ============================================================================
    ## ============================ OBSERVERS =====================================
    ## ============================================================================

    infotext <-
      '<center><iframe width="1120" height="630" src="https://www.youtube.com/embed/rRIRMW_RRS4"
        title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write;
        encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></center>'

    shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Multi-Omics WGCNA Board</strong>"),
        shiny::HTML(infotext),
        size = "xl",
        easyClose = TRUE
      ))
    })

    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Dendrograms" = list(disable = c("trait","module")),
      "Sample Clustering" = list(disable = c("trait","module")),
      "Module-Trait" = list(disable = c("trait","module")),
      "Preservation" = list(disable = c("trait","module")),      
      "Feature Table" = list(disable = c())
    )

    shiny::observeEvent( input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })

    shiny::observeEvent( input$useLLM, {
      if(input$useLLM) {
        shinyalert::shinyalert("",
          "Warning. Using LLM might expose some of your data to external LLM servers.")
      }
    })

    ## ============================================================================
    ## ============================ REACTIVES =====================================
    ## ============================================================================
   
    shiny::observe({
      cons <- r_consensusWGCNA()
      shiny::validate(shiny::need(!is.null(cons), "Please compute"))
    })

    shiny::observeEvent( pgx$X, {
      if(pgx$datatype == "multi-omics") {
        shinyjs::hide("splitpheno")
        shinyjs::show("splitdata")        
      } else {
        shinyjs::hide("splitdata")                
        shinyjs::show("splitpheno")
      }
    })

    shiny::observeEvent( list(pgx$X, pgx$samples), {

      dtypes <- names(playbase::mofa.split_data(pgx$X))
      sel.dtypes <- intersect(c("gx","px","me"), dtypes)
      shiny::updateSelectizeInput(session, "splitdata", choices = dtypes,
        selected = sel.dtypes)

      splitpheno <- colnames(pgx$samples)
      shiny::updateSelectInput(session, "splitpheno", choices = splitpheno,
        selected = splitpheno[1])
      
    })
    
    
    r_consensusWGCNA <- shiny::eventReactive( {
      list( input$compute, pgx$X ) 
    }, {

      shiny::req(pgx$X)
      shiny::req(input$splitpheno)      
      
      xx <- NULL
      if( pgx$datatype == "multi-omics" ) {        
        xx <- playbase::mofa.split_data(pgx$X)
        has.gxpx <- all(c("gx","px") %in% names(xx))
        has.gxpx
        shiny::validate(shiny::need(has.gxpx, "Your mulit-omics dataset is incompatible for Consensus WGCNA: You must have both transcriptomics (gx) and proteomics (px)"))

        ## Rename all tables to symbol
        xx <- xx[names(xx) %in% c("gx","px")]      
        xx <- lapply(xx, function(x) playbase::rename_by2(x, annot_table=pgx$genes))
        gg <- Reduce(intersect, lapply(xx, rownames))
        shiny::validate(shiny::need(length(gg)>0, "Your dataset is incompatible for consensus WGCNA: No overlapping features."))
        xx <- lapply(xx, function(x) x[gg,])
        
      } else if(!is.null(pgx$samples)) {

        pheno <- input$splitpheno
        if( is.null(pheno) || pheno == '') {
          pheno <- colnames(pgx$samples)[1]
        }
        shiny::req(pheno %in% colnames(pgx$samples))
        group <- pgx$samples[,pheno]
        if(is.numeric(group) && length(unique(group)) > 3) {
          group <- c("LO", "HI")[1 + (group >= median(group,na.rm=TRUE))]
        }
        group <- base::abbreviate(toupper(group),2L)
        xx <- tapply(1:ncol(pgx$X), group, function(ii) pgx$X[,ii,drop=FALSE])
      } else {
        ## should not come here???
        shiny::validate(shiny::need(has.gxpx, "Your dataset is incompatible for consensus WGCNA."))
      }
      
      progress <- shiny::Progress$new(session, min=0, max=1)
      on.exit(progress$close())
      progress$set(message = paste("computing consensus WGCNA..."), value = 0.33)
      pgx.showSmallModal("computing consensus WGCNA...")
      
      power <- input$power
      if(power == "<auto>") {
        power <- NULL
      } else {
        power <- as.numeric(power)
      }
      
      ## This runs consensus WGCNA on an expression list
      #ngenes=2000;minModuleSize=20;deepSplit=2
      ngenes = as.integer(input$ngenes)
      minModuleSize = as.integer(input$minmodsize)
      deepSplit = as.integer(input$deepsplit)
      
      cons <- playbase::wgcna.runConsensusWGCNA(
        exprList = xx,
        phenoData = pgx$samples,
        # GMT = pgx$GMT,
        annot = pgx$genes,
        power = power,
        ngenes = ngenes,
        minModuleSize = minModuleSize,
        maxBlockSize = 9999,
        minKME = 0.3,
        mergeCutHeight = 0.15,
        deepSplit = deepSplit,
        calcMethod = "fast",
        drop.ref = FALSE,
        addCombined = FALSE,
        compute.stats = TRUE,
        compute.enrichment = TRUE,
        ai_summary = input$useLLM,
        ai_model = opt$LLM_MODEL,
        ai_experiment = pgx$description,
        gsea.mingenes = 5,
        gsea.ntop = 1000,
        verbose = 1
      ) 
      
      shiny::removeModal()
            
      all_modules <- rownames(cons$modTraits)
      module1 <- all_modules[[1]][1]
      updateSelectInput(session, "module", choices = sort(all_modules),
        selected = module1)

      #traits <- colnames(cons$datTraits)
      traits <- colnames(cons$stats[[1]][['moduleTraitCor']])
      updateSelectInput(session, "trait", choices = sort(traits),
        selected = traits[1])
      
      return(cons)
    }, ignoreNULL=FALSE)


    ## ==========================================================================
    ## ========================== BOARD FUNCTIONS ===============================
    ## ==========================================================================

    
    ## ==========================================================================
    ## =========================== MODULES ======================================
    ## ==========================================================================
        
    consensusWGCNA_plot_dendrograms_server(
      id = "consensusWGCNADendro",
      mwgcna = r_consensusWGCNA
    )

    consensusWGCNA_plot_power_server(
      id = "consensusWGCNAPower",
      mwgcna = r_consensusWGCNA
    )

    consensusWGCNA_plot_moduletrait_server(
      "consensusWGCNATrait",
      mwgcna = r_consensusWGCNA
    )

    consensusWGCNA_plot_sampletree_server(
      "consensusWGCNASampleTree",
      mwgcna = r_consensusWGCNA
    )

    consensusWGCNA_table_modulegenes_server(
      id = "consensusWGCNATable",
      mwgcna = r_consensusWGCNA,
      r_annot = reactive(pgx$genes),
      r_trait = reactive(input$trait),
      r_module = reactive(input$module)      
    )

    consensusWGCNA_table_enrichment_server(
      id = "consensusWGCNAEnrichment",
      mwgcna = r_consensusWGCNA,      
      r_module = reactive(input$module)
    )

    consensusWGCNA_plot_preservation_server(
      id = "consensusWGCNAPreservation",
      mwgcna = r_consensusWGCNA
    )
    
    # Enrichment plot
    wgcna_html_module_summary_server(
      "consensusWGCNAmoduleSummary",
      wgcna = r_consensusWGCNA,
      multi = TRUE,
      pgx = pgx,
      r_module = shiny::reactive(input$module),      
      watermark = WATERMARK
    )

    return(NULL)
  })
} ## end of Board
