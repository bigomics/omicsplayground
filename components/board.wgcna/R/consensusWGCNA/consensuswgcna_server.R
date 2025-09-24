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
      "Dendrograms" = list(disable = c("phenotype","module")),
      "Sample Clustering" = list(disable = c("phenotype","module")),
      "Module-Trait" = list(disable = c("phenotype","module")),
      "Preservation" = list(disable = c("phenotype","module")),      
      "Feature Table" = list(disable = c())
    )

    shiny::observeEvent( input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })

    ## ============================================================================
    ## ============================ REACTIVES =====================================
    ## ============================================================================
   
    shiny::observe({
      cons <- r_consensusWGCNA()
      dbg("[ConsensusWGCNA_Board] isnull.cons=", is.null(cons))
      dbg("[ConsensusWGCNA_Board] len.cons=", length(cons))      
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

      phenotypes <- colnames(pgx$samples)
      shiny::updateSelectInput(session, "splitpheno", choices = phenotypes,
        selected = phenotypes[1])
      
    })
    
    
    r_consensusWGCNA <- shiny::eventReactive( {
      list( input$compute, pgx$X ) 
    }, {

      shiny::req(pgx$X)
      shiny::req(input$splitpheno)      

      dbg("[r_consensusWGCNA] dimX = ", dim(pgx$X))
      dbg("[r_consensusWGCNA] input.splitpheno = ", input$splitpheno)
      dbg("[r_consensusWGCNA] dimSamples = ", dim(pgx$samples))
      dbg("[r_consensusWGCNA] colnamesSamples = ", colnames(pgx$samples))
      
      xx <- NULL
      if( pgx$datatype == "multi-omics" ) {        
        if("mofa" %in% names(pgx)) {
          xx <- pgx$mofa$xx
        } else {
          xx <- playbase::mofa.split_data(pgx$X)
        }        
        has.gxpx <- all(c("gx","px") %in% names(xx))
        has.gxpx
        shiny::validate(shiny::need(has.gxpx, "Your mulit-omics dataset is incompatible for Consensus WGCNA: You must have both transcriptomics (gx) and proteomics (px)"))

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
      
      ## Rename all tables to symbol
      xx <- lapply(xx, playbase::rename_by2, annot_table=pgx$genes)
      gg <- Reduce(intersect, lapply(xx, rownames))
      shiny::validate(shiny::need(length(gg)>0, "Your dataset is incompatible for consensus WGCNA: No overlapping features."))
      xx <- lapply(xx, function(x) x[gg,])
      ##names(xx) <- substring(names(xx),1,2)
      names(xx) <- base::abbreviate(names(xx),2)
      
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
      cons <- playbase::wgcna.runConsensusWGCNA(
        exprList = xx,
        phenoData = pgx$samples,
        # GMT = pgx$GMT,
        annot = pgx$genes,
        power = power,
        ngenes = as.integer(input$ngenes),
        minModuleSize = as.integer(input$minmodsize),
        maxBlockSize = 9999,
        minKME = 0.3,
        mergeCutHeight = 0.15,
        deepSplit = as.integer(input$deepsplit),
        calcMethod = "fast",
        drop.ref = FALSE,
        addCombined = FALSE,
        gsea.mingenes = 10
      ) 
      
      shiny::removeModal()
      
      phenotypes <- colnames(cons$datTraits)
      updateSelectInput(session, "phenotype", choices = sort(phenotypes),
        selected = phenotypes[1])
      
      all_modules <- rownames(cons$modTraits)
      module1 <- all_modules[[1]][1]
      updateSelectInput(session, "module", choices = sort(all_modules),
        selected = module1)
      
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
      r_phenotype = reactive(input$phenotype),
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
    
    return(NULL)
  })
} ## end of Board
