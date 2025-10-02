##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

PreservationWGCNA_Board <- function(id, pgx) {
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
        title = shiny::HTML("<strong>Preservation WGCNA Board</strong>"),
        shiny::HTML(infotext),
        size = "xl",
        easyClose = TRUE
      ))
    })

    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Dendrograms" = list(disable = c("module","trait")),
      "Sample Clustering" = list(disable = c("module","trait")),
      "Module-Trait" = list(disable = c("module")),
      "Module Table" = list(disable = c(), enable = c("module","trait"))
    )

    shiny::observeEvent( input$tabs, {
      ##dbg("[PreservationWGCNA_Board] input$tabs = ", input$tabs)
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })

    ## ============================================================================
    ## ============================ REACTIVES =====================================
    ## ============================================================================
   

    shiny::observeEvent( list(pgx$X, pgx$samples), {
      splitpheno <- colnames(pgx$samples)
      shiny::updateSelectInput(session, "splitpheno", choices = splitpheno,
        selected = splitpheno[1])
    })
    
    
    r_wgcna <- shiny::eventReactive( {
      list( input$compute, pgx$X ) 
    }, {

      shiny::req(pgx$X)
      shiny::req(input$splitpheno)

      dbg("[PreservationWGCNA_Board] input.compute = ", input$compute)
      dbg("[PreservationWGCNA_Board] input.tabs = ", input$tabs)
      shiny::req(input$compute) ## refute first call     
      
      pheno="activated"
      phenoData <- pgx$samples

      pheno <- input$splitpheno
      if( is.null(pheno) || pheno == '') {
        pheno <- colnames(phenoData)[1]
      }
      shiny::req(pheno %in% colnames(phenoData))

      group <- phenoData[,pheno]
      if(is.numeric(group) && length(unique(group)) > 3) {
        group <- c("LO", "HI")[1 + (group >= median(group,na.rm=TRUE))]
      }
      group <- base::abbreviate(toupper(group),2L)
      exprList<- tapply(1:ncol(pgx$X), group, function(ii) pgx$X[,ii,drop=FALSE])
            
      progress <- shiny::Progress$new(session, min=0, max=1)
      on.exit(progress$close())
      progress$set(message = paste("computing preservation WGCNA..."), value = 0.33)
      pgx.showSmallModal("computing preservation WGCNA...")
      
      power <- input$power
      if(power == "<auto>") {
        power <- NULL
      } else {
        power <- as.numeric(power)
      }
      
      ## This runs preservation WGCNA on an expression list
      #ngenes=2000;minModuleSize=20;deepSplit=2
      ngenes = as.integer(input$ngenes)
      minModuleSize = as.integer(input$minmodsize)
      deepSplit = as.integer(input$deepsplit)
      
      dbg("[PreservationWGCNA_Board] power = ", power)
      dbg("[PreservationWGCNA_Board] ngenes = ", ngenes)
      dbg("[PreservationWGCNA_Board] minModuleSize = ", minModuleSize)
      dbg("[PreservationWGCNA_Board] deepSplit = ", deepSplit)      

      res <- playbase::wgcna.runPreservationWGCNA(
        exprList,
        phenoData,
        annot = pgx$genes,
        reference = 1,
        add.merged = FALSE,
        compute.stats = TRUE,
        compute.enrichment = TRUE
      ) 
      
      shiny::removeModal()
     
      all_modules <- rownames(res$modTraits[[1]])
      module1 <- all_modules[1]
      updateSelectInput(session, "module", choices = sort(all_modules),
        selected = module1)

      ##all_traits <- colnames(res$modTraits[[1]])
      all_traits <- colnames(res$zlist[[1]])
      trait1 <- all_traits[1]
      updateSelectInput(session, "trait", choices = sort(all_traits),
        selected = trait1)
      
      return(res)
    }, ignoreNULL=TRUE )


    ## ==========================================================================
    ## ========================== BOARD FUNCTIONS ===============================
    ## ==========================================================================

    
    ## ==========================================================================
    ## =========================== MODULES ======================================
    ## ==========================================================================
        
    preservationWGCNA_plot_dendrograms_server(
      id = "preservationWGCNADendro",
      rwgcna = r_wgcna
    )

    preservationWGCNA_plot_summaries_server(
      id = "preservationWGCNASummaries",
      rwgcna = r_wgcna
    )

    preservationWGCNA_plot_overlap_server(
      "preservationWGCNAOverlap",
      rwgcna = r_wgcna
    )

    preservationWGCNA_plot_eigenNetwork_server(
      "preservationWGCNAEigenNetwork",
      rwgcna = r_wgcna
    )
    
    preservationWGCNA_plot_moduletrait_server(
      "preservationWGCNAModuleTrait",
      rwgcna = r_wgcna
    )

    preservationWGCNA_plot_modulenetwork_server(
      id = "preservationWGCNAModuleNetwork",
      rwgcna = r_wgcna,
      rmodule = reactive(input$module)
    )

    preservationWGCNA_plot_traitsignificance_server(
      id = "preservationWGCNATraitSignificance",
      rwgcna = r_wgcna,
      rtrait = reactive(input$trait)
    )

    preservationWGCNA_table_modulegenes_server(
      id = "preservationWGCNATable",
      rwgcna = r_wgcna,
      rannot = reactive(pgx$genes),
      rtrait = reactive(input$trait),
      rmodule = reactive(input$module)      
    )

    preservationWGCNA_table_enrichment_server(
      id = "preservationWGCNAEnrichment",
      rwgcna = r_wgcna,      
      rmodule = reactive(input$module)
    )
    
    return(NULL)
  })
} ## end of Board
