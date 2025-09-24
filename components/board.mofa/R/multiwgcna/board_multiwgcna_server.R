##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

MultiWGCNA_Board <- function(id, pgx) {
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
        "Dendrograms" = list(disable = c("phenotype","module","lasagna_options")),
        "Module-Trait" = list(disable = c("phenotype","module","mwgcna_options",
          "lasagna_options")),
        "Module correlation" = list(disable = c("module", "mwgcna_options",
          "lasagna_options")),
        "WGCNA-Lasagna" = list(disable = c("module","mwgcna_options")),
        "Feature Table" = list(disable = c("layers","mwgcna_options",
          "lasagna_options"))
    )

    shiny::observeEvent( input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })

    ## ============================================================================
    ## ============================ REACTIVES =====================================
    ## ============================================================================

    shiny::observeEvent( pgx$mofa, {
      
      shiny::validate( shiny::need(!is.null(pgx$mofa), "missing MOFA slot"))
      
      datatypes <- setdiff( names(pgx$mofa$xx), c("gset"))
      datatypes2 <- unique(c(datatypes,"gset"))
      updateSelectInput(session, "layers", choices = datatypes2,
        selected = datatypes)
      
    }, ignoreNULL=FALSE)


    r_multiwgcna <- shiny::eventReactive( {
      list( input$updateplots, pgx$X ) 
    }, {
      shiny::validate( shiny::need( !is.null(pgx$mofa), "missing MOFA slot"))
      shiny::validate( shiny::need( pgx$datatype=="multi-omics", "not multi-omics data"))
      
      shiny::req(pgx$X)
      shiny::req(pgx$mofa)

      dataX = pgx$mofa$xx
      samples = pgx$mofa$samples
      names(dataX) <- substring(names(dataX),1,2)

      ## add geneset matrix
      dataX$gset <- pgx$gsetX

      if(input$power == "<auto>") {
        #power <- NULL
        power <- "iqr"
      } else {
        power <- as.numeric(input$power)
      }
      
      progress <- shiny::Progress$new(session, min=0, max=1)
      on.exit(progress$close())
      progress$set(message = paste("computing multi-omics WGCNA..."), value = 0.33)
      pgx.showSmallModal("computing multi-omics WGCNA...")
      
      wgcna <- playbase::wgcna.compute_multiomics(
        dataX = dataX,
        samples = samples,
        clustMethod = input$clust,
        cutMethod = "hybrid",
        deepsplit = as.integer(input$deepsplit),
        power = power,
        ngenes = as.integer(input$ngenes),
        minmodsize = 3,
        minKME = 0.3,
        compute.enrichment = TRUE,
        xtop = 100,
        annot = pgx$genes,
        GMT = pgx$GMT,  ##??
        gsetX = pgx$gsetX,
        progress = progress
      ) 

      shiny::removeModal()

      phenotypes <- colnames(wgcna[[1]]$datTraits)
      updateSelectInput(session, "phenotype", choices = phenotypes,
        selected = phenotypes[1])
      
      layers <- names(wgcna)
      sel.layers <- setdiff(layers, c("gset","gs","pheno","ph"))
      updateSelectInput(session, "layers", choices = layers,
        selected = sel.layers)


      all_modules <- lapply(wgcna, function(w) sort(names(w$me.genes)))
      module1 <- all_modules[[1]][1]
      updateSelectInput(session, "module", choices = all_modules,
        selected = module1)
      
      return(wgcna)
    }, ignoreNULL=FALSE)


    ## ==========================================================================
    ## ========================== BOARD FUNCTIONS ===============================
    ## ==========================================================================

    
    ## ==========================================================================
    ## =========================== MODULES ======================================
    ## ==========================================================================
        
    multiwgcna_plot_dendrograms_server(
      id = "multiwgcnaDendro",
      mwgcna = r_multiwgcna,
      r_layers = reactive(input$layers)
    )

    multiwgcna_plot_power_server(
      id = "multiwgcnaPower",
      mwgcna = r_multiwgcna,
      r_layers = reactive(input$layers)
    )

    multiwgcna_plot_moduletrait_server(
      "multiwgcnaTrait",
      mwgcna = r_multiwgcna,
      r_layers = reactive(input$layers)
    )
    
    multiwgcna_plot_modulecorr_server(
      "multiwgcnaCorr",
      mwgcna = r_multiwgcna,
      r_layers = reactive(input$layers),
      r_phenotype = reactive(input$phenotype)
    )

    multiwgcna_plot_lasagna_server(
      "multiwgcnaLasagna",
      mwgcna = r_multiwgcna,
      r_phenotype = reactive(input$phenotype),
      r_layers = reactive(input$layers)
    )

    multiwgcna_table_modulegenes_server(
      id = "multiwgcnaTable",
      mwgcna = r_multiwgcna,
      r_annot = reactive(pgx$genes),
      r_phenotype = reactive(input$phenotype),
      r_module = reactive(input$module)      
    )

    multiwgcna_table_enrichment_server(
      id = "multiwgcnaEnrichment",
      mwgcna = r_multiwgcna,      
      r_module = reactive(input$module)
    )

    multiwgcna_table_crossgenes_server(
      id = "multiwgcnaCrossgene",
      mwgcna = r_multiwgcna,
      r_annot = reactive(pgx$genes),
      r_module = reactive(input$module)      
    )

    
    return(NULL)
  })
} ## end of Board
