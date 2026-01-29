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
      "Dendrograms" = list(disable = c("phenotype", "module", "condition", "lasagna_options", "report_options")),
      "Module-Trait" = list(disable = c("phenotype", "module", "condition", "wgcna_options",
        "lasagna_options", "report_options")),
      "Module correlation" = list(disable = c("phenotype", "module", "wgcna_options", "lasagna_options",
        "report_options")),
      "WGCNA-Lasagna" = list(disable = c("module", "condition", "wgcna_options", "report_options")),
      "Feature Table" = list(disable = c("layers", "condition", "wgcna_options", "lasagna_options",
        "report_options")),
      "AI Report" = list(disable = c("phenotype","module", "condition", "layers", "lasagna_options",
        "wgcna_options"))
    )

    shiny::observeEvent(input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })

    ## ============================================================================
    ## ============================ REACTIVES =====================================
    ## ============================================================================
    
    r_multiwgcna <- shiny::eventReactive(
      {
        list(input$compute, pgx$X)
      },
      {
        shiny::req(pgx$X)
        shiny::validate(shiny::need(
          pgx$datatype %in% c("multi-omics","multiomics"),
          "ERROR: not multi-omics data"
        ))

        if (input$power == "<auto>") {
          power <- "iqr"
        } else {
          power <- as.numeric(input$power)
        }
        
        ## setup progress bars
        progress <- shiny::Progress$new(session, min = 0, max = 1)
        on.exit(progress$close())
        progress$set(message = paste("computing multi-omics WGCNA..."), value = 0.33)
        pgx.showSmallModal("computing multi-omics WGCNA...")
        
        dataX <- playbase::mofa.split_data(pgx$X)
        samples <- pgx$samples
        contrasts <- pgx$contrasts
        
        obj.wgcna <- playbase::wgcna.compute_multiomics(
          dataX = dataX,
          samples = samples,
          contrasts = contrasts,
          do.consensus = input$consensus,
          add.pheno = (ncol(samples) > 10),
          add.gsets = input$addgsets,
          cutMethod = "hybrid",
          deepsplit = as.integer(input$deepsplit),
          power = power,
          ngenes = as.integer(input$ngenes),
          minmodsize = as.integer(input$minmodsize),
          minKME = 0.3,
          compute.enrichment = TRUE,
          gset.xtop = 100,
          gset.ntop = 1000,        
          gset.methods = c("gsetcor","xcor","fisher"),        
          annot = pgx$genes,
          GMT = pgx$GMT,  
          report = TRUE,
          ai_model = NULL,
          experiment = pgx$description,
          progress = progress
        )         
        shiny::removeModal()

        wgcna <- obj.wgcna$layers
        phenotypes <- colnames(wgcna[[1]]$datTraits)

        updateSelectInput(session, "phenotype",
          choices = phenotypes,
          selected = phenotypes[1]
        )

        updateSelectInput(session, "condition",
          choices = c("None", phenotypes),
          selected = "None"
        )

        layers <- names(wgcna)
        sel.layers <- setdiff(layers, c("gset", "gs", "pheno", "ph"))
        updateSelectInput(session, "layers",
          choices = layers,
          selected = sel.layers
        )

        all_modules <- lapply(wgcna, function(w) sort(names(w$me.genes)))
        module1 <- all_modules[[1]][1]
        updateSelectInput(session, "module",
          choices = all_modules,
          selected = module1
        )

        return(obj.wgcna)
      },
      ignoreNULL = FALSE
    )


    ## ==========================================================================
    ## ========================== BOARD FUNCTIONS ===============================
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
      r_condition = reactive(input$condition)
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

    # Enrichment plot
    wgcna_html_module_summary_server(
      "multiwgcnaSummary",
      wgcna = r_multiwgcna,
      multi = TRUE,
      r_annot = reactive(pgx$genes),      
      r_module = shiny::reactive(input$module),      
      watermark = WATERMARK
    )

    wgcna_html_report_server(
      id = "multiwgcnaReport",
      wgcna = r_multiwgcna,
      multi = TRUE,
      r_annot = reactive(pgx$genes),      
      watermark = WATERMARK
    )

    return(NULL)
  })
} ## end of Board
