##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


WgcnaBoard <- function(id, pgx) {
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

    ## ================================================================================
    ## ========================== OBSERVE FUNCTIONS ===================================
    ## ================================================================================

    infotext <-
      '<center><iframe width="560" height="315" src="https://www.youtube.com/embed/BtMQ7Y0NoIA?si=WUBozFwNdZbwpT69" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe></center>'

    shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>WGCNA Analysis Board</strong>"),
        shiny::HTML(infotext),
        size = "xl",
        easyClose = TRUE
      ))
    })


    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "WGCNA" = list(disable = c("selected_module", "selected_trait")),
      "Eigengenes" = list(disable = c("selected_module", "selected_trait")),
      "Modules" = list(disable = c(NULL)),
      "Enrichment" = list(disable = c("selected_trait"))
    )

    shiny::observeEvent(input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })

    shiny::observeEvent(input$useLLM, {
      if (input$useLLM) {
        model <- getUserOption(session, "llm_model")
        dbg("[WgcnaBoard] input$useLLM => model = ", model)
        if (is.null(model) || model == "") {
          shinyalert::shinyalert(
            "ERROR",
            "No LLM server available. Please check your settings."
          )
          return(NULL)
        }
        shinyalert::shinyalert(
          "WARNING",
          "Using LLM might expose some of your data to external LLM servers."
        )
      }
    })

    ## ================================================================================
    ## ======================= PRECOMPUTE FUNCTION ====================================
    ## ================================================================================

    compute_wgcna <- function() {
      pgx.showSmallModal("Recalculating WGCNA with new parameters...")
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Calculating WGCNA...", value = 0)

      message("[WGCNA:compute_wgcna] >>> Calculating WGCNA...")

      # ai_model = opt$LLM_MODEL
      ai_model <- getUserOption(session, "llm_model")
      message("[WGCNA:compute_wgcna] ai_model = ", ai_model)

      out <- playbase::pgx.wgcna(
        pgx = pgx,
        ngenes = as.integer(input$ngenes),
        gset.filter = NULL,
        minmodsize = as.integer(input$minmodsize),
        power = as.numeric(input$power),
        numericlabels = FALSE,
        summary = TRUE,
        ai_model = ifelse(input$useLLM, ai_model, ""),
        progress = progress
      )
      shiny::removeModal()
      out
    }

    wgcna <- shiny::reactiveVal({
      require(WGCNA)
      all.req <- all(c("stats") %in% names(pgx$wgcna)) && any(c("TOM", "svTOM", "wTOM") %in% names(pgx$wgcna))
      # Use pre-computed results only if they exist, conditions are
      # met, AND we're not forcing recomputation
      if ("wgcna" %in% names(pgx) && all.req) {
        message("[wgcna] >>> using pre-computed WGCNA results...")
        out <- pgx$wgcna
        ## old style had these settings
        if (is.null(pgx$wgcna$networktype)) out$networktype <- "unsigned"
        if (is.null(pgx$wgcna$tomtype)) out$tomtype <- "signed"
        if (is.null(pgx$wgcna$power)) out$power <- 6
      } else {
        message("[wgcna] >>> COMPUTE1")
        out <- compute_wgcna()
      }
      out
    })

    shiny::observeEvent(input$compute,
      {
        message("[wgcna] >>> COMPUTE2")
        wgcna(compute_wgcna())
      },
      ignoreInit = TRUE
    )

    shiny::observeEvent(wgcna(), {
      ## update Inputs
      me <- sort(names(wgcna()$me.genes))
      shiny::updateSelectInput(session, "selected_module",
        choices = me,
        sel = me[1]
      )
      tt <- sort(colnames(wgcna()$datTraits))
      shiny::updateSelectInput(session, "selected_trait",
        choices = tt,
        selected = tt[1]
      )
    })


    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    # Gene dendrogram and gene modules
    wgcna_plot_gdendogram_server(
      "geneDendro",
      wgcna.compute = wgcna,
      watermark = WATERMARK
    )

    # Scale independence and mean connectivity
    wgcna_plot_s_independence_server(
      "topologyPlots",
      wgcna.compute = wgcna,
      watermark = WATERMARK
    )

    # TOM heatmap
    wgcna_plot_TOMheatmap_server(
      "TOMplot",
      wgcna.compute = wgcna,
      watermark = WATERMARK
    )

    # Gene clustering
    wgcna_plot_gclustering_server(
      "umap",
      wgcna = wgcna,
      pgx = pgx,
      watermark = WATERMARK
    )

    # Module graph
    wgcna_plot_module_graph_server(
      "moduleGraph",
      wgcna = wgcna,
      watermark = WATERMARK
    )

    # Module-Trait relationships
    wgcna_plot_MTrelationships_server(
      "moduleTrait",
      wgcna.compute = wgcna,
      watermark = WATERMARK
    )

    # Correlation network
    wgcna_plot_correlation_network_server(
      "corGraph",
      wgcna = wgcna,
      pgx = pgx,
      selected_module = shiny::reactive(input$selected_module),
      watermark = WATERMARK
    )

    # Module genes
    wgcna_table_genes_server(
      "geneTable",
      wgcna = wgcna,
      pgx = pgx,
      selected_module = shiny::reactive(input$selected_module),
      selected_trait = shiny::reactive(input$selected_trait)
    )

    # Eigengene clustering
    wgcna_plot_eigengene_clustering_server(
      "eigenClustering",
      wgcna.compute = wgcna,
      watermark = WATERMARK
    )

    wgcna_plot_eigengene_heatmap_server(
      "eigenHeatmap",
      wgcna = wgcna,
      watermark = WATERMARK
    )

    # Module membership (eigengene correlation)
    wgcna_plot_module_barplot_server(
      "moduleSize",
      wgcna = wgcna,
      watermark = WATERMARK
    )

    # Module membership (eigengene correlation)
    wgcna_plot_module_membership_server(
      "modulemembership",
      wgcna = wgcna,
      pgx = pgx,
      selected_module = shiny::reactive(input$selected_module),
      watermark = WATERMARK
    )

    # Module membership (eigengene correlation)
    wgcna_plot_module_heatmap_server(
      "moduleheatmap",
      wgcna = wgcna,
      pgx = pgx,
      selected_module = shiny::reactive(input$selected_module),
      selected_trait = shiny::reactive(input$selected_trait),
      watermark = WATERMARK
    )

    # Membership-trait heatmap
    ## wgcna_plot_heatmap_membership_server(
    ##   "intraHeatmap",
    ##   wgcna.compute = wgcna,
    ##   watermark = WATERMARK
    ## )

    # Membership vs. trait correlation
    wgcna_plot_membership_v_trait_server(
      "memberTrait",
      wgcna = wgcna,
      selected_module = shiny::reactive(input$selected_module),
      selected_trait = shiny::reactive(input$selected_trait),
      watermark = WATERMARK
    )

    wgcna_plot_module_significance_server(
      "moduleSignificance",
      wgcna.compute = wgcna,
      selected_module = shiny::reactive(input$selected_module),
      watermark = WATERMARK
    )

    ## wgcna_plot_MMvsGS_server(
    ##   "geneSignificance",
    ##   wgcna.compute = wgcna,
    ##   watermark = WATERMARK
    ## )

    wgcna_plot_sampledendrogram_server(
      "sampleDendrogram",
      wgcna = wgcna,
      what = "traits",
      watermark = WATERMARK
    )

    wgcna_plot_sampledendrogram_server(
      "sampleDendrogram2",
      wgcna = wgcna,
      what = "me",
      watermark = WATERMARK
    )

    wgcna_plot_geneset_heatmap_server(
      "genesetHeatmap",
      pgx = pgx,
      selected_module = shiny::reactive(input$selected_module),
      enrichTable = enrichTableModule,
      watermark = FALSE
    )

    wgcna_plot_gene_heatmap_server(
      "geneHeatmap",
      pgx = pgx,
      wgcna = wgcna,
      selected_module = shiny::reactive(input$selected_module),
      enrichTable = enrichTableModule,
      watermark = FALSE
    )

    # Enrichment plot
    wgcna_plot_enrichment_server(
      "enrichPlot",
      enrichTable = enrichTableModule,
      pgx = pgx,
      watermark = WATERMARK
    )

    # Module enrichment
    enrichTableModule <- wgcna_table_enrichment_server(
      "enrichTable",
      wgcna = wgcna,
      selected_module = shiny::reactive(input$selected_module)
    )

    # Enrichment plot
    wgcna_html_module_summary_server(
      "moduleSummary",
      wgcna = wgcna,
      multi = FALSE,
      r_module = shiny::reactive(input$selected_module),
      watermark = WATERMARK
    )

    return(NULL)
  })
} ## end of Board
