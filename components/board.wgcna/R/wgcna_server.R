##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


WgcnaBoard <- function(id, pgx, board_observers) {
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
    ## ======================= PRECOMPUTE FUNCTION ====================================
    ## ================================================================================

    ## wgcna.compute <- shiny::reactive({
    wgcna.compute <- shiny::eventReactive(
      {
        list(input$compute)
      },
      {
        require(WGCNA)

        if (input$compute == 0 && "wgcna" %in% names(pgx)) {
          message("[wgcna.compute] >>> using pre-computed WGCNA results...")
          me <- names(pgx$wgcna$me.genes)
          shiny::updateSelectInput(session, "selected_module", choices = me, selected = "ME1")
          return(pgx$wgcna)
        }

        pgx.showSmallModal("Calculating WGCNA using new parameters...<br>please wait")
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Calculating WGCNA...", value = 0)
        message("[wgcna.compute] >>> Calculating WGCNA...")

        WGCNA::enableWGCNAThreads()

        out <- playbase::pgx.wgcna(
          pgx = pgx,
          ngenes = as.integer(input$ngenes),
          minmodsize = as.integer(input$minmodsize),
          power = as.numeric(input$power),
          cutheight = as.numeric(input$cutheight),
          deepsplit = as.integer(input$deepsplit)
        )

        me <- names(out$me.genes)
        shiny::updateSelectInput(session, "selected_module", choices = me, sel = "ME1")

        beepr::beep(2)
        shiny::removeModal()

        out
      }
    )


    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================
    my_observers <- list()

    infotext <-
      '<center><iframe width="1120" height="630" src="https://www.youtube.com/embed/rRIRMW_RRS4"
        title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write;
        encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></center>'

    my_observers[[1]] <- shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>WGCNA Analysis Board</strong>"),
        shiny::HTML(infotext),
        size = "xl",
        easyClose = TRUE
      ))
    })

    ## add to list global of observers. suspend by default.
    my_observers <- my_observers[!sapply(my_observers,is.null)]
    lapply( my_observers, function(b) b$suspend() )
    if(!is.null(board_observers)) board_observers[[id]] <- my_observers
    
    ## ================================================================================
    ## ======================= PLOTTING FUNCTIONS =====================================
    ## ================================================================================

    ## ----------------------------------------
    ## --------- enrichment table -------------
    ## ----------------------------------------

    enrich_table <- shiny::reactive({
      out <- wgcna.compute()
      df <- out$gse
      k <- input$selected_module
      message("[enrichTable.RENDER] 1 : k = ", k)
      if (length(k) == 0 || k == "") k <- "<all>"
      message("[enrichTable.RENDER] 2 : k = ", k)
      if (k %in% df$module) {
        df <- df[df$module == k, , drop = FALSE]
      }
      df$score <- df$odd.ratio * -log10(df$p.value)
      df <- df[, c(
        "module", "geneset", "score", "p.value", "q.value",
        "odd.ratio", "overlap", "genes"
      )]
      df <- df[order(-df$score), ]
      df
    })

    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    # Gene dendrogram and gene modules

    wgcna_plot_gdendogram_server(
      "geneDendro",
      wgcna.compute = wgcna.compute,
      labels2rainbow = playbase::labels2rainbow,
      watermark = WATERMARK
    )

    # Scale independence and mean connectivity

    wgcna_plot_s_independence_server(
      "topologyPlots",
      wgcna.compute = wgcna.compute,
      watermark = WATERMARK
    )

    # TOM heatmap

    wgcna_plot_TOMheatmap_server(
      "TOMplot",
      wgcna.compute = wgcna.compute,
      labels2rainbow = playbase::labels2rainbow,
      power = shiny::reactive(input$power),
      watermark = WATERMARK
    )

    # Gene clustering

    wgcna_plot_gclustering_server(
      "umap",
      wgcna.compute = wgcna.compute,
      watermark = WATERMARK
    )

    # Module graph

    wgcna_plot_module_graph_server(
      "moduleGraph",
      wgcna.compute = wgcna.compute,
      labels2rainbow = playbase::labels2rainbow,
      watermark = WATERMARK
    )

    # Module-Trait relationships

    wgcna_plot_MTrelationships_server(
      "moduleTrait",
      wgcna.compute = wgcna.compute,
      labels2rainbow = playbase::labels2rainbow,
      watermark = WATERMARK
    )

    # Correlation network

    wgcna_plot_correlation_network_server(
      "corGraph",
      wgcna.compute = wgcna.compute,
      selected_module = shiny::reactive(input$selected_module),
      watermark = WATERMARK
    )

    # Enrichment plot

    wgcna_plot_enrichment_server(
      "enrichPlot",
      enrich_table = enrich_table,
      enrichTable_module = enrichTable_module,
      watermark = WATERMARK
    )

    # Module genes

    wgcna_table_genes_server(
      "geneTable",
      wgcna.compute = wgcna.compute,
      selected_module = shiny::reactive(input$selected_module)
    )

    # Module enrichment

    enrichTable_module <- wgcna_table_enrichment_server(
      "enrichTable",
      enrich_table
    )

    # Eigengene clustering

    wgcna_plot_eigengene_clustering_server(
      "eigenClustering",
      wgcna.compute = wgcna.compute,
      watermark = WATERMARK
    )

    # Module membership (eigengene correlation)

    wgcna_plot_module_membership_server(
      "eigenCorrelation",
      wgcna.compute = wgcna.compute,
      watermark = WATERMARK
    )

    # Membership-trait heatmap

    wgcna_plot_heatmap_membership_server(
      "intraHeatmap",
      wgcna.compute = wgcna.compute,
      watermark = WATERMARK
    )

    # Membership vs. trait correlation

    wgcna_plot_membership_v_trait_server(
      "intraScatter",
      wgcna.compute = wgcna.compute,
      selected_module = shiny::reactive(input$selected_module),
      watermark = WATERMARK
    )

    return(NULL)
  })
} ## end of Board
