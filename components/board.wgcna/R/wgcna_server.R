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

    infotext <- "<b>Weighted gene co-expression network analysis (WGCNA)</b> is a systems biology method for describing the correlation patterns among genes across microarray samples. Weighted correlation network analysis can be used for finding clusters (modules) of highly correlated genes, for summarizing such clusters using the module eigengene or an intramodular hub gene, for relating modules to one another and to external sample traits (using eigengene network methodology), and for calculating module membership measures. Correlation networks facilitate network based gene screening methods that can be used to identify candidate biomarkers or therapeutic targets.

<p>References:<br>
<ol>
<li>Langfelder, P. and Horvath, S., 2008. WGCNA: an R package for weighted correlation network analysis. BMC bioinformatics, 9(1), p.559.
<li>Zhang, B. and Horvath, S., 2005. A general framework for weighted gene co-expression network analysis. Statistical applications in genetics and molecular biology, 4(1).
</ol>
"

    intra_caption <-
      "<b>WGCNA intramodular analysis.</b> We quantify associations of individual genes with our trait of interest (weight) by defining Gene Significance GS as (the absolute value of) the correlation between the gene and the trait. For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules."

    output$intra_UI <- shiny::renderUI({
      shiny::fillCol(
        flex = c(NA, 0.04, 2, 1),
        height = fullH,
        shiny::div(shiny::HTML(intra_caption), class = "caption"),
        shiny::br(),
        shiny::fillRow(
          flex = c(1, 0.06, 2.5),
          ## plotWidget(ns('eigenHeatmap')),
          plotWidget(ns("intraHeatmap")),
          shiny::br(),
          plotWidget(ns("intraScatter"))
        )
      )
    })

    ## ================================================================================
    ## ======================= PRECOMPUTE FUNCTION ====================================
    ## ================================================================================

    ## wgcna.compute <- shiny::reactive({
    wgcna.compute <- shiny::eventReactive(
      {
        input$compute
        1
      },
      {
        require(WGCNA)

        if ("wgcna" %in% names(pgx)) {
          message("[wgcna.compute] >>> using pre-computed WGCNA results...")
          return(pgx$wgcna)
        }

        pgx.showSmallModal("Calculating WGCNA...<br>please wait")
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Calculating WGCNA...", value = 0)
        message("[wgcna.compute] >>> Calculating WGCNA...")
        if (0) {
          shinyalert::shinyalert(
            title = "",
            text = "No WGCNA data found in PGX object. Computing now.. "
          )
        }

        WGCNA::enableWGCNAThreads()

        if (0) {
          liv <- read.csv("~/Downloads/LiverFemale3600.csv")
          X <- as.matrix(liv[, 9:ncol(liv)])
          rownames(X) <- liv$gene_symbol
          X <- X[order(-apply(X, 1, sd)), ]
          X <- X[!duplicated(rownames(X)), ]
          dim(X)
        }

        out <- playbase::pgx.wgcna(
          pgx = pgx,
          ngenes = input$ngenes,
          minmodsize = as.integer(input$minmodsize),
          power = as.numeric(input$power),
          cutheight = as.numeric(input$cutheight),
          deepsplit = as.integer(input$deepsplit)
        )
        
        shiny::updateSelectInput(session, "selected_module", choices = names(out$me.genes), sel = "ME1")

        beepr::beep(2) ## short beep
        shiny::removeModal()

        out
      }
    )


    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>WGCNA Analysis Board</strong>"),
        shiny::HTML(infotext),
        size = 'l',
        easyClose = TRUE
      ))
    })


    ## ================================================================================
    ## ======================= PLOTTING FUNCTIONS =====================================
    ## ================================================================================

    ## ----------------------------------------
    ## ------------ samples dendro ------------
    ## ----------------------------------------

    sampleDendro.RENDER <- shiny::reactive({

      out <- wgcna.compute()
      datExpr <- out$datExpr
      pheno <- out$datTraits

      sampleTree2 <- hclust(dist(datExpr), method = "average")
      ipheno <- apply(pheno, 2, function(x) as.numeric(factor(x)))
      colnames(ipheno) <- colnames(pheno)
      rownames(ipheno) <- rownames(pheno)
      traitColors <- WGCNA::numbers2colors(ipheno, signed = FALSE)
      ## Plot the sample dendrogram and the colors underneath.
      par(mfrow = c(1, 1), mar = c(1, 1, 1, 1) * 0)
      WGCNA::plotDendroAndColors(
        sampleTree2, traitColors[, ],
        groupLabels = colnames(ipheno),
        cex.colorLabels = 0.8, cex.dendroLabels = 0.9, cex.rowText = 0.8,
        marAll = c(0.2, 5, 0.2, 0.2),
        ## main = "Sample dendrogram and trait heatmap"
        main = NULL
      )
    })

    sampleDendro_opts <- shiny::tagList()
    sampleDendro_info <- "<b>WGCNA sample dendrogram and trait heatmap.</b>"

    shiny::callModule(
      plotModule,
      id = "sampleDendro", ## ns=ns,
      title = "Sample dendrogram and trait heatmap", label = "b",
      func = sampleDendro.RENDER,
      func2 = sampleDendro.RENDER,
      download.fmt = c("png", "pdf"),
      ## options = sampleDendro_opts,
      info.text = sampleDendro_info,
      height = c(rowH1, 650), width = c("auto", 1000),
      pdf.width = 10, pdf.height = 5, res = c(72, 90),
      add.watermark = WATERMARK
    )

    ## ----------------------------------------
    ## --------- Eigengenes heatmap --------------
    ## ----------------------------------------

    eigenHeatmap.RENDER <- shiny::reactive({

      out <- wgcna.compute()
      net <- out$net
      datExpr <- out$datExpr
      datTraits <- out$datTraits

      MEs <- net$MEs
      rho <- cor(MEs, datExpr)
      sdx <- apply(datExpr, 2, sd)
      ## rho <- t(t(rho) * sdx**2)
      dim(rho)

      if (input$mask_markers) {
        ME <- as.character(net$colors)
        M <- t(model.matrix(~ 0 + ME))[rownames(rho), ]
        rho <- rho * pmax(M, 0.33)
      }

      ntop <- floor(200 / ncol(MEs))
      ntop
      ii <- apply(rho, 1, function(x) head(order(-x), ntop))
      ii <- unique(as.vector(t(ii)))
      ii <- head(ii, 70)

      playbase::gx.heatmap(t(rho[, ii]),
        keysize = 0.2, mar = c(4, 5), key = FALSE,
        cexRow = 0.85, cexCol = 1, scale = "none"
      )
    })

    eigenHeatmap_opts <- shiny::tagList(
      shiny::checkboxInput(ns("mask_markers"), "mask_markers", TRUE)
    )
    eigenHeatmap_info <- "<b>WGCNA Eigengene correlation heatmap.</b> The heatmap shows the correlation of genes to the module eigengenes."

    shiny::callModule(
      plotModule,
      id = "eigenHeatmap", ## ns=ns,
      title = "Eigengene correlation heatmap", label = "a",
      func = eigenHeatmap.RENDER,
      func2 = eigenHeatmap.RENDER,
      download.fmt = c("png", "pdf"),
      options = eigenHeatmap_opts,
      info.text = eigenHeatmap_info,
      height = c(fullH, 650), width = c("auto", 650),
      pdf.width = 6, pdf.height = 10, res = c(72, 90),
      add.watermark = WATERMARK
    )

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

    WATERMARK <- FALSE

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
