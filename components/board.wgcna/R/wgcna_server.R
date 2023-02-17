##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


WgcnaBoard <- function(id, inputData) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    fullH <- 700 ## full height of page
    rowH1 <- 250 ## row 1 height
    rowH2 <- 440 ## row 2 height

    infotext <- "Weighted gene co-expression network analysis (WGCNA) is a systems biology method for describing the correlation patterns among genes across microarray samples. Weighted correlation network analysis can be used for finding clusters (modules) of highly correlated genes, for summarizing such clusters using the module eigengene or an intramodular hub gene, for relating modules to one another and to external sample traits (using eigengene network methodology), and for calculating module membership measures. Correlation networks facilitate network based gene screening methods that can be used to identify candidate biomarkers or therapeutic targets."

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
        ngs <- inputData()
        1
      },
      {
        ngs <- inputData()
        require(WGCNA)

        if ("wgcna" %in% names(ngs)) {
          message("[wgcna.compute] >>> using pre-computed WGCNA results...")
          return(ngs$wgcna)
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

        X <- as.matrix(ngs$X)
        dim(X)
        X <- X[order(-apply(X, 1, sd, na.rm = TRUE)), ]
        X <- X[!duplicated(rownames(X)), ]

        minmodsize <- 30
        power <- 6
        cutheight <- 0.25
        deepsplit <- 2
        ngenes <- 1000

        ngenes <- input$ngenes
        minmodsize <- as.integer(input$minmodsize)
        power <- as.numeric(input$power)
        cutheight <- as.numeric(input$cutheight)
        deepsplit <- as.integer(input$deepsplit)

        datExpr <- t(head(X, ngenes))
        progress$inc(0.1, "Computing WGCNA modules...")
        require(WGCNA)
        net <- WGCNA::blockwiseModules(
          datExpr,
          power = power,
          TOMType = "unsigned", minModuleSize = minmodsize,
          reassignThreshold = 0, mergeCutHeight = cutheight,
          numericLabels = TRUE, pamRespectsDendro = FALSE,
          deepSplit = deepsplit,
          ## saveTOMs = TRUE, saveTOMFileBase = "WCNA.tom",
          verbose = 3
        )
        names(net)
        table(net$colors)

        ## clean up traits matrix
        datTraits <- ngs$samples
        ## no dates please...
        isdate <- apply(datTraits, 2, is.Date)
        datTraits <- datTraits[, !isdate, drop = FALSE]

        ## Expand multi-class discrete phenotypes into binary vectors
        ## datTraits1 <- datTraits
        tr.class <- sapply(type.convert(datTraits), class)
        sel1 <- which(tr.class %in% c("factor", "character"))
        sel2 <- which(tr.class %in% c("integer", "numeric"))

        tr1 <- datTraits[, 0]
        if (length(sel1)) {
          tr1 <- expandPhenoMatrix(datTraits[, sel1, drop = FALSE], drop.ref = FALSE)
        }
        ## keeping numeric phenotypes
        tr2 <- datTraits[, sel2, drop = FALSE]
        datTraits <- cbind(tr1, tr2)

        ## get colors of eigengene modules
        me.genes <- tapply(names(net$colors), net$colors, list)
        names(me.genes) <- paste0("ME", names(me.genes))
        color1 <- labels2rainbow(net)
        me.colors <- color1[!duplicated(color1)]
        names(me.colors) <- paste0("ME", names(me.colors))
        me.colors <- me.colors[names(me.genes)]
        progress$inc(0.4, "")

        if (1) {
          message("[wgcna.compute] >>> Calculating WGCNA clustering...")
          progress$inc(0.1, "Computing dim reductions...")

          X1 <- t(datExpr)
          X1 <- t(scale(datExpr))
          ## pos <- Rtsne::Rtsne(X1)$Y
          ## dissTOM <- 1 - abs(cor(datExpr))**6
          dissTOM <- 1 - WGCNA::TOMsimilarityFromExpr(datExpr, power = power)
          rownames(dissTOM) <- colnames(dissTOM) <- colnames(datExpr)
          clust <- pgx.clusterBigMatrix(dissTOM, methods = c("umap", "tsne", "pca"), dims = c(2))
          ## pos <- pgx.clusterBigMatrix(t(X1), methods="tsne", dims=2)[[1]]
          ## pos <- pgx.clusterBigMatrix(dissTOM, methods="pca", dims=2)[[1]]
          names(clust)
          if ("cluster.genes" %in% names(ngs)) {
            clust[["umap2d"]] <- ngs$cluster.genes$pos[["umap2d"]][colnames(datExpr), ]
          }
          progress$inc(0.2)
        }

        if (1) {
          message("[wgcna.compute] >>> Calculating WGCNA module enrichments...")
          progress$inc(0, "Calculating module enrichment...")

          gmt <- getGSETS(grep("HALLMARK|GOBP|^C[1-9]", names(iGSETS), value = TRUE))
          gse <- NULL
          ## bg <- unlist(me.genes)
          bg <- toupper(rownames(ngs$X))
          i <- 1
          for (i in 1:length(me.genes)) {
            gg <- toupper(me.genes[[i]])
            rr <- gset.fisher(gg, gmt, background = bg, fdr = 1)
            rr <- cbind(
              module = names(me.genes)[i],
              geneset = rownames(rr), rr
            )
            rr <- rr[order(rr$p.value), , drop = FALSE]
            if (i == 1) gse <- rr
            if (i > 1) gse <- rbind(gse, rr)
          }
          rownames(gse) <- NULL

          progress$inc(0.3)
        }

        ## construct results object
        out <- list(
          datExpr = datExpr,
          datTraits = datTraits,
          net = net,
          gse = gse,
          clust = clust,
          me.genes = me.genes,
          me.colors = me.colors
        )

        shiny::updateSelectInput(session, "selected_module", choices = names(me.genes), sel = "ME1")

        message("[wgcna.compute] >>> done!")
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
        easyClose = TRUE
      ))
    })


    ## ================================================================================
    ## ======================= PLOTTING FUNCTIONS =====================================
    ## ================================================================================

    ## ----------------------------------------
    ## ------------ samples dendro ------------
    ## ----------------------------------------

    labels2rainbow <- function(net) {
      hc <- net$dendrograms[[1]]
      nc <- length(unique(net$colors))
      n <- length(net$colors)
      ii <- hc$order
      col1 <- labels2colors(net$colors)
      col.rnk <- rank(tapply(1:n, col1[ii], mean))
      new.col <- rainbow(nc)[col.rnk]
      ## new.col <- heat.colors(nc)[col.rnk]
      names(new.col) <- names(col.rnk)
      new.col["grey"] <- "#AAAAAA"
      new.col
      new.col <- new.col[col1]
      names(new.col) <- net$colors
      new.col
    }

    ## ----------------------------------------
    ## ------------ samples dendro ------------
    ## ----------------------------------------

    sampleDendro.RENDER <- shiny::reactive({
      message("[sampleDendro.RENDER] reacted")

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
      message("[eigenHeatmap.RENDER] reacted")

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

      gx.heatmap(t(rho[, ii]),
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
      labels2rainbow = labels2rainbow,
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
      labels2rainbow = labels2rainbow,
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
      labels2rainbow = labels2rainbow,
      watermark = WATERMARK
    )

    # Module-Trait relationships

    wgcna_plot_MTrelationships_server(
      "moduleTrait",
      wgcna.compute = wgcna.compute,
      labels2rainbow = labels2rainbow,
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
