##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

LasagnaBoard <- function(id, pgx) {
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


    ## ============================================================================
    ## ============================ OBSERVERS =====================================
    ## ============================================================================

    infotext <-
      '<center><iframe width="1120" height="630" src="https://www.youtube.com/embed/rRIRMW_RRS4"
        title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write;
        encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></center>'

    shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>LASAGNA Analysis Board</strong>"),
        shiny::HTML(infotext),
        size = "xl",
        easyClose = TRUE
      ))
    })

    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Multi-layer model" = list(disable = c("mpartite_options", "gsfilter")),
      "Multi-partite graph" = list(disable = c("clust_options"))
      # "Multi-type network" = list(disable = c("clust_options"))
    )

    shiny::observeEvent(input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })

    ## ============================================================================
    ## ============================ REACTIVES =====================================
    ## ============================================================================

    shiny::observeEvent(pgx$mofa,
      {
        shiny::validate(shiny::need(!is.null(pgx$mofa), "missing MOFA slot"))

        ## update factors in selectInput
        ct1 <- colnames(pgx$mofa$contrasts)
        ct2 <- colnames(playbase::pgx.getMetaMatrix(pgx)$fc)
        contrasts <- intersect(ct1, ct2)
        contrasts <- contrasts[!grepl("^IA", contrasts)] ## no interaction contrasts
        updateSelectInput(session, "contrast",
          choices = contrasts,
          selected = contrasts[1]
        )

        datatypes <- setdiff(names(pgx$mofa$xx), c("gset"))
        datatypes2 <- unique(c(datatypes, "gset"))
        updateSelectInput(session, "layers",
          choices = datatypes2,
          selected = datatypes
        )
      },
      ignoreNULL = FALSE
    )

    lasagna_model <- shiny::eventReactive(
      {
        list(input$updateplots, pgx$X)
      },
      {
        shiny::validate(shiny::need(!is.null(pgx$mofa), "missing MOFA slot"))
        shiny::validate(shiny::need(pgx$datatype == "multi-omics", "not multi-omics data"))

        shiny::req(pgx$X)
        shiny::req(pgx$mofa)

        progress <- shiny::Progress$new(session, min = 0, max = 1)
        on.exit(progress$close())

        ## xdata <- playbase::mofa.split_data(pgx$X)
        layers <- input$layers
        add_gsets <- "gset" %in% layers
        gsetX <- NULL
        if (add_gsets) {
          gsetX <- pgx$gsetX
        }

        progress$set(message = paste("creating LASAGNA model"), value = 0.33)

        res <- playbase::lasagna.create_from_pgx(
          pgx,
          xdata = NULL,
          add_gsets = add_gsets,
          gsetX = gsetX,
          pheno = "contrasts",
          ntop = 2000,
          nc = 20,
          layers = layers,
          add.sink = TRUE,
          intra = TRUE,
          fully_connect = FALSE,
          add.revpheno = TRUE
        )

        progress$set(message = paste("computing positions..."), value = 0.33)

        layers <- setdiff(res$layers, c("SOURCE", "SINK"))
        posx <- list()
        posf <- list()
        k <- layers[1]
        for (k in layers) {
          ii <- which(igraph::V(res$graph)$layer == k)
          cx <- res$X[ii, , drop = FALSE]
          cx <- cx - rowMeans(cx, na.rm = TRUE)
          sv <- svd(cx, nv = 2, nu = 2)
          posx[[k]] <- sv$v[, 1:2]
          posf[[k]] <- sv$u[, 1:2]
          colnames(posf[[k]]) <- c("PC1", "PC2")
          if (input$clustmethod == "umap" && nrow(cx) > 10) {
            nb <- max(min(15, nrow(cx) / 3), 2)
            posf[[k]] <- uwot::umap(cx, n_neighbors = nb)
            colnames(posf[[k]]) <- c("UMAP-x", "UMAP-y")
          }
          if (input$clustmethod == "tsne" && nrow(cx) > 10) {
            nb <- max(min(15, nrow(cx) / 5), 2)
            posf[[k]] <- Rtsne::Rtsne(cx,
              check_duplicates = FALSE,
              perplexity = nb
            )$Y
            colnames(posf[[k]]) <- c("tSNE-x", "tSNE-y")
          }
          rownames(posx[[k]]) <- colnames(res$X)
          rownames(posf[[k]]) <- rownames(cx)
        }
        res$posf <- posf
        res$posx <- posx

        gs <- NULL
        ii <- which(igraph::V(res$graph)$layer == "gset")
        if (length(ii)) {
          gs <- gsub("gset:|:.*|_.*", "", igraph::V(res$graph)$name[ii])
          gs <- names(which(table(gs) >= 3))
          gs <- c("*", gs)
        }
        shiny::updateSelectInput(session, "gsfilter", choices = gs)

        return(res)
      },
      ignoreNULL = FALSE
    )

    solved_data <- reactive({
      shiny::req(input$contrast)

      res <- lasagna_model()
      pheno <- input$contrast
      value <- input$node_value

      progress <- shiny::Progress$new(session, min = 0, max = 1)
      on.exit(progress$close())
      progress$set(message = paste("computing connections..."), value = 0.66)

      solved <- playbase::lasagna.solve(
        res, pheno,
        value = value,
        min_rho = 0.1,
        max_edges = 1000,
        fc.weight = TRUE,
        sp.weight = input$sp_weight,
        prune = FALSE
      )
      res$graph <- solved

      return(res)
    })

    pruned_data <- reactive({
      res <- solved_data()
      ntop <- as.integer(input$ntop)
      gsfilter <- NULL

      if (!input$gsfilter %in% c("", "*")) {
        gsfilter <- list(gset = input$gsfilter)
      }

      edge.sign <- ifelse(input$consensus, "consensus", "both")

      ## prune graph for plotting
      pruned <- playbase::lasagna.prune_graph(
        res$graph,
        ntop = ntop,
        layers = NULL,
        filter = gsfilter,
        normalize.edges = TRUE,
        min.rho = input$minrho,
        edge.sign = edge.sign,
        edge.type = "both",
        prune = FALSE
      )

      res$graph <- pruned
      return(res)
    })


    ## ==========================================================================
    ## ========================== BOARD FUNCTIONS ===============================
    ## ==========================================================================


    ## ==========================================================================
    ## =========================== MODULES ======================================
    ## ==========================================================================

    mofa_plot_lasagna3D_server(
      "lasagna",
      data = solved_data,
      pgx = pgx,
      watermark = WATERMARK
    )

    mofa_plot_lasagna_clustering_server(
      "clusters",
      data = solved_data,
      ## pgx = pgx,
      input_contrast = reactive(input$contrast),
      watermark = WATERMARK
    )

    mofa_plot_lasagna_partite_server(
      "lasagnaPartite",
      data = pruned_data,
      pgx = pgx,
      input_nodevalue = reactive(input$node_value),
      watermark = WATERMARK
    )

    multipartite_nodes_table <- lasagna_multipartite_nodes_table_server(
      "multipartite_nodes_table",
      data = pruned_data,
      pgx = pgx,
      scrollY = "calc(100vh - (240px + 140px))"
    )

    multipartite_edges_table <- lasagna_multipartite_edges_table_server(
      "multipartite_edges_table",
      data = pruned_data,
      scrollY = "calc(100vh - (240px + 140px))"
    )

    return(NULL)
  })
} ## end of Board
