##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Importance plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
connectivity_plot_enrichmentGraph_ui <- function(
  id,
  title,
  info.text,
  caption,
  label = "",
  height,
  width
) {
  ns <- shiny::NS(id)

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::sliderInput(ns("enrichGraph_threshold"), "edge threshold:", 0, 1, 0.3, 0.01),
      "Threshold value for edges."
    ),
    hr(),
    withTooltip(
      shiny::radioButtons(ns("enrichGraph_ntop"), "N-neighbours:", c(5, 10, 25, 100),
        selected = 10, inline = TRUE
      ),
      "Number of simlar experiments to consider."
    ),
    hr(),
    withTooltip(
      shiny::checkboxInput(ns("enrichGraph_oddweighting"), "Odd ratio weighting", FALSE),
      "Odds ratio weighting."
    ),
    hr(),
    withTooltip(
      shiny::radioButtons(ns("enrichGraph_sizevar"), "Size:", c("FC", "cumFC", "centrality"),
        selected = "cumFC", inline = TRUE
      ),
      "Parameter for node size."
    )
  )

  PlotModuleUI(ns("plot"),
    title = title,
    label = label,
    plotlib = "visnetwork",
    info.text = info.text,
    caption = caption,
    options = plot_opts,
    width = width,
    height = height
  )
}

#' Importance plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
connectivity_plot_enrichmentGraph_server <- function(id,
                                                     getLeadingEdgeGraph,
                                                     getConnectivityScores,
                                                     connectivityScoreTable,
                                                     cumEnrichmentTable,
                                                     watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      getEnrichmentGraph <- shiny::reactive({
        ## get enrichment scores
        F <- cumEnrichmentTable()
        shiny::req(F)

        if (input$enrichGraph_oddweighting) {
          gr2 <- getLeadingEdgeGraph()
          le.genes <- igraph::V(gr2)$name
          gsets <- playdata::getGSETS(rownames(F))
          gsets <- gsets[sapply(gsets, length) >= 5]
          bg <- unique(unlist(gsets))
          ft <- playbase::gset.fisher(le.genes, gsets,
            fdr = 1.0,
            min.genes = 3, max.genes = 99999,
            background = bg
          )
          ft <- ft[match(rownames(F), rownames(ft)), ]
          F <- F * log(1 + ft$odd.ratio)
        }

        ## get top average-enriched
        ntop <- as.integer(input$enrichGraph_ntop)
        le.sets <- apply(F, 2, function(x) head(order(-abs(x)), ntop))
        le.sets <- apply(le.sets, 2, function(i) rownames(F)[i])
        A <- 1 * apply(le.sets, 2, function(g) rownames(F) %in% g)
        rownames(A) <- rownames(F)

        ## create graph
        A <- head(A[order(-rowMeans(A, na.rm = TRUE)), , drop = FALSE], 100)
        adjM <- (A %*% t(A))
        adjM <- adjM / max(adjM, na.rm = TRUE)

        gr <- igraph::graph_from_adjacency_matrix(
          adjM,
          mode = "undirected", weighted = TRUE, diag = FALSE
        )

        ## set graph threshold to some sensible value [0,1]
        #
        #

        return(gr)
      })

      plot_RENDER <- shiny::reactive({
        gr <- getEnrichmentGraph()
        if (is.null(gr)) {
          return(NULL)
        }

        gr2 <- getLeadingEdgeGraph()
        pw <- igraph::V(gr)$name
        le.genes <- igraph::V(gr2)$name

        gsets <- playdata::getGSETS(pw)
        pw.genes <- sapply(gsets, function(gs) intersect(gs, le.genes))
        pw.genes <- sapply(pw.genes, paste, collapse = " ")

        minwt <- 0.5
        minwt <- input$enrichGraph_threshold
        minwt <- min(c(minwt, 0.99 * max(abs(igraph::E(gr)$weight), na.rm = TRUE)))
        gr <- igraph::subgraph.edges(gr, which(abs(igraph::E(gr)$weight) >= minwt))
        if (length(igraph::V(gr)) == 0) {
          return(NULL)
        }

        fc <- cumFC <- NULL
        cumFC <- cumEnrichmentTable()
        cumFC <- cumFC[igraph::V(gr)$name, ]
        fc <- cumFC[, 1]
        fontsize <- 18
        fc <- fc / max(abs(fc))

        sizevar <- input$enrichGraph_sizevar
        vsize <- 15
        if (sizevar == "centrality") {
          vsize <- log(1 + igraph::betweenness(gr))
        } else if (sizevar == "cumFC") {
          fc1 <- rowMeans(cumFC, na.rm = TRUE)
          vsize <- abs(fc1[match(igraph::V(gr)$name, names(fc1))])**2
        } else if (sizevar == "FC") {
          vsize <- abs(fc[match(igraph::V(gr)$name, names(fc))])**2
        } else {
          vsize <- 1
        }
        vsize <- 3 + 12 * (abs(vsize) / max(abs(vsize), na.rm = TRUE))**0.5

        bluered.pal <- colorRampPalette(
          colors = c("royalblue4", "royalblue2", "grey90", "indianred3", "firebrick4")
        )
        vcolor <- bluered.pal(65)[33 + round(32 * fc)]
        vcolor <- paste0(vcolor, "AA") ## add transparency
        ## defaults graph parameters
        vname <- sub("H:HALLMARK_|C2:KEGG_", "", igraph::V(gr)$name)
        # Select the same names(pw.genes) as in `vname`, remove duplicates
        pw.genes.selector <- pw.genes[which(names(pw.genes) %in% igraph::V(gr)$name)] %>%
          names() %>%
          duplicated() %>%
          `!`()
        pw.genes <- pw.genes[which(names(pw.genes) %in% igraph::V(gr)$name)][pw.genes.selector]

        igraph::V(gr)$label <- vname
        igraph::V(gr)$title <- paste0("<b>", vname, "</b><br>", pw.genes)
        igraph::V(gr)$size <- vsize ## rather small
        igraph::V(gr)$color <- vcolor

        ew <- abs(igraph::E(gr)$weight)
        igraph::E(gr)$width <- 1.5 * (0.2 + 10 * (ew / max(ew, na.rm = TRUE))**2)
        igraph::E(gr)$color <- "#DDD" ## lightgrey

        visdata <- visNetwork::toVisNetworkData(gr, idToLabel = FALSE)

        ## ------------------ plot using visNetwork (zoomable) -----------------
        graph <- visNetwork::visNetwork(
          nodes = visdata$nodes,
          edges = visdata$edges
        ) %>%
          visNetwork::visNodes(font = list(size = fontsize)) %>%
          visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
          visNetwork::visIgraphLayout(layout = "layout_nicely")
        graph
      })


      PlotModuleServer(
        "plot",
        plotlib = "visnetwork",
        func = plot_RENDER,
        func2 = plot_RENDER,
        csvFunc = getEnrichmentGraph,
        pdf.height = 8, pdf.width = 8,
        res = c(90, 100),
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
