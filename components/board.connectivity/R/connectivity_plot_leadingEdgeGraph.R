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
connectivity_plot_leadingEdgeGraph_ui <- function(id,
                                                  label = "",
                                                  rowH = 600) {
  ns <- shiny::NS(id)

  info_text <- strwrap(
    "<b>Leading-edge graph.</b> Network of shared leading-edge genes between
    top-N most similar signatures. The edge width corresponds to the number of
    signatures that share that pair of genes in their top differentially expressed
    genes. In the plot options you can set the threshold of the edges."
  )

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::sliderInput(ns("LEgraph_threshold"), "edge threshold:", 0, 1, 0, 0.01),
      "Threshold value for edges."
    ),
    withTooltip(
      shiny::radioButtons(ns("LEgraph_ntop"), "N-neighbours:", c(5, 10, 25, 100),
        selected = 10, inline = TRUE
      ),
      "Number of simlar experiments to consider."
    ),
    withTooltip(
      shiny::radioButtons(ns("LEgraph_sizevar"), "Size:", c("FC", "cumFC", "centrality"),
        selected = "cumFC", inline = TRUE
      ),
      "Parameter for node size."
    )
  )

  PlotModuleUI(ns("plot"),
    title = "Leading-edge graph",
    label = label,
    plotlib = "visnetwork",
    info.text = info_text,
    options = plot_opts,
    height = c(680, 720),
    width = c("auto", 1300)
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
connectivity_plot_leadingEdgeGraph_server <- function(id,
                                                      getConnectivityScores,
                                                      connectivityScoreTable,
                                                      getCurrentContrast,
                                                      getTopProfiles,
                                                      watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      cumulativeFCtable <- shiny::reactive({
        F <- getTopProfiles()
        F[is.na(F)] <- 0

        ## maximum 10??
        MAXF <- 20

        ## multiply with sign of rho
        df <- getConnectivityScores()
        rho1 <- df$rho[match(colnames(F), df$pathway)]
        F <- t(t(F) * sign(rho1))

        ## add current contrast
        cc <- getCurrentContrast()
        shiny::req(cc)
        fc <- cc$fc[rownames(F)]
        fc[is.na(fc)] <- 0
        F <- cbind(fc[rownames(F)], F)
        colnames(F)[1] <- "thisFC"
        colnames(F)[1] <- cc$name

        F <- F[order(-rowMeans(F**2)), , drop = FALSE]
        F
      })

      getLeadingEdgeGraph <- shiny::reactive({
        df <- getConnectivityScores()
        if (is.null(df)) {
          return(NULL)
        }
        df$score[is.na(df$score)] <- 0
        df <- df[which(df$score > 0), ]

        ## always only top-N from selection
        ntop <- as.integer(input$LEgraph_ntop)
        ii <- connectivityScoreTable$rows_all()
        shiny::req(ii)
        ii <- head(order(-abs(df$score)), 25)
        ii <- head(ii, ntop)
        df <- df[ii, , drop = FALSE]

        le.genes <- sort(unique(unlist(df$leadingEdge)))
        A <- 1 * sapply(df$leadingEdge, function(g) le.genes %in% g)
        rownames(A) <- le.genes
        A <- head(A[order(-rowMeans(A)), , drop = FALSE], 100)
        adjM <- (A %*% t(A))
        adjM <- adjM / max(adjM, na.rm = TRUE)

        gr <- igraph::graph_from_adjacency_matrix(
          adjM,
          mode = "undirected", weighted = TRUE, diag = FALSE
        )

        ## set graph threshold to some sensible value [0,1]
        wt0 <- tail(sort(abs(igraph::E(gr)$weight)), 150)[1] ## about 150 edges
        shiny::updateSliderInput(session, "LEgraph_threshold", value = 0)

        return(gr)
      })

      plot_RENDER <- shiny::reactive({
        gr <- getLeadingEdgeGraph()
        if (is.null(gr)) {
          return(NULL)
        }
        minwt <- 0.5
        minwt <- input$LEgraph_threshold
        minwt <- min(c(minwt, 0.99 * max(abs(igraph::E(gr)$weight), na.rm = TRUE)))
        gr <- igraph::subgraph.edges(gr, which(abs(igraph::E(gr)$weight) >= minwt))
        if (length(igraph::V(gr)) == 0) {
          return(NULL)
        }

        fc <- cumFC <- NULL
        fc <- getCurrentContrast()$fc
        fc <- fc[igraph::V(gr)$name]
        cumFC <- cumulativeFCtable()

        cumFC <- cumFC[igraph::V(gr)$name, ]
        fontsize <- 22
        fc <- fc / max(abs(fc))

        sizevar <- input$LEgraph_sizevar
        vsize <- 15
        if (sizevar == "centrality") {
          vsize <- log(1 + igraph::betweenness(gr))
        } else if (sizevar == "cumFC") {
          fc1 <- rowMeans(cumFC)
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
        gene <- igraph::V(gr)$name
        igraph::V(gr)$label <- igraph::V(gr)$name
        igraph::V(gr)$title <- paste0("<b>", gene, "</b><br>", GENE.TITLE[toupper(gene)])
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


      plt <- PlotModuleServer(
        "plot",
        plotlib = "visnetwork",
        func = plot_RENDER,
        func2 = plot_RENDER,
        csvFunc = getLeadingEdgeGraph,
        pdf.height = 8, pdf.width = 8,
        res = c(90, 100),
        add.watermark = watermark
      )
      return(getLeadingEdgeGraph)
    } ## end of moduleServer
  )
}
