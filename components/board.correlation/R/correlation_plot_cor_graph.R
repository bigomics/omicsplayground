##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Expression plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
correlation_plot_cor_graph_ui <- function(
    id,
    title,
    caption,
    info.text,
    info.methods,
    info.references,
    info.extra_link,
    height,
    width) {
  ns <- shiny::NS(id)

  cor_graph.opts <- shiny::tagList(
    shiny::checkboxInput(ns("physics"), "physics", FALSE),
    shiny::checkboxInput(ns("as_mst"), "minimum spanning tree", FALSE),    
    shiny::sliderInput(ns("graph_threshold"), "correlation threshold:", 0, 0.5, 0.1),
    shiny::sliderInput(ns("graph_radius"), "maximum radius:", 1, 5, 3, 1)
    ##shiny::sliderInput(ns("edge_degree"), "edge connectivity:", 0.1, 2, 1, 0.1),    
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    plotlib = "visnetwork",
    label = "a",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
    options = cor_graph.opts,
    width = width,
    height = height,
    download.fmt = c("png", "pdf", "csv", "svg"),
  )
}

#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
correlation_plot_cor_graph_server <- function(
    id,
    gene,
    getPartialCorrelationMatrix,
    watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- shiny::reactive({
      shiny::req(gene())

      res <- getPartialCorrelationMatrix()
      gene <- "XIST"
      rho.min <- 0.3

      gene <- gene()
      rho.min <- input$graph_threshold
      numnodes <- nrow(res$cor)
      vsize <- ifelse(numnodes > 50, 15, 20)
      vsize <- ifelse(numnodes > 100, 10, vsize)
      radius <- as.integer(input$graph_radius)

      gr <- playbase::pgx.plotPartialCorrelationGraph(
        res,
        gene, ## what="graph", ## degree=deg,
        plot = FALSE,
        rho.min = rho.min,
        nb = Inf,
        nsize = -1,
        radius = radius,
        vsize = vsize,
        edge.width = 10
      )
      gr
    })

    cor_graph.VISNETWORK <- function() {
      gr <- plot_data()
      if (is.null(gr)) {
        return(NULL)
      }

      if(input$as_mst) {
        gr <- igraph::mst(gr)
      }
      
      visdata <- visNetwork::toVisNetworkData(gr, idToLabel = FALSE)
      visdata$edges$width <- 2 * visdata$edges$width

      if (nrow(visdata$edges) == 0 && nrow(visdata$nodes) == 1) {
        visdata$edges <- data.frame(
          from = rownames(visdata$nodes),
          to   = rownames(visdata$nodes),
          weight = 1,
          rho = 1,
          width = 1,
          color = "fff"
        )
      }
      
      graph <- visNetwork::visNetwork(
        nodes = visdata$nodes,
        edges = visdata$edges
      ) %>%
        visNetwork::visNodes(font = list(size = 12)) %>%
        visNetwork::visEdges(hidden = FALSE, width = 4, color = list(opacity = 0.9)) %>%
        visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
        visNetwork::visIgraphLayout(
          layout = "layout_with_kk",
          physics = input$physics
        )

      return(graph)
    }

    plot_csv_data <- function() {
      gr <- plot_data()
      if (is.null(gr)) {
        return(NULL)
      }
      visdata <- visNetwork::toVisNetworkData(gr, idToLabel = FALSE)
      visdata <- visdata$edges[, 1:5]
      return(visdata)
    }

    PlotModuleServer(
      "plot",
      plotlib = "visnetwork",
      func = cor_graph.VISNETWORK,
      csvFunc = plot_csv_data,
      res = c(72, 80), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
