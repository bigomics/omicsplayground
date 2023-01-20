##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics Sagl. All rights reserved.
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
correlation_plot_cor_graph_ui <- function(id,
                                          height = c(600, 800)) {
  ns <- shiny::NS(id)
  info_text <- "<b>Partial correlation network.</b> Partial correlation graph centered on selected gene with top most correlated features. Red edges correspond to negative correlation, grey edges to positive correlation. Width of the edges is proportional to the absolute partial correlation value of the gene pair."
  GRAPH.LAYOUTS <- c(
    "Fruchterman-Reingold" = "fr", "Kamada-Kawai" = "kk",
    "graphopt" = "graphopt", "tree layout" = "tree"
  )
  cor_graph.opts <- shiny::tagList(
    shiny::sliderInput(ns("cor_graph_radius"), "radius:", 1, 8, 2, 1),
    shiny::sliderInput(ns("cor_graph_threshold"), "pcor threshold:", 0, 1, 0.90),
    shiny::selectInput(ns("cor_graph_layout"), "layout:", choices = GRAPH.LAYOUTS)
  )

  PlotModuleUI(ns("plot"),
    title = "Partial correlation network",
    plotlib = "generic",
    outputFunc = uiOutput,
    outputFunc2 = uiOutput,
    label = "a",
    info.text = info_text,
    options = cor_graph.opts,
    download.fmt = c("png", "pdf", "csv"),
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
correlation_plot_cor_graph_server <- function(id,
                                              cor_gene,
                                              getPartialCorrelationMatrix,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    getCorGraph <- shiny::reactive({
      shiny::req(cor_gene)

      res <- getPartialCorrelationMatrix()
      gene <- "XIST"
      rho.min <- 0.3
      layout <- "kk"
      gene <- cor_gene

      rho.min <- input$cor_graph_threshold

      layout <- input$cor_graph_layout
      ## fixed <- input$cor_graph_fixed
      numnodes <- nrow(res$cor)
      vsize <- ifelse(numnodes > 50, 10, 12)
      vsize <- ifelse(numnodes > 100, 8, vsize)

      radius <- as.integer(input$cor_graph_radius)

      gr <- pgx.plotPartialCorrelationGraph(
        res, gene, ## what="graph", ## degree=deg,
        plot = FALSE,
        rho.min = rho.min, nsize = -1,
        layout = layout, radius = radius,
        vsize = vsize, edge.width = 10
      )
      gr
    })

    # reactive function listeninng for changes in input
    plot_data <- shiny::reactive({
      getCorGraph()
    })

    cor_graph.VISNETWORK <- shiny::reactive({
      gr <- plot_data()
      if (is.null(gr)) {
        return(NULL)
      }

      visdata <- visNetwork::toVisNetworkData(gr, idToLabel = FALSE)
      visdata$edges$width <- 2 * visdata$edges$width

      graph <- visNetwork::visNetwork(
        nodes = visdata$nodes,
        edges = visdata$edges
      ) %>%
        visNetwork::visNodes(font = list(size = 12)) %>%
        visNetwork::visEdges(hidden = FALSE, width = 4, color = list(opacity = 0.9)) %>%
        visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
        visNetwork::visIgraphLayout(layout = "layout_nicely")
      graph
    })

    func <- shiny::reactive({
      visNetwork::visNetworkOutput(
        outputId = ns("cor_graph"),
        width = 700,
        height = 700
      )
    })

    func2 <- shiny::reactive({
      visNetwork::visNetworkOutput(
        outputId = ns("cor_graph2"),
        width = 720,
        height = 1000
      )
    })

    output$cor_graph <- visNetwork::renderVisNetwork(cor_graph.VISNETWORK())
    output$cor_graph2 <- visNetwork::renderVisNetwork(cor_graph.VISNETWORK())

    PlotModuleServer(
      "plot",
      plotlib = "generic",
      renderFunc = renderUI,
      renderFunc2 = renderUI,
      func = func,
      func2 = func2,
      csvFunc = plot_data,
      res = c(72, 80), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
