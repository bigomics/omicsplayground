##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_modulegraph_ui <- function(
  id,
  title = "",
  info.text = "",
  caption = "",
  label = "",
  height = 400,
  width = 400
) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::checkboxInput(ns("mst"), "Min. spanning tree (MST)", TRUE),
    shiny::checkboxInput(ns("rm_singles"), "Remove singletons", TRUE),
    shiny::checkboxInput(ns("top20"), "Top 20", FALSE)
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    plotlib = "visnetwork",
    download.fmt = c("png", "pdf", "svg", "cx2")
  )
}

mofa_plot_modulegraph_server <- function(id,
                                         mofa,
                                         pgx,
                                         input_k = reactive(1),
                                         filter_types = reactive(NULL),
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    moduleGraphData <- shiny::reactive({
      graphs <- mofa()$graphs
      k <- input_k()
      shiny::req(k %in% names(graphs$features))
      shiny::req(filter_types())
      shiny::validate(shiny::need(
        length(filter_types()) > 0,
        "Please select at least one datatype"
      ))

      gr <- graphs$features[[k]]
      vtype <- sub(":.*", "", igraph::V(gr)$name)
      gr <- igraph::subgraph(gr, vids = vtype %in% filter_types())
      if (input$top20) {
        gr <- igraph::subgraph(gr, vids = head(order(-igraph::V(gr)$size), 20))
      }

      ## set labels
      vlabel <- pgx$genes[igraph::V(gr)$name, "gene_name"]
      igraph::V(gr)$label <- vlabel
      gr
    })

    plot.RENDER <- function(springLength = 25) {
      gr <- moduleGraphData()

      par(mar = c(0, 0, 0, 0))
      vis <- playbase::mofa.plot_module(
        gr,
        mst = input$mst,
        rm.single = input$rm_singles,
        nlabel = 999,
        cex = 0.5,
        plotlib = "visnet"
      )

      vis <- vis %>%
        visNetwork::visPhysics(
          barnesHut = list(
            gravitationalConstant = -1000,
            centralGravity = 0.3,
            springLength = springLength
          )
        )

      vis
    }

    plot.RENDER2 <- function() {
      plot.RENDER(springLength = 50)
    }

    ## CX2 export for Cytoscape (mirrors the visNetwork attributes)
    cx2_content <- function(file) {
      gr <- moduleGraphData()
      shiny::req(gr)
      ## mirror the displayed network's mst / singleton-removal options
      if (input$mst) {
        gr <- igraph::mst(gr)
      }
      if (input$rm_singles) {
        gr <- igraph::delete_vertices(gr, which(igraph::degree(gr) == 0))
      }
      xy <- igraph::layout_nicely(gr)
      igraph::V(gr)$x <- xy[, 1]
      igraph::V(gr)$y <- xy[, 2]
      write_cx2(gr, file)
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      func2 = plot.RENDER2,
      plotlib = "visnetwork",
      pdf.width = 8, pdf.height = 8,
      res = c(80, 100),
      add.watermark = watermark,
      download.handlers = list(
        cx2 = list(ext = "cx2", content = cx2_content)
      )
    )
  })
}
