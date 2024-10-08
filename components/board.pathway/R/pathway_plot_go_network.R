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
functional_plot_go_network_ui <- function(
    id,
    title,
    info.text,
    info.methods,
    info.references,
    info.extra_link,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(ns("GO_prunetree"), "Prune tree", TRUE),
      "Prune the tree with only significant branches."
    ),
    withTooltip(
      shiny::checkboxInput(ns("GO_colorclusters"), "Color clusters", FALSE),
      "Highlight clusters with different colors."
    )
  )

  PlotModuleUI(ns("plot"),
    title = title,
    label = label,
    plotlib = "visnetwork",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
    options = plot_opts,
    download.fmt = c("pdf", "png"),
    height = height,
    width = width
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
functional_plot_go_network_server <- function(id,
                                              pgx,
                                              fa_contrast,
                                              watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      plot_data <- shiny::reactive({
        fa_contrast <- fa_contrast()
        score <- pgx$meta.go$pathscore[, fa_contrast]
        res <- list(
          score = score,
          fa_contrast = fa_contrast
        )
        return(res)
      })

      plot_RENDER <- shiny::reactive({
        res <- plot_data()
        comparison <- res$fa_contrast
        require(igraph)

        methods <- c("fisher", "gsva", "camera")

        if (is.null(comparison)) {
          return(NULL)
        }

        sub2 <- go <- pgx$meta.go$graph
        if (is.null(go)) {
          shinyalert::shinyalert(
            title = "No GO graph in enrichment results",
            text = "",
            type = "warning"
          )
          return(NULL)
        }

        score <- pgx$meta.go$pathscore[, comparison]
        score[is.na(score) | is.infinite(score)] <- 0
        score <- (score / (1e-8 + max(abs(score), na.rm = TRUE)))
        igraph::V(sub2)$value <- score
        igraph::V(sub2)$color <- playdata::BLUERED(32)[16 + round(15 * score)]
        igraph::V(sub2)$label <- igraph::V(sub2)$Term
        igraph::V(sub2)$label[which(is.na(score) | score == 0)] <- ""
        pos <- sub2$layout

        all.zero <- all(score == 0)

        if (!all.zero && input$GO_prunetree) {
          vv <- igraph::V(sub2)[which(!is.na(score) & abs(score) > 0)]
          sp <- igraph::shortest_paths(sub2, from = "all", to = vv, mode = "all", output = "vpath")
          sp.vv <- unique(unlist(sp$vpath))
          sub2 <- igraph::induced.subgraph(sub2, sp.vv)
          pos <- igraph::layout_with_fr(sub2)
          score <- score[igraph::V(sub2)$name]
        }

        ## remove root?
        removeroot <- TRUE
        if (removeroot) {
          sub2 <- igraph::induced_subgraph(sub2, which(igraph::V(sub2)$name != "all"))
          if (input$GO_prunetree) pos <- igraph::layout_with_fr(sub2)
          score <- score[igraph::V(sub2)$name]
        }
        roots <- c("all", neighbors(go, igraph::V(go)["all"], mode = "all")$name)
        roots <- intersect(roots, igraph::V(sub2)$name)

        astree <- TRUE
        if (astree) {
          if ("all" %in% igraph::V(sub2)$name) {
            pos <- igraph::layout_as_tree(sub2, root = "all", mode = "all")
          } else {
            pos <- igraph::layout_as_tree(sub2, root = roots, mode = "all")
          }
          pos[, 2] <- -pos[, 2]
        }

        ## color clusters
        if (input$GO_colorclusters) {
          clust <- igraph::cluster_louvain(igraph::as.undirected(go))$membership
          names(clust) <- igraph::V(go)$name
          cc <- c(
            RColorBrewer::brewer.pal(12, "Set3"),
            RColorBrewer::brewer.pal(8, "Set2"),
            RColorBrewer::brewer.pal(8, "Set1")
          )
          igraph::V(sub2)$color <- rep(cc, 99)[clust[igraph::V(sub2)$name]]
          jj <- which(is.na(score) | score == 0)
          if (length(jj) > 0) igraph::V(sub2)$color[jj] <- NA
        }

        gr <- visNetwork::toVisNetworkData(sub2)
        gr$nodes$color[is.na(gr$nodes$color)] <- "#F9F9F9"
        gr$nodes$value <- pmax(abs(gr$nodes$value), 0.001)
        gr$nodes$x <- pos[, 1] * 60
        gr$nodes$y <- pos[, 2] * 90
        gr$nodes$label <- gr$nodes$Term
        no.score <- (is.na(score) | score == 0)
        gr$nodes$label[which(no.score)] <- NA

        gr$nodes$shape <- c("box", "circle")[1 + 1 * no.score]
        gr$nodes$label <- sapply(gr$nodes$label, playbase::breakstring, n = 25, nmax = 95, force = TRUE, brk = "\n")

        gr.def <- sapply(gr$nodes$Definition, playbase::breakstring, n = 50, brk = "<br>")
        gr$nodes$title <- paste0(
          gr$nodes$Term, "  (", gr$nodes$id, ")<br>",
          "<small>", gr.def, "</small>"
        )

        ## rendering
        font.size <- 20
        cex <- 1
        if (input$GO_prunetree) {
          font.size <- 20
          cex <- 0.6
        }

        visNetwork::visNetwork(gr$nodes, gr$edges) %>%
          visNetwork::visEdges(
            smooth = FALSE,
            hidden = FALSE,
            arrows = list(enabled = TRUE),
            scaling = list(min = 10 * cex, max = 30 * cex),
            width = 5 * cex
          ) %>%
          visNetwork::visNodes(
            font = list(size = font.size * cex, vadjust = 0),
            scaling = list(min = 1 * cex, max = 80 * cex)
          ) %>%
          visNetwork::visPhysics(enabled = FALSE, stabilization = FALSE) %>%
          visNetwork::visOptions(highlightNearest = list(enabled = T, degree = 1, hover = TRUE)) %>%
          visNetwork::visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE)
      })

      PlotModuleServer(
        "plot",
        plotlib = "visnetwork",
        func = plot_RENDER,
        csvFunc = plot_data,
        res = 72,
        pdf.width = 10, pdf.height = 8,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
