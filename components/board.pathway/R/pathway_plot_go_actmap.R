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
functional_plot_go_actmap_ui <- function(
    id,
    title,
    info.text,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("normalize"),
        "Normalize columns",
        FALSE
      ),
      "Click to normalize the columns of the activation matrices."
    ),
    withTooltip(
      shiny::checkboxInput(
        ns("rotate"),
        "Rotate",
        FALSE
      ),
      "Click to rotate the activation matrix."
    ),
    shiny::selectInput(
      ns("selected_contrasts"),
      "Select comparisons:",
      choices = NULL,
      multiple = TRUE
    )
  )

  PlotModuleUI(ns("plot"),
    title = title,
    label = label,
    caption = caption,
    plotlib = "plotly",
    info.text = info.text,
    options = plot_opts,
    height = height,
    width = width,
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
functional_plot_go_actmap_server <- function(id,
                                             pgx,
                                             watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      shiny::observe({
        shiny::req(pgx$X)
        ct <- colnames(pgx$model.parameters$contr.matrix)
        ct <- sort(ct)
        selected_ct <- head(ct, 7)
        shiny::updateSelectInput(
          session,
          "selected_contrasts",
          choices = ct,
          selected = selected_ct
        )
      })


      plotGOactmap <- function(score, go, normalize, rotate, maxterm, maxfc,
                               tl.cex = 0.85, row.nchar = 60, colorbar = FALSE) {
        rownames(score) <- igraph::V(go)[rownames(score)]$Term

        ## avoid errors!!!
        score[is.na(score) | is.infinite(score)] <- 0
        score[is.na(score)] <- 0

        shiny::validate(
          shiny::need(
            !is.null(input$selected_contrasts),
            "Please select at least one comparison."
          )
        )
        score <- score[, input$selected_contrasts, drop = FALSE]

        ## reduce score matrix
        score <- score[head(order(-rowSums(score**2, na.rm = TRUE)), maxterm), , drop = FALSE] ## max number terms
        score <- score[, head(order(-colSums(score**2, na.rm = TRUE)), maxfc), drop = FALSE] ## max comparisons/FC
        score <- score + 1e-3 * matrix(rnorm(length(score)), nrow(score), ncol(score))

        ## normalize colums
        if (normalize) {
          ## column scale???
          score <- t(t(score) / (1e-8 + sqrt(colMeans(score**2, na.rm = TRUE))))
        }
        score <- score / max(abs(score), na.rm = TRUE) ## global normalize
        score <- sign(score) * abs(score)**0.5 ## fudging for better colors

        d1 <- as.dist(1 - cor(t(score), use = "pairwise"))
        d2 <- as.dist(1 - cor(score, use = "pairwise"))
        d1 <- dist(score)
        d2 <- dist(t(score))
        d1[is.na(d1)] <- 1
        d2[is.na(d2)] <- 1
        ii <- 1:nrow(score)
        jj <- 1:ncol(score)
        if (NCOL(score) == 1) {
          score <- score[order(-score[, 1]), 1, drop = FALSE]
        } else {
          ii <- hclust(d1)$order
          jj <- hclust(d2)$order
          score <- score[ii, jj, drop = FALSE]
        }

        colnames(score) <- substring(colnames(score), 1, 30)
        rownames(score) <- substring(rownames(score), 1, row.nchar)
        colnames(score) <- paste0(colnames(score), " ")

        par(mfrow = c(1, 1), mar = c(1, 1, 1, 1), oma = c(0, 1.5, 0, 0.5))

        if (rotate) score <- t(score)

        bluered.pal <- colorRamp(colors = c("royalblue3", "#ebeffa", "white", "#faeeee", "indianred3"))
        score <- score[nrow(score):1, ]
        x_axis <- colnames(score)
        y_axis <- rownames(score)
        fig <- plotly::plot_ly(
            x = x_axis, y = y_axis,
            z = score, type = "heatmap",
            colors = bluered.pal,
            showscale = colorbar
        )
        return(fig)
      }

      plot_data <- shiny::reactive({
        shiny::req(pgx$meta.go)
        pathscore <- pgx$meta.go$pathscore
        graph <- pgx$meta.go$graph
        res <- list(
          pathscore = pathscore,
          graph = graph
        )
      })

      plot_RENDER <- function() {
        res <- plot_data()
        shiny::req(res)
        pathscore <- res$pathscore
        graph <- res$graph

        fig <- plotGOactmap(
          score = pathscore,
          go = graph,
          normalize = input$normalize,
          rotate = input$rotate,
          maxterm = 50,
          maxfc = 25,
          tl.cex = 1.05,
          row.nchar = 60
        )
      }

      plot_RENDER2 <- function() {
        res <- plot_data()
        shiny::req(res)

        pathscore <- res$pathscore
        graph <- res$graph
        rotate <- input$rotate

        plotGOactmap(
          score = pathscore,
          go = graph,
          normalize = input$normalize,
          rotate = rotate,
          maxterm = 50,
          maxfc = 100,
          tl.cex = 1.1,
          row.nchar = ifelse(rotate, 60, 200),
          colorbar = TRUE
        )
      }

      PlotModuleServer(
        "plot",
        plotlib = "plotly",
        func = plot_RENDER,
        func2 = plot_RENDER2,
        csvFunc = plot_data,
        res = c(75, 105),
        pdf.width = 9,
        remove_margins = FALSE,
        pdf.height = 9,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
