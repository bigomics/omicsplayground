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
        ct <- playbase::pgx.getContrasts(pgx)
        ct <- sort(ct[!grepl("^IA:", ct)])
        selected_ct <- head(ct, 8)
        shiny::updateSelectInput(
          session,
          "selected_contrasts",
          choices = ct,
          selected = selected_ct
        )
      })

      plot_data <- shiny::reactive({
        shiny::req(pgx$meta.go)
        score <- pgx$meta.go$pathscore
        score <- score[, input$selected_contrasts, drop = FALSE]
        graph <- pgx$meta.go$graph
        rownames(score) <- igraph::V(graph)[rownames(score)]$Term
        res <- list(score = score)
      })

      plot_RENDER <- function() {
        res <- plot_data()
        shiny::req(res)

        playbase::pgx.plotActivation(
          pgx,
          contrasts = input$selected_contrasts,
          what = "matrix",
          matrix = res$score,
          plotlib = "plotly",
          filter = NULL,
          cexCol = 1.4,
          cexRow = 1,
          normalize = input$normalize,
          rotate = input$rotate,
          maxterm = 30,
          maxfc = 20,
          mar = c(15, 30),
          tl.cex = 1.05,
          row.nchar = 60
        )
      }

      plot_RENDER2 <- function() {
        res <- plot_data()
        shiny::req(res)

        playbase::pgx.plotActivation(
          pgx,
          contrasts = input$selected_contrasts,
          what = "matrix",
          matrix = res$score,
          plotlib = "plotly",
          filter = NULL,
          cexCol = 1.4,
          cexRow = 1,
          normalize = input$normalize,
          rotate = input$rotate,
          maxterm = 40,
          maxfc = 100,
          mar = c(15, 30),
          tl.cex = 1.1,
          row.nchar = ifelse(input$rotate, 60, 200)
        )

        ## plotGOactmap(
        ##   score = pathscore,
        ##   go = graph,
        ##   normalize = input$normalize,
        ##   rotate = rotate,
        ##   maxterm = 50,
        ##   maxfc = 100,
        ##   tl.cex = 1.1,
        ##   row.nchar = ifelse(rotate, 60, 200),
        ##   colorbar = TRUE
        ## )
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
