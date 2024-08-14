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
functional_plot_reactome_actmap_ui <- function(
    id,
    title,
    caption,
    info.text,
    label = "",
    height) {
  ns <- shiny::NS(id)

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("normalize"),
        "normalize activation matrix",
        FALSE
      ),
      "Click to normalize the columns of the activation matrices."
    ),
    withTooltip(
      shiny::checkboxInput(
        ns("rotate"),
        "Rotate activation matrix.",
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
    width = c("100%", "100%")
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
functional_plot_reactome_actmap_server <- function(id,
                                                   r_meta,
                                                   pgx,
                                                   getReactomeTable,
                                                   watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      shiny::observe({
        shiny::req(pgx$X)
        ct <- colnames(pgx$model.parameters$contr.matrix)
        ct <- sort(ct)
        selected_ct <- head(ct, 8)
        shiny::updateSelectInput(
          session,
          "selected_contrasts",
          choices = ct,
          selected = selected_ct
        )
      })
      plot_data <- shiny::reactive({
        df <- getReactomeTable()
        meta <- r_meta()
        shiny::req(df, meta)

        shiny::validate(
          shiny::need(
            !is.null(input$selected_contrasts),
            message = "Please select at least one comparison."
          )
        )
        meta <- meta[input$selected_contrasts]
        res <- list(
          df = df,
          meta = meta
        )
      })

      plot_RENDER <- function() {
        res <- plot_data()
        df <- res$df
        meta <- res$meta
        rotate <- input$rotate

        playbase::pgx.plotActivation(
          pgx,
          features = df$pathway,
          contrasts = input$selected_contrasts,
          what = "geneset",
          plotlib = "plotly",
          filter = NULL,
          cexCol = 1.4,
          cexRow = 1,
          normalize = input$normalize,
          rotate = input$rotate,
          maxterm = 30,
          maxfc = 20,
          mar = c(15, 30),
          tl.cex = 0.85,
          row.nchar = 50,
          showscale = TRUE
        )
      }

      plot_RENDER2 <- function() {
        res <- plot_data()
        df <- res$df
        meta <- res$meta
        if (is.null(df) || nrow(df) == 0) {
          return(NULL)
        }

        playbase::pgx.plotActivation(
          pgx,
          features = df$pathway,
          contrasts = input$selected_contrasts,
          what = "geneset",
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
          row.nchar = ifelse(input$rotate, 60, 200),
          showscale = TRUE
        )
      }

      PlotModuleServer(
        "plot",
        plotlib = "plotly",
        func = plot_RENDER,
        func2 = plot_RENDER2,
        csvFunc = plot_data,
        res = c(90, 100),
        pdf.height = 10,
        pdf.width = 10,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
