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
    plotlib = "base",
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
                                                   plotActivationMatrix,
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
          selected = selected_ct)
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
        plotActivationMatrix(
          meta, df,
          normalize = input$normalize,
          rotate = rotate,
          nterms = 50,
          nfc = 20,
          tl.cex = 0.9,
          row.nchar = 60
        )
      }

      plot_RENDER2 <- function() {
        res <- plot_data()
        df <- res$df
        meta <- res$meta
        if (is.null(df) || nrow(df) == 0) {
          return(NULL)
        }
        rotate <- input$rotate
        plotActivationMatrix(
          meta, df,
          normalize = input$normalize,
          rotate = rotate,
          nterms = 50,
          nfc = 100,
          tl.cex = 1.1,
          row.nchar = ifelse(rotate, 60, 200)
        )
      }

      PlotModuleServer(
        "plot",
        plotlib = "base",
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
