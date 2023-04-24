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
functional_plot_wikipathway_actmap_ui <- function(
  id,
  title,
  caption,
  info.text,
  label = "",
  height
  ) {
  ns <- shiny::NS(id)
  
  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("normalize"),
        "normalize activation matrix",
        FALSE
      ),
      "Click to normalize the columns of the activation matrices."
    )
  )

  PlotModuleUI(
    id = ns("plot"),
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
functional_plot_wikipathway_actmap_server <- function(id,
                                               pgx,
                                               getWikiPathwayTable,
                                               plotActivationMatrix,
                                               watermark = FALSE) {
  moduleServer(
      id, function(input, output, session) {
          
      plot_data <- shiny::reactive({

        df <- getWikiPathwayTable()
        meta <- pgx$gset.meta$meta
        shiny::req(df,pgx$X,meta)
        
        res <- list(
          df = df,
          meta = meta
        )
      })

      plot_RENDER <- function() {
        res <- plot_data()
        df <- res$df
        meta <- res$meta

        if (is.null(df) || nrow(df) == 0) {
          return(NULL)
        }
        plotActivationMatrix(
            meta,
            df,
            normalize = input$normalize,
            nterms = 50,
            nfc = 20,
            tl.cex = 0.95,
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
        plotActivationMatrix(
            meta,
            df,
            normalize = input$normalize,
            nterms = 50,
            nfc = 100,
            tl.cex = 1.1,
            row.nchar = 200            
        )
      }

      PlotModuleServer(
        "plot",
        plotlib = "base",
        func = plot_RENDER,
        func2 = plot_RENDER2,
        csvFunc = plot_data,
        res = c(100,140),
        remove_margins = FALSE,
        pdf.height = 11,
        pdf.width = 6,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
