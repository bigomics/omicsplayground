##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
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
mofa_plot_pathbank_ui <- function(
    id,
    label = "",
    title = "",
    caption = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    id = ns("plotmodule"),
    title = title,
    label = label,
    plotlib = "svgPanZoom",
    options = NULL,
    download.fmt = c("png", "pdf"),
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
mofa_plot_pathbank_server <- function(id,
                                      pgx,
                                      ## pathbank_table
                                      ##fa_contrast,
                                      watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {

      ## reactive or function? that's the question...
      plot_data <- shiny::reactive({
        res <- list()
        return(res)
      })
      
      getPathwayImage <- shiny::reactive({

        # Get pathway image using WPID and fc values
        ## svg <- wikipathview(wp = "WP2446", val = fc)
        ##source("function_getpathbank.R")
        wp = "SMP0080852"        
        ##wp.features <- gmt[grep(wp, names(gmt))]

        svg <- get_pathbank_svg(wp = wp, val=NULL) 
        if (is.null(svg)) {
          return(NULL)
        }
        list(
          src = normalizePath(svg),
          contentType = "image/svg+xml",
          width = "100%", height = "100%", ## actual size: 1040x800
          alt = "pathway SVG"
        )
      })

      plot_RENDER <- function() {
        img <- getPathwayImage()
        validate(
          need(!is.null(img), "Could not retrieve pathway image")
        )
        shiny::req(img$width, img$height)
        filename <- img$src
        img.svg <- readChar(filename, nchars = file.info(filename)$size)
        pz <- svgPanZoom::svgPanZoom(
          img.svg,
          controlIconsEnabled = TRUE,
          zoomScaleSensitivity = 0.4,
          minZoom = 1,
          maxZoom = 5,
          viewBox = FALSE
        )
        return(pz)
      }

      PlotModuleServer(
        "plotmodule",
        plotlib = "svgPanZoom",
        func = plot_RENDER,
      )
    } ## end of moduleServer
  )
}
