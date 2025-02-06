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
functional_plot_reactome_graph_ui <- function(
    id,
    label = "",
    title,
    info.text,
    info.methods,
    info.references,
    info.extra_link,
    caption,
    info.width,
    height,
    width) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    id = ns("plotmodule"),
    title = title,
    label = label,
    plotlib = "svgPanZoom",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    info.width = info.width,
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
functional_plot_reactome_graph_server <- function(id,
                                                  pgx,
                                                  getFilteredReactomeTable,
                                                  reactome_table,
                                                  fa_contrast,
                                                  watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      ## preload... takes few seconds...
      #

      ## reactive or function? that's the question...
      plot_data <- function() {
        res <- list(
          df = getFilteredReactomeTable(),
          reactome_table = reactome_table,
          fa_contrast = fa_contrast()
        )
        return(res)
      }
      
      get_image <- shiny::reactive({
        res <- plot_data()
        shiny::req(res, res$df)

        df <- res$df
        comparison <- res$fa_contrast
        reactome_table <- res$reactome_table

        ###############
        shiny::req(df, comparison, pgx$X)

        NULL.IMG <- list(src = "", contentType = "image/png")

        ## get fold-change vector
        fc <- pgx$gx.meta$meta[[comparison]]$meta.fx
        pp <- rownames(pgx$gx.meta$meta[[comparison]])

        if (pgx$organism != "Human") {
          names(fc) <- pgx$genes[pp, "human_ortholog"]
        } else {
          names(fc) <- pgx$genes[pp, "symbol"]
        }
        fc <- fc[order(-abs(fc))]
        fc <- fc[which(!duplicated(names(fc)) & names(fc) != "")]

        sel.row <- reactome_table$rows_selected()
        if (is.null(sel.row) || length(sel.row) == 0 || nrow(df) == 0) {
          return(NULL.IMG)
        }

        sel.row <- as.integer(sel.row)
        pathway.id <- df[sel.row, "reactome.id"]
        pathway.name <- df[sel.row, "pathway"]
        pw.genes <- unlist(playdata::getGSETS(as.character(pathway.name)))
        
        ## folder with predownloaded SBGN files
        sbgn.dir <- pgx.system.file("sbgn/", package = "pathway")
        sbgn.dir <- normalizePath(sbgn.dir) ## absolute path
        imgfile <- playbase::getReactomeSVG(
          pathway.id, val=fc, sbgn.dir=sbgn.dir, as.img=TRUE)

        return(imgfile)
      })

      plot_RENDER <- function() {
        img <- get_image()
        shiny::req(img$src)
        filename <- img$src
        img.svg <- readChar(filename, nchars = file.info(filename)$size)
        pz <- svgPanZoom::svgPanZoom(
          img.svg,
          controlIconsEnabled = TRUE,
          zoomScaleSensitivity = 0.4,
          minZoom = 1,
          maxZoom = 25,
          viewBox = FALSE
        )
        return(pz)
      }

      PlotModuleServer(
        "plotmodule",
        func = plot_RENDER,
        plotlib = "svgPanZoom",
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
