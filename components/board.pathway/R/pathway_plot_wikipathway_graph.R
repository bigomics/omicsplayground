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
functional_plot_wikipathway_graph_ui <- function(
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
functional_plot_wikipathway_graph_server <- function(id,
                                                     pgx,
                                                     getFilteredWikiPathwayTable,
                                                     wikipathway_table,
                                                     fa_contrast,
                                                     watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      ## reactive or function? that's the question...
      plot_data <- shiny::reactive({
        res <- list(
          df = getFilteredWikiPathwayTable(),
          wikipathway_table = wikipathway_table,
          fa_contrast = fa_contrast()
        )
        return(res)
      })

      getPathwayImage <- shiny::reactive({
        res <- plot_data()
        shiny::req(res, res$df)

        df <- res$df
        comparison <- res$fa_contrast
        wikipathway_table <- res$wikipathway_table

        ###############

        NULL.IMG <- list(src = "", contentType = "image/png")
        if (is.null(pgx)) {
          return(NULL.IMG)
        }
        if (is.null(comparison) || length(comparison) == 0) {
          return(NULL.IMG)
        }
        if (comparison == "") {
          return(NULL.IMG)
        }

        if (is.null(df)) {
          return(NULL.IMG)
        }

        ## get fold-change vector
        fc <- pgx$gx.meta$meta[[comparison]]$meta.fx
        pp <- rownames(pgx$gx.meta$meta[[comparison]])

        # Rename to human orthologs for non-human species and sort

        if (pgx$organism != "Human") {
          names(fc) <- pgx$genes[pp, "human_ortholog"]
        } else {
          names(fc) <- pgx$genes[pp, "symbol"]
        }
        fc <- fc[order(-abs(fc))]
        fc <- fc[which(!duplicated(names(fc)) & names(fc) != "")]

        # Get selected row and path id
        sel.row <- wikipathway_table$rows_selected()
        if (is.null(sel.row) || length(sel.row) == 0) {
          return(NULL.IMG)
        }
        sel.row <- as.integer(sel.row)
        pathway.name <- pw.genes <- "x"
        if (is.null(sel.row) || length(sel.row) == 0) {
          return(NULL.IMG)
        }
        if (!is.null(sel.row) && length(sel.row) > 0) {
          pathway.id <- df[sel.row, "pathway.id"]
          pathway.name <- df[sel.row, "pathway"]
          pw.genes <- unlist(playdata::getGSETS(as.character(pathway.name)))
        }

        if (!interactive()) {
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          progress$set(message = "Rendering pathway...", value = 0.33)
        }

        # Get pathway image using WPID and fc values
        shiny::req(pathway.id)
        svg <- playbase::wikipathview(wp = pathway.id, val = fc)
        if (is.null(svg)) {
          return(NULL)
        }
        list(
          src = normalizePath(svg),
          contentType = "image/svg+xml",
          width = "100%", height = "100%", ## actual size: 1040x800
          alt = "wikipathway SVG"
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
