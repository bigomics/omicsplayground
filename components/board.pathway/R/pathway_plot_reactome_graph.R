##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics Sagl. All rights reserved.
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
  width
) {
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
    download.fmt = c("png", "pdf", "svg"),
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

        ## Fetch the native Reactome diagram SVG from our mirror. reactome.org's
        ## live exporter is Cloudflare-blocked for server-side fetches; the table
        ## above is filtered to pathways that resolve to a diagram, so this
        ## returns an image for every selectable row (see playbase::getReactomeSVG).
        imgfile <- playbase::getPathwayImage(
          pathway.id,
          val = fc, as.img = TRUE
        )

        return(imgfile)
      })

      plot_RENDER <- function() {
        img <- get_image()

        ## nothing selected yet -> keep the panel blank (get_image returns src = "")
        if (!is.null(img) && !is.null(img$src) && !nzchar(img$src)) {
          shiny::req(FALSE)
        }

        ## selected, but no local diagram could be rendered -> note instead of a crash
        filename <- img$src
        if (is.null(filename) || !file.exists(filename)) {
          note <- paste0(
            '<svg xmlns="http://www.w3.org/2000/svg" width="640" height="110">',
            '<text x="15" y="45" font-family="sans-serif" font-size="15">',
            "<tspan x=\"15\">Reactome diagram is not available for this pathway.</tspan>",
            "<tspan x=\"15\" dy=\"28\">Use the link in the table to open it on reactome.org.</tspan>",
            "</text></svg>"
          )
          return(svgPanZoom::svgPanZoom(note, controlIconsEnabled = FALSE, viewBox = FALSE))
        }

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
