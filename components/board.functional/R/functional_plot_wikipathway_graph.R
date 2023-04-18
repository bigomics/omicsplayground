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
    ##    plotlib = "image",
    plotlib = "generic",
    plotlib2 = "generic",
    outputFunc = svgPanZoom::svgPanZoomOutput,
    outputFunc2 = svgPanZoom::svgPanZoomOutput,
    info.text = info.text,
    info.width = info.width,
    options = NULL,
    download.fmt = NULL,
    # download.fmt = c("png","csv"),    
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

      ## preload... takes few seconds... 
      suppressMessages(require(SBGNview)) 

      ## reactive or function? that's the question...
      ## plot_data <- shiny::reactive({
      plot_data <- function() {
        ## folder with predownloaded SBGN files
        svg.dir <- file.path(FILES, "wikipathway-svg")
        svg.dir <- normalizePath(svg.dir) ## absolute path
        res <- list(
          df = getFilteredWikiPathwayTable(),
          wikipathway_table = wikipathway_table,
          fa_contrast = fa_contrast(),
          svg.dir = svg.dir
        )
        return(res)
      }#)

      #  getPathwayImage <- function() {
      getPathwayImage <- shiny::reactive({

        res <- plot_data()
        shiny::req(res, res$df)
        
        df <- res$df
        fa_contrast <- res$fa_contrast
        wikipathway_table <- res$wikipathway_table
        svg.dir <- res$svg.dir

        ###############

        NULL.IMG <- list(src = "", contentType = "image/png")
        if (is.null(pgx)) {
          return(NULL.IMG)
        }

        comparison <- fa_contrast
        if (is.null(comparison) || length(comparison) == 0) {
          return(NULL.IMG)
        }
        if (comparison == "") {
          return(NULL.IMG)
        }

        ## get fold-change vector
        fc <- pgx$gx.meta$meta[[comparison]]$meta.fx
        pp <- rownames(pgx$gx.meta$meta[[comparison]])

        if ("hgnc_symbol" %in% colnames(pgx$genes)) {
          names(fc) <- pgx$genes[pp, "hgnc_symbol"]
        } else {
          names(fc) <- toupper(pgx$genes[pp, "gene_name"])
        }
        fc <- fc[order(-abs(fc))]
        fc <- fc[which(!duplicated(names(fc)) & names(fc) != "")]
                
        if (is.null(df)) {
          return(NULL.IMG)
        }

        sel.row <- wikipathway_table$rows_selected()
        if (is.null(sel.row) || length(sel.row) == 0) {
          return(NULL.IMG)
        }
        sel.row <- as.integer(sel.row)
        
        pathway.id <- "WP179"
        pathway.name <- pw.genes <- "x"
        if (is.null(sel.row) || length(sel.row) == 0) {
          return(NULL.IMG)
        }
        dbg("[functional_plot_wikipathway_graph.R] sel.row = ",sel.row)
        
        if (!is.null(sel.row) && length(sel.row) > 0) {
          pathway.id   <- df[sel.row, "pathway.id"]
          pathway.name <- df[sel.row, "pathway"]
          pw.genes <- unlist(getGSETS(as.character(pathway.name)))
        }

        dbg("[functional_plot_wikipathway_graph.R] rendering WikiPathway",pathway.id)

        tmpfile <- paste0(tempfile(),".svg")
        svg <- try( wikipathview(wp=pathway.id, val=fc, dir=svg.dir) )
        fluctuator::write_svg(svg, file = tmpfile)

        list(
          src = normalizePath(tmpfile),
          contentType = "image/svg+xml",
          width = "100%", height = "100%", ## actual size: 1040x800
          alt = "wikipathway SVG"
        )
      })
      
      plot_RENDER <- function() {
        img <- getPathwayImage()
        shiny::req(img$width, img$height)
        filename <- img$src
        img.svg <-  readChar(filename, nchars = file.info(filename)$size)
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
        plotlib = "generic",
        func = plot_RENDER,
        func2 = plot_RENDER,
        renderFunc = svgPanZoom::renderSvgPanZoom,
        renderFunc2 = svgPanZoom::renderSvgPanZoom,
        ## csvFunc = plot_data,
        add.watermark = watermark
      )

    } ## end of moduleServer
  )
}
