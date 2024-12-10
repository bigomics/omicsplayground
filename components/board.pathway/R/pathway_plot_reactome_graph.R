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
        ## folder with predownloaded SBGN files
        sbgn.dir <- pgx.system.file("sbgn/", package = "pathway")
        sbgn.dir <- normalizePath(sbgn.dir) ## absolute path
        res <- list(
          df = getFilteredReactomeTable(),
          reactome_table = reactome_table,
          fa_contrast = fa_contrast(),
          sbgn.dir = sbgn.dir
        )
        return(res)
      }

      getPathwayImage <- shiny::reactive({
        res <- plot_data()
        shiny::req(res, res$df)

        df <- res$df
        comparison <- res$fa_contrast
        reactome_table <- res$reactome_table
        sbgn.dir <- res$sbgn.dir

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

        ## We temporarily switch the working directory to always readable
        ## TMP folder
        curwd <- getwd()
        tmpdir <- tempdir()
        setwd(tmpdir)

        suppressMessages(require(SBGNview)) ## slow!! but needed!!!


        ## this is a trick. the original object in SBGNview.data was 700MB!!
        sbgn.xmls <- dir(sbgn.dir, ".sbgn")
        names(sbgn.xmls) <- sbgn.xmls

        if (!interactive()) {
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          progress$set(message = "Rendering pathway...", value = 0.33)
        }

        obj <- tryCatch(
          {
            SBGNview::SBGNview(
              gene.data = fc,
              gene.id.type = "SYMBOL",
              sbgn.dir = sbgn.dir,
              input.sbgn = pathway.id,
              output.file = "reactome",
              output.formats = c("png")
            )
          },
          error = function(w) {
            SBGNview::SBGNview(
              gene.data = NULL,
              gene.id.type = "SYMBOL",
              sbgn.dir = sbgn.dir,
              input.sbgn = pathway.id,
              output.file = "reactome",
              output.formats = c("png")
            )
          }
        )
        if (class(obj) == "SBGNview") {
          try(print(obj))
        }
        Sys.sleep(0.2) ## wait for graph

        ## back to previous working folder
        setwd(curwd)

        imgfile <- "/tmp/hsa00010.png"
        imgfile <- file.path(tmpdir, paste0("reactome_", pathway.id, ".png"))
        svgfile <- file.path(tmpdir, paste0("reactome_", pathway.id, ".svg"))
        file.exists(imgfile)

        if (!file.exists(imgfile)) {
          return(NULL.IMG)
        }

        ## parse image dimensions from file
        img.dim <- NULL
        if (grepl("png|PNG", imgfile)) img.dim <- dim(png::readPNG(imgfile))[1:2]
        if (grepl("jpg|JPG", imgfile)) img.dim <- dim(jpeg::readJPEG(imgfile))[1:2]
        img.dim

        list(
          src = imgfile,
          svg = svgfile,
          contentType = "image/png",
          width = img.dim[2], height = img.dim[1], ## actual size
          alt = "reactome pathway"
        )
      })

      plot_RENDER <- function() {
        img <- getPathwayImage()
        shiny::req(img$width, img$height)
        filename <- img$svg
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
