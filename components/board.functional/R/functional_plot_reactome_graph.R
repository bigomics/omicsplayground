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
functional_plot_reactome_graph_server <- function(id,
                                              pgx,
                                              getFilteredReactomeTable,
                                              reactome_table,
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
        sbgn.dir <- file.path(FILES, "reactome-sbgn")
        sbgn.dir <- normalizePath(sbgn.dir) ## absolute path
        res <- list(
          df = getFilteredReactomeTable(),
          reactome_table = reactome_table,
          fa_contrast = fa_contrast(),
          sbgn.dir = sbgn.dir
        )
        return(res)
      }#)

      #  getPathwayImage <- function() {
      getPathwayImage <- shiny::reactive({

        res <- plot_data()
        shiny::req(res, res$df)
        
        df <- res$df
        fa_contrast <- res$fa_contrast
        reactome_table <- res$reactome_table
        sbgn.dir <- res$sbgn.dir

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

        sel.row <- reactome_table$rows_selected()
        if (is.null(sel.row) || length(sel.row) == 0) {
          return(NULL.IMG)
        }
        sel.row <- as.integer(sel.row)
        
        pathway.id <- "R-HSA-109704"
        pathway.name <- pw.genes <- "x"
        if (is.null(sel.row) || length(sel.row) == 0) {
          return(NULL.IMG)
        }

        if (!is.null(sel.row) && length(sel.row) > 0) {
          pathway.id   <- df[sel.row, "reactome.id"]
          pathway.name <- df[sel.row, "pathway"]
          pw.genes <- unlist(getGSETS(as.character(pathway.name)))
        }

        ## We temporarily switch the working directory to always readable
        ## TMP folder
        curwd <- getwd()
        tmpdir <- tempdir()
        setwd(tmpdir)

        ##require(SBGNview)
        ##data("mapped.ids","pathways.info", "sbgn.xmls")
        ##data("sbgn.xmls", package="SBGNview.data",verbose=1)  ### BIG 700MB!!! 
        ##object.size(sbgn.xmls)

        ## this is a trick. the original object in SBGNview.data was 700MB!!
        sbgn.xmls <- dir(sbgn.dir,".sbgn")
        names(sbgn.xmls) <- sbgn.xmls
        
        obj <- try( SBGNview::SBGNview(
          gene.data = fc, 
          ##gene.id.type = "ENTREZID",
          gene.id.type = "SYMBOL",
          ##sbgn.gene.id.type = "SYMBOL",
          sbgn.dir = sbgn.dir,
          # input.sbgn = pathway.id,
          input.sbgn = pathway.id,
          output.file = "reactome", 
          output.formats =  c("png")
        ))
        if(class(obj) == "SBGNview") {
          try(print(obj))
        }
        Sys.sleep(0.2) ## wait for graph
        
        ## back to previous working folder
        setwd(curwd)

        imgfile="/tmp/hsa00010.png"
        imgfile <- file.path(tmpdir, paste0("reactome_", pathway.id, ".png"))
        svgfile <- file.path(tmpdir, paste0("reactome_", pathway.id, ".svg"))
        file.exists(imgfile)
        
        if (!file.exists(imgfile)) {
          return(NULL.IMG)
        }

        ## parse image dimensions from file
        img.dim <- NULL
        if(grepl("png|PNG",imgfile)) img.dim <- dim(png::readPNG(imgfile))[1:2]
        if(grepl("jpg|JPG",imgfile)) img.dim <- dim(jpeg::readJPEG(imgfile))[1:2]
        img.dim

        list(
          src = imgfile,
          svg = svgfile,
          contentType = "image/png",
          #width = "100%", height = "100%", ## actual size: 1040x800
          width = img.dim[2], height = img.dim[1], ## actual size
          ##width = img.width, height = img.height, ## actual size
          alt = "reactome pathway"
        )
      })
      
      plot_RENDER <- function() {

        img <- getPathwayImage()
        shiny::req(img$width, img$height)
        filename <- img$svg
        img.svg <-  readChar(filename, nchars = file.info(filename)$size)
        pz <- svgPanZoom::svgPanZoom(
          img.svg,
          controlIconsEnabled = TRUE,
          zoomScaleSensitivity = 0.4,
          maxZoom = 25
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
