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
functional_plot_kegg_graph_ui <- function(
  id,
  title,
  info.text,
  caption,
  label = "",
  height,
  width,
  info.width
  ) {
  ns <- shiny::NS(id)
  
  PlotModuleUI(ns("plot"),
    title = title,
    label = label,
    plotlib = "image",
    info.text = info.text,
    caption = caption,
    info.width,
    options = NULL,
    download.fmt = c("png","csv"),
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
functional_plot_kegg_graph_server <- function(id,
                                              pgx,
                                              getFilteredKeggTable,
                                              kegg_table,
                                              fa_contrast,
                                              watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {


      ## plot_data <- shiny::reactive({
      plot_data <- function() {
        ## folder with predownloaded XML files
        xml.dir <- file.path(FILES, "kegg-xml")
        xml.dir <- normalizePath(xml.dir) ## absolute path
        res <- list(
          df = getFilteredKeggTable(),
          kegg_table = kegg_table,
          fa_contrast = fa_contrast(),
          xml.dir = xml.dir
        )
        return(res)
      }#)

      #      getPathwayImage <- function() {
      getPathwayImage <- shiny::reactive({
        res <- plot_data()
        shiny::req(res, res$df)
        
        df <- res$df
        fa_contrast <- res$fa_contrast
        kegg_table <- res$kegg_table
        xml.dir <- res$xml.dir

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

        ## get selected KEGG id
        # df <- getFilteredKeggTable()
        if (is.null(df)) {
          return(NULL.IMG)
        }

        sel.row <- kegg_table$rows_selected()
        if (is.null(sel.row) || length(sel.row) == 0) {
          return(NULL.IMG)
        }
        sel.row <- as.integer(sel.row)

        pathway.id <- "04110" ## CELL CYCLE
        pathway.name <- pw.genes <- "x"
        if (is.null(sel.row) || length(sel.row) == 0) {
          return(NULL.IMG)
        }

        if (!is.null(sel.row) && length(sel.row) > 0) {
          pathway.id <- df[sel.row, "kegg.id"]
          pathway.name <- df[sel.row, "pathway"]
          pw.genes <- unlist(getGSETS(as.character(pathway.name)))
        }
        
        ## We temporarily switch the working directory to always readable
        ## TMP folder
        curwd <- getwd()
        tmpdir <- tempdir()
        setwd(tmpdir)

        pv.out <- pathview::pathview(
          gene.data = fc, pathway.id = pathway.id, gene.idtype = "SYMBOL",
          gene.annotpkg = "org.Hs.eg.db", species = "hsa",
          out.suffix = "pathview", limit = list(gene = 2, cpd = 1),
          low = list(gene = "dodgerblue2", cpd = "purple"),
          high = list(gene = "firebrick2", cpd = "yellow"),
          kegg.dir = xml.dir, kegg.native = TRUE, same.layer = FALSE,
        )
        Sys.sleep(0.2) ## wait for graph

        ## back to previous working folder
        setwd(curwd)

        imgfile="/tmp/hsa00010.png"
        imgfile <- file.path(tmpdir, paste0("hsa", pathway.id, ".pathview.png"))
        if (!file.exists(imgfile)) {
          return(NULL.IMG)
        }

        img.dim <- NULL
        if(grepl("png|PNG",imgfile)) img.dim <- dim(png::readPNG(imgfile))[1:2]
        if(grepl("jpg|JPG",imgfile)) img.dim <- dim(jpeg::readJPEG(imgfile))[1:2]
        img.dim

        list(
          src = imgfile,
          contentType = "image/png",
          #width = "100%", height = "100%", ## actual size: 1040x800
          width = img.dim[2], height = img.dim[1], ## actual size
          ##width = img.width, height = img.height, ## actual size
          alt = "pathview image"
        )
      })

      calcImageSize <- function(img.height, img.width, client.height, client.width ) {

        client.aspectratio <- client.width / client.height
        img.aspectratio <- img.width / img.height
        dbg("[functional_plot_kegg_graph.R] img.height=", img.height)
        dbg("[functional_plot_kegg_graph.R] img.width=", img.width)        
        
        new.width=new.height=400
        if(client.aspectratio > img.aspectratio) {
          new.height <- "100%"
          new.width <- (client.height / img.height) * img.width
        } else {
          new.width <- "100%"
          new.height <- (client.width / img.width) * img.height
        }
        dbg("[functional_plot_kegg_graph.R] img.width=", new.width)
        dbg("[functional_plot_kegg_graph.R] img.height=", new.height)
        
        c(height = new.height, width = new.width)
      }
      
      plot_RENDER <- function() {

        width=height=NULL
        client.width  <- session$clientData$"output_pathway-kegg_graph-plot-renderfigure_width"
        client.height <- session$clientData$"output_pathway-kegg_graph-plot-renderfigure_height"
        client.pixelratio <- session$clientData$pixelratio
        client.pixelratio <- 1
        dbg("[functional_plot_kegg_graph.R] client.width=", client.width)
        dbg("[functional_plot_kegg_graph.R] client.height=", client.height)
        dbg("[functional_plot_kegg_graph.R] client.pixelratio=", client.pixelratio)        

        img <- getPathwayImage()
        shiny::req(img$width, img$height)

        if(0) {
          res <- calcImageSize(img$height, img$width, client.height, client.width)        
          dbg("[functional_plot_kegg_graph.R] res.width=", res["width"])
          dbg("[functional_plot_kegg_graph.R] res.height=", res["height"])
          img$width  <- res["width"]
          img$height <- res["height"]
        } else {
          img$width  <- "100%"
          img$height <- "100%"
        }
        
        return(img)
      }

      PlotModuleServer(
        "plot",
        plotlib = "image",
        func = plot_RENDER,
        func2 = plot_RENDER,
        csvFunc = plot_data,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
