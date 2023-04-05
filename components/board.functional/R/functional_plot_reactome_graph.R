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
functional_plot_reactome_graph_ui <- function(id,
                                          label = "",
                                          height = c(400,700),
                                          width = c("auto","100%")
                                          ) {
  ns <- shiny::NS(id)
  info_text <- strwrap("<strong>REACTOME pathway.</strong> Genes are colored according
    to their upregulation (red) or downregulation (blue) in the contrast profile.
    Each pathway is scored for the selected contrast profile and reported in
    the table below.")

  PlotModuleUI(ns("plot"),
    title = "Reactome pathway map",
    label = label,
    plotlib = "image",
    info.text = info_text,
    info.width = "350px",
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
functional_plot_reactome_graph_server <- function(id,
                                              pgx,
                                              getFilteredReactomeTable,
                                              reactome_table,
                                              fa_contrast,
                                              watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {


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

      #      getPathwayImage <- function() {
      getPathwayImage <- shiny::reactive({

        dbg("[functional_plot_reactome_graph.R] getPathwayImage  reacted!")
        res <- plot_data()
        shiny::req(res, res$df)
        
        df <- res$df
        fa_contrast <- res$fa_contrast
        reactome_table <- res$reactome_table
        sbgn.dir <- res$sbgn.dir

        dbg("[functional_plot_reactome_graph.R] 2:")        
        ###############

        NULL.IMG <- list(src = "", contentType = "image/png")
        if (is.null(pgx)) {
          return(NULL.IMG)
        }

        dbg("[functional_plot_reactome_graph.R] 3:")
        
        comparison <- fa_contrast
        if (is.null(comparison) || length(comparison) == 0) {
          return(NULL.IMG)
        }
        if (comparison == "") {
          return(NULL.IMG)
        }

        dbg("[functional_plot_reactome_graph.R] 4:")
        
        ## get fold-change vector
        fc <- pgx$gx.meta$meta[[comparison]]$meta.fx
        pp <- rownames(pgx$gx.meta$meta[[comparison]])

        dbg("[functional_plot_reactome_graph.R] 5:")
        
        if ("hgnc_symbol" %in% colnames(pgx$genes)) {
          names(fc) <- pgx$genes[pp, "hgnc_symbol"]
        } else {
          names(fc) <- toupper(pgx$genes[pp, "gene_name"])
        }
        fc <- fc[order(-abs(fc))]
        fc <- fc[which(!duplicated(names(fc)) & names(fc) != "")]

        dbg("[functional_plot_reactome_graph.R] 6:")
        
        ## get selected REACTOME id
        # df <- getFilteredReactomeTable()
        if (is.null(df)) {
          return(NULL.IMG)
        }

        dbg("[functional_plot_reactome_graph.R] 6b:")
        
        sel.row <- reactome_table$rows_selected()
        if (is.null(sel.row) || length(sel.row) == 0) {
          return(NULL.IMG)
        }
        sel.row <- as.integer(sel.row)
        
        dbg("[functional_plot_reactome_graph.R] 7:")
        
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

        dbg("[functional_plot_reactome_graph.R] 8:")
        
        ## We temporarily switch the working directory to always readable
        ## TMP folder
        curwd <- getwd()
        tmpdir <- tempdir()
        setwd(tmpdir)

        ## pv.out <- pathview::pathview(
        ##   gene.data = fc, pathway.id = pathway.id, gene.idtype = "SYMBOL",
        ##   gene.annotpkg = "org.Hs.eg.db", species = "hsa",
        ##   out.suffix = "pathview", limit = list(gene = 2, cpd = 1),
        ##   low = list(gene = "dodgerblue2", cpd = "purple"),
        ##   high = list(gene = "firebrick2", cpd = "yellow"),
        ##   reactome.dir = xml.dir, reactome.native = TRUE, same.layer = FALSE,
        ## )
        dbg("[functional_plot_reactome_graph.R] tmpdir = ",tmpdir)
        dbg("[functional_plot_reactome_graph.R] pathway.id = ",pathway.id)
        dbg("[functional_plot_reactome_graph.R] sbgn.dir = ",sbgn.dir)
        dbg("[functional_plot_reactome_graph.R] head.fc = ",head(fc))
        dbg("[functional_plot_reactome_graph.R] head.names.fc = ",head(names(fc)))

        require(SBGNview)
        ##data("mapped.ids","pathways.info", "sbgn.xmls")
        ##data("sbgn.xmls", package="SBGNview.data",verbose=1)  ### THIS IS BIG 700MB!!! 
        ##object.size(sbgn.xmls)
       
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
        file.exists(imgfile)
        dbg("[functional_plot_reactome_graph.R] img.file = ",imgfile)
        dbg("[functional_plot_reactome_graph.R] file.exists = ",file.exists(imgfile))
        
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
          contentType = "image/png",
          #width = "100%", height = "100%", ## actual size: 1040x800
          width = img.dim[2], height = img.dim[1], ## actual size
          ##width = img.width, height = img.height, ## actual size
          alt = "reactome pathway"
        )
      })

      calcImageSize <- function(img.height, img.width, client.height, client.width ) {

        client.aspectratio <- client.width / client.height
        img.aspectratio <- img.width / img.height
        dbg("[functional_plot_reactome_graph.R] img.height=", img.height)
        dbg("[functional_plot_reactome_graph.R] img.width=", img.width)        
        
        new.width=new.height=400
        if(client.aspectratio > img.aspectratio) {
          new.height <- "100%"
          new.width <- (client.height / img.height) * img.width
        } else {
          new.width <- "100%"
          new.height <- (client.width / img.width) * img.height
        }
        dbg("[functional_plot_reactome_graph.R] img.width=", new.width)
        dbg("[functional_plot_reactome_graph.R] img.height=", new.height)
        
        c(height = new.height, width = new.width)
      }
      
      plot_RENDER <- function() {

        width=height=NULL
        client.width  <- session$clientData$"output_pathway-reactome_graph-plot-renderfigure_width"
        client.height <- session$clientData$"output_pathway-reactome_graph-plot-renderfigure_height"
        client.pixelratio <- session$clientData$pixelratio
        client.pixelratio <- 1
        dbg("[functional_plot_reactome_graph.R] client.width=", client.width)
        dbg("[functional_plot_reactome_graph.R] client.height=", client.height)
        dbg("[functional_plot_reactome_graph.R] client.pixelratio=", client.pixelratio)        

        img <- getPathwayImage()
        shiny::req(img$width, img$height)

        if(0) {
          res <- calcImageSize(img$height, img$width, client.height, client.width)        
          dbg("[functional_plot_reactome_graph.R] res.width=", res["width"])
          dbg("[functional_plot_reactome_graph.R] res.height=", res["height"])
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
