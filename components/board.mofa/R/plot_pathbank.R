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
                                      sel_pathway = reactive(NULL),
                                      sel_contrast = reactive(NULL),
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
        wp.name <- sel_pathway()
        dbg("[mofa_plot_pathbank_server] selected wp.name = ", wp.name )
        if(length(wp.name)==0 || length(wp.name) > 1) return(NULL)
        
        ##wp.name = "TF_ARCHS4:NCOA3 human tf ARCHS4 coexpression [SMP0080852]"
        wp <- stringr::str_match(wp.name, "SMP[0-9]+|WP[0-9]+|R-HSA-[0-9]+")[,1]
        shiny::validate( shiny::need(!is.na(wp), "pathway diagram not available"))
        if(length(wp)>1) wp <- wp[1]
        
        k <- sel_contrast()
        val <- playbase::pgx.getMetaMatrix(pgx)$fc[,k]

        ## convert to UNIPROT and PATHBANK ID
        ##newnames <- convert2pathbankid(names(val))
        ##names(val) <- newnames

        val = NULL ## temporary...
        
        sbgn.dir <- pgx.system.file("sbgn/", package = "pathway")
        sbgn.dir <- normalizePath(sbgn.dir) ## absolute path
        ##wp = "SMP0080852"        
        img <- playbase::getPathwayImage(
          wp, val=val, sbgn.dir=sbgn.dir, as.img=TRUE)
        return(img)
      })

      plot_RENDER <- function() {
        img <- getPathwayImage()
        validate(
          need(!is.null(img), "Could not retrieve pathway image")
        )
        shiny::req(img$src)
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
