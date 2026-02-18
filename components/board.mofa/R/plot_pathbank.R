##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' Importance plot UI input function
#' @description A shiny Module for plotting (UI code).
#' @param id
#' @param label
#' @param height
#' @export
mofa_plot_pathbank_ui <- function(
  id,
  label = "",
  title = "",
  info.text = "",
  caption = "",
  height = 400,
  width = 400
) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    id = ns("plotmodule"),
    title = title,
    label = label,
    info.text = info.text,
    plotlib = "svgPanZoom",
    options = NULL,
    download.fmt = c("png", "pdf", "svg"),
    height = height,
    width = width
  )
}


#' Importance plot Server function
#' @description A shiny Module for plotting (server code).
#' @param id
#' @return
#' @export
mofa_plot_pathbank_server <- function(id,
                                      pgx,
                                      sel_pathway = reactive(NULL),
                                      sel_contrast = reactive(NULL),
                                      watermark = FALSE) {

  moduleServer(

    id, function(input, output, session) {

      getPathwayImage <- shiny::reactive({

        wp.name <- sel_pathway()
        if (length(wp.name) == 0 || length(wp.name) > 1) {
          return(NULL)
        }

        wp <- stringr::str_match(wp.name, "SMP[0-9]+|WP[0-9]+|R-HSA-[0-9]+")[, 1]
        shiny::validate(shiny::need(!is.na(wp), "pathway diagram not available"))
        if (length(wp) > 1) wp <- wp[1]

        val <- NULL
        k <- sel_contrast()
        if (!is.null(k)) {
          val <- playbase::pgx.getMetaMatrix(pgx)$fc[, k]
          val_names <- names(val)
          ids <- sapply(strsplit(val_names, ":"), function(x) x[2])
          id_mapping <- playdata::METABOLITE_ANNOTATION[, c("ID", "PATHBANK")]
          id_mapping <- id_mapping[!is.na(id_mapping$PATHBANK), ]
          id_mapping <- setNames(id_mapping$PATHBANK, id_mapping$ID)
          mapped_ids <- id_mapping[ids]
          mapped_ids <- mapped_ids[!is.na(mapped_ids)]
          names(val)[match(names(mapped_ids), ids)] <- mapped_ids
        }

        ## convert to UNIPROT and PATHBANK ID
        ## newnames <- convert2pathbankid(names(val))
        ## names(val) <- newnames
        sbgn.dir <- pgx.system.file("sbgn/", package = "pathway")
        sbgn.dir <- normalizePath(sbgn.dir) ## absolute path
        ## wp = "SMP0080852"
        img <- playbase::getPathwayImage(wp, val = val,
          sbgn.dir = sbgn.dir, as.img = TRUE)
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
    }
  )
}
