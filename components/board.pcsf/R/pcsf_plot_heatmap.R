##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Expression plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
pcsf_plot_heatmap_ui <- function(id, caption, info.text, height, width) {
  ns <- shiny::NS(id)

  plot_opts <- tagList(
    withTooltip(
      checkboxInput(ns("pcsf_option1"), "option1", FALSE),
      paste(
        "...."
      ),
      placement = "left",
      options = list(container = "body")
    ),
    withTooltip(
      radioButtons(ns("pcsf_option2"), "N cor genes:",
                   c(25, 100, 250, 1000), selected = 100, inline = TRUE),
      "...",
      placement = "left",
      options = list(container = "body")
    )
  )

  PlotModuleUI(ns("plotmodule"),
    title = "PCSF heatmap analysis",
    label = "a",
    plotlib = "base",
    caption = caption,
    info.text = info.text,
    height = height,
    width = width,
    options = plot_opts,
    download.fmt = c("png", "pdf"),
  )
}

#' PCSF network function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
pcsf_plot_heatmap_server <- function(id,
                                     pgx,
                                     pcsf_compute,
                                     colorby = reactive("cluster"),
                                     contrast = reactive(NULL),
                                     watermark = FALSE
) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    get_network <- reactive({
      colorby <- input$colorby
      res <-  pcsf_compute()
      net <- res$net
      fx <- res$meta[,1]
      vv <- igraph::V(net)$name
      igraph::V(net)$type <- c("down","up")[1 + 1*(sign(fx[vv])>0)]
      net
    })

    heatmap.RENDER <- function() {
      shiny::req(pgx$X)
      net <- get_network()
      dbg("[pcsf_plot_heatmap.R:heatmap.RENDER] reacted!")
      dbg("[pcsf_plot_heatmap.R:heatmap.RENDER] names.pgx=",names(pgx))      
      
      idx <- igraph::V(net)$group
      genes <- igraph::V(net)$name

      dbg("[pcsf_plot_heatmap.R:heatmap.RENDER] head.genes=",head(genes))
      dbg("[pcsf_plot_heatmap.R:heatmap.RENDER] head.idx=",head(idx))
      
      playbase::gx.splitmap(
        pgx$X[genes,],
        col.annot = pgx$samples,
        ## row.annot = (cc[genes,1:5]),  
        split = idx,
        scale = 'row',
        cexRow = 1,
        show_rownames = 0,
        ##keysize=0.5,
        mar=c(4,16),
        show_legend = FALSE,
        key.offset = c(0.05, 0.96)
      )
      
      dbg("[pcsf_plot_heatmap.R:heatmap.RENDER] done!")
    }
    
    PlotModuleServer(
      "plotmodule",
      func = heatmap.RENDER,
      plotlib = "base",
      pdf.width = 10,
      pdf.height = 10,
      add.watermark = watermark
    )
    
  }) ## end of moduleServer
}
