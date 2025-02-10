##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_weights_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::checkboxInput(ns("show_top"),"Show top features", TRUE)
  )

  
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

mofa_plot_weights_server <- function(id,
                                     mofa,
                                     pgx,
                                     input_factor = reactive(1),
                                     show_types = reactive(NULL),
                                     watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function(ntop=8) {
      res <- mofa()
      k <- input_factor()
      factors <- colnames(res$F)
      shiny::req(k %in% factors)

      show_types <- show_types()
      dtypes <- names(res$ww)
      show_types <- intersect(show_types, dtypes)
      if(is.null(show_types)) show_types <- dtypes
      shiny::validate(need(length(show_types)>0,
                           "Please select at least one datatype"))
      ww <- res$ww[show_types]
      ntypes <- length(ww)
      ntop <- ifelse(input$show_top, ntop, -1)

      ## switch label type
      ## dbg("[mofa_plot_weights_server:plot.RENDER] head.rownames.ww = ", head(rownames(ww[[1]])))
      ## dbg("[mofa_plot_weights_server:plot.RENDER] colnames(pgx.genes) = ", head(colnames(pgx$genes)))
      ## lab <- pgx$genes$label
      ## dbg("[mofa_plot_weights_server:plot.RENDER] head.lab = ", head(lab))
      ## dbg("[mofa_plot_weights_server:plot.RENDER] head.rownames.ww = ", head(rownames(ww[[1]])))
      ## ww <- playbase::mofa.prefix(ww)
      ## dbg("[mofa_plot_weights_server:plot.RENDER] head.rownames.ww = ", head(rownames(ww[[1]])))      
      ## dbg("[mofa_plot_weights_server:plot.RENDER] rownames(pgx.genes) = ", head(rownames(pgx$genes)))
      ## ww <- lapply(ww, function(w) playbase::rename_by(w, pgx$genes, "label"))
      
      mfrow=c(1,ntypes)
      par(mfrow=mfrow, mar=c(4,8,1.5,0))
      playbase::mofa.plot_weights(ww, k=k, ntop=ntop, maxchar=60)
    }

    plot.RENDER2 <- function() {
      plot.RENDER(n=16) 
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      func2 = plot.RENDER2,      
      pdf.width = 12, pdf.height = 5,
      res = c(75, 110),
      add.watermark = watermark
    )

    
  })
}
