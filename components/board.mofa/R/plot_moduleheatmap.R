##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_moduleheatmap_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(
      inputId = ns("split"),
      label = "Split heatmap by data types",
      value = FALSE
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    options = options,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

mofa_plot_moduleheatmap_server <- function(id,
                                           mofa,
                                           input_factor = reactive(1),
                                           show_types = reactive(NULL),
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function() {
      res <- mofa()      
      k <- input_factor()
      factors <- colnames(res$F)
      if(!is.null(k)) shiny::req(k %in% factors)

      show_types <- show_types()
      dtypes <- names(res$ww)
      if(is.null(show_types)) show_types <- dtypes
      shiny::validate(need(length(show_types)>0,
                           "Please select at least one datatype"))
      
      playbase::mofa.plot_loading_heatmap(
        res, k=k, main=k,
        ntop = 50, split = input$split,
        type="splitmap", annot = "pheno",
        maxchar = 40, show_types = show_types,
        mar = c(3,0,0,0), annot.ht = 3.5,
        cexRow = 0.9 )
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 8, pdf.height = 12,
      res = c(80, 90),
      add.watermark = watermark
    )

    
  })
}



