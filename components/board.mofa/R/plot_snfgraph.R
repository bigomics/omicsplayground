##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_snfgraph_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::selectInput(
      ns("labeltype"), "Label type",
      choices = "none"
    )
  )
  
  PlotModuleUI(
    ns("plot"),
    options = options,
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

mofa_plot_snfgraph_server <- function(id,
                                      mofa,
                                      watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    
    observeEvent( mofa()$snf, {
      snf <- mofa()$snf
      labeltypes <- c("<none>","<cluster>","<sample_id>",
                      colnames(snf$samples))      
      shiny::updateSelectInput(
        session, "labeltype", choices = labeltypes,
        selected="<sample_id>" )
                               
    })
    
    plot.RENDER <- function() {
      res <- mofa()
      snf <- res$snf
      validate(need(!is.null(res), "missing MOFA data."))              
      shiny::req(input$labeltype)
      labeltype <- input$labeltype
      label <- NULL
      if(labeltype == "<none>") {
        label <- rep(" ", nrow(snf$samples))
      }
      if(labeltype == "<cluster>") {
        label <- snf$cluster
      }
      if(labeltype == "<sample_id>") {
        label <- rownames(snf$samples)
      }
      if(labeltype %in% colnames(snf$samples)) {
        label <- snf$samples[,labeltype]
      }

      playbase::snf.plot_graph(snf, plot=TRUE, label=label)         
      
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(80, 100),
      add.watermark = watermark
    )

    
  })
}
