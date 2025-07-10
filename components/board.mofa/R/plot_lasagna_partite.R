##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_lasagna_partite_ui <- function(
    id,
    title = "",
    info.text = "",
    info.methods = "",    
    info.references = NULL,
    info.extra_link = NULL,    
    caption = info.text,
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)
  
  options = tagList(
    shiny::checkboxInput(ns("showintra"),"show intra-correlation"),
    shiny::checkboxInput(ns("prune"),"prune nodes"),
    shiny::hr(),
    shiny::sliderInput(ns("xdist"),"Layer spacing:",0.2,2,1,0.1),
        shiny::hr(),
    shiny::radioButtons(ns("labeltype"), "Label type:",
      c("feature","symbol","title"), selected="feature", inline=TRUE),
    shiny::checkboxInput(ns("strip_prefix"),"strip prefix", TRUE)        
  )
  
  PlotModuleUI(
    ns("plotmodule"),
    title = title,
    options = options,
    label = "",
    caption = caption,
    info.text = info.text,
    info.references = info.references,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg", "csv")
  )
  
}


mofa_plot_lasagna_partite_server <- function(id,
                                             data,
                                             pgx, 
                                             input_nodevalue = reactive(NULL),
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot_data <- function() {
      res <- data()
      data <- visNetwork::toVisNetworkData(res$graph)
      df <- data$edges
      return(df)
    }

    
    plot_render <- function() {

      res <- data()
      shiny::req(res)
      
      labtype <- input$labeltype
      if(labtype=="title") {
        labels <- paste0(pgx$genes[,"gene_title"],
          " (",pgx$genes[,"symbol"],")")
      } else {
        labels <- pgx$genes[,labtype]  ## user reactive
      }
      names(labels) <- rownames(pgx$genes)
      
      graph <- res$graph
      layers <- setdiff(res$layers,c("SOURCE","SINK"))
      value.name <- input_nodevalue()      
      edge.type <- ifelse(input$showintra, "both", "inter")      
      
      par(mar=c(0,0,0,0))
      playbase::plotMultiPartiteGraph2(
        graph,
        value.name = value.name,
        layers = layers,
        min.rho = 0,
        ntop = -1,
        # labpos = c(2,2,2,4,4),
        xdist = input$xdist,
        labels = labels,
        cex.label = 0.8,
        vx.cex = 1.1,
        edge.cex = 1.3,
        edge.alpha = 0.2,
        edge.type = edge.type,
        yheight = 3,
        normalize.edges = 1,
        strip.prefix2 = input$strip_prefix,
        prune = input$prune
      )
      
    }

    PlotModuleServer(
      "plotmodule",
      func = plot_render,
      csvFunc = plot_data,
      plotlib = "base",
      pdf.width = 12, pdf.height = 6,
      res = c(75, 90),
      add.watermark = watermark
    )
    
  })
}
