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
    shiny::radioButtons(ns("plottype"),"Plot type:",
      choices=c("parallel","hive plot"="hive"), inline=TRUE),
    shiny::radioButtons(ns("labeltype"), "Label type:",
      c("feature","symbol","title"), selected="feature", inline=TRUE),
    shiny::sliderInput(ns("xdist"),"layer spacing:",0.2,2,1,0.1)    
  )
  
  PlotModuleUI(
    ns("mpart"),
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
    download.fmt = c("png", "pdf", "svg")
  )
  
}


mofa_plot_lasagna_partite_server <- function(id,
                                             data,
                                             pgx, 
                                             #input_contrast = reactive(NULL),
                                             input_datatypes = reactive(NULL),
                                             input_minrho = reactive(NULL), 
                                             input_edgetype = reactive(input$edgetype),
                                             ##input_labeltype = reactive(NULL),
                                             input_nodevalue = reactive(NULL),
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {


    ##------------------------------------------------------------------
    ##------------------------------------------------------------------
    ##------------------------------------------------------------------
    
    plot_mpart.RENDER <- function() {

      res <- data()
      shiny::req(res)
      
      labvar <- input$labeltype
      if(labvar=="title") {
        labels <- paste0(pgx$genes[,"gene_title"],
          " (",pgx$genes[,"symbol"],")")
      } else {
        labels <- pgx$genes[,labvar]  ## user reactive
      }
      labels <- playbase::mofa.strip_prefix(labels)
      labels <- paste0(pgx$genes$data_type,":",labels)
      names(labels) <- rownames(pgx$genes)
      
      graph <- res$graph
      group <- igraph::V(graph)$layer
      layers <- setdiff(res$layers,c("SOURCE","SINK"))
      value.name <- input_nodevalue()

    fc <- igraph::V(graph)$value
      names(fc) <- igraph::V(graph)$name
      
      if(input$plottype == "parallel") {

        par(mar=c(0,0,0,0))
        playbase::plotMultiPartiteGraph2(
          graph,
          value.name = value.name,
          layers = layers,
          min.rho = input_minrho(),
          ntop = 50,
          # xpos = c(1,2,3,4,5)*2,
          # xlim = c(-0.5,6),
          # labpos = c(2,2,2,4,4),
          xdist = input$xdist,
          labels = labels,
          cex.label = 0.8,
          vx.cex = 1.1,
          edge.cex = 1.3,
          edge.alpha = 0.2,
          edge.type = input_edgetype(),  
          yheight = 3,
          normalize.edges = 1,
          strip.prefix = TRUE
        )
        
      }

      if(input$plottype == "hive") {

        par(mar=c(0,0,0,0))
        playbase::plotHivePlot(
          res$X,
          fc,
          group = group,
          groups = layers,
          labels = labels,
          cex.label = 0.9,
          vx.cex = 0.8,
          edge.alpha = 0.2,
          edge.type = input_edgetype(),
          min.rho = input_minrho(),
          #min.rho = 0.1,
          ntop = 20) 

      }
    }

    PlotModuleServer(
      "mpart",
      func = plot_mpart.RENDER,
      #csvFunc = plot_data,
      plotlib = "base",
      pdf.width = 12, pdf.height = 6,
      res = c(75, 90),
      add.watermark = watermark
    )
    
  })
}
