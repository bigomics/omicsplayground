##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

multiwgcna_plot_lasagna_ui <- function(
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
      inputId = ns("showintra"),
      label = "Show intra-layer edges",
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
    download.fmt = c("png", "pdf", "svg")
  )
}


multiwgcna_plot_lasagna_inputs <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::checkboxInput(ns("consensus"),"positive path",FALSE),              
    shiny::checkboxInput(ns("norm_edges"),"normalize edges",TRUE),
    shiny::checkboxInput(ns("sp_weight"),"SP weighting",FALSE),
    shiny::checkboxInput(ns("prune"),"Prune nodes ",FALSE),
    shiny::sliderInput(ns("minrho"),"Edge threshold:",0,0.95,0.5,0.05)
  )
}

multiwgcna_plot_lasagna_server <- function(id,
                                           mwgcna,
                                           r_phenotype = reactive(NULL),
                                           r_layers = reactive(NULL)
                                           ) {
  moduleServer(id, function(input, output, session) {

    lasagna_model <- reactive({
      
      wgcna <- mwgcna()

      ## Get eigengene matrices
      ww <- lapply(wgcna, function(w) t(w$net$MEs))
      ww <- lapply(ww, function(w) w[!grepl("[A-Z]{2}grey$",rownames(w)),,drop=FALSE])

      sel.layers <- r_layers()
      sel.layers <- intersect(sel.layers, names(ww))
      shiny::req(sel.layers)
      ww <- ww[sel.layers]

      data2 <- list(
        X = ww,
        samples = wgcna[[1]]$datTraits
      )
      progress <- shiny::Progress$new(session, min=0, max=1)
      on.exit(progress$close())
      progress$set(message = paste("creating Lasagna model..."), value = 0.33)
      
      ## Create lasagna model
      obj <- playbase::lasagna.create_model(
        data2,
        pheno = "expanded",
        ntop = 2000,
        nc = 20,
        add.sink = TRUE,
        intra = TRUE,
        fully_connect = FALSE,
        add.revpheno = TRUE
      )

      obj
    })

    
    plot.RENDER <- function() {

      obj <- lasagna_model()
      
      pheno = colnames(obj$Y)[1]
      pheno <- r_phenotype()
      shiny::req(!is.null(pheno))

      ## from my module inputs
      sp.weight <- input$sp_weight
      min.rho <- input$minrho
      consensus <- input$consensus
      normalize.edges <- input$norm_edges
      
      ## 'solve' model
      graph <- playbase::lasagna.solve(
        obj, pheno,
        min_rho = 0.01,
        max_edges = 1000,
        value = "logFC",
        sp.weight = sp.weight,
        prune = input$prune)
      
      ## prune graph for plotting
      ## subgraph <- playbase::lasagna.prune_graph(
      ##   graph,
      ##   ntop = 50,
      ##   layers = NULL,
      ##   normalize.edges = normalize.edges,
      ##   min.rho = min.rho,
      ##   edge.sign = "both",
      ##   edge.type = "both", ## inter or intra
      ##   prune = FALSE
      ## )
      ## subgraph
      
      par(mfrow=c(1,1), mar=c(1,1,1,1)*0)
      playbase::plotMultiPartiteGraph2(
        graph,
        min.rho = min.rho,
        ntop = -50,
        xdist = 1.25,
        cex.label = 0.9,
        vx.cex = 1.1,
        edge.cex = 1.1,
        edge.alpha = 0.3,
        edge.sign = ifelse(consensus, "consensus", "both"),
        edge.type = ifelse(input$showintra, "both2", "inter"),  
        yheight = 3.2,
        normalize.edges = normalize.edges,
        strip.prefix = TRUE,
        prune = input$prune
      ) 

    }
    
    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 8,
      pdf.height = 12,
      res = c(80, 100),
      add.watermark = FALSE
    )

    
  })
}



