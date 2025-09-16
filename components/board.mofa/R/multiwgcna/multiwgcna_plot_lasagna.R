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

multiwgcna_plot_lasagna_server <- function(id,
                                           mwgcna,
                                           r_phenotype = reactive(NULL),
                                           r_layers = reactive(NULL),
                                           r_edgenorm = reactive(NULL),
                                           r_edgepos = reactive(NULL),
                                           r_solvesp = reactive(NULL)
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
      shiny::req(!is.null(r_edgepos()))
      shiny::req(!is.null(r_edgenorm()))
      shiny::req(!is.null(r_solvesp()))
      
      graph <- playbase::lasagna.solve(
        obj, pheno,
        min_rho = 0.01,
        max_edges = 1000,
        value = "logFC",
        sp.weight = r_solvesp(),
        prune = FALSE) 

      ## prune graph for plotting
      subgraph <- playbase::lasagna.prune_graph(
        graph, ntop=50,
        layers = NULL,
        normalize.edges = r_edgenorm(),
        min.rho = 0.3,
        edge.sign = "both",
        edge.type = "both", ## inter or intra
        prune = FALSE
      )
      subgraph
      
      par(mfrow=c(1,1), mar=c(1,1,1,1)*0)
      playbase::plotMultiPartiteGraph2(
        subgraph,
        min.rho = 0.1,
        ntop = -50,
        xdist = 1.25,
        cex.label = 0.9,
        vx.cex = 1.1,
        edge.cex = 1.1,
        edge.alpha = 0.3,
        edge.sign = ifelse(r_edgepos(), "pos", "both"),
        edge.type = ifelse(input$showintra, "both", "inter"),  
        yheight = 3.2,
        normalize.edges = 1,
        strip.prefix = TRUE,
        prune = 0
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



