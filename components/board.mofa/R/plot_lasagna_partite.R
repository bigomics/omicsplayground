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

  PlotModuleUI(
    ns("mpart"),
    title = title,
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


mofa_plot_lasagna_partite_adjmat_ui <- function(
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

  PlotModuleUI(
    ns("adjmat"),
    title = title,
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
                                             input_contrast = reactive(NULL),
                                             input_datatypes = reactive(NULL),
                                             input_minrho = reactive(NULL), 
                                             input_edgetype = reactive(input$edgetype),
                                             input_labeltype = reactive(NULL),
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot_data <- function() {
      
      shiny::req(data())
      shiny::req(pgx$X)

      shiny::validate( shiny::need( !is.null(pgx$mofa), "Object has no mofa slot"))      
      shiny::validate( shiny::need( length(input_datatypes())>=2, "You need to select at least two datatypes"))
      
      X <- pgx$mofa$X
      F <- playbase::pgx.getMetaMatrix(pgx)$fc

      ct=1
      ct <- input_contrast()
      fc <- F[,ct]

      dtype <- sub(":.*","",names(fc))
      sel <- input_datatypes()
      sel <- intersect(sel, unique(dtype))
      dt.features <- tapply( fc, dtype, function(x) {
        head(names(sort(-abs(x))),50)
      })
      features <- lapply(sel, function(s) unlist(dt.features[names(dt.features) %in% s]))
      features  <- intersect(unlist(features),rownames(X))
      X <- X[features,,drop=FALSE]

      fc <- fc[rownames(X)]
      res <- list(fc=fc, X=X, groups=sel)
      res
    }

    ##------------------------------------------------------------------
    ##------------------------------------------------------------------
    ##------------------------------------------------------------------
    
    plot_mpart.RENDER <- function() {
      res <- plot_data()
      shiny::req(res)

      labvar <- input_labeltype()
      if(labvar=="title") {
        labels <- paste0(pgx$genes[,"gene_title"],
                         " (",pgx$genes[,"symbol"],")")
      } else {
        labels <- pgx$genes[,labvar]  ## user reactive
      }
      labels <- playbase::mofa.strip_prefix(labels)
      labels <- paste0(pgx$genes$data_type,":",labels)

      names(labels) <- rownames(pgx$genes)
      xt <- playbase::mofa.get_prefix(rownames(res$X))
      
      par(mar=c(0,0,0,0))
      playbase::plotMultiPartiteGraph(
        res$X,
        res$fc,
        group = xt,
        groups = res$groups,
        labels = labels,
        cex.label = 0.9,
        vx.cex = 1,
        edge.alpha = 0.2,
        edge.type = input_edgetype(),
        yheight = 2,
        min.rho = input_minrho(),
        ntop = 40) 
      
    }

    PlotModuleServer(
      "mpart",
      func = plot_mpart.RENDER,
      csvFunc = plot_data,
      plotlib = "base",
      pdf.width = 12, pdf.height = 6,
      res = c(75, 90),
      add.watermark = watermark
    )

    ##------------------------------------------------------------------
    ##------------------------------------------------------------------
    ##------------------------------------------------------------------
    
    plot_adjmat.RENDER <- function() {
      res <- plot_data()
      shiny::req(res)

      labvar <- input_labeltype()
      if(labvar=="title") {
        labels <- paste0(pgx$genes[,"gene_title"],
                         " (",pgx$genes[,"symbol"],")")
      } else {
        labels <- pgx$genes[,labvar]  ## user reactive
      }
      labels <- playbase::mofa.strip_prefix(labels)
      labels <- paste0(pgx$genes$data_type,":",labels)

      names(labels) <- rownames(pgx$genes)
      xt <- playbase::mofa.get_prefix(rownames(res$X))
      
      graph <- playbase::plotMultiPartiteGraph(
        res$X,
        res$fc,
        justgraph = TRUE,
        group = xt,
        groups = res$groups,
        labels = labels,
        min.rho = input_minrho(),
        ntop = 40) 

      par(mar=c(1,1,1,1)*2)
      par(oma=c(0,0,2,0))
      playbase::plotAdjacencyMatrixFromGraph(
        graph,
        nmax = 40,
        key=TRUE, keysize=0.85, 
        binary = FALSE,
        mar=c(10,8)*1.2,
        cexRow=0.9, cexCol=0.9
      ) 
      
    }

    PlotModuleServer(
      "adjmat",
      func = plot_adjmat.RENDER,
      csvFunc = plot_data,
      plotlib = "base",
      pdf.width = 12, pdf.height = 6,
      res = c(80, 85),
      add.watermark = watermark
    )

    
    
  })
}
