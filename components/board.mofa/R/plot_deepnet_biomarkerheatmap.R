##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

plot_deepnet_biomarkerheatmap_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = c("100%", TABLE_HEIGHT_MODAL),
    width = c("auto", "100%")
    ) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    selectInput(ns("ntop"), "Number of features:", c(20,30,50,100,200))
  )
  
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

plot_deepnet_biomarkerheatmap_server <- function(id,
                                                 net,
                                                 pgx,
                                                 update,
                                                 add_annot = c(0,1),
                                                 show_legend = c(0,1),
                                                 ntop = c(20,30),
                                                 rmar = c(0,20),
                                                 plot.res = c(90, 120),
                                                 datatypes = reactive(NULL),
                                                 watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    
    plot.RENDER <- function(n=12) {
      update()  ## react on updates
      net <- net()
      annot <- NULL
      if(add_annot[1]) annot <- pgx$samples[colnames(net$X[[1]]),]      

      # set labels
      gene.labels <- playbase::mofa.strip_prefix(pgx$genes$gene_name)
      gene.labels <- paste0(pgx$genes$data_type,":",gene.labels)
      names(gene.labels) <- pgx$genes$feature
      gset.labels <- paste0("GSET:",rownames(pgx$gsetX))
      names(gset.labels) <- paste0("GSET:",rownames(pgx$gsetX))
      labels <- c(gene.labels, gset.labels)
      
      playbase::deep.plotBiomarkerHeatmap(
        net, ntop = ntop[1],
        datatypes = datatypes(),
        labels = labels,
        balanced = TRUE,
        cexRow = 0.85,
        cexCol = 0.7,
        annot = annot,
        rowlab.maxlen = 25 + rmar[1],
        rownames_width = 45 + rmar[1],
        show_legend = show_legend[1],
        show_colnames = FALSE
      ) 
    }

    plot.RENDER2 <- function(n=12) {
      update()  ## react on updates
      net <- net()
      nsamples <- ncol(net$X[[1]])
      annot <- NULL
      if(add_annot[2]) annot <- pgx$samples[colnames(net$X[[1]]),]
      playbase::deep.plotBiomarkerHeatmap(
        net,
        ntop = ntop[2],
        datatypes = datatypes(),
        balanced = TRUE,        
        rowlab.maxlen = 35 + rmar[2],
        rownames_width = 60 + rmar[2],
        show_legend = show_legend[2],
        annot = annot,
        cexRow = 0.8,
        cexCol = 0.7,
        show_colnames = (nsamples<100)
      ) 
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      func2 = plot.RENDER2,      
      pdf.width = 12, pdf.height = 5,
      res = plot.res,
      add.watermark = watermark
    )


  })
}
