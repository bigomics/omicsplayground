##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_geneset_heatmap_ui <- function(
    id,
    label,
    title,
    info.text,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)
  
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    height = height,
    caption = caption,
    width = width,
    download.fmt = c("png", "pdf", "csv", "svg")
  )
}

wgcna_plot_geneset_heatmap_server <- function(id,
                                              pgx,
                                              wgcna,
                                              selected_module,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    
    plot.RENDER <- function() {
      wgcna <- wgcna()
      mod <- selected_module()
      df <- wgcna$gse[[mod]]
      ##sel <- head(rownames(df),20)
      sel <- head(df$geneset,20)
      shiny::req(mod, length(sel)>0)

      gsetX <- pgx$gsetX[sel,,drop=FALSE]
      playbase::gx.splitmap(
        gsetX, nmax=50,
        col.annot = pgx$samples,
        ##cexCol = 0.01, cexRow = 0.01,
        rowlab.maxlen = 120,
        show_legend = FALSE,
        show_colnames = FALSE,        
        split = 1,
        main = paste("Module",mod),
        verbose = 2
      )
      
    }

    plot.RENDER2 <- function() {
      wgcna <- wgcna()
      mod <- selected_module()
      df <- wgcna$gse[[mod]]
      sel <- head(rownames(df),40)

      playbase::gx.splitmap(
        pgx$gsetX[sel,], nmax=50,
        ##mar = c(1,1), # keysize=1,
        col.annot = pgx$samples,
        ##cexCol = 0.01, cexRow = 0.01,
        rowlab.maxlen = 200,
        show_legend = TRUE,
        split = 1,
        main = paste("Module",mod)
      )
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      func2 = plot.RENDER2,      
      ##csvFunc = csvFunc,
      pdf.width = 8, pdf.height = 6,
      res = c(80, 100),
      add.watermark = watermark
    )
  })
}
