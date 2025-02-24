##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_gene_heatmap_ui <- function(
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

wgcna_plot_gene_heatmap_server <- function(id,
                                           pgx,
                                           wgcna,
                                           enrich_table,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    
    plot.RENDER <- function() {
      wgcna <- wgcna()
      sel <- enrich_table$rows_selected()
      shiny::validate( shiny::need(length(sel)>0, "Please select a geneset"))      
      data <- enrich_table$data()
      gset <- rownames(data)[sel]
      gg <- strsplit(data$genes[sel],split="\\|")[[1]]
      gg <- intersect(gg, rownames(pgx$X))
      
      playbase::gx.splitmap(
        pgx$X[gg,], nmax=50,
        col.annot = pgx$samples,
        rowlab.maxlen = 40,
        show_legend = FALSE,
        split = 1,
        main = gset
      )
      
    }

    plot.RENDER2 <- function() {
      wgcna <- wgcna()
      sel <- enrich_table$rows_selected()
      shiny::validate( shiny::need(length(sel)>0, "Please select a geneset"))      
      data <- enrich_table$data()
      gset <- rownames(data)[sel]
      gg <- strsplit(data$genes[sel],split="\\|")[[1]]
      gg <- intersect(gg, rownames(pgx$X))      
      playbase::gx.splitmap(
        pgx$X[gg,], nmax=45,
        col.annot = pgx$samples,
        ##cexCol = 0.01, cexRow = 0.01,
        rowlab.maxlen = 80,
        show_legend = TRUE,
        split = 1,
        main = gset
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
