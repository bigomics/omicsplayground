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
  width
) {
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
                                           selected_module,
                                           enrichTable,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    render_plot <- function(nmax, maxlen, show_legend,
                            show_colnames) {
      wgcna <- wgcna()
      
      sel <- enrichTable$rows_selected()
      if(length(sel) > 0) {
        data <- enrichTable$data()
        maintxt <- data$geneset[sel]
        maintxt <- sub(".*:","",maintxt)
        gg <- strsplit(data$genes[sel], split = "\\|")[[1]]
        pp <- playbase::map_probes(pgx$genes, gg)
        pp <- intersect(pp, rownames(pgx$X))
      } else {
        mod <- selected_module()
        gg <- wgcna$me.genes[[mod]]
        pp <- playbase::map_probes(pgx$genes, gg)
        pp <- intersect(pp, rownames(pgx$X))
        maintxt <- mod
      }
      
      df <- pgx$X[pp, , drop = FALSE]
      shiny::validate(shiny::need(nrow(df) > 1, "Geneset should contain at least two genes to plot a heatmap."))
      rownames(df) <- playbase::probe2symbol(rownames(df), pgx$genes, "gene_name", fill_na = TRUE)

      playbase::gx.splitmap(
        df,
        nmax = nmax,
        col.annot = pgx$samples,
        rowlab.maxlen = maxlen,
        show_legend = show_legend,
        show_colnames = show_colnames,
        split = 1,
        main = maintxt
      )
    }

    plot.RENDER <- function() {
      render_plot(
        nmax = 30, maxlen = 40, show_legend = FALSE,
        show_colnames = FALSE
      )
    }

    plot.RENDER2 <- function() {
      render_plot(
        nmax = 40, maxlen = 80, show_legend = TRUE,
        show_colnames = TRUE
      )
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      func2 = plot.RENDER2,
      ## csvFunc = csvFunc,
      pdf.width = 8, pdf.height = 6,
      res = c(80, 100),
      add.watermark = watermark
    )
  })
}
