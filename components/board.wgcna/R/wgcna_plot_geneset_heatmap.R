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
                                              selected_module,
                                              enrichTable,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot_heatmap <- function(n=20, maxlen=120) {
      mod <- selected_module()

      df <- enrichTable$data()
      if (is.null(df) || nrow(df) == 0) {
        return(NULL)
      }
      ii <- enrichTable$rows_all()
      shiny::req(ii)
      sel <- head(df$geneset[ii], n)

      sel1 <- intersect(sel, rownames(pgx$gsetX))
      if(length(sel1)) {
        gsetX <- pgx$gsetX[sel1, , drop = FALSE]
      } else {
        sel1 <- intersect(sel, rownames(playdata::GSETxGENE))
        X <- playbase::rename_by( pgx$X, pgx$genes, "human_ortholog")
        G <- Matrix::t( playdata::GSETxGENE[sel1,] )
        gsetX <- plaid::plaid(X, G)
      }

      playbase::gx.splitmap(
        gsetX,
        nmax = 50,
        col.annot = pgx$samples,
        ## cexCol = 0.01, cexRow = 0.01,
        rowlab.maxlen = maxlen,
        show_legend = FALSE,
        show_colnames = FALSE,
        split = 1,
        main = mod,
        verbose = 2
      )

    }

    plot.RENDER <- function() {
      plot_heatmap(n=20, maxlen=80)
    }
    plot.RENDER2 <- function() {
      plot_heatmap(n=40, maxlen=240)
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
