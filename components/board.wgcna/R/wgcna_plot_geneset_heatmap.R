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
  width
) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxGroupInput(ns("pheno"), "Show phenotype:",
      choices = NULL, inline = TRUE
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    options = options,
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
    observeEvent(pgx$samples, {
      shiny::updateCheckboxGroupInput(session, "pheno",
        choices = colnames(pgx$samples),
        selected = head(colnames(pgx$samples), 8)
      )
    })

    get_data <- function(n = 20) {
      df <- enrichTable$data()
      if (is.null(df) || nrow(df) == 0) {
        return(NULL)
      }
      ii <- enrichTable$rows_all()
      shiny::req(ii)
      sel <- df$geneset[ii]
      sel1 <- intersect(sel, rownames(pgx$gsetX))
      sel1 <- head(sel1, n)
      shiny::validate(shiny::need(length(sel1), "Error getting data"))
      gsetX <- pgx$gsetX[sel1, , drop = FALSE]
      list(gsetX = gsetX)
    }

    plot_heatmap <- function(n = 20, maxlen = 120) {
      res <- get_data(n = n)
      gsetX <- res$gsetX
      mod <- selected_module()

      annot <- pgx$samples
      sel <- input$pheno
      shiny::req(sel)
      sel <- intersect(sel, colnames(annot))
      annot <- annot[, sel, drop = FALSE]

      playbase::gx.splitmap(
        gsetX,
        nmax = 9999,
        col.annot = annot,
        rowlab.maxlen = maxlen,
        show_legend = FALSE,
        show_colnames = FALSE,
        split = 1,
        main = mod,
        verbose = 2
      )
    }

    plot.RENDER <- function() {
      plot_heatmap(n = 20, maxlen = 80)
    }
    plot.RENDER2 <- function() {
      plot_heatmap(n = 40, maxlen = 240)
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      func2 = plot.RENDER2,
      csvFunc = get_data,
      pdf.width = 8, pdf.height = 6,
      res = c(80, 100),
      add.watermark = watermark
    )
  })
}
