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

  options <- shiny::tagList(
    shiny::checkboxGroupInput(ns("pheno"), "Show phenotype:",
      choices = NULL, inline = TRUE
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    options = options,
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
    observeEvent(pgx$samples, {
      shiny::updateCheckboxGroupInput(session, "pheno",
        choices = colnames(pgx$samples),
        selected = head(colnames(pgx$samples), 8)
      )
    })

    get_data <- function() {
      wgcna <- wgcna()

      sel <- enrichTable$rows_selected()
      if (length(sel) > 0) {
        data <- enrichTable$data()
        maintxt <- data$geneset[sel]
        maintxt <- sub(".*:", "", maintxt)
        gg <- strsplit(data$genes[sel], split = "\\|")[[1]]
        pp <- playbase::map_probes(pgx$genes, gg)
        pp <- intersect(pp, rownames(pgx$X))
      } else {
        mod <- selected_module()
        gg <- wgcna$me.genes[[mod]]
        pp <- playbase::map_probes(pgx$genes, gg)
        pp <- intersect(pp, rownames(pgx$X))
        maintxt <- paste(mod, " (top SD)")
      }

      df <- pgx$X[pp, , drop = FALSE]
      shiny::validate(shiny::need(nrow(df) > 1, "Geneset should contain at least two genes to plot a heatmap."))
      if (playbase::is.multiomics(rownames(pgx$X))) {
        rownames(df) <- playbase::probe2symbol(rownames(df), pgx$genes, "gene_name", fill_na = TRUE, add_datatype = TRUE)
      } else {
        rownames(df) <- playbase::probe2symbol(rownames(df), pgx$genes, "gene_name", fill_na = TRUE)
      }

      sel <- input$pheno
      shiny::req(sel)
      annot <- pgx$samples[, sel, drop = FALSE]

      list(df = df, annot = annot, main = maintxt)
    }

    render_plot <- function(nmax, maxlen, show_legend, show_colnames) {
      res <- get_data()
      df <- res$df
      annot <- res$annot

      playbase::gx.splitmap(
        df,
        nmax = nmax,
        col.annot = annot,
        rowlab.maxlen = maxlen,
        show_legend = show_legend,
        show_colnames = show_colnames,
        split = 1,
        main = res$main
      )
    }

    plot.RENDER <- function() {
      render_plot(
        nmax = 20, maxlen = 40, show_legend = FALSE,
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
      csvFunc = get_data,
      pdf.width = 8, pdf.height = 6,
      res = c(80, 100),
      add.watermark = watermark
    )
  })
}
