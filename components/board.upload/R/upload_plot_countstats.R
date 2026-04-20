##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

upload_plot_countstats_ui <- function(id,
                                      label = "",
                                      height,
                                      width,
                                      title,
                                      caption,
                                      info.text) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    plotlib = "base",
    info.text = info.text,
    caption = caption,
    options = NULL,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height
  )
}

upload_plot_countstats_server <- function(id, checkTables, countsRT, watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ## extract data from pgx object
    plot_data <- shiny::reactive({
      counts <- countsRT()
      has.counts <- !is.null(counts) && NCOL(counts) > 0
      check <- checkTables()
      shiny::req(check)

      status.ok <- check["counts.csv", "status"]
      status.ds <- tolower(check["counts.csv", "description"])
      error.msg <- paste(
        toupper(status.ok), "\nPlease upload 'counts.csv' (Required):",
        status.ds
      )
      shiny::validate(
        shiny::need(
          status.ok == "OK" && has.counts,
          error.msg
        )
      )
      return(counts)
    })

    plot.RENDER <- function() {
      counts <- plot_data()
      ## Subsample both genes and cells before densifying.
      ## Adding scalar 1 to a sparse matrix fills all structural zeros, breaking sparsity.
      ## For a density plot 500 cells is more than sufficient for visual purposes.
      n_genes <- nrow(counts)
      n_samples <- ncol(counts)
      if (nrow(counts) > 1000) counts <- counts[sample(nrow(counts), 1000), , drop = FALSE]
      if (ncol(counts) > 500) counts <- counts[, sample(ncol(counts), 500), drop = FALSE]
      if (inherits(counts, "sparseMatrix")) counts <- as.matrix(counts)
      xx <- log2(1 + counts)
      suppressWarnings(dc <- data.table::melt(xx))
      dc$value[dc$value == 0] <- NA
      tt2 <- paste(n_genes, tspan("genes x", js = FALSE), n_samples, "samples")
      ggplot2::ggplot(dc, ggplot2::aes(x = value, color = Var2)) +
        ggplot2::geom_density() +
        ggplot2::xlab("log2(1+counts)") +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ggtitle("COUNTS", subtitle = tt2)
    }

    modal_plot.RENDER <- function() {
      plot.RENDER()
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "base",
      func = plot.RENDER,
      func2 = modal_plot.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(90, 90), ## resolution of plots
      pdf.width = 4, pdf.height = 4,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
