##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_plot_freq_top_gsets_ui <- function(
    id,
    title,
    info.text,
    info.methods,
    info.extra_link,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)

  topEnrichedFreq.opts <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(ns("gs_enrichfreq_ntop"), "Number of top sets",
        c(5, 10, 15),
        inline = TRUE, selected = 15
      ),
      "Number of top genesets to consider for counting the gene frequency."
    ),
    withTooltip(
      shiny::checkboxInput(
        ns("gs_enrichfreq_gsetweight"),
        tspan("Weight by geneset size"), TRUE
      ),
      "Weight by (inverse) gene set size."
    ),
    withTooltip(
      shiny::checkboxInput(
        ns("gs_enrichfreq_fcweight"),
        "Weight by FC", TRUE
      ),
      "Weight by fold-change of current contrast."
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = "b",
    plotlib = "plotly",
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    caption = caption,
    options = topEnrichedFreq.opts,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "csv")
  )
}

enrichment_plot_freq_top_gsets_server <- function(id,
                                                  pgx,
                                                  getFilteredGeneSetTable,
                                                  gs_contrast,
                                                  gseatable,
                                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      rpt <- getFilteredGeneSetTable()
      shiny::req(pgx$X, rpt, gs_contrast())

      comp <- gs_contrast()
      if (is.null(comp)) {
        return(NULL)
      }
      if (!(comp %in% names(pgx$gx.meta$meta))) {
        return(NULL)
      }

      ## filter on active rows (using search)
      ii <- gseatable$rows_current()
      req(ii)
      rpt <- rpt[ii, , drop = FALSE]

      ntop <- as.integer(input$gs_enrichfreq_ntop)
      gset.weight <- input$gs_enrichfreq_gsetweight
      fcweight <- input$gs_enrichfreq_fcweight
      return(
        list(
          pgx,
          rpt,
          ntop,
          gset.weight,
          fcweight
        )
      )
    })

    topEnrichedFreq.RENDER <- function(return_csv = FALSE) {
      dt <- plot_data()
      shiny::req(dt)
      pgx <- dt[[1]]
      rpt <- dt[[2]]
      ntop <- dt[[3]]
      gset.weight <- dt[[4]]
      fcweight <- dt[[5]]
      ngenes <- 30
      fx.col <- grep("score|fx|fc|sign|NES|logFC", colnames(rpt))[1]
      fx <- rpt[, fx.col]
      names(fx) <- rownames(rpt)

      top <- rownames(rpt)
      top <- head(top, ntop)
      if (!all(top %in% colnames(pgx$GMT))) {
        return(NULL)
      }

      F <- 1 * (pgx$GMT[, top, drop = FALSE] > 0)
      F <- as.matrix(F)
      wt <- FALSE
      if (gset.weight) {
        F <- Matrix::t(Matrix::t(F) / Matrix::colSums(F, na.rm = TRUE))
        wt <- TRUE
      }
      F <- Matrix::t(Matrix::t(F) * sign(fx[top]))
      if (fcweight) {
        F <- Matrix::t(Matrix::t(F) * abs(fx[top]))
        wt <- TRUE
      }
      # sum duplicated rows
      F <- rowsum(F, row.names(F))

      F <- head(F[order(-Matrix::rowSums(abs(F), na.rm = TRUE)), , drop = FALSE], ngenes)
      F <- F[order(-Matrix::rowSums(F, na.rm = TRUE)), , drop = FALSE]

      sel.zero <- which(Matrix::rowSums(abs(F), na.rm = TRUE) < 1e-4)
      if (length(sel.zero)) F <- F[-sel.zero, , drop = FALSE]

      if (return_csv) {
        return(F)
      }

      playbase::pgx.stackedBarplot(
        x = F,
        ylab = ifelse(wt, "weighted frequency", "frequency"),
        xlab = tspan("genes", js = FALSE),
        showlegend = FALSE
      )
    }

    plot_data_csv <- function() {
      df <- topEnrichedFreq.RENDER(return_csv = TRUE)
      return(df)
    }

    PlotModuleServer(
      "plot",
      func = topEnrichedFreq.RENDER,
      plotlib = "plotly",
      pdf.width = 5,
      pdf.height = 5,
      res = c(68, 100),
      csvFunc = plot_data_csv,
      add.watermark = watermark
    )
  })
}
