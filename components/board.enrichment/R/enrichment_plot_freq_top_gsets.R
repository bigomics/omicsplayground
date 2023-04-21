##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_plot_freq_top_gsets_ui <- function(
  id,
  title,
  info.text,
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
        "Weight by geneset size", TRUE
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
    caption = caption,
    options = topEnrichedFreq.opts,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "csv")
  )
}

enrichment_plot_freq_top_gsets_server <- function(id,
                                                  pgx,
                                                  gseatable,
                                                  getFilteredGeneSetTable,
                                                  gs_contrast,
                                                  gseatable_rows_selected,
                                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx,  getFilteredGeneSetTable(), gs_contrast(),gseatable)
      rpt <- getFilteredGeneSetTable()

      comp <- gs_contrast()
      if (is.null(comp)) {
        shiny::validate(shiny::need(!is.null(comp)), "Please select geneset in the table below.")
      }
      if (!(comp %in% names(pgx$gx.meta$meta))) {
        shiny::validate(shiny::need((comp %in% names(pgx$gx.meta$meta)), "Please select geneset in the table below."))
      }

      ## filter on active rows (using search)
      ii <- gseatable_rows_selected()

      rpt <- rpt[ii, , drop = FALSE]
      
      if(nrow(rpt) == 0) {
        shiny::validate(shiny::need(nrow(rpt) > 0, "Please select geneset in the table below."))
      }
      
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

    topEnrichedFreq.RENDER <- function() {
      dt <- plot_data()
      shiny::req(dt)
      pgx <- dt[[1]]
      rpt <- dt[[2]]
      ntop <- dt[[3]]
      gset.weight <- dt[[4]]
      fcweight <- dt[[5]]
      ngenes <- 35
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
      F <- head(F[order(-Matrix::rowSums(abs(F))), , drop = FALSE], ngenes)
      F <- F[order(-Matrix::rowSums(F)), , drop = FALSE]

      sel.zero <- which(Matrix::rowSums(abs(F)) < 1e-4)
      if (length(sel.zero)) rownames(F)[sel.zero] <- ""

      playbase::pgx.stackedBarplot(
        x = F,
        ylab = ifelse(wt, "weighted frequency", "frequency"),
        showlegend = FALSE
      )
    }

    PlotModuleServer(
      "plot",
      func = topEnrichedFreq.RENDER,
      plotlib = "plotly",
      pdf.width = 5,
      pdf.height = 5,
      res = c(68, 100),
      csvFunc = plot_data,
      add.watermark = watermark
    )
  })
}
