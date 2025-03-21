##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Expression plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
tcga_plot_survival_ui <- function(id, caption, info.text, height, width) {
  ns <- shiny::NS(id)

  tcga_tcgasurv_opts <- tagList(
    withTooltip(
      checkboxInput(ns("tcga_surv_deceasedonly"), "deceased only", FALSE),
      paste(
        "Only include deceased cases in survival analysis,",
        "i.e. exclude censored cases (patients still alive at evaluation time).",
        "This compares strictly the deceased cases, early vs late."
      ),
      placement = "left",
      options = list(container = "body")
    ),
    withTooltip(
      radioButtons(ns("tcga_tcgasurv_ntop"), tspan("N cor genes:"), c(25, 100, 250, 1000), selected = 100, inline = TRUE),
      "Number of top genes for calculating the correlation.",
      placement = "left",
      options = list(container = "body")
    )
  )

  PlotModuleUI(ns("plot"),
    title = "TCGA survival analysis",
    label = "a",
    caption = caption,
    info.text = info.text,
    height = height,
    width = width,
    options = tcga_tcgasurv_opts,
    download.fmt = c("png", "pdf", "svg"),
  )
}

#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
tcga_plot_survival_server <- function(id,
                                      pgx,
                                      contrast,
                                      sigtype,
                                      genelist,
                                      watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    tcga_tcgasurv_test <- reactive({
      matrix_file <- playbase::search_path(c(FILES, FILESX), "tcga_matrix.h5")
      if (is.null(matrix_file)) {
        showNotification("FATAL ERROR: could not find tcga_matrix.h5")
        return(NULL)
      }

      req(pgx$X)

      if (sigtype() == "contrast") {
        req(contrast())

        res <- playbase::pgx.getMetaFoldChangeMatrix(pgx, what = "meta")
        sig <- res$fc[, contrast()]
      } else if (sigtype() == "genelist") {
        req(genelist())
        genes <- as.character(genelist())
        genes <- strsplit(genes, split = "[\t, \n]")[[1]]
        genes <- gsub("[ ]", "", genes)
        sig <- rownames(pgx$X) %in% genes
        names(sig) <- rownames(pgx$X)
      }

      showNotification("Computing survival probabilities...")
      playbase::pgx.testTCGAsurvival(
        sig,
        matrix_file,
        ntop = as.integer(input$tcga_tcgasurv_ntop),
        sortby.p = FALSE,
        deceased.only = input$tcga_surv_deceasedonly,
        min.cases = 10
      )
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = tcga_tcgasurv_test,
      res = c(80, 85),
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
