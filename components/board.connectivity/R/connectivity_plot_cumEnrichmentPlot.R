##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' Importance plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
connectivity_plot_cumEnrichmentPlot_ui <- function(id,
                                          label = ""
                                          ) {
  ns <- shiny::NS(id)
  info_text <- strwrap(
    "<b>Meta-enrichment.</b> The barplot visualizes the cumulative enrichment
    of the top-10 most similar profiles. Gene sets that are frequently shared
    with high enrichment will show a higher cumulative scores. You can choose
    between signed or absolute enrichment in the options."
  )

  plot_opts <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("cumgsea_absfc"), "Absolute foldchange", FALSE),
                "Take the absolute foldchange for calculating the cumulative sum.",
                placement = "right", options = list(container = "body")
    ),
    withTooltip(
      shiny::radioButtons(ns("cumgsea_order"), "Order:",
                          choiceValues = c("FC", "cumFC"),
                          choiceNames = c("this FC", "cumFC"),
                          selected = "cumFC", inline = TRUE
      ),
      "How to order the cumulative barplot.",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(ns("plot"),
               title = "Cumulative enrichment",
               label = label,
               plotlib = "plotly",
               info.text = info_text,
               options = plot_opts,
               height = c("auto", 720),
               width = c("auto", 1000)
  )
}

#' Importance plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
connectivity_plot_cumEnrichmentPlot_server <- function(id,
                                              inputData,
                                              cmap_sigdb,
                                              getConnectivityScores,
                                              connectivityScoreTable,
                                              getEnrichmentMatrix,
                                              getCurrentContrast,
                                              watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {

      cumEnrichmentTable <- shiny::reactive({
        cmap_sigdb <- cmap_sigdb()
        shiny::req(cmap_sigdb)
        if (!grepl(".h5$", cmap_sigdb)) {
          return(NULL)
        }

        df <- getConnectivityScores()
        if (is.null(df)) {
          return(NULL)
        }
        ii <- connectivityScoreTable$rows_all()
        shiny::req(ii)

        sel <- head(df$pathway[ii], 10)
        sigdb <- cmap_sigdb
        F <- getEnrichmentMatrix(sigdb, select = sel)
        if (is.null(F)) {
          return(NULL)
        }

        ## multiply with sign of enrichment
        rho1 <- df$rho[match(colnames(F), df$pathway)]
        F <- t(t(F) * sign(rho1))

        if (input$cumgsea_absfc) {
          F <- abs(F)
        }
        F <- F[order(-rowMeans(F**2)), , drop = FALSE]

        ## add current contrast
        ct <- getCurrentContrast()
        gx <- ct$gs[match(rownames(F), names(ct$gs))]
        names(gx) <- rownames(F)
        gx[is.na(gx)] <- 0
        F <- cbind(gx, F)
        colnames(F)[1] <- ct$name

        F
      })

      plot_RENDER <- shiny::reactive({
        ##
        F <- cumEnrichmentTable()
        if (is.null(F)) {
          frame()
          return(NULL)
        }

        NSETS <- 20
        if (input$cumgsea_order == "FC") {
          F <- F[order(-abs(F[, 1])), ]
          F <- head(F, NSETS)
          F <- F[order(F[, 1]), , drop = FALSE]
        } else {
          F <- F[order(-rowMeans(F**2)), ]
          F <- head(F, NSETS)
          F <- F[order(rowMeans(F)), , drop = FALSE]
        }

        rownames(F) <- gsub("H:HALLMARK_", "", rownames(F))
        rownames(F) <- gsub("C2:KEGG_", "", rownames(F))
        rownames(F) <- shortstring(rownames(F), 72)
        maxfc <- max(abs(rowSums(F, na.rm = TRUE)))
        xlim <- c(-1 * (min(F, na.rm = TRUE) < 0), 1.2) * maxfc

        pgx.stackedBarplot(x = data.frame(F),
                           ylab = "cumulative enrichment",showlegend = FALSE
                           )
      })

      plot_RENDER2 <- shiny::reactive({
        plot_RENDER() %>% plotly::layout(showlegend = TRUE)
      })

      PlotModuleServer(
        "plot",
        plotlib = "plotly",
        func = plot_RENDER,
        func2 = plot_RENDER2,
        csvFunc = cumEnrichmentTable,
        pdf.height = 8, pdf.width = 12,
        res = c(72, 90),
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
