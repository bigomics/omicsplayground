##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
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
connectivity_plot_cumEnrichmentPlot_ui <- function(
  id,
  label = "",
  title,
  info.text,
  caption,
  height,
  width
) {
  ns <- shiny::NS(id)

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
    title = title,
    caption = caption,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    options = plot_opts,
    height = height,
    width = width
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
                                                       pgx,
                                                       sigdb,
                                                       getConnectivityScores,
                                                       connectivityScoreTable,
                                                       getEnrichmentMatrix,
                                                       getCurrentContrast,
                                                       watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      getEnrichmentTable <- shiny::reactive({
        sigdb <- sigdb()
        shiny::req(sigdb)
        if (!grepl(".h5$", sigdb)) {
          return(NULL)
        }

        df <- getConnectivityScores()
        if (is.null(df)) {
          return(NULL)
        }
        ii <- connectivityScoreTable$rows_all()
        shiny::req(ii)
        sel <- head(df$pathway[ii], 10)
        F <- getEnrichmentMatrix(sigdb, select = sel)

        if (is.null(F)) {
          return(NULL)
        }

        ## multiply with sign of enrichment
        if (input$cumgsea_absfc) {
          F <- abs(F)
        } else {
          rho1 <- df$rho[match(colnames(F), df$pathway)]
          F <- t(t(F) * sign(rho1))
        }
        ## order by term/geneset
        F <- F[order(-rowMeans(F**2, na.rm = TRUE)), , drop = FALSE]

        ## add current contrast
        ct <- getCurrentContrast()
        gx <- ct$gs[match(rownames(F), names(ct$gs))]
        names(gx) <- rownames(F)
        gx[is.na(gx)] <- 0
        F <- cbind(gx, F)
        colnames(F)[1] <- ct$name

        ## normalize columns
        F <- t(t(F) / sqrt(colSums(F**2)))
        F
      })

      plot_stackedbar <- function(nsets) {
        ##
        F <- getEnrichmentTable()
        if (is.null(F)) {
          frame()
          return(NULL)
        }

        #
        if (input$cumgsea_order == "FC") {
          F <- F[order(-abs(F[, 1])), ]
          F <- head(F, nsets)
          F <- F[order(F[, 1]), , drop = FALSE]
        } else {
          F <- F[order(-rowMeans(F**2, na.rm = TRUE)), ]
          F <- head(F, nsets)
          F <- F[order(rowMeans(F, na.rm = TRUE)), , drop = FALSE]
        }

        ## strip comments/prefixes to shorten names
        rownames(F) <- gsub("H:HALLMARK_", "", rownames(F))
        rownames(F) <- gsub("C2:KEGG_", "", rownames(F))
        rownames(F) <- playbase::shortstring(rownames(F), 72)

        playbase::pgx.stackedBarplot(
          x = data.frame(F, check.names = FALSE),
          #
          ylab = "cumulative enrichment",
          xlab = "",
          showlegend = FALSE
        )
      }

      plot_RENDER <- function() {
        plot_stackedbar(30) %>%
          plotly::layout(showlegend = FALSE)
      }

      plot_RENDER2 <- function() {
        plot_stackedbar(50) %>%
          plotly::layout(showlegend = TRUE)
      }

      PlotModuleServer(
        "plot",
        plotlib = "plotly",
        func = plot_RENDER,
        func2 = plot_RENDER2,
        csvFunc = getEnrichmentTable,
        pdf.height = 8, pdf.width = 12,
        res = c(72, 90),
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
