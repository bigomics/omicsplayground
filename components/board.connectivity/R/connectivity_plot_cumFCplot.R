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
connectivity_plot_cumFCplot_ui <- function(id,
                                           label = "",
                                           height,
                                           width) {
  ns <- shiny::NS(id)
  info_text <- strwrap("<b>Meta-foldchange.</b> The barplot visualizes the
                       cumulative foldchange between the top-10 most similar
                       profiles. Genes that are frequently shared with high
                       foldchange will show a higher cumulative score. You can
                       choose between signed or absolute foldchange in the options.")

  plot_opts <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("cumFCplot_absfc"), "Absolute foldchange", FALSE),
      "Take the absolute foldchange for calculating the cumulative sum.",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(
      shiny::radioButtons(ns("cumFCplot_order"), "Order:",
        choiceValues = c("FC", "cumFC"),
        choiceNames = c("this FC", "cumFC"),
        selected = "cumFC", inline = TRUE
      ),
      "How to order the cumulative barplot.",
      placement = "right", options = list(container = "body")
    )
  )
  PlotModuleUI(ns("plot"),
    title = "Cumulative foldchange",
    label = label,
    plotlib = "plotly",
    info.text = info_text,
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
connectivity_plot_cumFCplot_server <- function(id,
                                               getTopProfiles,
                                               getConnectivityScores,
                                               getCurrentContrast,
                                               watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      cumulativeFCtable <- shiny::reactive({
        F <- getTopProfiles()
        F[is.na(F)] <- 0

        ## maximum 10??
        MAXF <- 20

        ## multiply with sign of rho
        df <- getConnectivityScores()
        rho1 <- df$rho[match(colnames(F), df$pathway)]
        F <- t(t(F) * sign(rho1))

        ## add current contrast
        cc <- getCurrentContrast()
        shiny::req(cc)
        fc <- cc$fc[rownames(F)]
        fc[is.na(fc)] <- 0
        F <- cbind(fc[rownames(F)], F)
        colnames(F)[1] <- "thisFC"
        colnames(F)[1] <- cc$name

        if (input$cumFCplot_absfc) {
          F <- abs(F)
        }
        F <- F[order(-rowMeans(F**2)), , drop = FALSE]
        F
      })

      plot_RENDER <- shiny::reactive({
        F <- cumulativeFCtable()
        shiny::req(F)

        MAXF <- 10
        NGENES <- 50

        F <- F[, 1:min(MAXF, ncol(F)), drop = FALSE]
        if (input$cumFCplot_order == "FC") {
          F <- F[order(-abs(F[, 1])), ]
          F1 <- head(F, NGENES)
          F1 <- F1[order(F1[, 1]), , drop = FALSE]
        } else {
          F1 <- head(F, NGENES)
          F1 <- F1[order(rowMeans(F1)), , drop = FALSE]
        }

        playbase::pgx.stackedBarplot(
          x = data.frame(F1, check.names=FALSE),
          ylab = "cumulative logFC", xlab="",
          showlegend = FALSE
        )
      })

      plot_RENDER2 <- shiny::reactive({
        plot_RENDER() %>% plotly::layout(showlegend = TRUE)
      })

      PlotModuleServer(
        "plot",
        func = plot_RENDER,
        func2 = plot_RENDER2,
        plotlib = "plotly",
        csvFunc = cumulativeFCtable,
        pdf.height = 6, pdf.width = 9,
        res = c(72, 90),
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
