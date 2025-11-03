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
connectivity_plot_cumFCplot_ui <- function(
  id,
  title,
  info.text,
  caption,
  label = "",
  height,
  width
) {
  ns <- shiny::NS(id)

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
    title = title,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    caption = caption,
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
                                               getProfiles,
                                               getConnectivityScores,
                                               getCurrentContrast,
                                               watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      cumulativeFCtable <- shiny::reactive({
        F <- getProfiles()
        F[is.na(F)] <- 0


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
        colnames(F)[1] <- paste("********", cc$name, "********")

        if (input$cumFCplot_absfc) {
          F <- abs(F)
        }
        F <- F[order(-rowMeans(F**2, na.rm = TRUE)), , drop = FALSE]
        F
      })

      render_stackedbar <- function(ngenes) {
        F <- cumulativeFCtable()
        shiny::req(F)

        MAXF <- 10 ## number of top signatures
        #

        F <- F[, 1:min(MAXF, ncol(F)), drop = FALSE]
        if (input$cumFCplot_order == "FC") {
          F <- F[order(-abs(F[, 1])), ]
          F1 <- head(F, ngenes)
          F1 <- F1[order(F1[, 1]), , drop = FALSE]
        } else {
          F1 <- head(F, ngenes)
          F1 <- F1[order(rowMeans(F1, na.rm = TRUE)), , drop = FALSE]
        }

        playbase::pgx.stackedBarplot(
          x = data.frame(F1, check.names = FALSE),
          ylab = "cumulative logFC", xlab = "",
          showlegend = FALSE
        )
      }

      plot_RENDER <- function() {
        render_stackedbar(40) %>%
          plotly::layout(showlegend = FALSE)
      }

      plot_RENDER2 <- function() {
        render_stackedbar(60) %>%
          plotly::layout(showlegend = TRUE)
      }

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
