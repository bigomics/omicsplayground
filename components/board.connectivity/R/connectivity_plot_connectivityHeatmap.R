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
connectivity_plot_connectivityHeatmap_ui <- function(id,
                                                     label = "",
                                                     rowH = 660) {
  ns <- shiny::NS(id)
  info_text <- strwrap(
    "<b>The Connectivity Heatmap</b> shows the most similar profiles as a heatmap.
    Contrasts that are similar will be clustered close together."
  )
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
    title = "Connectivity Heatmap",
    label = label,
    plotlib = "base",
    info.text = info_text,
    options = plot_opts,
    height = c(420, 650),
    width = c("auto", 1400)
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
connectivity_plot_connectivityHeatmap_server <- function(id,
                                                         getTopProfiles,
                                                         getConnectivityScores,
                                                         getCurrentContrast,
                                                         watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      cumulativeFCtable <- shiny::reactive({
        F <- getTopProfiles()
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
        colnames(F)[1] <- cc$name

        if (input$cumFCplot_absfc) {
          F <- abs(F)
        }
        F <- F[order(-rowMeans(F**2)), , drop = FALSE]
        F
      })

      plotFCheatmap <- function(F, maxfc, maxgenes=60) {
        F <- F[, 1:min(NCOL(F), maxfc), drop = FALSE]
        if (input$cumFCplot_order == "FC") {
          F <- F[order(-abs(F[, 1])), ]
        }
        F1 <- head(F, maxgenes)
        par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
        playbase::gx.splitmap(t(F1),
          split = 1,
          ## cluster_columns = FALSE,
          cluster_columns = TRUE,
          cluster_rows = TRUE,
          rowlab.maxlen = 80,
          ## zsym = TRUE,
          symm.scale = TRUE,
          mar = c(15, 0, 0, 110),
          key.offset = c(0.90, 0.2),
          cexRow = 0.9, cexCol = 0.75
        )
      }

      plot_RENDER <- shiny::reactive({
        ##
        F <- cumulativeFCtable()
        shiny::req(F)
        ## maximum rows
        plotFCheatmap(F, maxfc=20, maxgenes=60)
        p <- grDevices::recordPlot()
        p
      })

      plot_RENDER2 <- shiny::reactive({
        ##
        F <- cumulativeFCtable()
        shiny::req(F)
        ## maximum rows
        plotFCheatmap(F, maxfc=40, maxgenes=60)
        p <- grDevices::recordPlot()
        p
      })

      
      PlotModuleServer(
        "plot",
        plotlib = "base",
        func = plot_RENDER,
        func2 = plot_RENDER2,
        csvFunc = cumulativeFCtable,
        pdf.width = 14, pdf.height = 5.5,
        res = c(90, 90),
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
