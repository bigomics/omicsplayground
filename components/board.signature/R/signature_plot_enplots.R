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
signature_plot_enplots_ui <- function(
    id,
    title,
    info.text,
    info.methods,
    info.references,
    info.extra_link,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)
  info_text <- "<b>Enrichment plots.</b> Enrichment of the query signature in all constrasts. Positive enrichment means that this particular contrast shows similar expression changes as the query signature."

  PlotModuleUI(ns("plotmodule"),
    plotlib = "plotly",
    title = "Enrichment plots",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
    download.fmt = c("png", "pdf"),
    height = height,
    width = width
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
signature_plot_enplots_server <- function(id,
                                          pgx,
                                          sigCalculateGSEA,
                                          enrichmentContrastTable,
                                          watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      # get_plots <- function() {
      shiny::req(pgx$X)
      alertDataLoaded(session, pgx)
      if (is.null(pgx)) {
        return(NULL)
      }

      gsea <- sigCalculateGSEA()
      if (is.null(gsea)) {
        return(NULL)
      }

      ## filter with table selection/search
      ii <- enrichmentContrastTable$rows_selected()
      if (is.null(ii)) {
        ii <- enrichmentContrastTable$rows_all()
      }

      shiny::req(ii)

      ct <- rownames(gsea$output)[ii]
      F <- as.matrix(gsea$F[, ct, drop = FALSE])
      qv <- gsea$output[ct, "q"]
      gset <- gsea$gset

      return(list(
        gset = gset,
        F = F,
        qv = qv,
        ct = ct
      ))
    })

    get_plots <- function(cex = 1) {
      #
      pd <- plot_data()
      shiny::req(pd)

      F <- pd[["F"]]
      gset <- pd[["gset"]]
      qv <- pd[["qv"]]


      if (ncol(F) == 1) {
        nc <- 1
      }
      if (ncol(F) > 1) {
        nc <- 3
        par(mfrow = c(4, 3), mar = c(0.3, 3, 3, 0.5), mgp = c(1.9, 0.7, 0), oma = c(0, 1, 0, 0))
      }
      if (ncol(F) > 12) {
        par(mfrow = c(5, 4), mar = c(0.2, 2, 3, 0.6))
        nc <- 4
      }
      plt <- list()
      for (i in 1:min(20, ncol(F))) {
        f <- colnames(F)[i]
        tt <- sub(".*\\]", "", f)
        tt <- playbase::breakstring(substring(tt, 1, 50), 28, force = TRUE)
        ylab <- ""
        if (i %% nc == 1) ylab <- "rank metric"
        anntitle <- function(tt) {
          list(
            x = 0.5, y = 1.0,
            xref = "paper", yref = "paper",
            xanchor = "center", yanchor = "bottom",
            text = tt, font = list(size = 10 * 1.33),
            align = "center", showarrow = FALSE
          )
        }
        rownames(F)[is.na(rownames(F))] <- "NA"
        p <- playbase::gsea.enplotly(
          F[, i],
          gset,
          main = "",
          cex.text = 1,
          xlab = "",
          ylab = ylab
        ) %>%
          plotly::layout(
            showlegend = TRUE,
            plot_bgcolor = "#f8f8f8",
            margin = list(0, 0, 0, 0),
            annotations = anntitle(tt)
          )
        plt[[i]] <- p
      }
      return(plt)
      #
      ##      p
    } ## )


    plotly.RENDER <- function() {
      plt <- get_plots(cex = 0.9)
      shiny::req(plt)
      ## layout
      nc <- 2
      nr <- ceiling(length(plt) / nc)
      if (length(plt) > 4) nr <- 3
      if (length(plt) > 6) nr <- 4
      if (length(plt) > 12) nr <- 5
      fig <- plotly::subplot(
        plt,
        nrows = nr,
        shareX = TRUE,
        shareY = TRUE,
        titleX = TRUE, titleY = TRUE,
        margin = c(0.01, 0.01, 0.01, 0.045)
      ) %>%
        plotly_default() %>%
        plotly::layout(
          margin = list(l = 10, r = 10, b = 10, t = 20) # lrbt
        )
      return(fig)
    }


    plotly_modal.RENDER <- function() {
      plt <- get_plots(cex = 1.3)
      ## layout
      nc <- 2
      if (length(plt) > 4) nc <- 3
      if (length(plt) > 6) nc <- 4
      if (length(plt) > 12) nc <- 5
      nr <- ceiling(length(plt) / nc)
      fig <- plotly::subplot(
        plt,
        nrows = nr,
        titleX = TRUE,
        titleY = TRUE,
        margin = c(0.01, 0.01, 0.01, 0.045)
      ) %>%
        plotly_modal_default() %>%
        plotly::layout(
          margin = list(l = 20, r = 20, b = 20, t = 40) # lrbt
        )

      return(fig)
    }



    PlotModuleServer(
      "plotmodule",
      func = plotly.RENDER,
      func2 = plotly_modal.RENDER,
      plotlib = "plotly",
      res = c(90, 130), ## resolution of plots
      pdf.width = 8,
      pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
