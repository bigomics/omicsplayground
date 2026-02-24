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
  width
) {
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
    download.fmt = c("png", "pdf", "svg"),
    height = height,
    width = width,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "enrichment"
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
        rownames(F) <- playbase::probe2symbol(rownames(F), pgx$genes, "gene_name", fill_na = TRUE)
        gset <- playbase::probe2symbol(gset, pgx$genes, "gene_name", fill_na = TRUE)
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

        ## Editor: color overrides via plotly_build post-processing
        color_up <- if (!is.null(input$color_up)) input$color_up else "#f23451"
        color_down <- if (!is.null(input$color_down)) input$color_down else "#3181de"
        color_line <- if (!is.null(input$color_line)) input$color_line else "#00EE00"
        colors_changed <- (!is.null(input$color_up) && input$color_up != "#f23451") ||
          (!is.null(input$color_down) && input$color_down != "#3181de") ||
          (!is.null(input$color_line) && input$color_line != "#00EE00")
        if (colors_changed) {
          p <- plotly::plotly_build(p)

          ## Identify colorbar segment indices (wide segments, not the green line)
          cbar_indices <- which(sapply(p$x$data, function(t) {
            !is.null(t$line$width) && t$line$width >= 15 &&
              !is.null(t$line$color) && t$line$color != "#00EE00"
          }))

          ## Build custom colorbar palette
          if (length(cbar_indices) > 0) {
            suppressWarnings(
              custom_cc <- gplots::colorpanel(length(cbar_indices), color_down, "#CCCCCC", color_up)
            )
            for (ci in seq_along(cbar_indices)) {
              p$x$data[[cbar_indices[ci]]]$line$color <- custom_cc[ci]
            }
          }

          ## Override enrichment score line color
          for (j in seq_along(p$x$data)) {
            if (!is.null(p$x$data[[j]]$line$color) && p$x$data[[j]]$line$color == "#00EE00") {
              p$x$data[[j]]$line$color <- color_line
            }
          }
        }

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
      add.watermark = watermark,
      parent_session = session
    )
  }) ## end of moduleServer
}
