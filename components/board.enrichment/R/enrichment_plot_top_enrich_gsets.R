##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' Create UI for top enriched gene set plots
#'
#' @description
#' Generates the Shiny UI for plotting the top enriched gene sets.
#'
#' @param id Widget id to use for the output.
#' @param title Plot title.
#' @param info.text Info text to display.
#' @param caption Plot caption.
#' @param height Plot height.
#' @param width Plot width.
#'
#' @return Shiny UI for top enriched plots.
enrichment_plot_top_enrich_gsets_ui <- function(
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

  options <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("label_features"),
        "Display labels on the plot",
        TRUE
      ),
      "Display labels on the plot."
    ),
    withTooltip(
      shiny::checkboxInput(
        ns("full_yaxis"),
        "Show full y-axis",
        FALSE
      ),
      "Show full range on y-axis"
    )
  )

  PlotModuleUI(
    id = ns("plotmodule"),
    plotlib = "plotly",
    title = title,
    caption = caption,
    label = "a",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    height = height,
    width = width,
    options = options,
    download.fmt = c("png", "pdf", "svg")
  )
}


#' Plot top enriched gene sets
#'
#' @param id Module ID
#' @param pgx PGX object
#' @param getFilteredGeneSetTable Function to get filtered gene set table
#' @param gs_contrast Gene set contrast
#' @param gseatable Gene set enrichment analysis table
#' @param watermark Add watermark to plot
#'
#' @return Shiny module server function
enrichment_plot_top_enrich_gsets_server <- function(id,
                                                    pgx,
                                                    getFilteredGeneSetTable,
                                                    gs_contrast,
                                                    gseatable,
                                                    watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    get_TopEnriched <- reactive({
      shiny::req(pgx$X)
      rpt <- getFilteredGeneSetTable()

      shiny::req(rpt, gs_contrast())
      comp <- 1
      comp <- gs_contrast()
      if (!(comp %in% names(pgx$gx.meta$meta))) {
        return(NULL)
      }

      ## selected
      ii <- gseatable$rownames_selected()
      jj <- gseatable$rownames_current()
      shiny::req(jj)

      if (nrow(rpt) == 0) {
        return(NULL)
      }

      ## ENPLOT TYPE
      if (length(ii) > 0) {
        itop <- ii[1]
      } else {
        itop <- head(jj, 12)
      }

      # Retreive Enrichment info
      rpt <- rpt[itop, ]
      gx.meta <- pgx$gx.meta$meta[[comp]]
      rnk0 <- gx.meta$meta.fx
      names(rnk0) <- rownames(gx.meta)

      rnk0 <- rnk0[!duplicated(names(rnk0))]

      # if names of rnk0 does nto match rownames(pgx$GMT), then use gene symbols
      # this is needed when collapse by gene is not used
      if (!all(rownames(pgx$GMT) %in% names(rnk0))) {
        names(rnk0) <- pgx$genes[names(rnk0), "symbol"]

        # remove NAs
        rnk0 <- rnk0[!is.na(names(rnk0))]
      }

      fx.col <- grep("score|fx|fc|sign|NES|logFC", colnames(rpt))[1]
      qv.col <- grep("meta.q|q$", colnames(rpt))[1]
      fx <- rpt[, fx.col]
      qv <- rpt[, qv.col]
      names(qv) <- names(fx) <- rownames(rpt)
      top <- rownames(rpt)
      top <- setdiff(top, c(NA, "NA"))
      if (is.null(top) || is.na(top[1])) {
        return(NULL)
      }

      gmt.genes <- list()
      for (i in 1:length(top)) {
        gs <- top[i]
        genes <- names(which(pgx$GMT[, gs] != 0))
        gmt.genes[[gs]] <- genes
      }

      res <- list(
        rnk0 = rnk0,
        gmt.genes = gmt.genes,
        qv = qv,
        fx = fx
      )
      return(res)
    })

    get_plotly_plots <- function(cex.text) {
      res <- get_TopEnriched()
      shiny::req(res)

      rnk0 <- res$rnk0
      # gsea.enplotly cannot deal with duplicated names.
      rnk0 <- rnk0[!duplicated(names(rnk0))]
      gmt.genes <- res$gmt.genes
      fc <- res$fx
      qv <- res$qv

      x.title <- 0.01 * length(rnk0)
      y.title <- 0.96 * max(rnk0)

      ntop <- length(gmt.genes)
      if (ntop == 1) rowcol <- c(1, 1)
      if (ntop == 12) rowcol <- c(3, 4)

      plist <- list()

      for (i in 1:ntop) {
        gset.name <- names(gmt.genes)[i]
        genes <- gmt.genes[[i]]
        if (ntop == 1) {
          plt <- playbase::gsea.enplotly(
            rnk0,
            genes,
            main = gset.name,
            xlab = "Rank in ordered dataset",
            ylab = "Rank metric",
            ticklen = 0.25,
            yth = ifelse(input$label_features, 0.1, 999), ## threshold for which points get label
            yq = ifelse(input$full_yaxis, 0, 0.01), ## limits for y-range
            cbar.width = 32,
            tooltips = NULL,
            cex.text = cex.text,
            cex.title = 1.1,
            cex.axis = 1.3
          ) %>% plotly::layout(
            margin = list(l = 30, r = 10, t = 20, b = 40)
          )
        } else {
          plt <- playbase::gsea.enplotly(
            rnk0,
            genes,
            main = "",
            cex = 0.4,
            xlab = NULL,
            ylab = NULL,
            ticklen = 0.25,
            yth = 999, ## threshold for which points get label
            cbar.width = 15,
            tooltips = NULL,
            cex.text = cex.text,
            cex.axis = 0.8
          ) %>%
            plotly::add_text(
              x = x.title, y = y.title, text = gset.name,
              textfont = list(size = 12 * cex.text),
              textposition = "bottom right"
            ) %>%
            plotly::layout(
              xaxis = list(showticklabels = FALSE),
              margin = list(l = 0, r = 0, t = 80, b = 40)
            )
        }
        plist[[i]] <- plt
      }
      plist
    }

    plotly.RENDER <- function() {
      plist <- get_plotly_plots(cex.text = 0.7)
      ntop <- length(plist)
      if (ntop > 1) {
        plt <- plotly::subplot(plist,
          nrows = 3,
          shareX = TRUE, shareY = TRUE,
          titleX = TRUE, titleY = TRUE,
          margin = c(0.0, 0.0, 0.02, 0.02)
        ) %>%
          plotly::layout(
            margin = list(l = 20, r = 0, t = 20, b = 20),
            xaxis = list(showticklabels = FALSE)
          )
      } else {
        plt <- plist[[1]]
      }
      plt
    }

    plotly.RENDER2 <- function() {
      plist <- get_plotly_plots(cex.text = 1.1)
      ntop <- length(plist)
      if (ntop > 1) {
        plt <- plotly::subplot(plist,
          nrows = 3,
          shareX = TRUE, shareY = TRUE,
          titleX = TRUE, titleY = TRUE,
          margin = c(0.0, 0.0, 0.02, 0.02)
        ) %>%
          plotly::layout(
            margin = list(l = 20, r = 0, t = 20, b = 20),
            xaxis = list(showticklabels = FALSE)
          )
      } else {
        plt <- plist[[1]]
      }
      plt
    }


    PlotModuleServer(
      "plotmodule",
      func = plotly.RENDER,
      func2 = plotly.RENDER2,
      plotlib = "plotly",
      pdf.width = 5,
      pdf.height = 5,
      res = c(90, 120),
      add.watermark = watermark
    )
  })
}
