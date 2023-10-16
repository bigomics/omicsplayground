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
    caption,
    height,
    width) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    id = ns("plotmodule"),
    plotlib = "plotly",
    title = title,
    caption = caption,
    label = "a",
    info.text = info.text,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
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
      dbg("[enrichment_plot_top_enrich_gsets_server] reacted!")
      shiny::req(pgx$X)
      rpt <- getFilteredGeneSetTable()
      shiny::req(rpt, gs_contrast())

      comp <- 1
      comp <- gs_contrast()
      if (!(comp %in% names(pgx$gx.meta$meta))) {
        return(NULL)
      }

      ## selected
      sel <- as.integer(gseatable$rows_selected())
      sel.gs <- NULL
      if (!is.null(sel) && length(sel) > 0) sel.gs <- rownames(rpt)[sel]

      ii <- gseatable$rows_selected()
      jj <- gseatable$rows_current()
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
      rnk0 <- rnk0 - mean(rnk0, na.rm = TRUE) #

      # We currently use Human Homolog for GSET 
     
     if(!is.null(pgx$genes$human_ortholog)){
      genes_id <- ifelse(is.na(pgx$genes$human_ortholog), 
                               pgx$genes$gene_name, 
                               pgx$genes$human_ortholog)
     }else{
      genes_id <- pgx$genes$gene_name
     }

     
      # Temporary deal with NAs and duplicates
      rnk0 <- rnk0[!is.na(genes_id)]
      names(rnk0) <- genes_id[!is.na(genes_id)]
      rnk0 <- rnk0[!duplicated(names(rnk0))]

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
        gmt.genes[[gs]] <- toupper(genes)
      }

      dbg("[enrichment_plot_top_enrich_gsets_server] done!")

      res <- list(
        rnk0 = rnk0,
        gmt.genes = gmt.genes,
        qv = qv,
        fx = fx
      )
      return(res)
    })

    plot.RENDER <- function() {
      res <- get_TopEnriched()
      shiny::req(res)

      rnk0 <- res$rnk0
      gmt.genes <- res$gmt.genes
      fc <- res$fc
      qv <- res$qv

      ntop <- length(gmt.genes)
      if (ntop == 1) rowcol <- c(1, 1)
      if (ntop == 12) rowcol <- c(3, 4)

      par(mfrow = rowcol)
      if (ntop == 1) {
        par(mar = c(1, 6, 2, 6), mgp = c(1.6, 0.6, 0), oma = c(0.1, 1, 0, 0.1))
      } else {
        par(mar = c(0.2, 1.8, 2.3, 0.1), mgp = c(1.6, 0.6, 0), oma = c(0.1, 1, 0, 0.1))
      }

      for (i in 1:ntop) {
        gs <- names(gmt.genes)[i]
        genes <- gmt.genes[[i]]
        if (i > ntop || is.na(gs)) {
          frame()
        } else {
          genes <- toupper(genes)
          names(rnk0) <- toupper(names(rnk0))
          ylab <- ""
          if (i %% rowcol[2] == 1) ylab <- "Rank metric"
          xlab <- ""
          gs1 <- playbase::breakstring(gs, 28, 50, force = FALSE)
          if (ntop == 1) {
            gs1 <- playbase::breakstring(gs, 100, 200, force = FALSE)
            xlab <- "Rank in ordered dataset"
            ylab <- "Rank metric"
          }
          playbase::gsea.enplot(rnk0, genes,
            names = NULL, #
            main = gs1, xlab = xlab, ylab = ylab,
            lab.line = c(0, 1.8), cex.lab = 0.75,
            cex.main = 0.78, len.main = 200
          )
          qv1 <- formatC(qv[gs], format = "e", digits = 2)
          legend("topright", paste("q=", qv1), bty = "n", cex = 0.85)
        }
      }
    }

    get_plotly_plots <- function(cex.text) {
      dbg("[enrichment_plot_top_enrich_gsets_server] plotly.RENDER called!")

      res <- get_TopEnriched()
      shiny::req(res)

      rnk0 <- res$rnk0
      gmt.genes <- res$gmt.genes
      fc <- res$fc
      qv <- res$qv

      x.title <- 0.01 * length(rnk0)
      y.title <- max(rnk0)

      ntop <- length(gmt.genes)
      if (ntop == 1) rowcol <- c(1, 1)
      if (ntop == 12) rowcol <- c(3, 4)

      plist <- list()

      for (i in 1:ntop) {
        gset.name <- names(gmt.genes)[i]
        genes <- gmt.genes[[i]]
        genes <- toupper(genes)
        names(rnk0) <- toupper(names(rnk0))
        if (ntop == 1) {
          plt <- playbase::gsea.enplotly(
            rnk0,
            genes,
            main = gset.name,
            xlab = "Rank in ordered dataset",
            ylab = "Rank metric",
            ticklen = 0.25,
            yth = 1, ## threshold for which points get label
            cbar.width = 32,
            tooltips = NULL,
            cex.text = cex.text
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
            cex.text = cex.text
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
