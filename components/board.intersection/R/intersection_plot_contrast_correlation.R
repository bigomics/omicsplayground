##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

contrast_correlation_ui <- function(id, label = "", height = c(600, 800)) {
  ns <- shiny::NS(id)

  info_text <- "<strong>Constrast heatmap.</strong> Similarity of the contrasts visualized as a clustered heatmap. Contrasts that are similar will be clustered close together. The numeric value in the cells correspond to the Pearson correlation coefficient between contrast signatures. Red corresponds to positive correlation and blue to negative correlation."

  ctcorrplot.opts <- shiny::tagList(
    ## tipify( shiny::checkboxInput(ns('ctcorrplot_showrho'), "show correlation values", FALSE),
    ## "Show correlation values in cells."),
    withTooltip(
      shiny::checkboxInput(ns("ctcorrplot_allfc"), "show all contrasts", TRUE),
      "Show all contrasts or just the selected ones."
    ),
    ## tipify( shiny::checkboxInput('ctcorrplot_fixed', "fix heatmap", FALSE),
    ##       "Fix heatmap layout when changing number of top genes"),
    withTooltip(
      shiny::radioButtons(ns("ctcorrplot_ntop"), "number of top genes",
        c("100", "1000", "all"),
        selected = "1000", inline = TRUE
      ),
      "Number of top genes to compute correlation values."
    )
  )

  PlotModuleUI(
    ns("ctcorrplot"),
    title = "Contrast correlation",
    label = "b",
    plotlib = "plotly",
    info.text = info_text,
    options = ctcorrplot.opts,
    download.fmt = c("png", "pdf", "csv"),
    height = c(550, 720),
    width = c("auto", 1100)
  )
}


contrast_correlation_server <- function(id,
                                        getFoldChangeMatrix,
                                        pgx,
                                        input_comparisons,
                                        watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx)

      res <- getFoldChangeMatrix()
      if (is.null(res)) {
        return(NULL)
      }
      if (NCOL(res$fc) < 2) {
        return(NULL)
      }

      fc0 <- res$fc
      qv0 <- res$qv

      ntop <- 2000
      ntop <- input$ctcorrplot_ntop
      if (ntop == "all") ntop <- 999999
      ntop <- as.integer(ntop)

      allfc <- input$ctcorrplot_allfc
      if (!allfc) {
        comp <- input_comparisons()
        if (length(comp) < 2) {
          return(NULL)
        }
        kk <- match(comp, colnames(fc0))
        fc0 <- fc0[, kk, drop = FALSE]
      }

      ## R.full <- cor(fc0[,], use="pairwise", method="spearman")
      R.full <- cor(apply(fc0, 2, rank), use = "pairwise")
      jj <- head(order(-rowMeans(fc0**2)), ntop)
      ## R <- cor(fc0[jj,], use="pairwise", method="spearman")
      R <- cor(apply(fc0[jj, ], 2, rank), use = "pairwise")
      R <- round(R, digits = 2)
      R
    })

    ctcorrplot.PLOTLY <- function() {
      R <- plot_data()
      col <- BLUERED(16)
      col <- gplots::colorpanel(64, "royalblue3", "grey90", "indianred3")
      ## col <- tail(BLUERED(16),8)
      if (min(R, na.rm = TRUE) >= 0) col <- tail(col, 32)
      if (max(R, na.rm = TRUE) <= 0) col <- head(col, 32)

      bluered.pal <- colorRampPalette(colors = c("royalblue3", "grey90", "indianred3"))
      cellnote <- NULL
      ## if(input$ctcorrplot_showrho) cellnote <- R

      plt <- heatmaply::heatmaply(
        R,
        margins = c(250, 200, NA, 0),
        ## k_col = 5, k_row = 5,
        cellnote = cellnote, cellnote_size = 11,
        cellnote_textposition = "middle center",
        colors = bluered.pal,
        limits = c(-1, 1)
      )
      plt
    }

    PlotModuleServer(
      "ctcorrplot",
      func = ctcorrplot.PLOTLY,
      csvFunc = plot_data,
      plotlib = "plotly",
      res = c(80, 85),
      pdf.width = 5, pdf.height = 5,
      add.watermark = watermark
    )
  })
}

# OLD PLOTING FUNCTION

# ctcorrplot.PLOT <- shiny::reactive({
#
#     shiny::req(pgx)
#     shiny::req(input$comparisons)
#
#     ## res <- pgx.getMetaFoldChangeMatrix(pgx, what="meta")
#     res <- getFoldChangeMatrix()
#
#     if(is.null(res)) return(NULL)
#     ##validate(shiny::need(NCOL(res$fc)<2, "warning. need multiple comparisons."))
#     if(NCOL(res$fc)<2) return(NULL)
#
#     fc0 = res$fc
#     qv0 = res$qv
#
#     ntop = 9999
#     ntop <- input$ctcorrplot_ntop
#     if(ntop=="all") ntop <- 999999
#     ntop <- as.integer(ntop)
#
#     allfc <- input$ctcorrplot_allfc
#     if(!allfc) {
#         comp = input_comparisons()
#         if(length(comp)<2) return(NULL)
#         kk = match(comp, colnames(fc0))
#         fc0 <- fc0[,kk,drop=FALSE]
#     }
#
#     ##R.full <- cor(fc0[,], use="pairwise", method="spearman")
#     R.full <- cor(apply(fc0,2,rank), use="pairwise")
#     jj <- head(order(-rowMeans(fc0**2)),ntop)
#     ##R <- cor(fc0[jj,], use="pairwise", method="spearman")
#     R <- cor(apply(fc0[jj,],2,rank), use="pairwise")
#     R <- round(R,digits=2)
#     diag(R) <- 0
#
#     notecex=0.001
#     notecex=1.1; cex=1.3
#     if( nrow(R) > 8)  {notecex=0.95; cex=1.2}
#     if( nrow(R) > 20) {notecex=0.75; cex=1.05}
#     if( nrow(R) > 50) {notecex=0.65; cex=0.9}
#     if( nrow(R) > 80) {notecex=0.0001; cex=0.6}
#
#     mar1 <- c(16,18)*1.2
#     if(nrow(R) <= 8) { mar1=c(16,18)*2 }
#     if(nrow(R) > 30) { mar1=c(16,18)*0.9 }
#     if(nrow(R) > 80) { mar1=c(16,18)*0.6 }
#
#     col <- BLUERED(16)
#     col <- gplots::colorpanel(64,"royalblue3","grey90","indianred3")
#     ##col <- tail(BLUERED(16),8)
#     if(min(R, na.rm=TRUE)>=0) col <- tail(col,32)
#     if(max(R, na.rm=TRUE)<=0) col <- head(col,32)
#     cellnote <- NULL
#
#     col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
#                                "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
#                                "#4393C3", "#2166AC", "#053061"))
#
#
#     corrplot::corrplot(R, method = "circle", order="hclust",
#              is.corr = FALSE,
#              cl.lim = c(-1.02,1.02)*max(abs(R),na.rm=TRUE),
#              col = rev(col2(50)),
#              ## mar = c(1,0.2,0.2,1)*0.2*mean(mar1),
#              mar = c(0,0,0,0),
#              tl.cex = 0.8*cex,
#              tl.col = "black",
#              tl.srt = 90)
#     ##corrplot(R, method = "circle", order="AOE")
#     ##return(R)
# })
