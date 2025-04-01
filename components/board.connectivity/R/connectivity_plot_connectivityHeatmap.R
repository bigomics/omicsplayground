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
connectivity_plot_connectivityHeatmap_ui <- function(
    id,
    title,
    info.text,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(ns("nsig"), "Number of signatures",
        choices = c(10, 20, 40, 100), selected = 20, inline = TRUE
      ),
      "Number of nearest signatures to show."
    ),
    hr(),
    withTooltip(
      shiny::checkboxInput(ns("clusterx"), tspan("Cluster genes")),
      "Cluster genes or sort by expression.."
    ),
    withTooltip(shiny::checkboxInput(ns("cumFCplot_absfc"), "Use absolute foldchange", FALSE),
      "Take the absolute foldchange for calculating the cumulative sum.",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(shiny::checkboxInput(ns("reverse_neg"), "Reverse negative contrasts", TRUE),
      "Reverses negative correlated contrasts so they become positively correlated with respect to the query profile. We simply reverse the sign of the logFC profile. Contrast names (A vs B) are also correctly reversed (B vs A).",
      placement = "right", options = list(container = "body")
    )
  )
  PlotModuleUI(ns("plotmodule"),
    title = title,
    label = label,
    #
    plotlib = "plotly",
    #
    info.text = info.text,
    options = plot_opts,
    download.fmt = c("pdf", "png", "csv", "svg"),
    height = height,
    width = width,
    caption = caption
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
                                                         getProfiles,
                                                         getConnectivityScores,
                                                         getCurrentContrast,
                                                         watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      plot_data <- shiny::reactive({
        F <- getProfiles()
        F[is.na(F)] <- 0

        ## get correlation
        df <- getConnectivityScores()
        rho1 <- df$rho[match(colnames(F), df$pathway)]

        ## add current contrast
        cc <- getCurrentContrast()
        shiny::req(cc)

        fc <- cc$fc[match(rownames(F), names(cc$fc))]
        names(fc) <- rownames(F)

        F <- cbind(fc[rownames(F)], F)
        colnames(F)[1] <- "thisFC"
        colnames(F)[1] <- cc$name
        colnames(F)[1] <- paste("********", cc$name, "********")
        rho2 <- c(1, rho1)
        names(rho2) <- colnames(F)
        if (input$cumFCplot_absfc) {
          F <- abs(F)
        }
        F <- F[order(-rowMeans(F**2, na.rm = TRUE)), , drop = FALSE]

        list(
          F = F,
          score = rho2
        )
      })

      plot_heatmap <- function(F, maxfc, maxgenes = 60) {
        F <- F[, 1:min(NCOL(F), maxfc), drop = FALSE]
        F1 <- head(F, maxgenes)
        par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
        playbase::gx.splitmap(t(F1),
          split = 1,
          cluster_columns = TRUE,
          cluster_rows = TRUE,
          rowlab.maxlen = 80,
          symm.scale = TRUE,
          mar = c(15, 0, 0, 110),
          key.offset = c(0.85, 0.15),
          cexRow = 0.9,
          cexCol = 0.75
        )
      }

      reverse_negative <- function(F, k = 1) {
        ## reverse sign of negative profiles
        fsign <- sign(cor(F, F[, k])[, 1])
        F2 <- t(t(F) * fsign)
        ii <- which(fsign < 0)
        if (length(ii)) {
          reverse_contrast_name <- function(s) {
            prefix <- sub(":.*", "", s)
            if (prefix == s) {
              prefix <- NULL
            } else {
              prefix <- paste0(prefix, ":")
            }
            vs.name <- sub(".*[:]", "", s)
            vs.name2 <- strsplit(vs.name, split = "_vs_")[[1]]
            paste0(prefix, paste(rev(vs.name2), collapse = "_vs_"))
          }
          rev.name <- sapply(colnames(F2)[ii], reverse_contrast_name)
          colnames(F2)[ii] <- rev.name
        }
        F2
      }

      create_iheatmap <- function(F, score, maxfc = 20, maxgenes = 60) {
        sel <- 1:min(NCOL(F), maxfc)
        F <- F[, sel, drop = FALSE]
        score <- score[colnames(F)]
        F <- head(F, maxgenes)
        if (input$reverse_neg) {
          F <- reverse_negative(F, k = 1)
          score <- abs(score)
        }
        ##iheatmap does not like NA
        score[is.na(score)] <- 0
        F[is.na(F)] <- 0        

        ii <- order(rowMeans(F, na.rm = TRUE))
        F <- F[ii,,drop=FALSE ]
        
        plt <- iheatmapr::main_heatmap(
          data = t(F),
          layout = list(margin = list(r = 0))
        ) %>%
          iheatmapr::add_row_clustering() %>%
          iheatmapr::add_row_labels(size = 0.5) %>%
          iheatmapr::add_col_labels()

        if (input$clusterx) {
          plt <- plt %>%
            iheatmapr::add_col_clustering()
        }

        ## add average logFC barplot on top
        avgF <- rowMeans(F, na.rm = TRUE)
        plt <- plt %>%
          iheatmapr::add_col_barplot(
            y = avgF,
            layout = list(title = "avg logFC"),
            buffer = 0.10,
            size = 0.30
          ) %>%
          iheatmapr::add_row_barplot(
            x = score,
            layout = list(title = "similarity"),
            size = 0.10
          )

        plt <- plt %>% iheatmapr::to_plotly_list()
        plt <- plotly::as_widget(plt) %>%
          plotly::layout(
            margin = list(l = 0, r = 0, t = 0, b = 40)
          )
        plt
      }

      plot_RENDER <- function() {
        pd <- plot_data()
        score <- pd$score
        F <- pd$F
        shiny::req(F)
        nsig <- as.numeric(input$nsig)
        create_iheatmap(F, score, maxfc = nsig, maxgenes = 999)
      }

      plot_RENDER2 <- function() {
        pd <- plot_data()
        F <- pd$F
        score <- pd$score
        shiny::req(F)
        nsig <- as.numeric(input$nsig) * 2
        create_iheatmap(F, score, maxfc = nsig, maxgenes = 999)
      }

      PlotModuleServer(
        "plotmodule",
        plotlib = "plotly",
        func = plot_RENDER,
        func2 = plot_RENDER2,
        csvFunc = plot_data,
        pdf.width = 14,
        pdf.height = 5.5,
        res = c(90, 90),
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
