##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Clustering plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
clustering_plot_splitmap_ui <- function(
    id,
    label = "",
    title,
    caption,
    info.text,
    info.methods,
    info.references,
    info.extra_link,
    height,
    width) {
  ns <- shiny::NS(id)

  splitmap_opts <- shiny::tagList(
    shiny::fillRow(
      height = 50,
      withTooltip(shiny::numericInput(ns("hm_cexRow"), "cexRow:", 1, 0, 1.4, 0.1, width = "100%"),
        "Specify the row label size. Set to 0 to suppress row labels.",
        placement = "right", options = list(container = "body")
      ),
      withTooltip(shiny::numericInput(ns("hm_cexCol"), "cexCol:", 1, 0, 1.4, 0.1, width = "100%"),
        "Specify the column label size. Set to 0 to suppress column labels.",
        placement = "right", options = list(container = "body")
      )
    ),
    shiny::br(),
    withTooltip(
      shiny::checkboxInput(
        ns("hm_legend"), "show legend",
        value = TRUE
      ), "Show or hide the legend.",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    plotlib = c("plotly", "base"),
    info.text = info.text,
    info.methods,
    info.references,
    info.extra_link,
    caption = caption,
    options = splitmap_opts,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
    height = height,
    cards = TRUE,
    card_names = c("dynamic", "static")
  )
}

#' Clustering plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#' @param watermark
#'
#'
#'
#' @export
clustering_plot_splitmap_server <- function(id,
                                            pgx,
                                            getTopMatrix,
                                            selected_phenotypes,
                                            hm_level,
                                            hm_ntop,
                                            hm_scale,
                                            hm_topmode,
                                            hm_clustk,
                                            watermark = FALSE,
                                            labeltype) {
  moduleServer(id, function(input, output, session) {
    fullH <- 850

    ns <- session$ns

    shiny::observeEvent(pgx$Y, {
      if (nrow(pgx$Y) > 100) { # Put cexCol (heatmap) to 0 if more than 100 samples
        shiny::updateNumericInput(session, "hm_cexCol", value = 0)
      } else {
        shiny::updateNumericInput(session, "hm_cexCol", value = 1)
      }
    })

    plot_data <- shiny::reactive({
      ## ComplexHeatmap based splitted heatmap ##########

      filt <- getTopMatrix()
      shiny::req(filt)

      zx <- filt$mat
      annot <- filt$annot
      zx.idx <- filt$idx

      shiny::validate(shiny::need(
        ncol(zx) > 1, "Filtering too restrictive. Please change 'Filter samples' settings."
      ))

      return(list(
        zx = zx,
        annot = annot,
        zx.idx = zx.idx,
        filt = filt
      ))
    })

    base_splitmap.RENDER <- function() {
      ## extract from plot data
      pd <- plot_data()
      zx <- pd[["zx"]]
      annot <- pd[["annot"]]
      zx.idx <- pd[["zx.idx"]]
      filt <- pd[["filt"]]

      if (nrow(zx) <= 1) {
        return(NULL)
      }

      show_rownames <- TRUE
      if (nrow(zx) > 100) show_rownames <- FALSE

      cex1 <- ifelse(ncol(zx) > 50, 0.75, 1)
      cex1 <- ifelse(ncol(zx) > 100, 0.5, cex1)
      cex1 <- ifelse(ncol(zx) > 200, 0, cex1)

      scale.mode <- "none"
      if (hm_scale() == "relative") scale.mode <- "row.center"
      if (hm_scale() == "BMC") scale.mode <- "row.bmc"
      scale.mode

      ## split genes dimension in 5 groups
      splity <- 5
      splity <- 6
      if (!is.null(zx.idx)) splity <- zx.idx

      ## split samples
      splitx <- NULL
      splitx <- filt$grp

      show_legend <- show_colnames <- TRUE
      show_legend <- input$hm_legend
      if (hm_level() == "geneset" || !is.null(splitx)) show_legend <- FALSE

      show_colnames <- (input$hm_cexCol != 0)

      rownames(zx) <- sub("HALLMARK:HALLMARK_", "HALLMARK:", rownames(zx))
      rownames(zx) <- gsub(playdata::GSET_PREFIX_REGEX, "", rownames(zx))
      rownames(zx) <- substring(rownames(zx), 1, 50) ## cut long names...
      if (hm_level() == "gene") {
        ## strip any prefix
        rownames(zx) <- sub(".*:", "", rownames(zx))

        rownames(zx) <- playbase::probe2symbol(rownames(zx), pgx$genes, labeltype(), fill_na = TRUE)
      }

      if (hm_level() == "geneset") rownames(zx) <- tolower(rownames(zx))

      cex2 <- ifelse(nrow(zx) > 60, 0.8, 0.9)
      cex1 <- as.numeric(input$hm_cexCol) * 0.85
      cex2 <- as.numeric(input$hm_cexRow) * 0.75
      shiny::validate(shiny::need(
        all(is.numeric(input$hm_cexCol), is.numeric(input$hm_cexRow)),
        "cexRow and cexCol must not be empty."
      ))
      cex0 <- ifelse(!is.null(splitx) && length(splitx) <= 10, 1.05, 0.85) ## title

      crot <- 0
      totnchar <- nchar(paste0(unique(splitx), collapse = ""))
      totnchar
      nx <- length(unique(splitx))
      if (!is.null(splitx) & (totnchar > 44 || nx >= 6)) crot <- 90

      nrownames <- 60
      nrownames <- 9999
      if (input$hm_cexRow == 0) nrownames <- 0

      shiny::showNotification("Rendering heatmap...")

      playbase::gx.splitmap(
        zx,
        split = splity, splitx = splitx,
        scale = scale.mode, show_legend = show_legend,
        show_colnames = show_colnames, column_title_rot = crot,
        column_names_rot = 45,
        show_rownames = nrownames, rownames_width = 40,
        softmax = 0,
        na_col = "green", na_text = "x",
        title_cex = cex0, cexCol = cex1, cexRow = cex2,
        col.annot = annot, row.annot = NULL, annot.ht = 2.3,
        key.offset = c(0.89, 1.01),
        main = " ", nmax = -1, mar = c(8, 16)
      )
      p <- grDevices::recordPlot()
      p

      # plt
    }

    plotly_splitmap.RENDER_get <- function() {
      ## iHeatmap based splitted heatmap #########

      shiny::req(pgx$genes)

      ## -------------- variable to split samples
      scale <- "none"
      if (hm_scale() == "relative") scale <- "row.center"
      if (hm_scale() == "BMC") scale <- "row.bmc"

      plt <- NULL

      ## extract from plot data
      pd <- plot_data()

      filt <- pd[["filt"]]
      X <- pd[["zx"]]
      annot <- pd[["annot"]]
      splity <- pd[["zx.idx"]]

      ## sample clustering index
      splitx <- NULL
      splitx <- filt$grp

      ## iheatmapr needs factors for sharing between groups
      annotF <- data.frame(as.list(annot), stringsAsFactors = TRUE, check.names = FALSE)
      rownames(annotF) <- rownames(annot)

      sel <- selected_phenotypes()
      sel <- intersect(sel, colnames(annotF))
      if (length(sel) == 0) {
        annotF <- NULL
      } else {
        annotF <- annotF[, sel, drop = FALSE]
      }

      colcex <- as.numeric(input$hm_cexCol)
      rowcex <- as.numeric(input$hm_cexRow)
      shiny::validate(shiny::need(
        all(is.numeric(input$hm_cexCol), is.numeric(input$hm_cexRow)),
        "cexRow and cexCol must not be empty."
      ))

      tooltips <- NULL
      if (hm_level() == "gene") {
        getInfo <- function(g) {
          aa <- paste0(
            "<b>", pgx$genes[g, "gene_name"], "</b>. ",
            pgx$genes[g, "gene_title"], "."
          )
          playbase::breakstring2(aa, 50, brk = "<br>")
        }
        tooltips <- sapply(rownames(X), getInfo)
        labeled_features <- NULL

        rownames(X) <- playbase::probe2symbol(rownames(X), pgx$genes, labeltype(), fill_na = TRUE)
      } else {
        aa <- gsub("_", " ", rownames(X)) ## just geneset names
        tooltips <- sapply(aa, function(x) {
          playbase::breakstring2(x, 50, brk = "<br>")
        })
        names(tooltips) <- rownames(X)
      }
      shiny::showNotification("Rendering iHeatmap...")

      plt <- playbase::pgx.splitHeatmapFromMatrix(
        X = X, annot = annotF, ytips = tooltips,
        idx = splity, splitx = splitx, scale = scale,
        row_annot_width = 0.025, rowcex = rowcex,
        colcex = colcex, show_legend = input$hm_legend,
        ## na_text = 'x',
        return_x_matrix = TRUE
      )
      return(plt)
    }

    plotly_splitmap.RENDER <- function() {
      plt <- plotly_splitmap.RENDER_get()$plt
      if (any(grepl("Iheatmap", class(plt)))) {
        plt <- plt %>% iheatmapr::to_plotly_list()
        plt <- plotly::as_widget(plt)
      }
      plt <- plt %>%
        plotly::layout(
          margin = list(l = 10, r = 5, t = 5, b = 5)
        )
      return(plt)
    }

    plot_grid <- list(
      list(plotlib = "plotly", func = plotly_splitmap.RENDER, card = 1),
      list(plotlib = "base", func = base_splitmap.RENDER, card = 2)
    )

    plot_data_csv <- function() {
      plotly_splitmap.RENDER_get()$X
    }

    lapply(plot_grid, function(x) {
      PlotModuleServer(
        "pltmod",
        plotlib = x$plotlib,
        func = x$func,
        csvFunc = plot_data_csv,
        res = c(80, 95), # resolution of plots
        pdf.width = 10,
        pdf.height = 8,
        add.watermark = watermark,
        card = x$card
      )
    })

    # return(list(
    #   hm_ntop = shiny::reactive(input$hm_ntop),
    #   hm_scale = shiny::reactive(input$hm_scale),
    #   hm_topmode = shiny::reactive(input$hm_topmode),
    #   hm_clustk = shiny::reactive(input$hm_clustk)
    # ))
  }) ## end of moduleServer
}
