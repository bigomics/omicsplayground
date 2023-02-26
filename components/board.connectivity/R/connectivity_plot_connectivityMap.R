##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
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
connectivity_plot_connectivityMap_ui <- function(id,
                                                 label = "",
                                                 fullH = 750) {
  ns <- shiny::NS(id)
  info_text <- strwrap(
    "<b>The Connectivity Map</b> shows the similarity of the contrasts profiles
    as a t-SNE plot. Contrasts that are similar will be clustered close together,
    contrasts that are different are placed farther away."
  )

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(
        ## "Choose the plot layout: t-SNE, PCA or UMAP",
        ns("cmap_layout"), "Layout:", c("pca", "tsne", "volcano"),
        selected = "tsne", inline = TRUE
      ),
      "Choose the plot layout: t-SNE, PCA, or volcano-type",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(shiny::sliderInput(ns("cmap_scorethreshold"), "Score threshold:", 0, 1, 0, step = 0.01),
      "Threshold the points by minimum score",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(
      shiny::radioButtons(
        ns("cmap_cmapcolorby"), "Color by:", c("score", "dataset", "hallmark"),
        inline = TRUE
      ), "Color the points by score, dataset or hallmark",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(shiny::sliderInput(ns("cmap_scoregamma"), "Color gamma:", 0.1, 2, 0.5, step = 0.1),
      "Gamma for color adjustments",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(
      shiny::checkboxGroupInput(
        ns("cmap_plotoptions"), "Other options:",
        choiceValues = c("label", "grouped", "3D", "dark", "large"),
        choiceNames = c(
          "show label", "group by dataset", "3D plot", "dark mode",
          "larger points"
        ),
        selected = c("label", "3D")
      ),
      "Show labels, group by dataset, show 3D plot, dark mode.",
      placement = "top", options = list(container = "body")
    )
  )

  PlotModuleUI(ns("plot"),
    title = "Connectivity map",
    label = label,
    plotlib = "plotly",
    info.text = info_text,
    options = plot_opts,
    download.fmt = c("pdf", "png", "html"),
    height = c(fullH - 100, 750), width = c("auto", 1000)
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
connectivity_plot_connectivityMap_server <- function(id,
                                                     inputData,
                                                     cmap_sigdb,
                                                     getConnectivityScores,
                                                     getEnrichmentMatrix,
                                                     watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      getConnectivityPositions <- shiny::reactive({
        ngs <- inputData()
        cmap_sigdb <- cmap_sigdb()
        shiny::req(ngs)

        ## get the foldchanges of selected comparison and neighbourhood
        dims <- 2
        method <- input$cmap_layout

        if (method == "volcano") {
          res <- getConnectivityScores()
          pos <- data.frame(rho = res$rho, NES = res$NES)
          rownames(pos) <- res$pathway
          colnames(pos) <- c("rho", "NES")
          return(pos)
        }

        do3d <- "3D" %in% input$cmap_plotoptions
        dims <- 2 + 1 * do3d

        sigdb <- cmap_sigdb
        shiny::req(sigdb)
        h5.ref <- grepl("h5$", sigdb)
        h5.file <- NULL
        pos <- NULL

        if (!h5.ref) {
          stop("sigdb must be H5 format!")
        }

        h5.file <- "/home/kwee/Playground/omicsplayground/libx/sigdb-archs4.h5"
        db.exists <- sapply(SIGDB.DIR, function(d) file.exists(file.path(d, sigdb)))
        if (!any(db.exists)) {
          warning("*** WARNING *** cannot locate signature matrix file")
          return(NULL)
        }
        db.dir <- names(which(db.exists))[1]
        h5.file <- file.path(db.dir, sigdb)

        rhdf5::h5closeAll()
        rhdf5::h5ls(h5.file)

        pos <- NULL
        if (method == "pca" && dims == 2) try(pos <- rhdf5::h5read(h5.file, "clustering/pca2d"))
        if (method == "pca" && dims == 3) try(pos <- rhdf5::h5read(h5.file, "clustering/pca3d"))
        if (method == "tsne" && dims == 2) try(pos <- rhdf5::h5read(h5.file, "clustering/tsne2d"))
        if (method == "tsne" && dims == 3) try(pos <- rhdf5::h5read(h5.file, "clustering/tsne3d"))
        if (method == "umap" && dims == 2) try(pos <- rhdf5::h5read(h5.file, "clustering/umap2d"))
        if (method == "umap" && dims == 3) try(pos <- rhdf5::h5read(h5.file, "clustering/umap3d"))

        ## normalize
        pos <- scale(pos)

        ## set row/colnames
        cn <- rhdf5::h5read(h5.file, "data/colnames")
        rownames(pos) <- cn
        xyz <- c("x", "y", "z")
        colnames(pos) <- paste0(toupper(method), "-", xyz[1:ncol(pos)])

        return(pos)
      })

      getThresholdedConnectivityScores <- shiny::reactive({
        res <- getConnectivityScores()
        if (is.null(res)) {
          return(NULL)
        }

        ## add all points (not within score table)
        pos0 <- getConnectivityPositions() ## all points
        pp <- rownames(pos0)
        res <- res[match(pp, rownames(res)), , drop = FALSE]
        rownames(res) <- rownames(pos0)
        res$pathway <- rownames(pos0)
        res$score[is.na(res$score)] <- 0
        res$rho[is.na(res$rho)] <- 0

        ## select on minimum score
        minscore <- input$cmap_scorethreshold
        shiny::req(input$cmap_scorethreshold)
        minscore <- min(minscore, 0.999 * max(abs(res$score), na.rm = TRUE))
        sel <- which(abs(res$score) >= minscore)
        if (length(sel) == 0) {
          return(NULL)
        }
        res <- res[sel, , drop = FALSE]

        ## rownames(res) <- res$pathway
        cmap_grouped <- "grouped" %in% input$cmap_plotoptions
        if (cmap_grouped) {
          res <- res[order(-res$score), ]
          dset <- sub("].*", "]", res$pathway)
          names(dset) <- res$pathway
          idx <- which(!duplicated(dset))
          resx <- res[idx, , drop = FALSE]
          tt <- table(dset)[dset[resx$pathway]]
          resx$pathway <- paste0(resx$pathway, " (+", tt, " contrasts)")
          resx$datasets <- sub("].*", "]", resx$pathway)
          rownames(resx) <- resx$datasets
          res <- resx
        }

        return(res)
      })

      plot_RENDER <- shiny::reactive({
        pgx <- inputData()
        cmap_sigdb <- cmap_sigdb()
        shiny::req(pgx)

        ## get positions
        res0 <- getThresholdedConnectivityScores()
        if (is.null(res0)) {
          return(NULL)
        }
        res0 <- res0[, c("pathway", "score", "rho")]
        pos0 <- getConnectivityPositions()
        if (is.null(pos0)) {
          return(NULL)
        }
        if (is.null(res0)) {
          return(NULL)
        }

        cmap_grouped <- "grouped" %in% input$cmap_plotoptions
        if (cmap_grouped) {
          pset <- sub("].*", "]", rownames(pos0))
          pos0 <- apply(pos0, 2, function(x) tapply(x, pset, mean))
        }
        if (is.null(res0) || nrow(res0) == 0) {
          return(NULL)
        }
        if (is.null(pos0) || nrow(pos0) == 0) {
          return(NULL)
        }

        pp <- intersect(rownames(pos0), rownames(res0))
        res <- res0[match(pp, rownames(res0)), , drop = FALSE]
        pos <- pos0[match(pp, rownames(pos0)), , drop = FALSE]

        ## draw high score last
        jj <- order(res$score)
        res <- res[jj, , drop = FALSE]
        pos <- pos[jj, , drop = FALSE]

        ## parse groups/dataset from name
        jj <- grep("\\]", rownames(pos))
        dset <- rep("[this_data]", nrow(pos))
        pw <- res$pathway ## full name (may be modified by collapsing)
        dset[jj] <- sub("].*", "]", pw[jj])
        ct.name <- sub(".*\\]", "", pw)
        table(dset)

        ## make dataframe
        score <- rho <- NULL
        df <- data.frame(pos,
          dataset = dset, contrast = ct.name, name = res$pathway,
          score = res$score, rho = res$rho, check.names = FALSE
        )
        rownames(df) <- rownames(pos)

        do3d <- ncol(pos) == 3
        mode <- "markers"
        ann.text <- rep(NA, nrow(df))
        showlabels <- "label" %in% input$cmap_plotoptions
        if (showlabels) ann.text <- df$name

        tt.info <- paste(
          "Contrast:", df$contrast,
          "</br>Dataset:", df$dataset,
          "</br>Similarity score:", round(df$score, 3),
          "</br>Correlation:", round(df$rho, 3)
        )
        cex1 <- c(1.0, 0.8, 0.6, 0.4)[1 + 1 * (nrow(pos) > 30) + 1 * (nrow(pos) > 200) +
          1 * (nrow(pos) > 500)]
        this.data <- 1 * grepl("this_data", dset)
        shapevar <- 1 + 1 * this.data
        symbols <- c("circle", "x")

        bluered64 <- colorRampPalette(
          colors = c("navyblue", "royalblue4", "grey90", "indianred3", "firebrick4")
        )(64)
        greyred64 <- colorRampPalette(colors = c("grey85", "grey70", "indianred3", "firebrick4"))(64)
        if ("dark" %in% input$cmap_plotoptions) {
          greyred64 <- colorRampPalette(colors = c("grey15", "grey30", "indianred3", "firebrick4"))(64)
        }

        ## add transparency
        greyred64 <- sapply(1:64, function(i) paste0(greyred64[i], sprintf("%02X", 0 + 4 * i - 1)))
        bluered64 <- sapply(1:64, function(i) paste0(bluered64[i], sprintf("%02X", abs(-259 + 8 * i - 1))))
        colorby <- "score"
        colorby <- input$cmap_cmapcolorby
        sizevar <- (df$score)**2
        colorpal <- NULL
        marker.col <- NULL

        if (colorby == "dataset") {
          ## dataset
          colorvar <- dset ## [GSE1234-xxx]
          marker.col <- list()
          colorpal <- rep(RColorBrewer::brewer.pal(8, "Set2"), 99)
        } else if (colorby == "hallmark") {
          Y <- t(getEnrichmentMatrix(cmap_sigdb, nc = 15))
          Y <- Y[match(rownames(df), rownames(Y)), ]
          colorvar <- colnames(Y)[max.col(Y)]
          colorvar <- shortstring(colorvar, 40) ## too long!!!
          marker.col <- list()
          colorpal <- rep(RColorBrewer::brewer.pal(8, "Set2"), 99)
        } else if (colorby == "score") {
          gamma1 <- input$cmap_scoregamma
          score1 <- sign(df$score) * (abs(df$score) / max(abs(df$score), na.rm = TRUE))**gamma1
          colorvar <- score1
          colorpal <- greyred64
          cmin <- 0
          if (min(score1) < 0) {
            colorpal <- bluered64
            cmin <- -1
          }
          marker.col <- list(
            colorscale = colorpal,
            cmin = cmin, cmax = 1, colorbar = list(title = "score", len = 0.33),
            opacity = 0.966, reversescale = FALSE
          )
        } else {
          stop("[connectivityMap.RENDER] FATAL:: invalid colorby") ## should not come here
        }

        if (do3d) {
          ## 3D plot
          sizeref <- 0.06 * max(1, nrow(df) / 1000)**0.33
          sizeref
          if ("large" %in% input$cmap_plotoptions) sizeref <- 0.25 * sizeref

          plt <- plotly::plot_ly(
            df,
            x = df[, 1], y = df[, 2], z = df[, 3],
            mode = "markers", type = "scatter3d",
            color = colorvar, colors = colorpal,
            size = sizevar, sizes = c(5, 35) * cex1,
            marker = c(
              marker.col,
              list(
                sizeref = sizeref,
                line = list(color = "grey50", width = 0, opacity = 0.5)
              )
            ),
            text = tt.info
          )
          plt <- plt %>%
            plotly::layout(
              xaxis = list(title = colnames(pos)[1]),
              yaxis = list(title = colnames(pos)[2]),
              zaxis = list(title = colnames(pos)[3])
            )
        } else {
          ## 2D plot
          ##
          sizeref <- 0.08 * max(1, nrow(df) / 1000)**0.33
          sizeref
          if ("large" %in% input$cmap_plotoptions) sizeref <- 0.85 * sizeref

          plt <- plotly::plot_ly(
            df,
            x = df[, 1], y = df[, 2],
            source = "cmap2d", key = rownames(df),
            mode = "markers", type = "scattergl",
            color = colorvar, colors = colorpal,
            marker = c(
              marker.col,
              list(
                sizeref = sizeref,
                line = list(color = "grey50", width = 0, opacity = 0.5)
              )
            ),
            size = sizevar, sizes = c(5, 30) * cex1,
            text = tt.info
          )

          if (showlabels && nrow(pos) < 100) {
            ## this is slow if not careful. too many annotation labels slows down
            ##
            plt <- plt %>%
              plotly::add_annotations(
                x = pos[, 1], y = pos[, 2],
                text = ann.text,
                ## xref = "x", yref = "y",
                yanchor = "bottom",
                yshift = 2,
                textposition = "top",
                showarrow = FALSE
              )
          }
          plt <- plt %>%
            plotly::layout(
              xaxis = list(title = colnames(pos)[1]),
              yaxis = list(title = colnames(pos)[2])
            )
        }

        ## scale range to cover 95%
        xrng <- quantile(pos[, 1], probs = c(0.025, 0.975))
        yrng <- quantile(pos[, 2], probs = c(0.025, 0.975))
        xrng <- xrng + 0.15 * c(-1, 1) * diff(xrng)
        yrng <- yrng + 0.15 * c(-1, 1) * diff(yrng)

        plt <- plt %>%
          plotly::config(
            toImageButtonOptions =
              list(format = "svg", height = 800, width = 800, scale = 1.1)
          ) %>%
          plotly::event_register("plotly_selected")

        if ("dark" %in% input$cmap_plotoptions) {
          plt <- darkmode(plt, dim = ncol(pos))
        }

        return(plt)
      })

      plot_RENDER1 <- shiny::reactive({
        ## Hide colorbar
        ##
        plt <- plot_RENDER()
        if (is.null(plt)) {
          return(NULL)
        }
        plt <- plt %>%
          plotly::layout(showlegend = FALSE) %>%
          plotly::hide_colorbar()
        return(plt)
      })

      PlotModuleServer(
        "plot",
        plotlib = "plotly",
        func = plot_RENDER1,
        func2 = plot_RENDER,
        csvFunc = getThresholdedConnectivityScores,
        pdf.width = 8, pdf.height = 8,
        res = 90,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
