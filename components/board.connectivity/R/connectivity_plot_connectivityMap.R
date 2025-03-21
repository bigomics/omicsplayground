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
connectivity_plot_connectivityMap_ui <- function(
    id,
    label = "",
    title,
    info.text,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(
        ## "Choose the plot layout: t-SNE, PCA or UMAP",
        ns("layout"), "Layout:", c("pca", "tsne", "volcano"),
        selected = "tsne", inline = TRUE
      ),
      "Choose the plot layout: t-SNE, PCA, or volcano-type",
      placement = "right", options = list(container = "body")
    ),
    hr(),
    withTooltip(shiny::sliderInput(ns("scorethreshold"), "Score threshold:", 0, 1, 0, step = 0.01),
      "Threshold the points by minimum score",
      placement = "right", options = list(container = "body")
    ),
    hr(),
    withTooltip(
      shiny::radioButtons(
        ns("cmapcolorby"), "Color by:", c("score", "dataset"),
        inline = TRUE
      ),
      "Color the points by score or dataset",
      placement = "right", options = list(container = "body")
    ),
    hr(),
    withTooltip(shiny::sliderInput(ns("scoregamma"), "Color gamma:", 0.1, 2, 0.5, step = 0.1),
      "Gamma for color adjustments",
      placement = "right", options = list(container = "body")
    ),
    hr(),
    withTooltip(
      shiny::checkboxGroupInput(
        ns("plotoptions"), "Other options:",
        choiceValues = c("label", "3D", "dark", "large"),
        choiceNames = c(
          "show label", "3D plot", "dark mode",
          "larger points"
        ),
        selected = c()
      ),
      "Show labels, group by dataset, show 3D plot, dark mode.",
      placement = "top", options = list(container = "body")
    )
  )

  PlotModuleUI(ns("plot"),
    title = title,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    options = plot_opts,
    caption = caption,
    download.fmt = c("pdf", "png", "html", "svg"),
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
connectivity_plot_connectivityMap_server <- function(id,
                                                     pgx,
                                                     sigdb,
                                                     getConnectivityScores,
                                                     getEnrichmentMatrix,
                                                     watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      #' Get the XY positions of signature points
      getConnectivityPositions <- shiny::reactive({
        sigdb <- sigdb()
        shiny::req(pgx$X)
        shiny::req(sigdb)

        ## get the foldchanges of selected comparison and neighbourhood
        dims <- 2
        method <- input$layout

        if (method == "volcano") {
          res <- getConnectivityScores()
          pos <- data.frame(rho = res$rho, NES = res$NES)
          rownames(pos) <- res$pathway
          colnames(pos) <- c("rho", "NES")
          return(pos)
        }

        do3d <- "3D" %in% input$plotoptions
        dims <- 2 + 1 * do3d

        if (!grepl("h5$", sigdb)) {
          stop("sigdb must be H5 format!")
        }
        if (!file.exists(sigdb)) {
          warning("*** WARNING *** cannot locate signature matrix file")
          return(NULL)
        }

        h5.file <- sigdb
        rhdf5::h5closeAll()
        rhdf5::h5ls(h5.file)

        pos <- NULL
        if (method == "pca" && dims == 2) try(pos <- rhdf5::h5read(h5.file, "clustering/pca2d"))
        if (method == "pca" && dims == 3) try(pos <- rhdf5::h5read(h5.file, "clustering/pca3d"))
        if (method == "tsne" && dims == 2) try(pos <- rhdf5::h5read(h5.file, "clustering/tsne2d"))
        if (method == "tsne" && dims == 3) try(pos <- rhdf5::h5read(h5.file, "clustering/tsne3d"))
        if (method == "umap" && dims == 2) try(pos <- rhdf5::h5read(h5.file, "clustering/umap2d"))
        if (method == "umap" && dims == 3) try(pos <- rhdf5::h5read(h5.file, "clustering/umap3d"))

        if (is.null(pos)) {
          return(NULL) ## oops...
        }

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
        minscore <- input$scorethreshold
        shiny::req(input$scorethreshold)
        minscore <- min(minscore, 0.999 * max(abs(res$score), na.rm = TRUE))
        sel <- which(abs(res$score) >= minscore)
        if (length(sel) == 0) {
          return(NULL)
        }
        res <- res[sel, , drop = FALSE]

        return(res)
      })

      plot_RENDER <- shiny::reactive({
        #
        shiny::req(pgx$X)

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
        showlabels <- "label" %in% input$plotoptions
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
        if ("dark" %in% input$plotoptions) {
          greyred64 <- colorRampPalette(colors = c("grey15", "grey30", "indianred3", "firebrick4"))(64)
        }

        ## add transparency
        greyred64 <- sapply(1:64, function(i) paste0(greyred64[i], sprintf("%02X", 0 + 4 * i - 1)))
        bluered64 <- sapply(1:64, function(i) paste0(bluered64[i], sprintf("%02X", abs(-259 + 8 * i - 1))))
        colorby <- "score"
        colorby <- input$cmapcolorby
        sizevar <- (df$score)**2
        colorpal <- NULL
        marker.col <- NULL

        if (colorby == "dataset") {
          ## dataset
          colorvar <- dset ## [GSE1234-xxx]
          marker.col <- list()
          colorpal <- rep(RColorBrewer::brewer.pal(8, "Set2"), 99)
        } else if (colorby == "score") {
          gamma1 <- input$scoregamma
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
          if ("large" %in% input$plotoptions) sizeref <- 0.25 * sizeref

          plt <- plotly::plot_ly(
            df,
            x = df[, 1],
            y = df[, 2],
            z = df[, 3],
            mode = "markers",
            type = "scatter3d",
            color = colorvar,
            colors = colorpal,
            marker = list(
              #              size = sizevar,
              #              sizes = c(5, 35) * cex1,
              #              color = marker.col,
              sizeref = sizeref,
              line = list(color = "grey50", width = 0, opacity = 0.5)
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
          if ("large" %in% input$plotoptions) sizeref <- 0.45 * sizeref

          plt <- plotly::plot_ly(
            df,
            x = df[, 1],
            y = df[, 2],
            source = "cmap2d",
            key = rownames(df),
            mode = "markers",
            type = "scattergl",
            color = colorvar,
            colors = colorpal,
            marker = list(
              #              size = sizevar,
              #              sizes = c(5, 30) * cex1,
              #              color = marker.col,
              sizeref = sizeref,
              line = list(color = "grey50", width = 0, opacity = 0.5)
            ),
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

        plt <- plt %>%
          plotly::layout(
            margin = list(l = 0, r = 0, b = 0, t = 0)
          )

        if ("dark" %in% input$plotoptions) {
          plt <- playbase::darkmode(plt, dim = ncol(pos))
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
          plotly::layout(
            showlegend = FALSE,
            margin = list(l = 0, r = 0, b = 0, t = 0)
          ) %>%
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
