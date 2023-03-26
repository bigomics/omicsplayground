##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##



## Annotate clusters ############

plot_clustpca_ui <- function(id,
                             label = "",
                             height = c(600, 800),
                             parent) {
  ns <- shiny::NS(id)

  info_text <- tagsub(paste0(" The <b>PCA/tSNE</b> panel visualizes unsupervised clustering obtained by the principal components analysis (", a_PCA, ") or t-distributed stochastic embedding (", a_tSNE, ") algorithms. This plot shows the relationship (or similarity) between the samples for visual analytics, where similarity is visualized as proximity of the points. Samples that are ‘similar’ will be placed close to each other.
<br><br>Users can customise the PCA/tSNE plot in the plot settings, including the {color} and {shape} of points using a phenotype class, choose t-SNE or PCA layout, label the points, or display 2D and 3D visualisation of the PCA/tSNE plot."))

  caption <- "<b>PCA/tSNE plot.</b> The plot visualizes the similarity in expression of samples as a scatterplot in reduced dimension (2D or 3D). Samples that are similar are clustered near to each other, while samples with different expression are positioned farther away. Groups of samples with similar profiles will appear as <i>clusters</i> in the plot."


  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::selectInput(parent("hmpca.colvar"), "Color/label:", choices = NULL, width = "100%"),
      "Set colors/labels according to a given phenotype."
    ),
    withTooltip(
      shiny::selectInput(parent("hmpca.shapevar"), "Shape:", choices = NULL, width = "100%"),
      "Set shapes according to a given phenotype."
    ),
    withTooltip(
      shiny::radioButtons(
        ns("hmpca_legend"),
        label = "Legend:",
        choices = c("group label", "bottom"), inline = TRUE
      ),
      "Normalize matrix before calculating distances."
    ),
    withTooltip(
      shiny::checkboxGroupInput(ns("hmpca_options"), "Other:",
        choices = c("sample label", "3D", "normalize"), inline = TRUE
      ),
      "Normalize matrix before calculating distances."
    ),
    withTooltip(
      shiny::radioButtons(ns("hm_clustmethod"), "Layout:",
        c("default", "tsne", "pca", "umap"),
        inline = TRUE
      ),
      "Choose the layout method for clustering to visualise.",
    )
  )

  PlotModuleUI(
    ns("pltmod"),
    title = "PCA/tSNE plot",
    label = label,
    plotlib = "plotly",
    info.text = info_text,
    caption = caption,
    options = plot_opts,
    download.fmt = c("png", "pdf", "csv"),
    width = c("auto", "100%"),
    height = height
  )
}

plot_clustpca_server <- function(id,
                                 pgx,
                                 r.samples = reactive(""),
                                 hmpca.colvar,
                                 hmpca.shapevar,
                                 watermark = FALSE,
                                 parent) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Functions ############

    hm_getClusterPositions <- shiny::reactive({

      ## shiny::req(pgx$tsne2d,pgx$tsne3d,pgx$cluster)

      ## take full matrix
      # flt <- getFilteredMatrix()
      # zx <- flt$zx
      sel.samples <- r.samples()

      clustmethod <- "tsne"
      pdim <- 2
      do3d <- ("3D" %in% input$hmpca_options)
      pdim <- c(2, 3)[1 + 1 * do3d]

      pos <- NULL
      force.compute <- FALSE
      clustmethod <- input$hm_clustmethod
      clustmethod0 <- paste0(clustmethod, pdim, "d")

      if (clustmethod == "default" && !force.compute) {
        if (pdim == 2 && !is.null(pgx$tsne2d)) {
          pos <- pgx$tsne2d[sel.samples, ]
        } else if (pdim == 3 && !is.null(pgx$tsne3d)) {
          pos <- pgx$tsne3d[sel.samples, ]
        }
      } else if (clustmethod0 %in% names(pgx$cluster$pos)) {
        shiny::showNotification(paste("switching to ", clustmethod0, " layout...\n"))
        pos <- pgx$cluster$pos[[clustmethod0]]
        if (pdim == 2) pos <- pos[sel.samples, 1:2]
        if (pdim == 3) pos <- pos[sel.samples, 1:3]
      } else {
        ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ## This should not be necessary anymore as we prefer to
        ## precompute all clusterings.
        shiny::showNotification(paste("computing ", clustmethod, "...\n"))

        ntop <- 1000
        ## ntop = as.integer(input$hm_ntop2)
        zx <- pgx$X
        zx <- zx[order(-apply(zx, 1, sd)), , drop = FALSE] ## OK?
        if (nrow(zx) > ntop) {
          ## zx = head(zx,ntop)  ## OK?
          zx <- zx[1:ntop, , drop = FALSE] ## OK?
        }
        if ("normalize" %in% input$hmpca_options) {
          zx <- scale(t(scale(t(zx))))
        }
        perplexity <- max(1, min((ncol(zx) - 1) / 3, 30))
        perplexity
        res <- pgx.clusterMatrix(
          zx,
          dims = pdim, perplexity = perplexity,
          ntop = 999999, prefix = "C",
          find.clusters = FALSE, kclust = 1,
          row.center = TRUE, row.scale = FALSE,
          method = clustmethod
        )
        if (pdim == 2) pos <- res$pos2d
        if (pdim == 3) pos <- res$pos3d
      }

      pos <- pos[sel.samples, ]
      pos <- scale(pos) ## scale
      ## colnames(pos) = paste0("dim",1:ncol(pos))
      ## rownames(pos) = colnames(zx)

      idx <- NULL

      clust <- list(pos = pos, clust = idx)

      return(clust)
    })

    ## Plot ############

    plot_data <- shiny::reactive({
      clust <- hm_getClusterPositions()
      ## data.frame( x=clust$pos[,1], y=clust$pos[,2], clust=clust$clust )
      return(
        list(
          hmpca_options = input$hmpca_options,
          hmpca.colvar = hmpca.colvar(),
          hmpca.shapevar = hmpca.shapevar(),
          df = data.frame(x = clust$pos[, 1], y = clust$pos[, 2])
        )
      )
    })

    plot.RENDER <- function() {
      pd <- plot_data()
      hmpca_options <- pd[["hmpca_options"]]
      hmpca.colvar <- pd[["hmpca.colvar"]]
      hmpca.shapevar <- pd[["hmpca.shapevar"]]
      pos <- pd[["df"]]


      shiny::req(pgx$Y)
      shiny::req(df)

      do3d <- ("3D" %in% hmpca_options)
      ## clust <- hm_getClusterPositions()
      sel <- rownames(pos)
      df <- cbind(pos, pgx$Y[sel, ])
      # if(!is.null(clust$clust)) df[["<cluster>"]] <- clust$clust

      colvar <- shapevar <- linevar <- textvar <- NULL
      if (hmpca.colvar %in% colnames(df)) colvar <- factor(df[, hmpca.colvar])
      if (hmpca.shapevar %in% colnames(df)) shapevar <- factor(df[, hmpca.shapevar])
      ## if(input$hmpca.line %in% colnames(df)) linevar = factor(df[,input$hmpca.line])
      ## if(input$hmpca.text %in% colnames(df)) textvar = factor(df[,input$hmpca.text])
      mode <- "markers"
      ann.text <- rep(" ", nrow(df))
      if (!do3d && "sample label" %in% hmpca_options) ann.text <- rownames(df)
      if (!is.null(colvar)) {
        colvar <- factor(colvar)
        textvar <- factor(df[, hmpca.colvar])
      }
      symbols <- c(
        "circle", "square", "star", "triangle-up", "triangle-down", "pentagon",
        "bowtie", "hexagon", "asterisk", "hash", "cross", "triangle-left",
        "triangle-right", "+", c(15:0)
      )

      Y <- cbind("sample" = rownames(pos), pgx$Y[sel, ])
      ## tt.info <- paste('Sample:', rownames(df),'</br>Group:', df$group)
      tt.info <- apply(Y, 1, function(y) paste0(colnames(Y), ": ", y, "</br>", collapse = ""))
      tt.info <- as.character(tt.info)
      cex1 <- c(1.0, 0.8, 0.6)[1 + 1 * (nrow(pos) > 30) + 1 * (nrow(pos) > 200)]

      if (do3d) {
        ## 3D plot
        j0 <- 1:nrow(df)
        j1 <- NULL
        if (!is.null(linevar)) {
          linevar <- factor(linevar)
          j0 <- which(linevar == levels(linevar)[1])
          j1 <- which(linevar != levels(linevar)[1])
        }
        plt <- plotly::plot_ly(df, mode = mode) %>%
          plotly::add_markers(
            x = df[j0, 1], y = df[j0, 2], z = df[j0, 3], type = "scatter3d",
            color = colvar[j0], ## size = sizevar, sizes=c(80,140),
            ## marker = list(size = 5*cex1),
            marker = list(size = 5 * cex1, line = list(color = "grey10", width = 0.1)),
            symbol = shapevar[j0], symbols = symbols,
            text = tt.info[j0]
          ) %>%
          plotly::add_annotations(
            x = pos[, 1], y = pos[, 2], z = pos[, 3],
            text = ann.text,
            ## xref = "x", yref = "y",
            showarrow = FALSE
          )
        if (!is.null(j1) & length(j1) > 0) {
          plt <- plt %>% plotly::add_markers(
            x = df[j1, 1], y = df[j1, 2], z = df[j1, 3], type = "scatter3d",
            color = colvar[j1], ## size = sizevar, sizes=c(80,140),
            ## marker = list(size=5*cex1, line=list(color="grey10", width=2)),
            symbol = shapevar[j1], symbols = symbols,
            text = tt.info[j1]
          )
        }
        ## add cluster annotation labels
        if (0 && length(unique(colvar)) > 1) {
          ## add cluster annotation labels
          grp.pos <- apply(pos, 2, function(x) tapply(x, colvar, median))
          ## grp.pos <- matrix(grp.pos, ncol=3)
          cex2 <- ifelse(length(grp.pos) > 20, 0.8, 1)
          plt <- plt %>% plotly::add_annotations(
            x = grp.pos[, 1], y = grp.pos[, 2], z = grp.pos[, 3],
            text = rownames(grp.pos),
            font = list(size = 24 * cex2, color = "#555"),
            showarrow = FALSE
          )
        }
      } else {
        ## 2D plot
        j0 <- 1:nrow(df)
        j1 <- NULL
        if (!is.null(linevar)) {
          linevar <- factor(linevar)
          j0 <- which(linevar == levels(linevar)[1])
          j1 <- which(linevar != levels(linevar)[1])
        }
        plt <- plotly::plot_ly(df, mode = mode) %>%
          plotly::add_markers(
            x = df[j0, 1], y = df[j0, 2], type = "scatter",
            color = colvar[j0], ## size = sizevar, sizes=c(80,140),
            marker = list(size = 16 * cex1, line = list(color = "grey20", width = 0.6)),
            symbol = shapevar[j0], symbols = symbols,
            text = tt.info[j0]
          ) %>%
          plotly::add_annotations(
            x = pos[, 1], y = pos[, 2],
            text = ann.text,
            ## xref = "x", yref = "y",
            showarrow = FALSE
          )

        ## add node labels
        if (!is.null(j1) & length(j1) > 0) {
          plt <- plt %>% plotly::add_markers(
            x = df[j1, 1], y = df[j1, 2], type = "scatter",
            color = colvar[j1], ## size = sizevar, sizes=c(80,140),
            marker = list(size = 16 * cex1, line = list(color = "grey20", width = 1.8)),
            symbol = shapevar[j1], symbols = symbols,
            text = tt.info[j1]
          )
        }

        ## add group/cluster annotation labels
        req(input$hmpca_legend)
        if (input$hmpca_legend == "inside") {
          plt <- plt %>%
            plotly::layout(legend = list(x = 0.05, y = 0.95))
        } else if (input$hmpca_legend == "bottom") {
          plt <- plt %>%
            plotly::layout(legend = list(orientation = "h"))
        } else {
          if (!is.null(textvar) && length(unique(textvar)) > 1) {
            grp.pos <- apply(pos, 2, function(x) tapply(x, as.character(textvar), median))
            cex2 <- 1
            if (length(grp.pos) > 20) cex2 <- 0.8
            if (length(grp.pos) > 50) cex2 <- 0.6
            plt <- plt %>% plotly::add_annotations(
              x = grp.pos[, 1], y = grp.pos[, 2],
              text = paste0("<b>", rownames(grp.pos), "</b>"),
              font = list(size = 24 * cex2, color = "#555"),
              showarrow = FALSE
            )
          }
          plt <- plt %>%
            plotly::layout(showlegend = FALSE)
        }
      }
      title <- paste0("<b>PCA</b>  (", nrow(pos), " samples)")
      if (input$hm_clustmethod == "tsne") title <- paste0("<b>tSNE</b>  (", nrow(pos), " samples)")
      ## plt <- plt %>% plotly::layout(title=title) %>%
      ##     plotly::config(displayModeBar = FALSE)
      plt <- plt %>%
        ## config(displayModeBar = FALSE) %>%
        plotly::config(displayModeBar = TRUE) %>%
        ## config(modeBarButtonsToRemove = all.plotly.buttons ) %>%
        plotly::config(displaylogo = FALSE) %>%
        plotly::config(toImageButtonOptions = list(format = "svg", height = 800, width = 800))
      ## print(plt)
      return(plt)
    }

    modal_plot.RENDER <- function() {
      plot.RENDER()
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      ## plotlib2 = "plotly",
      func = plot.RENDER,
      func2 = modal_plot.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      ## renderFunc = plotly::renderPlotly,
      ## renderFunc2 = plotly::renderPlotly,
      res = c(90, 170), ## resolution of plots
      pdf.width = 8, pdf.height = 8,
      add.watermark = watermark
    )
  })
}
