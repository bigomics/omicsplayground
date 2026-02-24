##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

clustering_plot_clustpca_ui <- function(
  id,
  label = "",
  height,
  width,
  title,
  info.text,
  info.methods,
  info.references,
  info.extra_link,
  caption,
  parent
) {
  ns <- shiny::NS(id)

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::selectInput(ns("hmpca.colvar"), "Color/label:", choices = NULL, width = "100%"),
      "Set colors/labels according to a given phenotype."
    ),
    withTooltip(
      shiny::selectInput(ns("hmpca.shapevar"), "Shape:", choices = NULL, width = "100%"),
      "Set shapes according to a given phenotype."
    ),
    withTooltip(
      shiny::selectInput(
        ns("pca_label"),
        label = "Label:",
        choices = list("group", "bottom", "sample", "<none>")
      ),
      "Place group labels as legend at the bottom or in plot as group or sample labels."
    ),
    withTooltip(
      shiny::checkboxInput(ns("all_clustmethods"), "show all methods"),
      "Show an overview of all dimensionality reduction methods."
    ),
    withTooltip(
      shiny::checkboxInput(ns("plot3d"), "plot 3D"),
      "Show 3D plot."
    )
  )

  quick_buttons <- tagList(
    div(shiny::checkboxInput(ns("plot3d"), "3D"), class = "header-btn")
  )

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
    options = plot_opts,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "clustering"
  )
}

clustering_plot_clustpca_server <- function(id,
                                            pgx,
                                            selected_samples,
                                            clustmethod,
                                            watermark = FALSE,
                                            parent) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Editor: dynamic color pickers for custom palette
    output$custom_palette_ui <- shiny::renderUI({
      shiny::req(input$palette == "custom")
      colvar <- input$hmpca.colvar
      shiny::req(colvar)
      samples <- selected_samples()
      groups <- sort(unique(as.character(pgx$samples[samples, colvar])))
      ## use muted_light as default colors for the pickers
      default_clrs <- rep(omics_pal_d(palette = "muted_light")(8), ceiling(length(groups) / 8))
      pickers <- lapply(seq_along(groups), function(i) {
        colourpicker::colourInput(
          ns(paste0("custom_color_", i)),
          label = groups[i],
          value = default_clrs[i]
        )
      })
      shiny::tagList(pickers)
    })

    plot_data <- shiny::reactive({
      samples <- selected_samples()
      cluster.pos <- pgx$cluster$pos
      for (m in names(cluster.pos)) {
        colnames(cluster.pos[[m]]) <- paste0(m, ".", colnames(cluster.pos[[m]]))
      }
      all.pos <- do.call(cbind, cluster.pos)
      all.pos <- all.pos[samples, ]
      pd <- list(pos = all.pos)
      return(pd)
    })


    create_plot <- function(pgx, pos, pca2d.varexp, method, colvar, shapevar, label, cex, palette = "muted_light", custom_colors = NULL) {
      do3d <- (ncol(pos) == 3)
      sel <- rownames(pos)
      df <- cbind(pos, pgx$samples[sel, , drop = FALSE])

      textvar <- NULL
      if (colvar %in% colnames(df)) {
        colvar <- factor(df[, colvar])
      }

      if (shapevar %in% colnames(df)) {
        shapevar <- factor(df[, shapevar])
      }

      ann.text <- rep(" ", nrow(df))

      label.samples <- (label == "sample")

      if (!do3d && label.samples) ann.text <- rownames(df)
      if (!is.null(colvar)) textvar <- factor(colvar)

      symbols <- c(
        "circle", "square", "star", "triangle-up", "triangle-down", "pentagon",
        "bowtie", "hexagon", "asterisk", "hash", "cross", "triangle-left",
        "triangle-right", "+", c(15:0)
      )

      Y <- cbind("sample" = rownames(pos), pgx$Y[sel, ])
      tt.info <- apply(Y, 1, function(y) paste0(colnames(Y), ": <b>", y, "</b></br>", collapse = ""))
      tt.info <- I(as.character(tt.info))
      cex1 <- c(1.0, 0.8, 0.6)[1 + 1 * (nrow(pos) > 30) + 1 * (nrow(pos) > 200)]
      clrs.length <- length(unique(colvar))
      if (!is.null(custom_colors) && length(custom_colors) >= clrs.length) {
        clrs <- custom_colors[1:clrs.length]
      } else {
        clrs <- rep(omics_pal_d(palette = palette)(8), ceiling(clrs.length / 8))[1:clrs.length]
      }

      if (do3d) {
        plt <- plotly::plot_ly(df, mode = "markers") %>%
          plotly::add_markers(
            x = df[, 1],
            y = df[, 2],
            z = df[, 3],
            type = "scatter3d",
            color = colvar,
            colors = clrs,
            marker = list(
              size = 6 * cex1 * cex,
              line = list(color = "grey10", width = 0.1)
            ),
            symbol = shapevar,
            symbols = symbols,
            text = tt.info
          ) %>%
          plotly::add_annotations(
            x = pos[, 1],
            y = pos[, 2],
            z = pos[, 3],
            text = ann.text,
            showarrow = FALSE
          )
        ## add cluster annotation labels
        if (0 && length(unique(colvar)) > 1) {
          ## add cluster annotation labels
          grp.pos <- apply(pos, 2, function(x) tapply(x, colvar, median))
          cex2 <- ifelse(length(grp.pos) > 20, 0.8, 1)
          plt <- plt %>% plotly::add_annotations(
            x = grp.pos[, 1], y = grp.pos[, 2], z = grp.pos[, 3],
            text = rownames(grp.pos),
            font = list(size = 24 * cex2 * cex, color = "#555"),
            showarrow = FALSE
          )
        }
        if (label == "<none>") {
          plt <- plt %>%
            plotly::layout(showlegend = FALSE)
        }
      } else {
        plt <- plotly::plot_ly(
          df,
          mode = "markers",
          hovertemplate = "</br>%{text}<extra></extra>"
        ) %>%
          plotly::add_markers(
            x = df[, 1],
            y = df[, 2],
            type = "scattergl",
            color = colvar,
            colors = clrs,
            marker = list(
              size = 16 * cex1 * cex,
              line = list(color = "grey20", width = 0.6)
            ),
            symbol = shapevar,
            symbols = symbols,
            text = tt.info
          ) %>%
          plotly::add_annotations(
            x = pos[, 1],
            y = pos[, 2],
            text = ann.text,
            showarrow = FALSE
          )

        if (method == "pca") {
          plt <- plt %>% plotly::layout(
            xaxis = list(title = paste0(toupper(method), "1 (", pca2d.varexp[1], "%)")),
            yaxis = list(title = paste0(toupper(method), "2 (", pca2d.varexp[2], "%)"))
          )
        } else {
          plt <- plt %>% plotly::layout(
            xaxis = list(title = paste0(toupper(method), "1")),
            yaxis = list(title = paste0(toupper(method), "2"))
          )
        }

        ## add group/cluster annotation labels
        if (label == "inside") {
          plt <- plt %>%
            plotly::layout(legend = list(x = 0.05, y = 0.95))
        } else if (label == "bottom") {
          plt <- plt %>%
            plotly::layout(legend = list(orientation = "h"))
        } else if (label == "group") {
          if (!is.null(textvar) && length(unique(textvar)) > 1) {
            grp.pos <- apply(pos, 2, function(x) tapply(x, as.character(textvar), median))
            cex2 <- 1
            if (length(grp.pos) > 20) cex2 <- 0.8
            if (length(grp.pos) > 50) cex2 <- 0.6
            plt <- plt %>% plotly::add_annotations(
              x = grp.pos[, 1],
              y = grp.pos[, 2],
              text = paste0("<b>", rownames(grp.pos), "</b>"),
              font = list(size = 24 * cex2 * cex, color = "#555"),
              showarrow = FALSE
            )
          }
          plt <- plt %>%
            plotly::layout(showlegend = FALSE)
        } else if (label == "sample") {
          plt <- plt %>%
            plotly::layout(showlegend = FALSE)
        } else if (label == "<none>") {
          plt <- plt %>%
            plotly::layout(showlegend = FALSE)
        }
      }
      return(plt)
    }

    create_plotlist <- function() {
      samples <- selected_samples()
      options <- input$hmpca_options
      colvar <- input$hmpca.colvar
      shapevar <- input$hmpca.shapevar
      clustmethod <- clustmethod()
      label <- input$pca_label

      shiny::validate(shiny::need(
        length(samples) > 1,
        "Filtering too restrictive. Please change 'Filter samples' settings."
      ))
      shiny::req(samples, colvar, shapevar, clustmethod, legend)

      methods <- clustmethod()
      if (input$all_clustmethods) {
        cluster.names <- names(pgx$cluster$pos)
        methods <- sub("2d", "", grep("2d", cluster.names, value = TRUE))
      }
      do3d <- (input$plot3d)
      multiplot <- length(methods) > 1

      ## Editor: palette and custom colors
      palette <- if (!is.null(input$palette)) input$palette else "muted_light"
      if (palette %in% c("original", "default")) palette <- "muted_light"
      custom_colors <- NULL
      if (palette == "custom") {
        groups <- sort(unique(as.character(pgx$samples[samples, colvar])))
        custom_colors <- sapply(seq_along(groups), function(j) {
          val <- input[[paste0("custom_color_", j)]]
          if (is.null(val)) omics_pal_d(palette = "muted_light")(8)[(j - 1) %% 8 + 1] else val
        })
      }

      plist <- list()
      for (i in 1:length(methods)) {
        m <- methods[i]
        m1 <- paste0(m, "2d")
        if (do3d) m1 <- paste0(m, "3d")
        pos <- pgx$cluster$pos[[m1]]
        pos <- pos[samples, ]
        pca2d.varexp <- pgx$cluster$pos$pca2d.varexp
        plist[[i]] <- create_plot(
          pgx = pgx,
          pos = pos,
          pca2d.varexp = pca2d.varexp,
          method = m,
          colvar = colvar,
          shapevar = shapevar,
          label = label,
          cex = ifelse(length(methods) > 1, 0.6, 1),
          palette = palette,
          custom_colors = custom_colors
        )
      }
      plist
    }

    plot.RENDER <- reactive({
      if (length(create_plotlist()) == 1) {
        # this is necessary to show axis titles (subplot errases them)
        create_plotlist()[[1]]
      } else {
        plist <- create_plotlist()
        nc <- ceiling(sqrt(length(plist)))
        plotly::subplot(plist, nrows = nc, margin = 0.04)
      }
    })

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plot.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(90, 170), ## resolution of plots
      pdf.width = 8,
      pdf.height = 8,
      add.watermark = watermark,
      parent_session = session
    )
  })
}
