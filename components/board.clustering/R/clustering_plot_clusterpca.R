##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
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
    ## arrows only make sense for PCA: loadings and phenotype correlations are
    ## defined against ordered components, which t-SNE/UMAP do not have
    shiny::conditionalPanel(
      "input.hm_clustmethod == 'pca'",
      ns = parent,
      withTooltip(
        shiny::selectInput(ns("pca_arrows"), "Arrows:",
          choices = c("<none>", "features", "phenotypes"), width = "100%"
        ),
        "Overlay a biplot: arrows for the top feature loadings, or for the direction in which each phenotype increases."
      ),
      shiny::conditionalPanel(
        "input.pca_arrows != '<none>'",
        ns = ns,
        withTooltip(
          shiny::numericInput(ns("pca_numarrows"), "Number of arrows:",
            value = 3, min = 1, max = 50, step = 1, width = "100%"
          ),
          "How many arrows to label and draw, ranked by length. Lower numbers keep the plot readable."
        ),
        withTooltip(
          shiny::checkboxInput(ns("pca_repel"), "repel labels", value = TRUE),
          "Move overlapping arrow labels apart. Switch off to pin each label to its own arrow tip."
        )
      ),
      withTooltip(
        shiny::checkboxInput(ns("pca_varplot"), "variance explained"),
        "Show how much variance each component explains, instead of the samples."
      )
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
    plot_type = "clustering_prism"
  )
}

clustering_plot_clustpca_server <- function(id,
                                            pgx,
                                            selected_samples,
                                            clustmethod,
                                            pca_components,
                                            pca_dims,
                                            labeltype = shiny::reactive("feature"),
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
      custom_palette_pickers(groups, ns)
    })

    plot_data <- shiny::reactive({
      samples <- selected_samples()
      ## drop pca2d.varexp and friends: only the position matrices are tabular
      cluster.pos <- Filter(is.matrix, pgx$cluster$pos)
      ## the board recomputes PCA on the fly, so export those instead of the
      ## stored PC1/PC2 only
      cluster.pos <- cluster.pos[grep("^pca", names(cluster.pos), invert = TRUE)]
      cluster.pos$pca <- pca_components()$pos
      for (m in names(cluster.pos)) {
        colnames(cluster.pos[[m]]) <- paste0(m, ".", colnames(cluster.pos[[m]]))
      }
      all.pos <- do.call(cbind, cluster.pos)
      all.pos <- all.pos[samples, ]
      pd <- list(pos = all.pos)
      return(pd)
    })


    ## arrows scaled to fill the score cloud (loadings and correlations live on
    ## a much smaller scale than the scores), with the label just past the tip
    scale_arrows <- function(arrows, pos) {
      s <- 0.65 * max(abs(pos)) / max(abs(arrows))
      ax <- arrows[, 1] * s
      ay <- arrows[, 2] * s
      r <- sqrt(ax^2 + ay^2)
      r[r == 0] <- 1
      nudge <- 0.05 * max(abs(pos))
      list(
        x = ax, y = ay,
        lx = ax + nudge * ax / r, ly = ay + nudge * ay / r,
        text = rownames(arrows)
      )
    }

    create_plot <- function(pgx, pos, xlab, ylab, method, colvar, shapevar, label, cex, palette = "muted_light", custom_colors = NULL, arrows = NULL) {
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

        plt <- plt %>% plotly::layout(
          xaxis = list(title = xlab),
          yaxis = list(title = ylab)
        )

        if (!is.null(arrows) && nrow(arrows) && max(abs(arrows)) > 0) {
          aa <- scale_arrows(arrows, pos)
          ## Arrows are drawn as traces, NOT as annotations: plotly.repel takes
          ## ownership of layout.annotations at render time (plotly-repel.js
          ## relayouts the whole array), which wipes anything annotation-drawn.
          ## Shafts are one NA-separated line trace, heads are markers rotated
          ## along each arrow (angle is the bearing, clockwise from up).
          plt <- plt %>%
            plotly::add_trace(
              x = as.vector(rbind(0, aa$x, NA)),
              y = as.vector(rbind(0, aa$y, NA)),
              type = "scatter", mode = "lines",
              line = list(color = "grey30", width = 1.8),
              showlegend = FALSE, hoverinfo = "skip", inherit = FALSE
            ) %>%
            plotly::add_trace(
              x = aa$x, y = aa$y,
              type = "scatter", mode = "markers",
              marker = list(
                symbol = "arrow", size = 12, color = "grey30",
                angle = atan2(aa$x, aa$y) * 180 / pi
              ),
              text = aa$text, hoverinfo = "text",
              showlegend = FALSE, inherit = FALSE
            )
          ## Arrow tips bunch up, so fixed-offset labels collide. plotly.repel
          ## places them in the browser and re-solves on zoom/pan. Off by
          ## request, or when the optional package is not installed: labels then
          ## sit nudged just past their own tip.
          if (isTRUE(input$pca_repel) && requireNamespace("plotly.repel", quietly = TRUE)) {
            plt <- plotly.repel::add_text_repel(
              plt,
              data = data.frame(x = aa$x, y = aa$y, text = aa$text),
              x = ~x, y = ~y, text = ~text,
              box_padding = 0.4, point_padding = 0.3, max_overlaps = 20,
              font = list(color = "grey20", size = 13 * cex),
              segment = list(color = "grey60", width = 1)
            )
          } else {
            plt <- plt %>% plotly::add_annotations(
              x = aa$lx, y = aa$ly, text = aa$text,
              showarrow = FALSE,
              font = list(color = "grey20", size = 13 * cex),
              bgcolor = "rgba(255,255,255,0.65)", borderpad = 2
            )
          }
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

    create_ggplot <- function(pgx, pos, xlab, ylab, method, colvar, label, cex, palette = "muted_light", custom_colors = NULL, arrows = NULL) {
      sel <- rownames(pos)
      df <- cbind(pos, pgx$samples[sel, , drop = FALSE])

      if (colvar %in% colnames(df)) {
        colvar_fct <- factor(df[, colvar])
      } else {
        colvar_fct <- factor(rep("all", nrow(df)))
      }

      clrs.length <- length(unique(colvar_fct))
      if (!is.null(custom_colors) && length(custom_colors) >= clrs.length) {
        clrs <- custom_colors[1:clrs.length]
      } else {
        clrs <- rep(omics_pal_d(palette = palette)(8), ceiling(clrs.length / 8))[1:clrs.length]
      }

      showlabels <- (label %in% c("group", "inside"))

      p <- playbase::pgx.scatterPlotXY.GGPLOT(
        pos,
        var = colvar_fct,
        col = clrs,
        cex = cex,
        xlab = xlab,
        ylab = ylab,
        title = toupper(method),
        cex.title = 1.2,
        cex.clust = 1.1,
        label.clusters = showlabels,
        legend = !(label %in% c("<none>", "group", "sample"))
      )

      if (!is.null(arrows) && nrow(arrows) && max(abs(arrows)) > 0) {
        aa <- scale_arrows(arrows, pos)
        df.a <- data.frame(x = aa$x, y = aa$y, lx = aa$lx, ly = aa$ly, text = aa$text)
        p <- p +
          ggplot2::geom_segment(
            data = df.a,
            ggplot2::aes(x = 0, y = 0, xend = x, yend = y),
            arrow = ggplot2::arrow(length = ggplot2::unit(0.18, "cm"), type = "closed"),
            colour = "grey30", linewidth = 0.5, inherit.aes = FALSE
          ) +
          if (isTRUE(input$pca_repel)) {
            ggrepel::geom_text_repel(
              data = df.a,
              ggplot2::aes(x = x, y = y, label = text),
              colour = "grey20", size = 3 * cex, inherit.aes = FALSE,
              box.padding = 0.4, point.padding = 0.3,
              segment.colour = "grey60", max.overlaps = 20
            )
          } else {
            ggplot2::geom_text(
              data = df.a,
              ggplot2::aes(x = lx, y = ly, label = text),
              colour = "grey20", size = 3 * cex, inherit.aes = FALSE
            )
          }
      }
      p
    }

    ## Scree data: percentage of variance per component plus its running total.
    ## Both series are percentages, so they share one axis. The percentages come
    ## from pca_components() and are relative to the total variance, so they do
    ## not sum to 100 over the five components we compute -- that is correct.
    varexp_data <- function() {
      v <- pca_components()$varexp
      data.frame(
        pc = factor(paste0("PC", seq_along(v)), levels = paste0("PC", seq_along(v))),
        varexp = as.numeric(v),
        cumulative = cumsum(as.numeric(v))
      )
    }

    create_varplot <- function(palette = "muted_light") {
      df <- varexp_data()
      clr <- omics_pal_d(palette = palette)(8)[1]
      plotly::plot_ly(df) %>%
        plotly::add_bars(
          x = ~pc, y = ~varexp, name = "per component",
          marker = list(color = clr),
          text = ~ paste0(varexp, "%"), textposition = "outside",
          hovertemplate = "%{x}: %{y}%<extra></extra>"
        ) %>%
        plotly::add_trace(
          x = ~pc, y = ~cumulative, name = "cumulative",
          type = "scatter", mode = "lines+markers",
          line = list(color = "grey40"), marker = list(color = "grey40"),
          hovertemplate = "up to %{x}: %{y}%<extra></extra>"
        ) %>%
        plotly::layout(
          xaxis = list(title = ""),
          yaxis = list(title = "variance explained (%)", rangemode = "tozero"),
          showlegend = FALSE
        )
    }

    create_gg_varplot <- function(palette = "muted_light") {
      df <- varexp_data()
      clr <- omics_pal_d(palette = palette)(8)[1]
      ggplot2::ggplot(df, ggplot2::aes(x = pc)) +
        ggplot2::geom_col(ggplot2::aes(y = varexp), fill = clr) +
        ggplot2::geom_text(
          ggplot2::aes(y = varexp, label = paste0(varexp, "%")),
          vjust = -0.5, size = 3
        ) +
        ggplot2::geom_line(ggplot2::aes(y = cumulative, group = 1), colour = "grey40") +
        ggplot2::geom_point(ggplot2::aes(y = cumulative), colour = "grey40") +
        ggplot2::labs(x = NULL, y = "variance explained (%)") +
        ggplot2::expand_limits(y = 0)
    }

    ## Biplot arrow endpoints on the two selected components, in loading or
    ## correlation units. Ranked by length and cut to the top N; the caller
    ## scales them to the score cloud.
    arrow_ends <- function(pc, k) {
      mode <- input$pca_arrows
      if (is.null(mode) || mode == "<none>") {
        return(NULL)
      }
      A <- if (mode == "features") pc$loadings else pc$pheno.cor
      if (is.null(A) || max(k) > ncol(A)) {
        return(NULL)
      }
      A <- A[stats::complete.cases(A[, k, drop = FALSE]), k, drop = FALSE]
      if (!nrow(A)) {
        return(NULL)
      }
      ## free-text numeric: empty or nonsense while typing must not blank the plot
      n <- suppressWarnings(as.integer(input$pca_numarrows))
      if (length(n) != 1 || is.na(n)) n <- 3
      n <- max(1, n)
      A <- A[head(order(-(A[, 1]^2 + A[, 2]^2)), n), , drop = FALSE]
      if (mode == "features") {
        rownames(A) <- playbase::probe2symbol(
          rownames(A), pgx$genes, labeltype(),
          fill_na = TRUE
        )
      }
      A
    }

    ## Positions and axis labels for one layout method. PCA is recomputed on
    ## the fly (pgx only stores PC1-PC3) so any pair of PC1-PC5 can be plotted.
    method_pos <- function(m, samples, do3d = FALSE) {
      arrows <- NULL
      if (m == "pca") {
        pc <- pca_components()
        k <- if (do3d) seq_len(min(3, ncol(pc$pos))) else pca_dims()
        pos <- pc$pos[samples, k, drop = FALSE]
        labs <- paste0("PC", k, " (", pc$varexp[k], "%)")
        if (!do3d) arrows <- arrow_ends(pc, k)
      } else {
        m1 <- paste0(m, ifelse(do3d, "3d", "2d"))
        pos <- pgx$cluster$pos[[m1]][samples, , drop = FALSE]
        labs <- paste0(toupper(m), seq_len(ncol(pos)))
      }
      list(pos = pos, xlab = labs[1], ylab = labs[2], arrows = arrows)
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
      groups <- sort(unique(as.character(pgx$samples[samples, colvar])))
      n_groups <- length(groups)
      custom_colors <- resolve_palette_colors(input, n_groups, fallback_colors = omics_pal_d("muted_light")(n_groups))
      palette <- if (!is.null(input$palette)) input$palette else "muted_light"
      if (palette %in% c("original", "default", "custom", "")) palette <- "muted_light"

      plist <- list()
      for (i in 1:length(methods)) {
        m <- methods[i]
        mp <- method_pos(m, samples, do3d = do3d)
        plist[[i]] <- create_plot(
          pgx = pgx,
          pos = mp$pos,
          xlab = mp$xlab,
          ylab = mp$ylab,
          method = m,
          colvar = colvar,
          shapevar = shapevar,
          label = label,
          cex = ifelse(length(methods) > 1, 0.6, 1),
          palette = palette,
          custom_colors = custom_colors,
          arrows = mp$arrows
        )
      }
      plist
    }

    create_gg_plotlist <- function() {
      samples <- selected_samples()
      colvar <- input$hmpca.colvar
      clustmethod <- clustmethod()
      label <- input$pca_label

      shiny::validate(shiny::need(
        length(samples) > 1,
        "Filtering too restrictive. Please change 'Filter samples' settings."
      ))
      shiny::req(samples, colvar, clustmethod)

      methods <- clustmethod()
      if (input$all_clustmethods) {
        cluster.names <- names(pgx$cluster$pos)
        methods <- sub("2d", "", grep("2d", cluster.names, value = TRUE))
      }

      groups <- sort(unique(as.character(pgx$samples[samples, colvar])))
      n_groups <- length(groups)
      custom_colors <- resolve_palette_colors(input, n_groups, fallback_colors = omics_pal_d("muted_light")(n_groups))
      palette <- if (!is.null(input$palette)) input$palette else "muted_light"
      if (palette %in% c("original", "default", "custom", "")) palette <- "muted_light"

      gp <- extract_ggprism_params(input)

      plts <- list()
      for (i in seq_along(methods)) {
        m <- methods[i]
        mp <- method_pos(m, samples)

        p <- create_ggplot(
          pgx = pgx,
          pos = mp$pos,
          xlab = mp$xlab,
          ylab = mp$ylab,
          method = m,
          colvar = colvar,
          label = label,
          cex = ifelse(length(methods) > 1, 0.6, 1),
          palette = palette,
          custom_colors = custom_colors,
          arrows = mp$arrows
        )
        p <- apply_ggprism_theme(p, gp)
        p <- apply_editor_theme(p, input)
        plts[[length(plts) + 1]] <- p
      }
      plts
    }

    plot.RENDER <- reactive({
      gp <- extract_ggprism_params(input)
      do3d <- isTRUE(input$plot3d)

      ## the variance toggle replaces the scatter in this same panel. Check the
      ## layout too: the checkbox keeps its value while hidden, so switching to
      ## tsne/umap must not leave a PCA scree plot on screen.
      if (isTRUE(input$pca_varplot) && clustmethod() == "pca") {
        palette <- if (!is.null(input$palette)) input$palette else "muted_light"
        if (palette %in% c("original", "default", "custom", "")) palette <- "muted_light"
        if (gp$use_ggprism) {
          p <- apply_editor_theme(apply_ggprism_theme(create_gg_varplot(palette), gp), input)
          return(ggplot_as_plotly_image(p, width = 6, height = 6))
        }
        return(apply_plotly_editor_theme(create_varplot(palette), input))
      }

      if (gp$use_ggprism && !do3d) {
        plts <- create_gg_plotlist()
        shiny::req(length(plts) > 0)
        if (length(plts) == 1) {
          ggplot_as_plotly_image(plts[[1]], width = 6, height = 6)
        } else {
          nc <- ceiling(sqrt(length(plts)))
          combined <- patchwork::wrap_plots(plts, ncol = nc)
          nr <- ceiling(length(plts) / nc)
          ggplot_as_plotly_image(combined, width = nc * 4, height = nr * 4)
        }
      } else {
        fig <- if (length(create_plotlist()) == 1) {
          create_plotlist()[[1]]
        } else {
          plist <- create_plotlist()
          nc <- ceiling(sqrt(length(plist)))
          plotly::subplot(plist, nrows = nc, margin = 0.04)
        }
        apply_plotly_editor_theme(fig, input)
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
