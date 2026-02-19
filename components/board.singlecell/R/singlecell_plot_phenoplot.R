##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Single cell plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
singlecell_plot_phenoplot_ui <- function(
  id,
  label = "",
  title,
  info.text,
  caption,
  height,
  width
) {
  ns <- shiny::NS(id)

  phenoplot.opts <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(ns("labelmode"), "Label:", c("groups", "legend"), inline = TRUE),
      "Select whether you want the group labels to be plotted inside the plots or in a seperate legend."
    )
  )

  PlotModuleUI(
    id = ns("plotmodule"),
    plotlib = "plotly",
    label = label,
    info.text = info.text,
    title = title,
    caption = caption,
    options = phenoplot.opts,
    download.fmt = c("png", "pdf", "svg"),
    height = height,
    width = width,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "clustering"
  )
}

#' Single cell plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @export
singlecell_plot_phenoplot_server <- function(id,
                                             pgx,
                                             pfGetClusterPositions,
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Editor: dynamic color pickers for custom palette
    output$custom_palette_ui <- shiny::renderUI({
      shiny::req(input$palette == "custom")
      pd <- plot_data()
      shiny::req(pd)
      Y <- pd$pgx$samples
      sel <- pd$sel
      pheno <- pd$pheno

      is.num.check <- function(y) {
        suppressWarnings(numy <- as.numeric(as.character(y)))
        !all(is.na(numy)) && (length(unique(y)) / length(y)) > 0.1
      }

      max_levels <- 0
      level_names <- vector("list", 8)
      for (ph in pheno) {
        y <- Y[sel, ph]
        y <- y[!y %in% c(NA, "", " ", "NA", "na")]
        if (length(y) == 0 || is.num.check(y)) next
        lvls <- sort(unique(as.character(y)))
        max_levels <- max(max_levels, length(lvls))
        for (j in seq_along(lvls)) {
          if (j <= 8) level_names[[j]] <- c(level_names[[j]], lvls[j])
        }
      }
      max_levels <- min(max(max_levels, 1), 8)

      default_clrs <- omics_pal_d(palette = "muted_light")(8)
      pickers <- lapply(seq_len(max_levels), function(i) {
        nms <- level_names[[i]]
        group_label <- if (length(nms) > 0) paste(unique(nms), collapse = ", ") else paste("Color", i)
        colourpicker::colourInput(
          ns(paste0("custom_color_", i)),
          label = group_label,
          value = default_clrs[i]
        )
      })
      shiny::tagList(pickers)
    })

    plot_data <- shiny::reactive({
      #
      shiny::req(pgx$X)
      clust.pos <- pfGetClusterPositions()
      if (is.null(clust.pos)) {
        return(NULL)
      }

      pos <- pgx$tsne2d
      pos <- clust.pos
      sel <- rownames(pos)
      pheno <- colnames(pgx$Y)
      return(list(
        pos = pos,
        pgx = pgx,
        sel = sel,
        pheno = pheno
      ))
    })

    get_plots <- function(cex = 1) {
      #
      pd <- plot_data()
      shiny::req(pd)

      sel <- pd[["sel"]]
      pheno <- pd[["pheno"]]
      Y <- pd[["pgx"]]$samples
      pos <- pd[["pos"]]

      cex1 <- 1.2 * c(1.8, 1.3, 0.8, 0.5)[cut(nrow(pos), breaks = c(-1, 40, 200, 1000, 1e10))]
      cex1 <- cex1 * ifelse(length(pheno) > 6, 0.8, 1)
      cex1 <- cex1 * ifelse(length(pheno) > 12, 0.8, 1)

      ## Editor: palette
      palette <- if (!is.null(input$palette)) input$palette else "muted_light"
      if (palette %in% c("original", "default")) palette <- "muted_light"
      if (palette == "custom") {
        base_clrs <- sapply(1:8, function(j) {
          val <- input[[paste0("custom_color_", j)]]
          if (is.null(val)) omics_pal_d(palette = "muted_light")(8)[j] else val
        })
      } else {
        base_clrs <- omics_pal_d(palette = palette)(8)
      }

      ## is it a float/number???
      is.num <- function(y, fmin = 0.1) {
        suppressWarnings(numy <- as.numeric(as.character(y)))
        t1 <- !all(is.na(numy)) && is.numeric(numy)
        t2 <- (length(unique(y)) / length(y)) > fmin
        (t1 && t2)
      }

      plt <- list()
      i <- 6
      for (i in 1:min(20, length(pheno))) {
        px <- 4
        px <- pheno[i]
        y <- Y[sel, px]
        y[which(y %in% c(NA, "", " ", "NA", "na"))] <- NA
        if (sum(!is.na(y)) == 0) next

        if (playbase::is.num(y)) {
          klrpal <- colorRampPalette(c("grey90", "grey50", "red3"))(16)
          y <- rank(as.numeric(y))
          ny <- round(1 + 15 * (y - min(y)) / (max(y) - min(y)))
          klr0 <- klrpal[ny]
        } else {
          y <- factor(as.character(y))
          klrpal <- rep(paste0(gplots::col2hex(base_clrs), "99"), length.out = 100)
          klr0 <- klrpal[y]
        }

        p <- playbase::pgx.scatterPlotXY.PLOTLY(
          pos,
          var = y,
          type = "factor", ## always factor?
          col = klrpal,
          cex = 0.5 * cex1 * cex,
          xlab = NA,
          ylab = NA,
          xlim = 1.2 * range(pos[, 1]),
          ylim = 1.2 * range(pos[, 2]),
          axis = FALSE,
          title = tolower(pheno[i]),
          cex.title = 1.2 * cex,
          title.y = 1,
          cex.clust = 1.2 * cex,
          label.clusters = TRUE,
          legend = FALSE,
          gridcolor = "fff"
        ) %>%
          plotly::layout(
            #
            plot_bgcolor = "#f8f8f8",
            margin = list(0, 0, 0, 0)
          )

        plt[[i]] <- p
      }
      return(plt)
    }

    plotly.RENDER <- function() {
      pd <- plot_data()
      plt <- get_plots(cex = 0.9)
      shiny::req(plt)
      ## layout
      nr <- 2
      if (length(plt) > 4) nr <- 3
      if (length(plt) > 6) nr <- 4
      if (length(plt) > 12) nr <- 5
      if (length(plt) > 1) {
        fig <- plotly::subplot(
          plt,
          nrows = nr,
          #
          margin = c(0.01, 0.01, 0.01, 0.045)
        ) %>%
          plotly_default() %>%
          plotly::layout(
            margin = list(l = 10, r = 10, b = 10, t = 20) # lrbt
          )
      } else {
        fig <- plt[[1]] %>%
          plotly_default()
      }
      return(fig)
    }

    plotly_modal.RENDER <- function() {
      pd <- plot_data()
      plt <- get_plots(cex = 1.3)
      ## layout
      nc <- 2
      if (length(plt) > 4) nc <- 3
      if (length(plt) > 6) nc <- 4
      if (length(plt) > 12) nc <- 5
      nr <- ceiling(length(plt) / nc)
      fig <- plotly::subplot(
        plt,
        nrows = nr,
        margin = c(0.01, 0.01, 0.01, 0.045)
      ) %>%
        plotly_modal_default() %>%
        plotly::layout(
          margin = list(l = 20, r = 20, b = 20, t = 40) # lrbt
        )

      return(fig)
    }

    PlotModuleServer(
      "plotmodule",
      func = plotly.RENDER,
      func2 = plotly_modal.RENDER,
      plotlib = "plotly",
      res = c(85, 95),
      pdf.width = 6,
      pdf.height = 10,
      add.watermark = watermark,
      parent_session = session
    )
  }) ## end of moduleServer
}
