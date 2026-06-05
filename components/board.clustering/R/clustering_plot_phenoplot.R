##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


clustering_plot_phenoplot_ui <- function(
  id,
  title,
  info.text,
  info.methods,
  info.extra_link,
  caption,
  label = "",
  height,
  width
) {
  ns <- shiny::NS(id)

  phenoplot.opts <- shiny::tagList(
    shiny::checkboxInput(ns("showlabels"), "Show group labels", TRUE)
  )

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    caption = caption,
    options = phenoplot.opts,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "clustering_prism"
  )
}

clustering_plot_phenoplot_server <- function(id,
                                             pgx,
                                             selected_phenotypes,
                                             clustmethod,
                                             selected_samples,
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Editor: dynamic color pickers for custom palette
    output$custom_palette_ui <- shiny::renderUI({
      shiny::req(input$palette == "custom")
      jj <- selected_samples()
      kk <- selected_phenotypes()
      shiny::req(length(jj) > 0, length(kk) > 0)
      kk <- kk[which(kk %in% colnames(pgx$Y))]
      shiny::req(length(kk) > 0)
      Y <- pgx$Y[jj, kk, drop = FALSE]

      ## collect unique group names per color position across all phenotypes
      max_levels <- 0
      level_names <- vector("list", 8)
      for (ph in kk) {
        vals <- as.character(Y[, ph])
        vals[vals %in% c("", " ", "NA", "na")] <- NA
        lvls <- sort(unique(vals[!is.na(vals)]))
        if (length(lvls) == 0) next
        max_levels <- max(max_levels, length(lvls))
        for (j in seq_along(lvls)) {
          if (j <= 8) level_names[[j]] <- c(level_names[[j]], lvls[j])
        }
      }
      max_levels <- min(max(max_levels, 1), 8)

      picker_labels <- sapply(seq_len(max_levels), function(i) {
        nms <- level_names[[i]]
        if (length(nms) > 0) paste(unique(nms), collapse = ", ") else paste("Color", i)
      })
      custom_palette_pickers(picker_labels, ns)
    })

    plot_data <- reactive({
      pgx <- pgx
      shiny::req(pgx$Y)
      Y <- pgx$Y

      ## get t-SNE positions
      clustmethod1 <- paste0(clustmethod(), "2d")
      pos <- pgx$cluster$pos[[clustmethod1]]
      colnames(pos) <- c("x", "y")
      jj <- selected_samples()
      kk <- selected_phenotypes()
      kk <- kk[which(kk %in% colnames(Y))]
      pos <- pos[jj, , drop = FALSE]
      shiny::validate(shiny::need(
        nrow(pos) > 1,
        "Filtering too restrictive. Please change 'Filter samples' settings."
      ))
      Y <- Y[jj, kk, drop = FALSE]
      ## complete dataframe for downloading
      df <- data.frame(pos, Y, check.names = FALSE)
      return(df)
    })

    create_plots <- function(cex = 1) {
      pd <- plot_data()
      showlabels <- input$showlabels
      pheno <- selected_phenotypes()
      shiny::validate(shiny::need(
        length(pheno) > 0,
        "Please select at least one phenotype."
      ))
      pheno.ex <- setdiff(pheno, colnames(pd))
      shiny::validate(shiny::need(
        length(pheno.ex) == 0,
        paste0(pheno.ex, " appear(s) not valid as phenotype(s). Please remove from selection on the menu on the right.")
      ))
      pheno <- pheno[which(pheno %in% colnames(pd))]
      Y <- pd[, pheno, drop = FALSE]
      pos <- pd[, c("x", "y")]

      ## points size depending on how many points we have
      ncex <- cut(nrow(pos), breaks = c(-1, 40, 200, 1000, 1e10))
      cex1 <- 0.8 * cex * c(1.8, 1.3, 0.8, 0.5)[ncex]
      cex1 <- cex1 * ifelse(length(pheno) > 6, 0.8, 1)
      cex1 <- cex1 * ifelse(length(pheno) > 12, 0.8, 1)

      ## Editor: palette
      base_clrs <- resolve_palette_colors(input, 8, fallback_colors = omics_pal_d("muted_light")(8))

      plt <- list()
      for (i in 1:min(20, length(pheno))) {
        ## ------- set colors
        colvar <- factor(Y[, pheno[i]])
        colvar[which(colvar %in% c(NA, "", " ", "NA", "na"))] <- NA
        colvar <- factor(as.character(colvar))
        klrpal <- rep(base_clrs, 10)
        klr1 <- klrpal[colvar]
        jj <- which(is.na(klr1))
        if (length(jj)) klr1[jj] <- "#AAAAAA22"
        tt <- tolower(pheno[i])

        ## ------- start plot
        p <- playbase::pgx.scatterPlotXY.PLOTLY(
          pos,
          var = colvar,
          col = klrpal,
          cex = cex1,
          xlab = "",
          ylab = "",
          title = tt,
          cex.title = cex * 1.2,
          cex.clust = cex * 1.1,
          label.clusters = showlabels
        )
        plt[[i]] <- p
      }
      return(plt)
    }

    create_gg_plots <- function(cex = 1) {
      pd <- plot_data()
      showlabels <- input$showlabels
      pheno <- selected_phenotypes()
      shiny::validate(shiny::need(length(pheno) > 0, "Please select at least one phenotype."))
      pheno <- pheno[which(pheno %in% colnames(pd))]
      shiny::req(length(pheno) > 0)
      Y <- pd[, pheno, drop = FALSE]
      pos <- pd[, c("x", "y")]

      cex1 <- 0.8 * cex
      cex1 <- cex1 * ifelse(length(pheno) > 6, 0.8, 1)
      cex1 <- cex1 * ifelse(length(pheno) > 12, 0.8, 1)

      base_clrs <- resolve_palette_colors(input, 8, fallback_colors = omics_pal_d("muted_light")(8))
      gp <- extract_ggprism_params(input)

      plts <- list()
      for (i in 1:min(20, length(pheno))) {
        colvar <- factor(Y[, pheno[i]])
        colvar[which(colvar %in% c(NA, "", " ", "NA", "na"))] <- NA
        colvar <- factor(as.character(colvar))
        klrpal <- rep(base_clrs, 10)
        tt <- tolower(pheno[i])

        p <- playbase::pgx.scatterPlotXY.GGPLOT(
          pos,
          var = colvar,
          col = klrpal,
          cex = cex1,
          xlab = "",
          ylab = "",
          title = tt,
          cex.title = cex * 1.2,
          cex.clust = cex * 1.1,
          label.clusters = showlabels
        )
        p <- apply_ggprism_theme(p, gp)
        p <- apply_editor_theme(p, input)
        plts[[length(plts) + 1]] <- p
      }
      plts
    }

    plotly.RENDER <- function() {
      gp <- extract_ggprism_params(input)

      if (gp$use_ggprism) {
        plts <- create_gg_plots(cex = 0.85)
        shiny::req(length(plts) > 0)
        nc <- floor(sqrt(length(plts)))
        if (nc < 1) nc <- 1
        combined <- patchwork::wrap_plots(plts, ncol = nc)
        nr <- ceiling(length(plts) / nc)
        ggplot_as_plotly_image(combined, width = nc * 3, height = nr * 3)
      } else {
        plt <- create_plots(cex = 0.85)
        nc <- floor(sqrt(length(plt)))
        combined_plots <- plotly::subplot(
          plt,
          nrows = nc,
          margin = 0.04
        )
        combined_plots <- plotly::layout(combined_plots,
          margin = list(t = 40)
        )
        return(apply_plotly_editor_theme(combined_plots, input))
      }
    }

    plotly_modal.RENDER <- function() {
      gp <- extract_ggprism_params(input)

      if (gp$use_ggprism) {
        plts <- create_gg_plots(cex = 1.3)
        shiny::req(length(plts) > 0)
        nc <- ceiling(sqrt(length(plts)))
        combined <- patchwork::wrap_plots(plts, ncol = nc)
        nr <- ceiling(length(plts) / nc)
        ggplot_as_plotly_image(combined, width = nc * 4, height = nr * 4)
      } else {
        plt <- create_plots(cex = 1.3)
        nc <- ceiling(sqrt(length(plt)))
        combined_plots <- plotly::subplot(
          plt,
          nrows = nc,
          margin = 0.04
        )
        combined_plots <- plotly::layout(combined_plots,
          margin = list(t = 40)
        )
        return(apply_plotly_editor_theme(combined_plots, input))
      }
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plotly.RENDER,
      func2 = plotly_modal.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(85), ## resolution of plots
      pdf.width = 6,
      pdf.height = 9,
      add.watermark = watermark,
      parent_session = session
    )
  })
}
