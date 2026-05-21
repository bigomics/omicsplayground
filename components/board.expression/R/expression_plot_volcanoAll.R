##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Expression plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
expression_plot_volcanoAll_ui <- function(id,
                                          title,
                                          caption,
                                          info.text,
                                          info.methods,
                                          info.references,
                                          info.extra_link,
                                          ## labeltype,
                                          label = "",
                                          height,
                                          width) {
  ns <- shiny::NS(id)

  plot_options <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("scale_per_plot"), "scale per plot", FALSE),
      "Scale each volcano plots individually.",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(
    id = ns("pltmod"),
    title = title,
    label = label,
    plotlib = c("plotly", "ggplot"),
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
    options = plot_options,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height,
    width = width,
    cards = TRUE,
    card_names = c("dynamic", "static"),
    editor = TRUE,
    ns_parent = ns,
    plot_type = "volcano"
  )
}

#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
expression_plot_volcanoAll_server <- function(id,
                                              pgx,
                                              getAllContrasts,
                                              fdr,
                                              lfc,
                                              show_pv,
                                              genes_selected,
                                              labeltype = reactive("symbol"),
                                              pval_cap,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ## reactive function listening for changes in input
    plot_data <- shiny::reactive({
      shiny::req(pgx$X)

      # Input variables
      ct <- getAllContrasts()
      F <- ct$F
      Q <- ct$Q
      P <- ct$P

      fdr <- as.numeric(fdr())
      lfc <- as.numeric(lfc())
      comp <- colnames(F)
      shiny::req(length(comp) > 0)

      ## combined matrix for output
      FQ <- data.frame(fc = F, q = Q, p = P)
      features <- rownames(FQ)
      symbols <- pgx$genes[rownames(FQ), "symbol"]
      names <- pgx$genes[rownames(FQ), "gene_title"]

      label.names <- playbase::probe2symbol(rownames(FQ), pgx$genes, labeltype(), fill_na = TRUE)

      # Input vars
      if (show_pv()) {
        ## P <- P
        title_y <- "Significance (-log10p)"
      } else {
        P <- Q
        title_y <- "Significance (-log10q)"
      }

      ## ps: FQ contains log2FC+q-value or log2FC+p-value. Depends on show_pv option.
      ## Position data for click-to-label (all facets, so nearest-neighbor
      ## matches the actual plotted point regardless of which facet is clicked)
      click_df <- data.frame(
        x = as.vector(F),
        y = as.vector(-log10(P)),
        feature_name = rep(features, ncol(F))
      )

      pd <- list(
        FQ = FQ, ## Remember: the first element is returned as downloadable CSV
        comp = comp,
        fdr = fdr,
        lfc = lfc,
        F = F,
        P = P,
        title_y = title_y,
        symbols = symbols,
        features = features,
        names = names,
        label.names = label.names,
        sel.genes = genes_selected()$sel.genes,
        lab.genes = genes_selected()$lab.genes,
        df = click_df
      )

      return(pd)
    })

    plotly_plots <- function(cex = 2, yrange = 0.5, n_rows = 2,
                             margin_l = 50, margin_b = 50) {
      pd <- plot_data()
      shiny::req(pd)

      pval_cap <- pval_cap()

      # Call volcano plots
      all_plts <- playbase::plotlyVolcano_multi(
        FC = pd$F,
        Q = pd$P,
        fdr = pd$fdr,
        lfc = pd$lfc,
        cex = cex,
        names = pd$features,
        label.names = pd$label.names,
        share_axis = !input$scale_per_plot,
        title_y = pd$title_y,
        title_x = "Effect size (log2FC)",
        yrange = yrange,
        n_rows = n_rows,
        margin_l = margin_l,
        margin_b = margin_b,
        color_up_down = TRUE,
        colors = extract_volcano_colors(input),
        highlight = pd$sel.genes,
        label = pd$lab.genes,
        by_sig = FALSE,
        pval_cap = pval_cap
      )

      return(all_plts)
    }

    plotly.RENDER <- function() {
      fig <- plotly_plots(
        cex = 3, yrange = 0.05, n_rows = 2, margin_b = 40, margin_l = 50
      ) %>%
        playbase::plotly_build_light(.)
      return(fig)
    }

    big_plotly.RENDER <- function() {
      fig <- plotly_plots(
        yrange = 0.02, n_rows = 3, margin_b = 70, margin_l = 70
      ) %>%
        plotly::style(
          marker.size = 6
        ) %>%
        playbase::plotly_build_light(.)

      return(fig)
    }

    base.plots <- function(label.cex = 4) {
      pd <- plot_data()
      shiny::req(pd)

      fc <- pd$F
      qv <- pd$P
      feature_names <- rep(rownames(fc), each = ncol(fc))
      label.names <- rep(pd$label.names, each = ncol(fc))
      pivot.fc <- data.frame(fc) %>%
        tidyr::pivot_longer(
          cols = everything(), # Select all columns to pivot
          names_to = "facet", # Name of the new column for timepoints
          values_to = "fc"
        )
      pivot.qv <- data.frame(qv) %>%
        tidyr::pivot_longer(
          cols = everything(), # Select all columns to pivot
          names_to = "facet", # Name of the new column for timepoints
          values_to = "qv"
        )
      facet <- pivot.fc$facet
      x <- pivot.fc$fc

      pval_cap <- pval_cap()
      y <- -log10(pivot.qv$qv + pval_cap)

      ## Editor: custom labels
      label_features <- get_custom_labels(input, pd[["features"]], defaults = pd[["lab.genes"]])

      highlight <- if (isTRUE(input$color_selection)) {
        label_features
      } else {
        pd[["sel.genes"]]
      }
      if (!is.null(input$cutoff_type) && input$cutoff_type == "hyperbolic") {
        highlight <- label_features
      }

      ## Editor: extract settings via helpers
      plot_colors <- extract_volcano_colors(input)
      ls <- extract_label_settings(input, defaults = list(label_size = label.cex))
      gp <- extract_ggprism_params(input)

      ## Editor: hyperbolic cutoff settings
      use_hyperbola <- !is.null(input$cutoff_type) && input$cutoff_type == "hyperbolic"

      p <- playbase::ggVolcano(
        x,
        y,
        names = feature_names,
        facet = facet,
        label = label_features,
        highlight = highlight,
        label.names = label.names,
        label.cex = ls$label_size,
        xlab = "Effect size (log2FC)",
        ylab = pd$title_y,
        psig = pd[["fdr"]],
        lfc = pd[["lfc"]],
        marker.size = ls$marker_size,
        showlegend = FALSE,
        title = NULL,
        axis.text.size = ls$axis_text_size,
        colors = plot_colors,
        box.padding = ls$box_padding,
        min.segment.length = ls$min_segment_length,
        label.box = ls$label_box,
        segment.linetype = ls$segment_linetype,
        use_hyperbola = use_hyperbola,
        hyperbola_k = ls$hyperbola_k,
        use_ggprism = gp$use_ggprism,
        ggprism_border = gp$ggprism_border,
        ggprism_axis_guide = gp$ggprism_axis_guide,
        ggprism_show_legend = gp$ggprism_show_legend,
        ggprism_legend_x = gp$ggprism_legend_x,
        ggprism_legend_y = gp$ggprism_legend_y,
        ggprism_legend_border = gp$ggprism_legend_border
      )

      ## Editor: margins & aspect ratio
      p <- apply_editor_theme(p, input)

      p
    }

    big_base.plots <- function() {
      base.plots(label.cex = 5)
    }

    plot_grid <- list(
      list(plotlib = "plotly", func = plotly.RENDER, func2 = big_plotly.RENDER, card = 1),
      list(plotlib = "ggplot", func = base.plots, func2 = big_base.plots, card = 2)
    )

    lapply(plot_grid, function(x) {
      PlotModuleServer(
        "pltmod",
        plotlib = x$plotlib,
        func = x$func,
        func2 = x$func2,
        csvFunc = plot_data,
        res = c(70, 90), # resolution of plots
        pdf.width = 12,
        pdf.height = 5,
        add.watermark = watermark,
        card = x$card,
        parent_session = session
      )
    })
  }) ## end of moduleServer
}
