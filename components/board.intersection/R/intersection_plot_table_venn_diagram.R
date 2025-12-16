##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

intersection_plot_venn_diagram_ui <- function(id,
                                              title,
                                              caption,
                                              info.text,
                                              label = "",
                                              height,
                                              width) {
  ns <- shiny::NS(id)

  FDR.VALUES2 <- c(1e-9, 1e-6, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.5, 1)

  venndiagram.opts <- shiny::tagList(
    shiny::fillRow(
      flex = c(1, 1),
      withTooltip(shiny::selectInput(ns("fdr"), "FDR", choices = FDR.VALUES2, selected = 0.20),
        "Threshold for false discovery rate",
        placement = "right", options = list(container = "body")
      ),
      withTooltip(
        shiny::selectInput(ns("lfc"), "logFC",
          choices = c(0, 0.1, 0.2, 0.5, 1, 2, 5),
          selected = 0.2
        ),
        "Threshold for fold-change (log2 scale)",
        placement = "right", options = list(container = "body")
      )
    ),
    shiny::br(), br(), br(), br(),
    shiny::radioButtons(ns("include"), "Counting:", choices = c("both", "up/down"), inline = TRUE)
  )

  PlotModuleUI(
    ns("vennplot"),
    title = title,
    label = "b",
    info.text = info.text,
    options = venndiagram.opts,
    caption = caption,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height,
    width = width
  )
}

intersection_table_venn_diagram_ui <- function(id,
                                               title,
                                               caption,
                                               info.text,
                                               label = "",
                                               height,
                                               width) {
  ns <- shiny::NS(id)

  venntable_opts <- shiny::tagList(
    shiny::selectInput(ns("venntable_intersection"), "Filter intersection:", choices = NULL)
  )

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    options = venntable_opts,
    caption = caption,
    height = height,
    width = width,
    title = title,
    label = "e"
  )
}


intersection_plot_venn_diagram_server <- function(id,
                                                  pgx,
                                                  level,
                                                  input_comparisons,
                                                  getFoldChangeMatrix,
                                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- shiny::reactive({
      dt <- getSignificanceCalls()
      shiny::req(dt)
      if (is.null(dt) || nrow(dt) == 0) {
        return(NULL)
      }

      dt1 <- dt[, 2:ncol(dt), drop = FALSE]
      label <- LETTERS[1:ncol(dt1)]
      colnames(dt1) <- label
      return(list(dt1, colnames(dt)[-1]))
    })

    venndiagram.RENDER <- function() {
      data <- plot_data()
      dt1 <- data[[1]]
      shiny::req(dt1)
      plot_names <- data[[2]]

      label <- colnames(dt1)
      include <- "both"
      if (input$include == "up/down") {
        include <- c("up/down")
      }

      par(mfrow = c(1, 1), mar = c(1, 1, 3, 1) * 0, bty = "n")
      par(oma = c(0.0, 0, 0, 0))
      if (ncol(dt1) == 1) {
        frame()
        text(0.5, 0.5, "Error: Venn diagram needs at least two groups")
        return(NULL)
      }
      if (ncol(dt1) > 7) {
        frame()
        text(0.5, 0.5, "Error: too many groups for Venn diagram")
        return(NULL)
      }

      p <- NULL
      if (0) {
        limma::vennDiagram(
          dt1,
          main = NULL, cex.main = 0.2, cex = 1.2, mar = c(0, 0, 2, 0),
          include = include, bty = "n", fg = grey(0.7),
          circle.col = c("turquoise", "salmon", "lightgreen", "orange")
        )
        tt <- paste(label, "=", plot_names)
        legend("topleft",
          legend = tt, bty = "n", cex = 0.9, y.intersp = 0.95,
          inset = c(0.04, -0.01), xpd = TRUE
        )
      } else {
        # convert matrix to list, simplify = FALSE will avoid apply function to simplify the result to a matrix
        x <- apply(dt1, 2, function(x) rownames(dt1)[which(x != 0)], simplify = FALSE)

        xlen <- sapply(x, length)

        if (length(x) == 0 || all(xlen == 0)) {
          frame()
          text(0.5, 0.5, tspan("Error: no valid genes. Please adjust thresholds.", js = FALSE),
            col = "red"
          )
          return(NULL)
        }

        if (include == "up/down") {
          x_up <- apply(dt1, 2, function(x) paste(rownames(dt1)[which(x > 0)], x[which(x > 0)]))
          x_down <- apply(dt1, 2, function(x) paste(rownames(dt1)[which(x < 0)], x[which(x < 0)]))

          if(length(x_up) != 0) {
            count_up <- ggVennDiagram::process_region_data(ggVennDiagram::Venn(x_up))$count
          } else {
            count_up <- 0
          }

          if(length(x_down) != 0) {
            count_down <- ggVennDiagram::process_region_data(ggVennDiagram::Venn(x_down))$count
          } else {
            count_down <- 0
          }

          count_both <- paste0(count_up, "\n", count_down)

          p <- ggVennDiagram::ggVennDiagram(
            x,
            label = "both",
            edge_size = 0.4,
            label_alpha = 0
          ) +
            ggplot2::scale_fill_gradient(low = "grey90", high = "red") +
            ggplot2::theme(
              legend.position = "none",
              plot.margin = ggplot2::unit(c(1, 1, 1, 1) * 0.3, "cm")
            )

          p$layers[[4]]$data$both <- count_both
        } else {
          # convert matrix x to a list of vectors
          p <- ggVennDiagram::ggVennDiagram(
            x,
            label = "count",
            edge_size = 0.4,
            label_alpha = 0
          ) +
            ggplot2::scale_fill_gradient(low = "grey90", high = "red") +
            ggplot2::theme(
              legend.position = "none",
              plot.margin = ggplot2::unit(c(1, 1, 1, 1) * 0.3, "cm")
            )
        }

        ## legend
        tt <- paste(label, "=", plot_names)
        n1 <- ceiling(length(tt) / 2)
        tt1 <- tt[1:n1]
        tt2 <- tt[(n1 + 1):length(tt)]
        if (length(tt2) < length(tt1)) tt2 <- c(tt2, "   ")
        tt1 <- paste(tt1, collapse = "\n")
        tt2 <- paste(tt2, collapse = "\n")

        xlim <- ggplot2::ggplot_build(p)$layout$panel_scales_x[[1]]$range$range
        ylim <- ggplot2::ggplot_build(p)$layout$panel_scales_y[[1]]$range$range
        x1 <- xlim[1] - 0.1 * diff(xlim)
        x2 <- xlim[1] + 0.6 * diff(xlim)
        y1 <- ylim[2] + 0.12 * diff(xlim)

        p <- p +
          ggplot2::annotate("text",
            x = x1, y = y1, hjust = "left",
            label = tt1, size = 4, lineheight = 0.83
          ) +
          ggplot2::annotate("text",
            x = x2, y = y1, hjust = "left",
            label = tt2, size = 4, lineheight = 0.83
          ) +
          ggplot2::coord_sf(clip = "off")
      }
      p
    }

    PlotModuleServer(
      "vennplot",
      func = venndiagram.RENDER,
      csvFunc = plot_data,
      res = c(72, 90), ## resolution of plots
      pdf.width = 8, pdf.height = 5,
      add.watermark = watermark
    )

    # Table

    getSignificanceCalls <- shiny::reactive({
      ## Gets the matrix of significance calls.

      sel <- head(names(pgx$gset.meta$meta), 7)
      sel <- input_comparisons()
      sel <- intersect(sel, names(pgx$gset.meta$meta))
      if (length(sel) == 0) {
        return(NULL)
      }

      res <- getFoldChangeMatrix()
      fc <- res$fc[, sel, drop = FALSE]
      qv <- res$qv[, sel, drop = FALSE]

      fdr <- 0.05
      lfc <- 0.2
      fdr <- as.numeric(input$fdr)
      lfc <- as.numeric(input$lfc)
      dt <- sign(fc) * (qv <= fdr & abs(fc) >= lfc)
      dt[is.na(dt)] <- 0
      ## add label of venn intersection region
      dt.labels <- LETTERS[1:ncol(dt)]

      venn.intersection <- apply(1 * (dt != 0), 1, function(x) {
        paste(dt.labels[which(x == 1)], collapse = "")
      })
      dt <- data.frame(intersection = venn.intersection, dt, check.names = FALSE)

      ## update filter
      choices <- c("<all>", sort(unique(venn.intersection)))
      selected <- venn.intersection[which.max(nchar(venn.intersection))]
      shiny::updateSelectInput(
        session, "venntable_intersection",
        choices = choices,
        selected = selected
      )

      return(dt)
    })

    getSignificantFoldChangeMatrix <- shiny::reactive({
      ## Filters FC matrix with significance and user-defined
      ## intersection region.
      dt <- getSignificanceCalls()
      shiny::req(dt)

      isect <- input$intersection
      fc0 <- getFoldChangeMatrix()$fc
      if (length(isect) == 0) {
        fc1 <- fc0
      } else {
        ## only genes at least significant in one group
        jj <- which(rowSums(dt[, 2:ncol(dt), drop = FALSE] != 0, na.rm = TRUE) > 0)
        if (length(jj) == 0) {
          return(NULL)
        }
        dt <- dt[jj, , drop = FALSE]

        ## check same sign
        if (input$include == "up/down") {
          kk <- 1 + match(c("B", "C"), LETTERS[1:10])
          kk <- 1 + match(isect, LETTERS[1:10])
          kk <- intersect(kk, 1:ncol(dt))

          dt1 <- dt[, kk, drop = FALSE]
          jj <- which(rowMeans(sign(dt1) == +1, na.rm = TRUE) == 1 |
            (rowMeans(sign(dt1) == -1, na.rm = TRUE) == 1))
          dt <- dt[jj, , drop = FALSE]
          remove(dt1)
        }

        ## only genes in the selected intersection
        intersection <- "ABC"
        intersection <- paste0(input$intersection, collapse = "")
        dt <- dt[which(dt$intersection == intersection), , drop = FALSE]
      }

      ## filtered by family/collection
      fc1 <- fc0[intersect(rownames(dt), rownames(fc0)), , drop = FALSE]
      if (nrow(dt) == 1) {
        fc1 <- matrix(fc1, nrow = 1)
        rownames(fc1) <- rownames(dt)
        colnames(fc1) <- colnames(fc0)
      }

      ## filtered by SPLOM selection
      splom.sel <- plotly::event_data("plotly_selected", source = "splom")
      sel.keys <- as.character(splom.sel$key)
      if (1 && length(sel.keys) > 0) {
        sel <- intersect(sel.keys, rownames(fc1))
        fc1 <- fc1[sel, , drop = FALSE]
      }

      ## only active/selected comparisons
      sel <- colnames(dt)[-1]
      kk <- match(sel, gsub(" \\(-\\)", "", colnames(fc1)))
      fc1 <- fc1[, kk, drop = FALSE]

      ## order
      fc1 <- fc1[order(-rowMeans(fc1, na.rm = TRUE)), , drop = FALSE]
      fc1 <- round(fc1, digits = 3)
      colnames(fc1) <- LETTERS[1:ncol(fc1)]

      ## add intersection code
      sel <- match(rownames(fc1), rownames(dt))
      fc1 <- data.frame(intersection = dt$intersection[sel], fc = fc1)

      ## filter on user selection
      vv <- input$venntable_intersection
      if (vv != "<all>") {
        sel <- which(fc1$intersection == vv)
        fc1 <- fc1[sel, , drop = FALSE]
      }
      return(fc1)
    })

    venntable.RENDER <- shiny::reactive({
      shiny::req(pgx$X)

      ## get foldchanges
      fc0 <- getSignificantFoldChangeMatrix() ## isolate??
      if (is.null(fc0) || nrow(fc0) == 0) {
        return(NULL)
      }

      ## add gene name/title
      if (level() == "gene") {
        gene <- as.character(pgx$genes[rownames(fc0), "gene_name"])
        gene.tt <- substring(playdata::GENE_TITLE[gene], 1, 50)
        gene.tt <- as.character(gene.tt)
        fc0 <- data.frame(name = gene, fc0, check.names = FALSE)
      } else {
        name <- substring(rownames(fc0), 1, 50)
        name[is.na(name)] <- "NA"
        fc0 <- data.frame(name = name, fc0, check.names = FALSE)
      }

      df <- data.frame(fc0, check.names = FALSE)
      nsc <- setdiff(1:ncol(df), 2)
      DT::datatable(df,
        class = "compact cell-border stripe",
        rownames = FALSE,
        extensions = c("Scroller"),
        plugins = "scrollResize",
        selection = "none",
        options = list(
          dom = "tip",
          scrollX = TRUE,
          scrollY = 215,
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    })

    venntable.RENDER2 <- shiny::reactive({
      dt <- venntable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    TableModuleServer(
      "datasets",
      func = venntable.RENDER,
      func2 = venntable.RENDER2,
      selector = "none"
    )
  })
}
