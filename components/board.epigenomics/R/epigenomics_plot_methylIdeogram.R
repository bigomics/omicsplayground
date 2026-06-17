## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.

epigenomics_plot_methylIdeogram_ui <- function(id,
                                               label = "",
                                               title,
                                               height,
                                               width,
                                               caption,
                                               info.text,
                                               info.methods,
                                               info.references,
                                               info.extra_link) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::div(
      id = ns("show_pheno_lines_wrap"),
      withTooltip(
        shiny::checkboxInput(
          ns("show_pheno_lines"),
          "Show phenotypes' average (if phenotype is selected)",
          FALSE
        ),
        "Compute and display average beta values per phenotype",
        placement = "top"
      )
    ),
    shiny::div(
      id = ns("remove_probecounts_bar_wrap"),
      withTooltip(
        shiny::checkboxInput(
          ns("remove_probecounts_bar"),
          "Remove barplots of feature counts",
          FALSE
        ),
        "Remove barplots of feature counts in the lower panel of each ideogram.",
        placement = "top"
      )
    )
  )

  PlotModuleUI(
    ns("pltmod"),
    plotlib = "base",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height,
    options = options,
    label = label,
    caption = caption,
    title = title,
    ns_parent = ns,
    editor = TRUE,
    plot_type = "scatterplot"
  )
}

epigenomics_plot_methylIdeogram_server <- function(id,
                                                   pgx,
                                                   r.chromosome = reactive(""),
                                                   r.samples = reactive(""),
                                                   r.pheno = reactive(""),
                                                   watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    shiny::observe({
      pheno_val <- r.pheno()
      if (is.null(pheno_val) || pheno_val == "<ungrouped>") {
        shinyjs::hide("show_pheno_lines_wrap")
        shinyjs::show("remove_probecounts_bar_wrap")
      } else {
        shinyjs::show("show_pheno_lines_wrap")
        shinyjs::hide("remove_probecounts_bar_wrap")
      }
    })

    plot_data <- shiny::reactive({
      shiny::req(pgx$X, pgx$genes, pgx$samples)
      X <- playbase::mToBeta(pgx$X)
      genes <- pgx$genes
      rownames(X) <- sub("_.*", "", rownames(X))
      rownames(genes) <- sub("_.*", "", rownames(genes))

      kk <- intersect(rownames(X), rownames(genes))
      if (length(kk) == 0) {
        return(NULL)
      }
      samples <- r.samples()
      if (!all(samples %in% colnames(X))) {
        return(NULL)
      }

      annot <- genes[kk, , drop = FALSE]
      pos_col <- colnames(annot)[grep("^pos$|position|location", tolower(colnames(annot)))[1]]
      if (!is.na(pos_col)) {
        annot[[pos_col]] <- as.numeric(sub(";.*", "", as.character(annot[[pos_col]])))
      }

      return(list(
        X = X[kk, samples, drop = FALSE],
        samples = pgx$samples[samples, , drop = FALSE],
        annot = annot,
        chromosomes = r.chromosome(),
        pheno = r.pheno(),
        dma = pgx$dma
      ))
    })

    plot.RENDER <- function() {
      res <- plot_data()
      shiny::req(res)

      X <- res[["X"]]
      samples <- res[["samples"]]
      pheno <- res[["pheno"]]
      annot <- res[["annot"]]
      chromosomes <- res[["chromosomes"]]
      dma <- res[["dma"]]

      if (!is.null(pheno) && pheno != "<ungrouped>") {
        kk <- which(colnames(samples) == pheno)
        if (length(kk) > 0) {
          pheno <- as.character(samples[, kk[1]])
          nas <- which(is.na(pheno) | pheno == "")
          if (length(nas) > 0) {
            pheno <- pheno[-nas]
            X <- X[, -nas, drop = FALSE]
          }
          if (length(unique(pheno)) != 2) pheno <- NULL
        } else {
          pheno <- NULL
        }
      } else {
        pheno <- NULL
      }

      if (!is.null(dma)) {
        loess_bins <- 5L
        bin_size <- 1e6
        if (dma == "Differentially methylated regions") {
          loess_bins <- 20L
          bin_size <- 2e6
        }
      }

      playbase::plotMethylIdeogram(X, annot, pheno, chromosomes,
        bin_size = bin_size,
        probe_count_bars = !input$remove_probecounts_bar,
        pheno_lines = input$show_pheno_lines,
        loess_bins = loess_bins
      )
    }

    modal_plot.RENDER <- function() {
      plot.RENDER()
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "base",
      plotlib2 = "base",
      func = plot.RENDER,
      func2 = modal_plot.RENDER,
      csvFunc = plot_data,
      renderFunc = shiny::renderPlot,
      renderFunc2 = shiny::renderPlot,
      res = c(100, 170) * 0.85,
      pdf.width = 6,
      pdf.height = 6,
      add.watermark = watermark,
      parent_session = session
    )
  })
}
