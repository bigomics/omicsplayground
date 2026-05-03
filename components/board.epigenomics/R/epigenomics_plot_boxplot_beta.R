## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.

epigenomics_plot_boxplot_beta_ui <- function(id,
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
    withTooltip(
      shiny::checkboxInput(
        ns("remove_sex_chr"),
        "Remove sex chromosomes",
        FALSE
      ),
      "Remove sex chromosomes from boxplot",
      placement = "top"
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

epigenomics_plot_boxplot_beta_server <- function(id,
                                                 pgx,
                                                 r.chromosome = reactive(""),
                                                 r.samples = reactive(""),
                                                 r.pheno = reactive(""),
                                                 watermark = FALSE) {

  moduleServer(id, function(input, output, session) {

    plot_data <- shiny::reactive({

      shiny::req(pgx$X, pgx$genes, pgx$samples)
      X <- playbase::mToBeta(pgx$X)
      Y <- pgx$samples
      annot <- pgx$genes
      rownames(X) <- sub("_.*", "", rownames(X))
      rownames(annot) <- sub("_.*", "", rownames(annot))
      kk <- intersect(rownames(X), rownames(annot))
      if (length(kk) == 0) return(NULL)
      samples <- r.samples()
      if (!all(samples %in% colnames(X))) return(NULL)
      X <- X[kk, samples, drop = FALSE]
      annot <- annot[kk, , drop = FALSE]
      Y <- Y[samples, , drop = FALSE]
      
      if (input$remove_sex_chr) {
        kk <- grep("chr|chromosome|chrom", tolower(colnames(annot)))[1]
        chroms <- paste0("chr", sub("^chr", "", sub("[pq].*", "", annot[, kk])))
        kk <- which(chroms %in% c("chrX","chrY"))
        if (length(kk) > 0) {
          X <- X[-kk, , drop = FALSE]
          annot <- annot[-kk, , drop = FALSE]
        }
      }

      return(list(X = X, samples = Y, annot = annot, pheno = r.pheno()))

    })

    plot.RENDER <- function() {

      res <- plot_data()
      shiny::req(res)

      X <- res[["X"]]
      samples <- res[["samples"]]
      pheno <- res[["pheno"]]
      annot <- res[["annot"]]
      
      if (!is.null(pheno) & pheno != "<ungrouped>") {
        kk <- which(colnames(samples) == pheno)
        if (length(kk) > 0) {
          pheno <- as.character(samples[, kk[1]])
          kk <- which(!is.na(pheno) & pheno != "")
          if (length(unique(pheno[kk])) == 2) {
            pheno <- pheno[kk]            
            samples <- samples[kk, , drop = FALSE]
            X <- X[, kk, drop = FALSE]
          } else {
            pheno <- "<ungrouped>"
          }
        }
      }
      playbase::plotMethylOverview(X, annot, pheno,
        plot.beta.dist = FALSE, plot.beta.boxplots = TRUE)

    }

    modal_plot.RENDER <- function() { plot.RENDER() }

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
