## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.

epigenomics_plot_beta_dist_ui <- function(id,
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
    label = label,
    caption = caption,
    title = title,
    ns_parent = ns,
    editor = TRUE,
    plot_type = "scatterplot"
  )

}

epigenomics_plot_beta_dist_server <- function(id,
                                              pgx,
                                              r.samples = reactive(""),
                                              r.pheno = reactive(""),
                                              watermark = FALSE) {

  moduleServer(id, function(input, output, session) {
    
    plot_data <- shiny::reactive({

      shiny::req(pgx$X, pgx$genes, pgx$samples)
      X <- playbase::mToBeta(pgx$X)
      genes <- pgx$genes
      rownames(X) <- sub("_.*", "", rownames(X))
      rownames(genes) <- sub("_.*", "", rownames(genes))
      kk <- intersect(rownames(X), rownames(genes))
      if (length(kk) == 0) return(NULL)
      samples <- r.samples()
      if (!all(samples %in% colnames(X))) return(NULL)

      return(list(
        X = X[kk, samples, drop = FALSE],
        samples = pgx$samples[samples, , drop = FALSE],
        annot = genes[kk, , drop = FALSE],
        pheno = r.pheno()
      ))

    })

    plot.RENDER <- function() {

      res <- plot_data()
      shiny::req(res)

      X <- res[["X"]]
      samples <- res[["samples"]]
      pheno <- res[["pheno"]]
      annot <- res[["annot"]]

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

      playbase::plotMethylOverview(X, annot, pheno, NULL, plot.beta.dist = TRUE)
        
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
