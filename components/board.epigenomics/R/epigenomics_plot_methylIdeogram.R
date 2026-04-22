## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.

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

epigenomics_plot_methylIdeogram_server <- function(id,
                                                   pgx,
                                                   r.chromosome = reactive(""),
                                                   r.samples = reactive(""),
                                                   r.groupby = reactive(""),
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
      chromosomes <- r.chromosome()
      samples <- r.samples()
      if (!all(samples %in% colnames(X))) return(NULL)
      return(list(
        X = X[kk, samples, drop = FALSE],
        annot = genes[kk, , drop = FALSE],
        chromosomes = chromosomes
      ))    

    })

    plot.RENDER <- function() {
      shiny::req(plot_data())
      X <- plot_data()[["X"]]
      annot <- plot_data()[["annot"]]
      chromosomes <- plot_data()[["chromosomes"]]
      shiny::req(X, annot)
      playbase::plotMethylIdeogram(X, annot, chromosomes = chromosomes)
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
