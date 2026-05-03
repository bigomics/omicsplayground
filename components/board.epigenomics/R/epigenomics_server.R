## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.

EpigenomicsBoard <- function(id, pgx) {

  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    shiny::observeEvent(input$board_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Epigenomics Board</strong>"),
        shiny::HTML("Epigenomics visualizations and analyses for methylomics data."),
        easyClose = TRUE,
        size = "l"
      ))
    })

    shiny::observe({
      shiny::req(pgx$X, pgx$samples, pgx$genes)
      kk <- grep("chr|chromosome|chromosomes|chrom|chroms", tolower(colnames(pgx$genes)))
      chroms <- unique(na.omit(sub("(p|q|cen).*", "", as.character(pgx$genes[, kk[1]]))))
      chroms <- unique(sub("^chr", "", chroms))
      sex.chr <- intersect(c("X", "Y"), chroms)
      chroms <- paste0("chr", c(sort(as.numeric(setdiff(chroms, sex.chr))), sex.chr))
      shiny::updateSelectizeInput(session, "select_chromosome",
        choices = chroms, selected = chroms[1:4], server = T)
      shiny::updateSelectizeInput(session, "search_gene", choices = rownames(pgx$X), server = T)
      shiny::updateSelectInput(session, "data_samplefilter", choices = colnames(pgx$X))
      Y <- pgx$samples
      pheno <- colnames(Y)[sapply(Y, function(v) any(!is.na(v) & v != ""))]
      pheno <- pheno[!grepl("cell_cycle", pheno, ignore.case = TRUE)]
      grps <- c("<ungrouped>", pheno)
      shiny::updateSelectInput(session, "select_pheno", choices = grps, selected = "<ungrouped>")
    })

    chromosomes <- shiny::reactive({
      chroms <- input$select_chromosome
      if (!is.null(chroms)) chroms <- sort(as.numeric(unique(sub("^chr", "", chroms))))
      validate(
        need(length(chroms) > 0, "Select at least 1 chromosome to be plotted."),
        need(length(chroms) < 7, "Select max 6 chroms at a time for better graphics.")
      )
      return(chroms)
    })

    chromosomes_all <- shiny::reactive({
      chroms <- input$select_chromosome
      if (!is.null(chroms)) chroms <- sort(as.numeric(unique(sub("^chr", "", chroms))))
      validate(need(length(chroms) > 0, "Select at least 1 chromosome to be plotted."))
      return(chroms)
    })

    samples <- shiny::reactive({
      samples <- rownames(pgx$samples)
      if (!is.null(input$data_samplefilter)) {
        kk <- input$data_samplefilter
        kk <- kk[which(!is.na(kk) & kk != "")]
        kk <- intersect(samples, kk)
        if (length(kk) > 0) samples <- samples[!samples %in% kk]
      }
      validate(need(length(samples) > 0, "No samples remaining after filtering."))
      return(samples)
    })

    epigenomics_plot_methylIdeogram_server(
      "methylIdeogram",
      pgx,
      r.chromosome = chromosomes,
      r.samples = samples,
      r.pheno = shiny::reactive({input$select_pheno}),
      watermark = WATERMARK
    )

    dataview_table_beta_server(
      "methyltable",
      pgx,
      r.samples = samples,
      r.pheno = shiny::reactive({input$select_pheno}),
      scrollY = "30vh"
    )
    
    epigenomics_plot_beta_dist_server(
      "betaDist",
      pgx,
      r.samples = samples,
      r.pheno = shiny::reactive({input$select_pheno}),
      watermark = WATERMARK
    )

    epigenomics_plot_boxplot_beta_server(
      "boxplotBeta",
      pgx,
      r.chromosome = chromosomes_all,
      r.samples = samples,
      r.pheno = shiny::reactive({input$select_pheno}),
      watermark = WATERMARK
    )

  })

}
