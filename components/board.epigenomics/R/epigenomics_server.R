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
      grps <- c("<ungrouped>", pheno)
      shiny::updateSelectInput(session, "select_pheno", choices = grps, selected = "<ungrouped>")
      ## grps <- colnames(Y)
      ## grps <- c(
      ##   grep("^[.]", grps, value = TRUE, invert = TRUE),
      ##   grep("^[.]", grps, value = TRUE)
      ## )
      ## grouped_choices <- setNames(
      ##   lapply(grps, function(p) {
      ##     kk <- sort(unique(as.character(na.omit(Y[[p]]))))
      ##     setNames(paste0(p, "::", kk), kk)
      ##   }),
      ##   grps
      ## )
      ## all_choices <- c(list("<ungrouped>" = "<ungrouped>"), grouped_choices)
      ## selgrp <- "<ungrouped>"
      ## if (nrow(Y) > 20) {
      ##   if (length(grps) > 0) {
      ##     kk <- sort(unique(as.character(na.omit(Y[[grps[1]]])))[1])
      ##     selgrp <- paste0(grps[1], "::", kk)
      ##   }
      ##   if ("group" %in% grps) {
      ##     kk <- sort(unique(as.character(na.omit(Y[["group"]]))))[1]
      ##     selgrp <- paste0("group::", kk)
      ##   }
      ##   if ("condition" %in% grps) {
      ##     kk <- sort(unique(as.character(na.omit(Y[["condition"]]))))[1]
      ##     selgrp <- paste0("condition::", kk)
      ##   }
      ## }    
      ## shiny::updateSelectInput(session, "select_pheno", choices = all_choices, selected = selgrp)
    })

    chromosomes <- shiny::reactive({
      chroms <- input$select_chromosome
      if (!is.null(chroms)) chroms <- sort(as.numeric(unique(sub("^chr", "", chroms))))
      validate(
        need(length(chroms) > 0, "Select at least 1 chromosome to be plotted."),
        need(length(chroms) < 5, "Select max 4 chroms at a time for better graphics.")
      ) 
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
      r.groupby = shiny::reactive({input$select_pheno)},
      watermark = WATERMARK
    )

  })

}
