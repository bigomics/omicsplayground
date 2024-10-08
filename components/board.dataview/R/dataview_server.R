##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' DataView module server function
#'
#' @description A shiny Module (server code).
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#' @param pgx Reactive expression that provides the input pgx data object
#'
#' @export
DataViewBoard <- function(id, pgx, labeltype = shiny::reactive("feature")) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    rowH <- 355 ## row height of panels
    imgH <- 315 ## height of images
    fullH <- 750 ## full height of panel
    tabH <- 600 ## height of tables

    ## ----------------------------------------------------------------------
    ## More Info (pop up window)
    ## ----------------------------------------------------------------------

    data_infotext <- HTML('
        <center><iframe width="1120" height="630" src="https://www.youtube.com/embed/S32SPINqO8E"
        title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write;
        encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></center>')


    ## ------- observe functions -----------
    shiny::observeEvent(input$board_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>DataView Board</strong>"),
        shiny::HTML(data_infotext),
        easyClose = TRUE, size = "xl"
      ))
    })

    ## update filter choices upon change of data set
    shiny::observe({
      shiny::req(pgx$Y, pgx$samples)

      ## levels for sample filter
      levels <- playbase::getLevels(pgx$Y)
      shiny::updateSelectInput(session, "data_samplefilter", choices = levels)

      grps <- playbase::pgx.getCategoricalPhenotypes(pgx$samples, min.ncat = 2, max.ncat = 999)
      grps <- sort(grps)
      grps <- c(grep("^[.]", grps, value = TRUE, invert = TRUE), grep("^[.]", grps, value = TRUE))
      selgrp <- grps[1]
      grps <- c("<ungrouped>", grps)
      if ("group" %in% grps) selgrp <- "group"
      if ("condition" %in% grps) selgrp <- "condition"
      if (nrow(pgx$samples) <= 20) selgrp <- "<ungrouped>"
      shiny::updateSelectInput(session, "data_groupby", choices = grps, selected = selgrp)
    })

    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Overview" = list(disable = NULL),
      "Sample QC" = list(disable = c("search_gene")),
      "Data table" = list(disable = NULL),
      "Sample information" = list(disable = c("search_gene", "data_groupby", "data_type")),
      "Contrasts" = list(disable = c("search_gene", "data_groupby", "data_type"))
    )
    shiny::observeEvent(input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })


    shiny::observeEvent(
      {
        list(
          input$data_type,
          pgx$X,
          pgx$counts,
          labeltype()
        )
      },
      {
        shiny::req(input$data_type)


        # X should be labelled as features, so rownames(counts) and rownames(x) shoud match (???)
        features <- rownames(pgx$X)
        # if (input$data_type %in% c("counts", "abundance")) {
        #   features <- rownames(pgx$counts)
        # } else {
        #   ## log2CPM
        #   features <- rownames(pgx$X)
        # }

        ## gene filter.
        fc2 <- rowMeans(playbase::pgx.getMetaFoldChangeMatrix(pgx)$fc**2, na.rm = TRUE)
        features <- intersect(names(sort(-fc2)), features) ## most var gene??
        sel.feature <- features[1]
        features <- sort(features)
        p1 <- head(rownames(pgx$genes), 1000)
        p2 <- head(pgx$genes$symbol, 1000)
        by.symbol <- mean(p1 == p2, na.rm = TRUE) > 0.8
        if (0 && !by.symbol) {
          gene <- pgx$genes[match(features, rownames(pgx$genes)), "symbol"]
          feature_gene <- paste0(gene, "_", features)
          names(features) <- feature_gene
          features <- features[order(names(features))]
        }
        i <- match(sel.feature, features)
        features <- c(features[i], features[-i])
        if (length(features) > 1000) {
          features <- c(
            features[1:1000], tspan("(type for more genes...)", js = FALSE),
            features[1001:length(features)]
          )
        }

        names(features) <- playbase::probe2symbol(features, pgx$genes, labeltype(), fill_na = TRUE)

        shiny::updateSelectizeInput(
          session, "search_gene",
          choices = features,
          selected = sel.feature,
          options = list(maxOptions = 1001),
          server = TRUE
        )
      }
    )

    last_search_gene <- reactiveVal()

    input_search_gene <- reactive({
      if (input$search_gene == "" || grepl("type for more", input$search_gene)) {
        gene1 <- last_search_gene()
        return(gene1)
      }
      last_search_gene(input$search_gene)
      return(input$search_gene)
    })


    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    ## get selected samples after sample filtering
    selected_samples <- reactive({
      samples <- colnames(pgx$X)

      if (!is.null(input$data_samplefilter)) {
        samples <- playbase::selectSamplesFromSelectedLevels(
          pgx$Y, input$data_samplefilter
        )
      }
      # validate samples
      validate(need(length(samples) > 0, "No samples remaining after filtering."))
      samples
    })

    #
    dataview_module_geneinfo_server(
      "geneinfo",
      pgx,
      r.gene = reactive(input$search_gene),
      watermark = WATERMARK
    )

    ## first tab ---------------------------------------
    dataview_plot_tsne_server(
      "tsneplot",
      pgx,
      r.gene = reactive(input$search_gene),
      r.samples = selected_samples,
      r.data_type = reactive(input$data_type),
      r.groupby = reactive(input$data_groupby),
      watermark = WATERMARK,
      info = DATAVIEW_TSNE_INFO
    )

    dataview_plot_averagerank_server(
      "averagerankplot",
      pgx,
      r.gene = reactive(input$search_gene),
      r.samples = selected_samples,
      r.data_type = reactive(input$data_type),
      info = DATAVIEW_AVERAGERANK_INFO,
      watermark = WATERMARK
    )

    dataview_plot_correlation_server(
      "correlationplot",
      pgx,
      r.gene = reactive(input$search_gene),
      r.samples = selected_samples,
      info = DATAVIEW_CORRELATION_INFO,
      watermark = WATERMARK
    )

    dataview_plot_tissue_server(
      "tissueplot",
      pgx,
      r.gene = reactive(input$search_gene),
      r.data_type = reactive(input$data_type),
      info = DATAVIEW_TISSUE_INFO,
      watermark = WATERMARK
    )

    dataview_plot_expression_server(
      "expressionplot",
      pgx,
      r.gene = reactive(input$search_gene),
      r.samples = selected_samples,
      r.data_type = reactive(input$data_type),
      r.data_groupby = reactive(input$data_groupby),
      info = DATAVIEW_EXPRESSION_INFO,
      watermark = WATERMARK
    )

    ## second tab -----------------------------------
    dataview_plot_totalcounts_server(
      "counts_total",
      getCountStatistics,
      r.data_type = reactive(input$data_type),
      watermark = WATERMARK
    )

    dataview_plot_boxplot_server(
      "counts_boxplot",
      input,
      getCountStatistics,
      r.data_type = reactive(input$data_type),
      watermark = WATERMARK
    )

    dataview_plot_histogram_server(
      "counts_histplot",
      getCountStatistics,
      watermark = WATERMARK
    )

    dataview_plot_abundance_server(
      "counts_abundance",
      getCountStatistics,
      watermark = WATERMARK
    )

    dataview_plot_genetypes_server(
      "counts_genetypes",
      getCountStatistics,
      watermark = WATERMARK
    )

    ## fourth tab
    dataview_plot_phenoheatmap_server(
      "phenoheatmap",
      pgx,
      r.samples = selected_samples,
      watermark = WATERMARK
    )

    dataview_plot_phenoassociation_server(
      "phenoassociation",
      pgx,
      r.samples = selected_samples,
      watermark = WATERMARK
    )

    ## ================================================================================
    ## ===============================  TABLES ========================================
    ## ================================================================================

    dataview_table_rawdata_server(
      "rawdatatable", pgx,
      r.gene = reactive(input$search_gene),
      r.data_type = reactive(input$data_type),
      r.samples = selected_samples,
      r.groupby = reactive(input$data_groupby),
      scrollY = "calc(100vh - (240px + 140px))"
    )

    dataview_table_samples_server(
      "sampletable", pgx,
      r.samples = selected_samples,
      scrollY = "calc(35vh - 140px)"
    )


    dataview_table_contrasts_server(
      "contrastTable", pgx,
      r.samples = selected_samples,
      scrollY = "calc(100vh - (240px + 140px))"
    )

    ## ================================================================================
    ## ========================= FUNCTIONS ============================================
    ## ================================================================================

    ##    getCountStatistics <- reactiveVal()
    ##    observeEvent(
    getCountStatistics <- eventReactive(
      {
        list(
          pgx$X,
          input$data_groupby,
          input$data_samplefilter
        )
      },
      {
        shiny::req(pgx$X, pgx$Y, pgx$samples)
        shiny::req(input$data_groupby)
        shiny::validate(shiny::need("counts" %in% names(pgx), "no 'counts' in object."))

        subtt <- NULL
        samples <- colnames(pgx$X)
        samples <- playbase::selectSamplesFromSelectedLevels(pgx$Y, input$data_samplefilter)
        nsamples <- length(samples)
        counts <- pgx$counts[, samples, drop = FALSE]

        grpvar <- input$data_groupby
        if (grpvar != "<ungrouped>") {
          gr <- pgx$Y[samples, grpvar]
          grps <- sort(unique(gr))
          ## check if there are any samples remaining in the group after filtering
          if (length(grps) > 1) {
            mx <- tapply(samples, gr, function(ii) {
              rowMeans(counts[, ii, drop = FALSE], na.rm = TRUE)
            })
            counts <- do.call(cbind, mx)
          }
        }

        ## if too many samples (like scRNA-seq do subsampling...)
        if (ncol(counts) > 500) {
          kk <- sample(ncol(counts), 400, replace = TRUE)
          counts <- counts[, kk, drop = FALSE]
          subtt <- c(subtt, "random subset")
        }
        colnames(counts) <- substring(colnames(counts), 1, 24)

        gset <- list()
        gg <- pgx$genes[rownames(counts), ]$symbol
        tt <- pgx$genes[rownames(counts), ]$gene_title
        g1 <- gg[grep("^rpl|^rps", gg, ignore.case = TRUE)]
        g2 <- gg[grep("^mrpl|^mrps", gg, ignore.case = TRUE)]
        g3 <- gg[grep("^MT-", gg, ignore.case = TRUE)]
        g4 <- gg[grep("mitochondr", tt, ignore.case = TRUE)]
        gset[["Ribosomal (RPL/RPS)"]] <- g1
        gset[["Mitochondrial ribosomal (MRPL/MRPS)"]] <- g2
        gset[["Mitochondrial (MT)"]] <- g3
        gset[["Other mitochondrial"]] <- setdiff(g4, g3)
        jj <- grep("mitochondr|ribosom", names(playdata::FAMILIES), invert = TRUE, ignore.case = TRUE)
        gset.other <- lapply(playdata::FAMILIES[jj], function(x) setdiff(x, c(g1, g2, g3, g4)))
        gset <- c(gset, gset.other)
        gset <- gset[grep("<all>", names(gset), invert = TRUE)]
        gset <- gset[sapply(gset, length) > 10]

        ## Counts per sample or group
        total.counts <- Matrix::colSums(counts, na.rm = TRUE)
        summed.counts <- t(sapply(gset, function(f) {
          Matrix::colSums(counts[which(gg %in% f), , drop = FALSE], na.rm = TRUE)
        }))
        prop.counts <- 100 * t(t(summed.counts) / total.counts)

        ## N. of detected features per sample or group AZ
        n.detected.features <- apply(counts, 2, function(x) length(which(x > 0)))

        ## get variation per group
        log2counts <- log2(1 + counts)
        varx <- apply(log2counts, 1, var)
        gset.var <- sapply(gset, function(s) mean(varx[s], na.rm = TRUE))
        gset.var
        tail(sort(gset.var), 10)

        ## sort get top 20 gene families
        jj <- head(order(-rowSums(prop.counts, na.rm = TRUE)), 20)
        prop.counts <- prop.counts[jj, , drop = FALSE]
        gset <- gset[rownames(prop.counts)]

        gset.genes <- sapply(gset, function(gg) {
          gg <- strwrap(paste(c(head(gg, 20), "+ ..."), collapse = " "), 40)
          paste(gg, collapse = "<br>")
        })

        ## align
        ss <- names(total.counts)
        prop.counts <- prop.counts[, ss, drop = FALSE]
        counts <- counts[, ss, drop = FALSE]
        if (any(pgx$X[, samples, drop = FALSE] < 0, na.rm = TRUE)) {
          offset <- 1e-6
        } else {
          offset <- 1
        }
        log2counts <- log2(offset + counts)

        names(total.counts) <- substring(names(total.counts), 1, 30)
        colnames(log2counts) <- substring(colnames(log2counts), 1, 30)
        colnames(prop.counts) <- substring(colnames(prop.counts), 1, 30)

        res <- list(
          total.counts = total.counts,
          subtt = subtt,
          log2counts = log2counts,
          prop.counts = prop.counts,
          n.detected.features = n.detected.features, ## AZ
          gset.genes = gset.genes
        )
        ## getCountStatistics(res)
        return(res)
      },
      ignoreNULL = TRUE
    )

    ## ================================================================================
    ## ================================= END ====================================
    ## ================================================================================
  })
}
