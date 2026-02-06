##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


SingleCellBoard <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    fullH <- 750 ## full height of panel
    imgH <- 680 ## row height of panel
    tabH <- 200 ## row height of panel

    infotext <- paste0(
      tspan("The <strong>Cell Profiling Board</strong> infers the type of cells using computational deconvolution methods and reference datasets from the literature. Currently, we have implemented a total of 8 methods and 9 reference datasets to predict immune cell types (4 datasets), tissue types (2 datasets), cell lines (2 datasets) and cancer types (1 dataset). However, we plan to expand the collection of methods and databases and to infer other cell types.

<br><br>The <strong>Proportions tab</strong> visualizes the interrelationships between two categorical variables (so-called cross tabulation). Although this feature is very suitable for a single-cell sequencing data, it provides useful information about the proportion of different cell types in samples obtained by the bulk sequencing method.

<br><br>For each combination of gene pairs, the platform can generate a cytometry-like plot of samples under the <strong>Cytoplot</strong> tab. The aim of this feature is to observe the distribution of samples in relation to the selected gene pairs. For instance, when applied to single-cell sequencing data from immunological cells, it can mimic flow cytometry analysis and distinguish T helper cells from the other T cells by selecting the CD4 and CD8 gene combination.

<br><br>The <strong>Markers</strong> section provides potential marker genes, which are the top N=36 genes with the highest standard deviation within the expression data across the samples. For every gene, it produces a t-SNE plot of samples, with samples colored in red when the gene is overexpressed in corresponding samples. Users can also restrict the marker analysis by selecting a particular functional group in which genes are divided into 89 groups, such as chemokines, transcription factors, genes involved in immune checkpoint inhibition, and so on.

<br><br>It is also possible to perform a copy number variation analysis under the <strong>CNV tab</strong>. The copy number is estimated from gene expression data by computing a moving average of the relative expression along the chromosomes. CNV creates a heatmap of samples versus chromosomes, where samples can be annotated further with a phenotype class provided in the data.", js = FALSE),
      '<center><iframe width="560" height="315" src="https://www.youtube.com/embed/BtMQ7Y0NoIA?si=pebNlzthvdZF7h5o&amp;start=35" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe></center>'
    )

    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    shiny::observeEvent(input$infotext, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Single Cell Board</strong>"),
        shiny::HTML(infotext),
        easyClose = TRUE, size = "l"
      ))
    })

    ## update filter choices upon change of data set
    shiny::observe({
      shiny::req(pgx$X, pgx$Y)
      ## levels for sample filter
      levels <- playbase::getLevels(pgx$Y)
      shiny::updateSelectInput(session, "samplefilter", choices = levels)

      ## update cluster methods if available in object
      if ("cluster" %in% names(pgx)) {
        clustmethods <- names(pgx$cluster$pos)
        clustmethods <- grep("3d", clustmethods, value = TRUE, invert = TRUE, ignore.case = TRUE) ## no 3D
        clustmethods <- c("default", clustmethods)
        shiny::updateSelectInput(session, "clustmethod",
          choices = clustmethods
        )
      }
    })

    shiny::observe({
      shiny::req(pgx$X)
      refsets <- "LM22"
      refsets <- sort(names(pgx$deconv))
      refsel <- unique(c(grep("LM22", refsets, value = TRUE), refsets))[1]
      shiny::updateSelectInput(session, "refset", choices = refsets, selected = refsel)
      shiny::updateSelectInput(session, "refset2", choices = refsets, selected = refsel)
      grpvars <- c("<ungrouped>", colnames(pgx$samples))
      sel <- grpvars[1]
      if (ncol(pgx$X) > 30) sel <- grpvars[2]
      shiny::updateSelectInput(session, "group2", choices = grpvars, selected = sel)
    })

    shiny::observeEvent(input$refset, {
      shiny::req(input$refset)
      shiny::req(pgx$X)

      dcmethods <- names(pgx$deconv[[input$refset]])
      dcsel <- intersect(c("meta.prod", "meta"), dcmethods)[1]
      shiny::updateSelectInput(session, "dcmethod", choices = dcmethods, selected = dcsel)
    })

    shiny::observeEvent(input$refset2, {
      shiny::req(input$refset2)
      dcmethods <- names(pgx$deconv[[input$refset2]])
      dcsel <- intersect(c("meta.prod", "meta"), dcmethods)[1]
      shiny::updateSelectInput(session, "dcmethod2", choices = dcmethods, selected = dcsel)
    })

    shiny::observe({
      shiny::req(pgx$X)

      pheno0 <- grep("group|sample|donor|id|batch", colnames(pgx$samples), invert = TRUE, value = TRUE)
      pheno0 <- grep("sample|donor|id|batch", colnames(pgx$samples), invert = TRUE, value = TRUE)
      kk <- tryCatch(playbase::selectSamplesFromSelectedLevels(pgx$Y, input$samplefilter),
        error = function(w) {
          0
        }
      )
      nphenolevel <- apply(pgx$samples[kk, pheno0, drop = FALSE], 2, function(v) length(unique(v)))
      pheno0 <- pheno0[which(nphenolevel > 1)]
      genes <- sort(as.character(rownames(pgx$X)))
      pheno1 <- c("<cell type>", pheno0)
      genes1 <- c("<none>", genes)
      shiny::updateSelectInput(session, "crosstabvar", choices = pheno1)
      shiny::updateSelectInput(session, "crosstabpheno", choices = pheno1, , selected = pheno1[1])
      names(genes1) <- playbase::probe2symbol(genes1, pgx$genes, "gene_name", fill_na = TRUE)
      shiny::updateSelectizeInput(session, "crosstabgene", choices = genes1, server = TRUE, selected = genes1[2])
    })

    shiny::observe({
      shiny::req(pgx$X, input$mrk_level)

      choices <- names(pgx$families)
      selected <- grep("^CD", choices, ignore.case = TRUE, value = TRUE)[1]
      if (input$mrk_level == "geneset") {
        gset_collections <- playbase::pgx.getGeneSetCollections(gsets = rownames(pgx$gsetX))
        nn <- sapply(gset_collections, function(k) sum(k %in% rownames(pgx$gsetX)))
        choices <- names(gset_collections)[nn >= 5]
        selected <- grep("HALLMARK", names(gset_collections), ignore.case = TRUE, value = TRUE)
      }
      shiny::updateSelectInput(session, "features", choices = choices, selected = selected)
      shiny::updateSelectInput(session, "mrk_features", choices = choices, selected = selected)
    })

    shiny::observe({
      shiny::req(pgx$X)
      ## just at new data load
      genes <- NULL
      g1 <- g2 <- NULL

      F <- playbase::pgx.getMetaFoldChangeMatrix(pgx)$fc
      F <- F[order(-apply(F, 1, sd)), , drop = FALSE]
      genes <- rownames(F)
      g1 <- rownames(F)[1]
      g2 <- rownames(F)[2]

      if (length(g1) == 0) g1 <- genes[1]
      if (length(g2) == 0) g2 <- genes[2]

      ## NOTE: server=TRUE sometime not renders plot. please check. (XM)
      shiny::updateSelectizeInput(session, "cytovar1", choices = genes, selected = g1, server = TRUE)
      shiny::updateSelectizeInput(session, "cytovar2", choices = genes, selected = g2, server = TRUE)
    })

    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Cell type" = list(disable = NULL),
      "Mapping" = list(disable = c("clustmethod")),
      "Markers" = list(disable = NULL),
      "AI Summary" = list(disable = NULL)
    )

    shiny::observeEvent(input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })

    ## ========================================================================
    ## ============================ REACTIVE ==================================
    ## ========================================================================

    pfGetClusterPositions <- shiny::reactive({ # used by many plots
      shiny::req(pgx$X)

      zx <- pgx$X
      kk <- colnames(zx)
      kk <- playbase::selectSamplesFromSelectedLevels(pgx$Y, input$samplefilter)
      if (length(kk) == 0) {
        return(NULL)
      }
      zx <- zx[, kk, drop = FALSE]
      zx <- head(zx[order(-apply(zx, 1, sd)), ], 1000)
      zx <- t(scale(t(zx))) ## scale??
      pos <- NULL
      m <- "tsne"
      m <- input$clustmethod
      has.clust <- ("cluster" %in% names(pgx) && m %in% names(pgx$cluster$pos))
      if (!has.clust && m == "pca") {
        pos <- irlba::irlba(zx, nv = 3)$v
        rownames(pos) <- colnames(zx)
      } else if (has.clust) {
        pos <- pgx$cluster$pos[[m]][, 1:2]
      } else {
        pos <- pgx$tsne2d
      }
      pos <- pos[colnames(zx), ]
      pos <- scale(pos) ## scale
      colnames(pos) <- paste0("dim", 1:ncol(pos))
      rownames(pos) <- colnames(zx)

      return(pos)
    }) %>% bindEvent(
      input$samplefilter,
      pgx$X
    )

    # Type mapping (heatmap) reactivity ##########

    getDeconvResults2 <- shiny::reactive({ # used by many functions
      shiny::req(pgx$X)
      shiny::req(input$dcmethod2)
      shiny::req(input$refset2)
      method <- "meta"
      method <- input$dcmethod2

      refset <- "LM22"
      refset <- input$refset2
      if (!("deconv" %in% names(pgx))) {
        return(NULL)
      }
      if (length(pgx) == 0) {
        return(NULL)
      }

      results <- pgx$deconv[[refset]][[method]]
      ## threshold everything (because DCQ can be negative!!!)
      results <- pmax(results, 0)

      return(results)
    })


    # plots -------------------------------------------------------------------

    singlecell_plot_icpplot_server(
      id = "icpplot",
      pgx = pgx,
      pfGetClusterPositions = pfGetClusterPositions,
      method = shiny::reactive(input$dcmethod),
      refset = shiny::reactive(input$refset),
      layout = shiny::reactive(input$layout),
      sortby = shiny::reactive(input$sortby),
      watermark = WATERMARK
    )

    singlecell_plot_phenoplot_server(
      id = "phenoplot",
      pgx = pgx,
      pfGetClusterPositions = pfGetClusterPositions,
      watermark = WATERMARK
    )

    singlecell_plot_mappingplot_server(
      id = "mappingplot",
      pgx = pgx,
      getDeconvResults2 = getDeconvResults2,
      pfGetClusterPositions = pfGetClusterPositions,
      grpvar = shiny::reactive(input$group2),
      refset = shiny::reactive(input$refset2),
      group = shiny::reactive(input$group2),
      view = shiny::reactive(input$view2),
      watermark = WATERMARK
    )

    singlecell_plot_crosstabPlot_server(
      id = "crosstabPlot",
      pgx = pgx,
      samplefilter = shiny::reactive(input$samplefilter),
      crosstabvar = shiny::reactive(input$crosstabvar),
      pheno = shiny::reactive(input$crosstabpheno),
      gene = shiny::reactive(input$crosstabgene),
      getDeconvResults2 = getDeconvResults2,
      watermark = WATERMARK
    )

    singlecell_plot_markersplot_server(
      id = "markersplot",
      pgx = pgx,
      pfGetClusterPositions = pfGetClusterPositions,
      mrk_level = shiny::reactive(input$mrk_level),
      mrk_features = shiny::reactive(input$mrk_features),
      mrk_search = shiny::reactive(input$mrk_search),
      mrk_sortby = shiny::reactive(input$mrk_sortby),
      watermark = WATERMARK
    )

    singlecell_plot_cytoplot_server(
      id = "cytoplot",
      pgx = pgx,
      pfGetClusterPositions = pfGetClusterPositions,
      samplefilter = shiny::reactive(input$samplefilter),
      cytovar1 = shiny::reactive(input$cytovar1),
      cytovar2 = shiny::reactive(input$cytovar2),
      nbins = shiny::reactive(input$nbins),
      selectSamplesFromSelectedLevels = selectSamplesFromSelectedLevels,
      watermark = WATERMARK
    )

    # AI cell profiling summary

    singlecell_ai_summary_server(
      "singlecellAISummary",
      pgx = pgx,
      refset = shiny::reactive(input$refset),
      dcmethod = shiny::reactive(input$dcmethod),
      session = session,
      watermark = WATERMARK
    )

    return(NULL)
  })
}
