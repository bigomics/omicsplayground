##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ConsensusWGCNA_Board <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    fullH <- 700 ## full height of page
    rowH1 <- 250 ## row 1 height
    rowH2 <- 440 ## row 2 height

    infotext <- tspan("<b>Multi-Omics WGCNA</b> ...
<p>References:<br>
<ol>
<li>Langfelder, P. and Horvath, S., 2008. WGCNA: an R package for weighted correlation network analysis. BMC bioinformatics, 9(1), p.559.
<li>Zhang, B. and Horvath, S., 2005. A general framework for weighted gene co-expression network analysis. Statistical applications in genetics and molecular biology, 4(1).
</ol>
", js = FALSE)


    ## ============================================================================
    ## ============================ OBSERVERS =====================================
    ## ============================================================================

    infotext <-
      '<center><iframe width="1120" height="630" src="https://www.youtube.com/embed/rRIRMW_RRS4"
        title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write;
        encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></center>'

    shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Multi-Omics WGCNA Board</strong>"),
        shiny::HTML(infotext),
        size = "xl",
        easyClose = TRUE
      ))
    })

    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Dendrograms" = list(disable = c("trait", "module")),
      "Sample Clustering" = list(disable = c("trait", "module")),
      "Module-Trait" = list(disable = c("module")),
      "Feature Table" = list(disable = c())
    )

    shiny::observeEvent(input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })

    ## ============================================================================
    ## ============================ REACTIVES =====================================
    ## ============================================================================

    shiny::observe({
      cons <- r_wgcna()
      shiny::validate(shiny::need(!is.null(cons), "Please compute"))
    })

    shiny::observeEvent(list(pgx$X, pgx$samples), {
      splitby <- colnames(pgx$samples)
      if (pgx$datatype == "multi-omics") {
        dtypes <- names(playbase::mofa.split_data(pgx$X))
        has.pxgx <- all(c("gx", "px") %in% dtypes)
        if (has.pxgx) {
          splitby <- c("<multi-omics>", splitby)
        }
      }

      shiny::updateSelectInput(session, "splitby",
        choices = splitby,
        selected = splitby[1]
      )
    })


    r_wgcna <- shiny::eventReactive(
      {
        list(input$compute, pgx$X)
      },
      {
        shiny::req(pgx$X)
        shiny::req(input$splitby)

        xx <- NULL
        splitby <- input$splitby
        if (pgx$datatype == "multi-omics" && splitby == "<multi-omics>") {
          xx <- playbase::mofa.split_data(pgx$X)
          has.gxpx <- all(c("gx", "px") %in% names(xx))
          has.gxpx
          shiny::validate(shiny::need(
            has.gxpx,
            "Your multi-omics dataset is incompatible for consensus WGCNA: both transcriptomics & proteomics data are needed."
          ))

          xx <- xx[names(xx) %in% c("gx", "px")]
          xx <- lapply(xx, function(x) playbase::rename_by2(x, annot_table = pgx$genes))
          gg <- Reduce(intersect, lapply(xx, rownames))
          shiny::validate(
            shiny::need(
              length(gg) > 0,
              "Your dataset is incompatible for consensus WGCNA: no shared features between transcriptomics and proteomics."
            ),
            shiny::need(
              length(gg) >= 50,
              paste0(
                "Only ", length(gg), " shared features found between transcriptomics and proteomics layers ",
                "(minimum required: 50). ",
                "Ensure both layers use compatible feature identifiers (gene symbols, UniProt, or Ensembl IDs)."
              )
            )
          )
          xx <- lapply(xx, function(x) x[gg, , drop = FALSE])
        } else if (!is.null(pgx$samples)) {
          splitby <- input$splitby
          if (is.null(splitby) || splitby == "") {
            splitby <- colnames(pgx$samples)[1]
          }
          shiny::req(splitby %in% colnames(pgx$samples))
          group <- pgx$samples[, splitby]
          if (is.numeric(group) && length(unique(group)) > 3) {
            group <- c("LO", "HI")[1 + (group >= median(group, na.rm = TRUE))]
          }
          if (max(nchar(group)) > 10) {
            group <- base::abbreviate(toupper(group), 4L)
          }
          xx <- tapply(1:ncol(pgx$X), group, function(ii) pgx$X[, ii, drop = FALSE])
        } else {
          shiny::validate(shiny::need(!is.null(xx), "Your dataset is incompatible for consensus WGCNA."))
        }

        ## exclude sample matrices with less than 4 samples.
        ## WGCNA::blockwiseConsensusModules fails when nsamples<4
        xx <- xx[sapply(xx, function(x) ncol(x) >= 4)]
        shiny::validate(shiny::need(
          length(xx) > 1,
          "Your selected phenotype is incompatible for consensus WGCNA: less than 4 samples available for any given phenotype level. Please select another trait from the 'Consensus by' menu."
        ))
        samples <- unique(unlist(lapply(xx, colnames)))
        phenoData <- pgx$samples[samples, , drop = FALSE]
        contrasts <- pgx$contrasts[samples, , drop = FALSE]

        ## random noise: avoids ME with NaNs; increase robustness
        for (i in 1:length(xx)) {
          mat <- xx[[i]]
          sdx0 <- matrixStats::rowSds(mat, na.rm = TRUE)
          sdx1 <- 0.05 * sdx0 + 0.5 * mean(sdx0, na.rm = TRUE)
          xx[[i]] <- mat + 0.1 * sdx1 * matrix(rnorm(length(mat)), nrow(mat), ncol(mat))
        }

        progress <- shiny::Progress$new(session, min = 0, max = 1)
        on.exit(progress$close())
        progress$set(message = paste("computing consensus WGCNA..."), value = 0.33)
        pgx.showSmallModal("computing consensus WGCNA...")

        power <- input$power
        if (power == "<auto>") {
          power <- NULL
        } else {
          power <- as.numeric(power)
        }

        ## This runs consensus WGCNA on an expression list
        ngenes <- as.integer(input$ngenes)
        minModuleSize <- as.integer(input$minmodsize)
        deepSplit <- as.integer(input$deepsplit)

        cons <- playbase::wgcna.runConsensusWGCNA(
          exprList = xx,
          phenoData = phenoData,
          contrasts = contrasts,
          GMT = pgx$GMT,
          annot = pgx$genes,
          power = power,
          ngenes = ngenes,
          minModuleSize = minModuleSize,
          maxBlockSize = 9999,
          minKME = 0.3,
          mergeCutHeight = 0.15,
          deepSplit = deepSplit,
          calcMethod = "fast",
          drop.ref = FALSE,
          addCombined = FALSE,
          compute.stats = TRUE,
          compute.enrichment = TRUE,
          summary = TRUE,
          ai_model = NULL,
          experiment = pgx$description,
          gsea.mingenes = 5,
          gsea.ntop = 1000,
          progress = progress,
          verbose = 1
        )

        shiny::removeModal()

        all_modules <- rownames(cons$modTraits)
        module1 <- all_modules[[1]][1]
        updateSelectInput(session, "module",
          choices = sort(all_modules),
          selected = module1
        )

        # traits <- colnames(cons$datTraits)
        traits <- colnames(cons$stats[[1]][["moduleTraitCor"]])
        updateSelectInput(session, "trait",
          choices = sort(traits),
          selected = traits[1]
        )

        return(cons)
      },
      ignoreNULL = FALSE

    )


    ## ==========================================================================
    ## ========================== BOARD FUNCTIONS ===============================
    ## ==========================================================================


    ## ==========================================================================
    ## =========================== MODULES ======================================
    ## ==========================================================================

    consensusWGCNA_plot_dendrograms_server(
      id = "consensusWGCNADendro",
      mwgcna = r_wgcna
    )

    consensusWGCNA_plot_power_server(
      id = "consensusWGCNAPower",
      mwgcna = r_wgcna
    )

    consensusWGCNA_plot_moduletrait_server(
      "consensusWGCNATrait",
      mwgcna = r_wgcna
    )

    consensusWGCNA_plot_sampletree_server(
      "consensusWGCNASampleTree",
      mwgcna = r_wgcna
    )

    consensusWGCNA_table_modulegenes_server(
      id = "consensusWGCNATable",
      mwgcna = r_wgcna,
      r_annot = reactive(pgx$genes),
      r_trait = reactive(input$trait),
      r_module = reactive(input$module)
    )

    consensusWGCNA_table_enrichment_server(
      id = "consensusWGCNAEnrichment",
      mwgcna = r_wgcna,
      r_module = reactive(input$module)
    )

    consensusWGCNA_plot_preservation_server(
      id = "consensusWGCNAPreservation",
      mwgcna = r_wgcna
    )

    # Enrichment plot
    wgcna_html_module_summary_server(
      "consensusWGCNAmoduleSummary",
      wgcna = r_wgcna,
      multi = FALSE,
      r_module = shiny::reactive(input$module),
      watermark = WATERMARK
    )

    consensusWGCNA_plot_traitsignificance_server(
      id = "consensusWGCNATraitSignificance",
      rwgcna = r_wgcna,
      rtrait = reactive(input$trait)
    )

    return(NULL)
  })
} ## end of Board
