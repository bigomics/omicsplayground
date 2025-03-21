##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

FeatureMapBoard <- function(id, pgx, labeltype = shiny::reactive("feature"),
                            board_observers = NULL ) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    fullH <- 800 ## full height of page
    rowH1 <- 220 ## row 1 height
    rowH2 <- 460 ## row 2 height

    infotext <- tspan("Visually explore and compare expression signatures on UMAP plots. Feature-level clustering is based on pairwise co-expression between genes (or genesets). This is in contrast to sample-level clustering which clusters samples by similarity of their expression profile. Feature-level clustering allows one to detect gene modules, explore gene neighbourhoods, and identify potential drivers, to study the relationships between features.
<br><br>The tabs present Gene and Geneset UMAP dimensionality reduction plots and are computed for gene and geneset features, respectively. The clustering of features is computed using UMAP from either the normalized log-expression matrix (logCPM) or the log-foldchange matrix (logFC), with the covariance as distance metric. The UMAP from the logCPM is the default, but in cases of strong batch/tissue effects the UMAP from the logFC matrix is a better choice. We prefer the covariance distance metric instead of the correlation because it takes the size of the foldchange into account. Doing so, genes that are close together in corners in the outer rim are those with high pairwise covariance, i.e. have high correlation and high FC.
<br><br>The maps can be colored according to the foldchange signature of the group contrasts (i.e. comparisons), or colored by the average relative log-expression according to some phenotype condition. Multiple signatures can then be easily compared by visually inspection of the colors.
", js = FALSE)

    ## ========================================================================
    ## ======================= OBSERVE FUNCTIONS ==============================
    ## ========================================================================

    my_observers <- list()
    
    # Observer (1):
    my_observers[[1]] <- shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Feature Map Analysis</strong>"),
        shiny::HTML(infotext),
        easyClose = TRUE
      ))
    })

    # Observer (2): tabPanel change to update Settings visibility
    tab_elements <- list(
      "Gene" = list(
        enable = c("filter_genes"),
        disable = c("filter_gsets")
      ),
      "Geneset" = list(
        enable = c("filter_gsets"),
        disable = c("filter_genes")
      )
    )
    my_observers[[2]] <- shiny::observeEvent(input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })

    # Observer (3):
    my_observers[[3]] <- shiny::observeEvent(
      {
        list(pgx$name, pgx$X, pgx$gsetX)
      },
      {
        shiny::req(pgx$X, pgx$gsetX)
        dbg("[FeatureMapBoard] set families and geneset db..")

        ## set gene families
        families <- names(playdata::FAMILIES)
        families <- c("<custom>", families)
        shiny::updateSelectInput(session, "filter_genes",
          choices = families,
          selected = "<all>"
        )

        ## set geneset categories
        gsetcats <- sort(unique(gsub(":.*", "", rownames(pgx$gsetX))))
        gsetcats <- c("<all>", gsetcats)
        sel0 <- grep("^H$|hallmark", gsetcats, ignore.case = TRUE, value = TRUE)
        sel0 <- "<all>"

        if (length(sel0) == 0) sel0 <- 1
        shiny::updateSelectInput(session, "filter_gsets",
          choices = gsetcats, selected = sel0
        )
        shiny::updateTextAreaInput(session, "customlist", placeholder = tspan("Paste your custom gene list", js = FALSE))
      }
    )

    my_observers[[4]] <- observeEvent(
      {
        list(input$sigvar, pgx$samples)
      },
      {
        shiny::req(pgx$samples, input$sigvar)
        if (input$sigvar %in% colnames(pgx$samples)) {
          y <- setdiff(pgx$samples[, input$sigvar], c(NA))
          y <- c("<average>", sort(unique(y)))
          shiny::updateSelectInput(session, "ref_group", choices = y)
        }
      }
    )
    my_observers[[5]] <- observeEvent(
      {
        list(pgx$samples, pgx$X, input$showvar)
      },
      {
        shiny::req(pgx$samples)
        shiny::req(pgx$X)
        shiny::req(input$showvar)

        if (input$showvar == "phenotype") {
          cvar <- playbase::pgx.getCategoricalPhenotypes(pgx$samples, max.ncat = 99)
          cvar0 <- head(grep("^[.]", cvar, invert = TRUE, value = TRUE), 1)
          shiny::updateSelectInput(session, "sigvar", choices = cvar, selected = cvar0)
        }
        if (input$showvar == "comparisons") {
          cvar <- colnames(pgx$model.parameters$contr.matrix)
          sel.cvar <- head(cvar, 8)
          shiny::updateSelectizeInput(session, "selcomp",
            choices = cvar,
            selected = sel.cvar
          )
          shiny::updateSelectInput(session, "ref_group", choices = "  ")
        }
      }
    )
    my_observers[[6]] <- observeEvent(
      {
        list(pgx$samples, input$showvar, input$sigvar)
      },
      {
        shiny::req(pgx$samples, input$sigvar, input$showvar)
        if (input$sigvar %in% colnames(pgx$samples)) {
          y <- setdiff(pgx$samples[, input$sigvar], c(NA))
          y <- c("<average>", sort(unique(y)))
          shiny::updateSelectInput(session, "ref_group", choices = y)
        }
      }
    )
    my_observers[[7]] <- observeEvent(
      {
        list(input$selcomp)
      },
      {
        ## shiny::req(pgx$samples, input$sigvar, input$showvar)
        shiny::updateSelectInput(session, "ref_group", choices = " ")
      }
    )

    ## add to list global of observers. suspend by default.
    my_observers <- my_observers[!sapply(my_observers,is.null)]
    # lapply( my_observers, function(b) b$suspend() )
    if(!is.null(board_observers)) board_observers[[id]] <- my_observers
    
    ## =========================================================================
    ## ============================= FUNCTIONS =================================
    ## =========================================================================

    plotUMAP <- function(pos, var, hilight = NULL, nlabel = 20, title = "",
                         labels = NULL, zlim = NULL, cex = 0.9, cex.label = 1,
                         source = "", plotlib = "base", ...) {
      opc.low <- 1
      if (!is.null(hilight) && !all(rownames(pos) %in% hilight)) {
        opc.low <- 0.2
      }

      hilight2 <- NULL
      if (!is.null(hilight)) {
        hilight <- hilight[hilight %in% names(var)]
        hilight <- hilight[order(-abs(var[hilight]))]
        if (min(var, na.rm = TRUE) < 0) {
          hilight2 <- c(head(hilight, nlabel / 2), tail(hilight, nlabel / 2))
          hilight2 <- hilight2[!duplicated(hilight2)]
        } else {
          hilight2 <- head(hilight, nlabel)
        }
      }

      if (length(hilight) > 0.33 * length(var)) hilight <- hilight2 ## ??? IK
      if (length(hilight) == 0) {
        hilight <- NULL
        hilight2 <- NULL
      }
      cexlab <- ifelse(length(hilight2) <= 20, 1, 0.85)
      cexlab <- ifelse(length(hilight2) <= 8, 1.15, cexlab)
      opacity <- ifelse(length(hilight2) > 0, 0.4, 0.90)
      if (plotlib == "plotly") opacity <- sqrt(opacity) ## less opacity..

      p <- playbase::pgx.scatterPlotXY(
        pos,
        var = var,
        plotlib = plotlib,
        softmax = TRUE,
        cex.lab = 1.2 * cex.label * cexlab,
        opacity = opacity,
        cex = cex,
        zsym = (min(var, na.rm = TRUE) < 0),
        zlim = zlim,
        hilight.cex = cex,
        hilight.col = NULL,
        hilight.lwd = 0.8,
        hilight = hilight,
        hilight2 = hilight2,
        labels = labels,
        title = title,
        source = source,
        key = rownames(pos),
        ...
      )

      p
    }

    plotFeaturesPanel <- function(pos, F, ntop, nr, nc, sel, progress) {
      par(mar = c(1.6, 1.5, 0.5, 0), oma = c(1, 1, 0, 0) * 2)
      par(mar = c(1.1, 1.0, 0.5, 0), oma = c(1, 1, 0, 0) * 2)
      par(mgp = c(1.35, 0.5, 0), las = 0, cex.axis = 0.85, cex.lab = 0.9, xpd = TRUE)
      cex <- ifelse(nc > 3, 0.5, 0.7)
      jj <- 1:nrow(F)
      if (ncol(F) > 4 && nrow(F) > 8000) jj <- sample(1:nrow(F), 8000) ## subsample for speed
      if (ncol(F) > 9 && nrow(F) > 4000) jj <- sample(1:nrow(F), 4000) ## subsample for speed
      if (ncol(F) > 16 && nrow(F) > 2000) jj <- sample(1:nrow(F), 2000) ## subsample for speed

      ## global zlim
      qq <- quantile(F, probs = c(0.05, 0.95), na.rm = TRUE)
      qq <- quantile(F, probs = c(0.01, 0.99), na.rm = TRUE)
      qq <- quantile(F, probs = c(0.002, 0.998), na.rm = TRUE)
      c(min(qq, na.rm = TRUE), max(qq, na.rm = TRUE))
      qq
      zlim <- qq

      i <- 1
      for (i in 1:ncol(F)) {
        if (!interactive()) progress$inc(1 / ncol(F))
        var <- F[, i]
        var <- var[match(rownames(pos), names(var))]
        zsym <- ifelse(min(var, na.rm = TRUE) >= 0, FALSE, TRUE)
        hmarks <- NULL
        if (!is.null(sel)) {
          hmarks <- intersect(sel, names(var))
          hmarks <- head(hmarks[order(var[hmarks])], ntop)
        }
        opacity <- ifelse(is.null(hmarks), 0.9, 0.4)
        xlab <- "UMAP-x"
        ylab <- "UMAP-y"
        xaxs <- yaxs <- TRUE
        if (i %% nc != 1) {
          ylab <- ""
          yaxs <- FALSE
        }
        if ((i - 1) %/% nc < (nr - 1)) {
          xlab <- ""
          xaxs <- FALSE
        }

        playbase::pgx.scatterPlotXY.BASE(
          pos[jj, ],
          var = var[jj],
          zsym = zsym,
          zlim = zlim,
          set.par = FALSE,
          softmax = 1,
          cex = cex,
          cex.legend = 0.9,
          cex.lab = 1.2,
          bty = "n",
          dlim = c(0.05, 0.05),
          hilight = hmarks,
          hilight2 = NULL,
          hilight.col = NULL,
          opacity = opacity,
          xlab = "", ylab = "",
          xaxs = xaxs,
          yaxs = yaxs,
          hilight.lwd = 0.5,
          hilight.cex = 1.3
        )

        cex1 <- ifelse(ncol(F) <= 16, 1.2, 1)
        title(colnames(F)[i], cex.main = cex1, line = -0.75)
        mtext(xlab, 1, line = 1.5, cex = 0.6)
        mtext(ylab, 2, line = 1.6, cex = 0.6)
      }
    }

    filteredGenes <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::validate(need(input$filter_genes, tspan("Please input at least one value in Annotate genes!", js = FALSE)))
      sel <- input$filter_genes
      filtgenes <- c()
      if (is.null(pgx$version) | pgx$organism == "Human") {
        filtgenes <- unlist(lapply(sel, function(genes) playdata::FAMILIES[[genes]]))
      } else {
        filtgenes <- unlist(lapply(sel, function(genes) {
          if (genes == "<all>") {
            x <- pgx$genes$symbol
          } else {
            x <- playdata::FAMILIES[[genes]]
            x <- pgx$genes$symbol[match(x, pgx$genes$human_ortholog, nomatch = 0)]
          }
          return(x)
        }))
      }
      if ("<custom>" %in% sel) {
        genes <- strsplit(input$customlist, split = "[, ;]")[[1]]
        if (length(genes) > 0) {
          filtgenes <- c(playbase::filterProbes(pgx$genes, genes), filtgenes)
        }
      }
      filtgenes
    })

    ## =========================================================================
    ## =========================== MODULES =====================================
    ## =========================================================================

    sigvar2 <- shiny::reactive({
      if (input$showvar == "phenotype") {
        return(input$sigvar)
      }
      if (input$showvar == "comparisons") {
        return(input$selcomp)
      }
    })

    # Gene Map
    featuremap_plot_gene_map_server(
      "geneUMAP",
      pgx = pgx,
      plotUMAP = plotUMAP,
      sigvar = sigvar2,
      filteredGenes = filteredGenes,
      watermark = WATERMARK,
      labeltype = labeltype
    )

    # Gene Signatures
    featuremap_plot_gene_sig_server(
      "geneSigPlots",
      pgx = pgx,
      sigvar = sigvar2,
      ref_group = shiny::reactive(input$ref_group),
      plotFeaturesPanel = plotFeaturesPanel,
      watermark = WATERMARK
    )

    # Geneset map
    featuremap_plot_table_geneset_map_server(
      "gsetUMAP",
      pgx = pgx,
      plotUMAP = plotUMAP,
      filter_gsets = shiny::reactive(input$filter_gsets),
      sigvar = sigvar2,
      watermark = WATERMARK
    )

    # Geneset signatures
    featuremap_plot_gset_sig_server(
      "gsetSigPlots",
      pgx = pgx,
      sigvar = sigvar2,
      ref_group = shiny::reactive(input$ref_group),
      plotFeaturesPanel = plotFeaturesPanel,
      watermark = WATERMARK
    )
  }) ## end of serverModule
} ## end of Board

## ========================================================================
## ========================= END OF FILE ==================================
## ========================================================================
