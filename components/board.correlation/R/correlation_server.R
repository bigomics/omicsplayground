##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

CorrelationBoard <- function(id, pgx, labeltype = shiny::reactive("feature"),
                             board_observers = board_observers) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    fullH <- 800 ## full height of page
    rowH <- 340 ## full height of page

    cor_infotext <- tspan("The <strong>Correlation Analysis Board</strong> provides statistical correlation analysis on gene level with visualisations. During the visual analysis, users can filter out some samples or collapse the samples by predetermined groups. The dark shaded area in the barplot estimates the partial correlation.", js = FALSE)

    COL <- RColorBrewer::brewer.pal(12, "Paired")[seq(1, 12, 2)]
    COL <- RColorBrewer::brewer.pal(9, "Set1")[c(2, 1, 3:9)]
    COL2 <- rev(grey.colors(2))
    COL2 <- RColorBrewer::brewer.pal(2, "Paired")[1:2]
    COL2 <- COL[1:2]

    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    my_observers <- list()

    my_observers[[1]] <- shiny::observeEvent(input$data_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Expression: Correlation analysis board</strong>"),
        shiny::HTML(cor_infotext),
        easyClose = TRUE, size = "l"
      ))
    })

    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Correlation" = list(
        enable = NULL,
        disable = c("pcor_ntop")
      ),
      "Graph" = list(
        enable = NULL,
        disable = NULL
      )
    )
    my_observers[[2]] <- shiny::observeEvent(input$tabs1, {
      bigdash::update_tab_elements(input$tabs1, tab_elements)
    })
    
    my_observers[[3]] <- shiny::observeEvent(input$cor_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Correlation Analysis Board</strong>"),
        shiny::HTML(cor_infotext),
        easyClose = TRUE
      ))
    })

    ## update filter choices upon change of data set
    my_observers[[4]] <- shiny::observe({
      req(pgx$X)
      genes <- rownames(pgx$X[complete.cases(pgx$X), ])
      ## genes <- sort(pgx$genes[rownames(pgx$X), ]$gene_name)
      ## genes <- sort(pgx$genes[genes, ]$gene_name)
      sel <- genes[1] ## most var gene
      ## sel <- names(head(sort(-rowMeans(playbase::pgx.getMetaMatrix(pgx)$fc**2)), 1))
      shiny::updateSelectizeInput(
        session, "gene",
        choices = genes, selected = sel, server = TRUE
      )

      fam <- playbase::pgx.getFamilies(pgx, nmin = 10, extended = FALSE)
      fam <- sort(c("<custom>", fam))
      names(fam) <- sub(".*:", "", fam)
      shiny::updateSelectInput(session, "cor_filter", choices = fam)

      px <- colnames(pgx$Y)
      s1 <- grep("^[.]", px, value = TRUE, invert = TRUE)[1]
      shiny::updateSelectInput(session, "cor_group", choices = px, selected = s1)
    })

    my_observers[[5]] <- shiny::observeEvent(pgx$X, {
      shiny::updateTextAreaInput(session, "cor_customfeatures", placeholder = tspan("Paste your custom gene list", js = FALSE))
    })

    ## add to list global of observers. suspend by default.
    my_observers <- my_observers[!sapply(my_observers,is.null)]
    lapply( my_observers, function(b) b$suspend() )
    if(!is.null(board_observers)) board_observers[[id]] <- my_observers
    

    ## =========================================================================
    ## ============================= FUNCTIONS =================================
    ## =========================================================================

    
    filterExpression <- function(X, gene) {
      ## filter genes
      ft <- input$cor_filter
      if (ft == "<custom>" && input$cor_customfeatures != "") {
        genes <- toupper(pgx$genes$symbol)
        gg1 <- strsplit(input$cor_customfeatures, split = "[, ;\n\t]")[[1]]
        if (length(gg1) == 1) gg1 <- paste0(gg1, "*")
        gg1 <- gsub("[ \n\t]", "", gg1)
        starred <- grep("[*]", gg1)
        if (length(starred) > 0) {
          gg2 <- lapply(gg1[starred], function(a) {
            genes[grep(paste0("^", sub("[*]", "", a)), genes, ignore.case = TRUE)]
          })
          gg1 <- unique(c(gg1, unlist(gg2)))
        }
        gg1 <- gg1[which(toupper(gg1) %in% toupper(pgx$genes$symbol))]
        shiny::validate(shiny::need(
          length(gg1) > 1,
          tspan(
            "Custom filtering does not match any gene. Please make sure the genes on the <custom> filter are present on your data.",
            js = FALSE
          )
        ))
        psel <- playbase::filterProbes(pgx$genes, c(gg1, gene))
        psel <- intersect(psel, rownames(X))
        X <- X[psel, , drop = FALSE]
      } else if (ft != "<all>" && ft %in% names(playdata::iGSETS)) {
        ft <- input$cor_filter
        psel <- playbase::filterProbes(
          pgx$genes, c(gene, unlist(playdata::getGSETS(ft)))
        )
        #
        psel <- intersect(psel, rownames(X))
        X <- X[psel, , drop = FALSE]
      }

      X <- X + 1e-3 * matrix(rnorm(length(X)), nrow(X), ncol(X))
      X
    }

    getFilteredExpression <- shiny::reactive({
      shiny::req(pgx$X, input$gene)
      X <- pgx$X
      gene <- rownames(X)[1]
      gene <- input$gene
      zx <- filterExpression(X, gene)
      zx
    })

    getFilteredImpExpression <- shiny::reactive({
      shiny::req(pgx$impX, input$gene)
      X <- pgx$impX
      gene <- input$gene
      filterExpression(X, gene)
    })

    getPartialCorrelationMatrix <- shiny::reactive({
      shiny::req(pgx$X, input$gene)

      gene <- rownames(pgx$X)[1]
      gene <- input$gene

      ## filter gene expression matrix
      X <- getFilteredExpression()
      nmissing <- sum(is.na(X))
      if (nmissing > 0) {
        X <- getFilteredImpExpression()
      }

      NTOP <- 50
      NTOP <- as.integer(input$pcor_ntop)
      res <- playbase::pgx.computeGlassoAroundGene(X, gene, nmax = NTOP)
      res$meta.pcor <- res$pcor

      j <- which(rownames(res$pcor) == gene)
      P <- res$pcor
      diag(P) <- 0
      rho1 <- min(head(sort(P, decreasing = TRUE), 200))
      max1 <- round(max(P), digits = 3)
      shiny::updateSliderInput(session, "cor_graph-cor_graph_threshold", value = rho1, max = max1)
      shiny::updateSliderInput(session, "dcga_graph_threshold", value = rho1, max = max1)

      res
    })

    getPartialCorrelation <- shiny::reactive({
      res <- getPartialCorrelationMatrix()
      gene <- input$gene
      rho <- res$cor[gene, ]
      prho <- res$pcor[gene, ]
      df <- data.frame(cor = rho, pcor = prho)
      df
    })

    getGeneCorr <- shiny::reactive({
      shiny::req(pgx$X)
      gene <- input$gene
      if (is.null(gene)) {
        return(NULL)
      }

      ## corr always in log.scale and restricted to selected samples subset
      zx <- pgx$X
      zx <- getFilteredExpression()

      xref <- list(
        "cor" = 2**zx
        ## "cor.HPA" = as.matrix(playdata::TISSUE),
        ## "cor.ImmProt" = as.matrix(playdata::IMMPROT)
      )
      R <- playbase::pgx.getGeneCorrelation(gene, xref = xref)
      if (is.null(R)) {
        return(NULL)
      }
      cm <- intersect(rownames(R), rownames(zx))
      R <- R[cm, , drop = FALSE]
      zx <- zx[cm, , drop = FALSE]
      zx <- zx - rowMeans(zx, na.rm = TRUE)
      sdx <- sqrt(rowMeans(zx**2, na.rm = TRUE))
      R <- cbind(R, cov = R[, "cor"] * sdx * sdx[gene])
      R <- R[order(R[, "cor"], decreasing = TRUE), , drop = FALSE]
      R
    })

    getFullGeneCorr <- shiny::reactive({
      shiny::req(pgx$X)

      gene <- input$gene
      if (is.null(gene)) {
        return(NULL)
      }

      ## corr always in log.scale and restricted to selected samples subset
      zx <- pgx$X
      xref <- list(
        "cor" = 2**zx
        ##        "cor.HPA" = as.matrix(playdata::TISSUE),
        ##        "cor.ImmProt" = as.matrix(playdata::IMMPROT)
      )
      R <- playbase::pgx.getGeneCorrelation(gene, xref = xref)
      if (is.null(R)) {
        return(NULL)
      }
      # Conform datasets
      gg <- intersect(rownames(zx), rownames(R))
      R <- R[gg, , drop = FALSE]
      zx <- zx[gg, , drop = FALSE]
      zx <- zx - rowMeans(zx, na.rm = TRUE)
      sdx <- sqrt(rowMeans(zx**2, na.rm = TRUE))
      R <- cbind(R, cov = R[, "cor"] * sdx * sdx[gene])
      R <- R[order(R[, "cor"], decreasing = TRUE), , drop = FALSE]
      R
    })

    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    cor_table <- correlation_table_corr_server(
      "cor_table",
      getPartialCorrelation = getPartialCorrelation,
      getGeneCorr           = getGeneCorr,
      pgx                   = pgx,
      watermark             = WATERMARK
    )

    correlation_plot_barplot_server(
      "cor_barplot",
      getPartialCorrelation = getPartialCorrelation,
      getGeneCorr           = getGeneCorr,
      cor_table             = cor_table,
      watermark             = WATERMARK,
      pgx                   = pgx,
      labeltype             = labeltype
    )

    correlation_plot_scattercorr_server(
      "cor_scatter",
      getFilteredExpression = getFilteredExpression,
      pgx = pgx,
      cor_table = cor_table,
      getPartialCorrelationMatrix = getPartialCorrelationMatrix,
      getGeneCorr = getGeneCorr,
      sel_gene = reactive(input$gene),
      COL = COL,
      watermark = WATERMARK,
      labeltype = labeltype
    )

    correlation_plot_correlation_UMAP_server(
      "cor_umap",
      pgx             = pgx,
      gene            = reactive(input$gene),
      getFullGeneCorr = getFullGeneCorr,
      getGeneCorr     = getGeneCorr,
      watermark       = WATERMARK
    )

    correlation_plot_cor_graph_server(
      "cor_graph",
      gene = reactive(input$gene),
      getPartialCorrelationMatrix = getPartialCorrelationMatrix,
      watermark = WATERMARK
    )
  })
} ## end of Board
