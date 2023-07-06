##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

CorrelationBoard <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    fullH <- 800 ## full height of page
    rowH <- 340 ## full height of page

    cor_infotext <- "The <strong>Correlation Analysis Board</strong> provides statistical correlation analysis on gene level with visualisations. During the visual analysis, users can filter out some samples or collapse the samples by predetermined groups. The dark shaded area in the barplot estimates the partial correlation."

    ## ------- observe functions -----------
    shiny::observeEvent(input$data_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Expression: Correlation analysis board</strong>"),
        shiny::HTML(cor_infotext),
        easyClose = TRUE, size = "l"
      ))
    })

    COL <- RColorBrewer::brewer.pal(12, "Paired")[seq(1, 12, 2)]
    COL <- RColorBrewer::brewer.pal(9, "Set1")[c(2, 1, 3:9)]
    COL2 <- rev(grey.colors(2))
    COL2 <- RColorBrewer::brewer.pal(2, "Paired")[1:2]
    COL2 <- COL[1:2]

    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    shiny::observeEvent(input$cor_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Correlation Analysis Board</strong>"),
        shiny::HTML(cor_infotext),
        easyClose = TRUE
      ))
    })

    ## update filter choices upon change of data set
    shiny::observe({
      req(pgx$X)

      genes <- sort(pgx$genes[rownames(pgx$X), ]$gene_name)
      sel <- genes[1] ## most var gene
      sel <- names(head(sort(-rowMeans(playbase::pgx.getMetaMatrix(pgx)$fc**2)), 1))
      shiny::updateSelectizeInput(session, "cor_gene", choices = genes, selected = sel, server = TRUE)

      fam <- playbase::pgx.getFamilies(pgx, nmin = 10, extended = FALSE)
      fam <- sort(c("<custom>", fam))
      names(fam) <- sub(".*:","",fam)
      shiny::updateSelectInput(session, "cor_features", choices = fam)

      px <- colnames(pgx$Y)
      s1 <- grep("^[.]", px, value = TRUE, invert = TRUE)[1]
      shiny::updateSelectInput(session, "cor_group", choices = px, selected = s1)
    })


    getFilteredExpression <- shiny::reactive({
      shiny::req(pgx, input$cor_gene)
      X <- pgx$X
      gene <- rownames(X)[1]
      gene <- input$cor_gene

      ## filter genes
      ft <- input$cor_features
      if (ft == "<custom>" && input$cor_customfeatures != "") {
        genes <- toupper(pgx$genes$gene_name)
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
        gg1 <- gg1[which(toupper(gg1) %in% toupper(pgx$genes$gene_name))]
        psel <- playbase::filterProbes(pgx$genes, c(gg1, gene))
        psel <- intersect(psel, rownames(X))
        X <- X[psel, , drop = FALSE]
      } else if (ft != "<all>" && ft %in% names(playdata::iGSETS)) {
        ft <- input$cor_features
        psel <- playbase::filterProbes(pgx$genes, c(gene, unlist(playdata::getGSETS(ft))))
        #
        psel <- intersect(psel, rownames(X))
        X <- X[psel, , drop = FALSE]
      }

      X <- X + 1e-3 * matrix(rnorm(length(X)), nrow(X), ncol(X))
      X
    })

    getPartialCorrelationMatrix <- shiny::reactive({
      shiny::req(pgx, input$cor_gene)

      gene <- rownames(pgx$X)[1]
      gene <- input$cor_gene

      ## filter gene expression matrix
      X <- getFilteredExpression()

      
      NTOP <- 50
      NTOP <- as.integer(input$pcor_ntop)
      ## res <- playbase::pgx.computePartialCorrelationAroundGene(
      ##    X, gene, method=methods, nmax=NTOP, fast=FALSE)
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
      gene <- rownames(res$cor)[1]
      gene <- input$cor_gene
      rho <- res$cor[gene, ]
      prho <- res$pcor[gene, ]
      df <- data.frame(cor = rho, pcor = prho)
      df
    })

    getGeneCorr <- shiny::reactive({
      shiny::req(pgx)
      gene <- input$cor_gene
      if (is.null(gene)) {
        return(NULL)
      }
      
      ## corr always in log.scale and restricted to selected samples subset
      zx <- pgx$X
      zx <- getFilteredExpression()
      dim(zx)
      
      zx.genes0 <- rownames(zx)
      zx.genes <- as.character(pgx$genes[rownames(zx), ]$gene_name)
      rownames(zx) <- toupper(zx.genes)
      
      xref <- list(
        "cor" = 2**zx,
        "cor.HPA" = as.matrix(playdata::TISSUE),
        "cor.ImmProt" = as.matrix(playdata::IMMPROT)
      )
      gene0 <- toupper(gene) ## uppercase mouse
      R <- playbase::pgx.getGeneCorrelation(gene0, xref = xref)
      if (is.null(R)) {
        return(NULL)
      }
      R <- R[rownames(zx), , drop = FALSE]
      
      zx <- zx - rowMeans(zx, na.rm = TRUE)
      sdx <- sqrt(rowMeans(zx**2))      
      R <- cbind(R, cov = R[, "cor"] * sdx * sdx[gene0])
      
      rho.genes <- rownames(zx)
      if ("hgnc_symbol" %in% colnames(pgx$genes)) {
        rho.genes <- as.character(pgx$genes[zx.genes0, ]$hgnc_symbol)
      }
      
      R <- R[match(rho.genes, rownames(R)), , drop = FALSE]
      rownames(R) <- zx.genes0
      R <- R[order(R[, "cor"], decreasing = TRUE), , drop = FALSE]
      
      R
    })

    getFullGeneCorr <- shiny::reactive({
      shiny::req(pgx)

      gene <- input$cor_gene
      if (is.null(gene)) {
        return(NULL)
      }

      ## corr always in log.scale and restricted to selected samples subset
      zx <- pgx$X
      zx.genes0 <- rownames(zx)
      zx.genes <- as.character(pgx$genes[rownames(zx), ]$gene_name)
      rownames(zx) <- toupper(zx.genes)
      xref <- list(
        "cor" = 2**zx,
        "cor.HPA" = as.matrix(playdata::TISSUE),
        "cor.ImmProt" = as.matrix(playdata::IMMPROT)
      )
      gene0 <- toupper(gene) ## uppercase mouse
      R <- playbase::pgx.getGeneCorrelation(gene0, xref = xref)
      if (is.null(R)) {
        return(NULL)
      }
      R <- R[rownames(zx), , drop = FALSE]

      zx <- zx - rowMeans(zx, na.rm = TRUE)
      sdx <- sqrt(rowMeans(zx**2))
      R <- cbind(R, cov = R[, "cor"] * sdx * sdx[gene0])

      rho.genes <- rownames(zx)
      if ("hgnc_symbol" %in% colnames(pgx$genes)) {
        rho.genes <- as.character(pgx$genes[zx.genes0, ]$hgnc_symbol)
      }
      R <- R[match(rho.genes, rownames(R)), , drop = FALSE]
      rownames(R) <- zx.genes0
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
      watermark             = WATERMARK
    )

    correlation_plot_scattercorr_server(
      "cor_scatter",
      getFilteredExpression = getFilteredExpression,
      pgx = pgx,
      cor_table = cor_table,
      getPartialCorrelationMatrix = getPartialCorrelationMatrix,
      getGeneCorr = getGeneCorr,
      cor_gene = reactive(input$cor_gene),
      COL = COL,
      watermark = WATERMARK
    )

    correlation_plot_correlation_UMAP_server(
      "cor_umap",
      pgx             = pgx,
      cor_gene        = reactive(input$cor_gene),
      getFullGeneCorr = getFullGeneCorr,
      getGeneCorr     = getGeneCorr,
      watermark       = WATERMARK
    )

    correlation_plot_cor_graph_server(
      "cor_graph",
      cor_gene = reactive(input$cor_gene),
      getPartialCorrelationMatrix = getPartialCorrelationMatrix,
      watermark = WATERMARK
    )

    
  })
} ## end of Board
