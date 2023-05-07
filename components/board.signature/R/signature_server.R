##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

SignatureBoard <- function(id, pgx, selected_gxmethods) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    fullH <- 800 ## full height of page
    tabH <- "70vh"

    infotext <-
      "In the <strong>Signature Analysis module</strong>, users can test their gene signature by calculating an enrichment score. They can use a sample list provided on the platform or upload their own gene list. Instead of a short list, a profile can also be selected, which is a complete gene list resulted from one of the contrasts in the analysis.

<br><br>After uploading a gene list, the <strong>Markers</strong> section produces a t-SNE plot of samples for each gene, where the samples are colored with respect to the upregulation (in red) or downregulation (in blue) of that particular gene.

<br><br>The <strong>Enrichment tab</strong> performs the enrichment analysis of the gene list against all contrasts by running the GSEA algorithm and plots enrichment outputs. The enrichment statistics can be found in the corresponding table

<br><br>Under the <strong>Overlap/similarity tab</strong>, users can find the similarity of their gene list with all the gene sets and pathways in the platform, including statistics such as the total number of genes in the gene set (K), the number of intersecting genes between the list and the gene set (k), the overlapping ratio of k/K, as well as the p and q values by the Fisher’s test for the overlap test.

<br><br><br><br>
<center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=7' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>
"

    ## ================================================================================
    ## ========================= INPUTS UI ============================================
    ## ================================================================================

    IMMCHECK.GENES <- "ADORA2A ARHGEF5 BTLA CD160 CD244 CD27 CD274 CD276 CD47 CD80 CEACAM1 CTLA4 GEM HAVCR2 ICOS IDO1 LAG3 PDCD1 TNFSF4 VISTA VTCN1 TIGIT PVR CD28 CD40 CD40LG ICOSLG TNFRSF9 TNFSF9 CD70 TNFRSF4 TNFRSF18 TNFSF18 SIRPA LGALS9 ARG1 CD86 IDO2 PDCD1LG2 KIR2DL3"
    APOPTOSIS.GENES <- "BAD CRADD AGT FAS BCL2 PPIF S100A9 S100A8 BBC3 BCL2L11 FADD CTSH MLLT11 TRAF7 BCL2L1 HTRA2 BNIP3 BAK1 PMAIP1 LGALS9 BID"
    CELLCYCLE.GENES = "MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8 HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA"

    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Signature Analysis Board</strong>"),
        shiny::HTML(infotext),
        easyClose = TRUE, size = "xl"
      ))
    })

    ## ------------------------ observe/reactive function  -----------------------------

    shiny::observeEvent(input$example1, {
      shiny::updateTextAreaInput(session, "genelistUP", value = IMMCHECK.GENES)
    })
    shiny::observeEvent(input$example2, {
      shiny::updateTextAreaInput(session, "genelistUP", value = APOPTOSIS.GENES)
    })
    shiny::observeEvent(input$example3, {
      shiny::updateTextAreaInput(session, "genelistUP", value = CELLCYCLE.GENES)
    })

    shiny::observe({
      if (is.null(pgx)) {
        return(NULL)
      }
      type <- "contrast"
      type <- input$type
      if (is.null(type)) type <- "<custom>"

      if (type == "contrast") {
        contr <- sort(names(pgx$gx.meta$meta))
        shiny::updateSelectInput(session, "feature", choices = contr, selected = contr[1])
      } else if (type == "hallmark") {
        ## collection
        gsets <- sort(grep("HALLMARK", names(iGSETS), value = TRUE))
        shiny::updateSelectInput(session, "feature", choices = gsets, selected = gsets[1])
      } else if (type == "KEGG") {
        ## collection
        gsets <- sort(grep("KEGG", names(iGSETS), value = TRUE))
        shiny::updateSelectInput(session, "feature", choices = gsets, selected = gsets[1])
      } else if (type == "geneset") {
        ## all genesets... this is a bit too much for selectInput (DO NOT USE!!)
        gsets <- sort(names(iGSETS))
        shiny::updateSelectizeInput(session, "feature", choices = gsets, selected = gsets[1], server = TRUE)
      } else {
        ## custom
        shiny::updateSelectInput(session, "feature", choices = "<custom>", selected = "<custom>")
      }
    })


    ## ================================================================================
    ## ======================= REACTIVE FUNCTIONS =====================================
    ## ================================================================================

    input_genelistUP <- shiny::reactive({
      gg <- input$genelistUP
      if (is.null(gg)) {
        return(NULL)
      }
      gg <- strsplit(as.character(gg), split = "[, \n\t]")[[1]]
      if (length(gg) == 1 && gg[1] != "") gg <- c(gg, gg) ## hack to allow single gene....
      return(gg)
    }) %>% shiny::debounce(1000)


    getCurrentMarkers <- shiny::reactive({
      ##
      ## Get current selection of markers/genes
      ##
      ##

      if (is.null(pgx)) {
        return(NULL)
      }

      type <- input$type
      shiny::req(input$type, input$feature)

      dbg("<signature:getCurrentMarkers> called\n")

      level <- "gene"
      features <- toupper(pgx$genes$gene_name)
      xfeatures <- toupper(pgx$genes[rownames(pgx$X), "gene_name"])
      gset <- NULL
      if (input$feature == "<custom>") {
        gset <- input_genelistUP()
        if (is.null(gset) || length(gset) == 0 || gset[1] == "") {
          return(NULL)
        }
        if (length(gset) == 1) {
          gene <- sub("^[@#]", "", gset[1])
          if (grepl("^@", gset[1]) && gene %in% xfeatures) {
            ## most correlated with this genes
            jj <- match(gene, xfeatures) ## single gene
            rho <- cor(t(pgx$X), pgx$X[jj, ])[, 1]
            gset <- head(names(sort(abs(rho), decreasing = TRUE)), 36) ## how many?
          } else {
            ## grep-like match
            rx <- toupper(gset[1])
            rx <- grep(rx, xfeatures, value = TRUE, ignore.case = TRUE)
            gset <- rownames(pgx$X)[which(xfeatures %in% rx)] ## all probes matching gene
          }
        }
      } else if (type == "contrast" &&
        input$feature %in% names(pgx$gx.meta$meta)) {
        contr <- input$feature
        fx <- pgx$gx.meta$meta[[contr]]$meta.fx
        probes <- rownames(pgx$gx.meta$meta[[contr]])
        genes <- toupper(pgx$genes[probes, "gene_name"])
        top.genes <- genes[order(-fx)]
        top.genes <- head(top.genes, 100)
        top.genes0 <- paste(top.genes, collapse = " ")
        shiny::updateTextAreaInput(session, "genelistUP", value = top.genes0)
        gset <- top.genes
      } else if (input$feature %in% names(iGSETS)) {
        gset <- toupper(unlist(getGSETS(input$feature)))
        gset0 <- paste(gset, collapse = " ")
        shiny::updateTextAreaInput(session, "genelistUP", value = gset0)
      } else {
        return(NULL)
      }

      return(gset)
    })

    sigCalculateGSEA <- shiny::reactive({
      ##
      ## Calculate fgsea for current marker selection and active
      ## datasets.
      ##
      if (is.null(pgx)) {
        return(NULL)
      }


      ## observe input list
      gset <- head(rownames(pgx$X), 100)
      gset <- getCurrentMarkers()
      if (is.null(gset)) {
        return(NULL)
      }

      ## get all logFC of this dataset
      meta <- playbase::pgx.getMetaFoldChangeMatrix(pgx, what = "meta")
      F <- meta$fc
      rownames(F) <- toupper(rownames(F))

      ## cleanup matrix
      F <- as.matrix(F)
      F <- F[, which(!duplicated(colnames(F))), drop = FALSE]

      ## cleanup names and uppercase for mouse genes
      rownames(F) <- toupper(sub(".*:", "", rownames(F)))
      gset <- toupper(sub(".*:", "", gset))
      gset <- intersect(toupper(gset), rownames(F))

      if (length(gset) == 0) {
        cat("FATAL:: sigCalculateGSEA : gset empty!\n")
        return(NULL)
      }

      ## ------------ prioritize with quick correlation

      y <- 1 * (toupper(rownames(F)) %in% toupper(gset))
      ss.rank <- function(x) scale(sign(x) * rank(abs(x)), center = FALSE)[, 1]
      rho <- cor(apply(F, 2, ss.rank), y, use = "pairwise")[, 1]
      rho[is.na(rho)] <- 0
      names(rho) <- colnames(F)

      ## ------- restrict to top 100 comparisons (fgsea is otherwise to
      ## ------- slow) but we do not expect many with so many
      ## ------- comparisons
      ntop <- 100
      jj <- head(order(-abs(rho)), ntop)
      F <- F[, jj, drop = FALSE]
      F <- F[!duplicated(rownames(F)), , drop = FALSE]
      F <- F + 1e-4 * matrix(rnorm(length(F)), nrow(F), ncol(F))

      ## ------------- do fast GSEA
      gmt <- list("gset" = unique(gset))
      res <- NULL
      enrich_method <- "fgsea"
      ## enrich_method <- input$rankmethod

      if (enrich_method == "fgsea") {
        i <- 1
        shiny::withProgress(message = "Computing GSEA ...", value = 0.8, {
          res <- lapply(1:ncol(F), function(i) {
            suppressWarnings(suppressMessages(
              res <- fgsea::fgsea(gmt, stats = F[, i], nperm = 1000)
            ))
            res <- as.data.frame(res[, c("pval", "padj", "ES", "NES")])
            rownames(res)[1] <- colnames(F)[i]
            return(res)
          })
        })
        res1 <- data.frame(do.call(rbind, res))
        res1$ES <- NULL
      } else {
        i <- 1
        fx <- 1 * (rownames(F) %in% gmt[[1]])
        rho <- cor(apply(F, 2, rank, na.last = "keep"), fx, use = "pairwise")[, 1]
        pv <- playbase::cor.pvalue(rho, nrow(F))
        qv <- p.adjust(pv, method = "fdr")
        res1 <- data.frame(pval = pv, padj = qv, rho = rho, NES = NA)
        rownames(res1) <- names(pv)
      }

      ## columns are: NES, pval, fdr, contrast
      res1 <- as.matrix(res1)
      res1 <- res1[match(colnames(F), rownames(res1)), , drop = FALSE]

      if (nrow(res1) != ncol(F)) {
        cat("WARNING sigCalculateGSEA:: fgsea results are corrupted?\n")
        cat("WARNING sigCalculateGSEA:: got contrasts: ", res$contrast, "\n")
        cat("WARNING sigCalculateGSEA:: colnames.F= ", colnames(F), "\n")
      }

      ## make nice table
      nes <- res1[, "NES"]
      pval <- res1[, "pval"]
      qval <- p.adjust(pval, method = "fdr")
      rho <- rho[colnames(F)]

      output <- as.matrix(cbind(NES = nes, p = pval, q = qval, rho = rho))
      rownames(output) <- colnames(F)
      output <- output[order(-abs(output[, "NES"])), , drop = FALSE]
      F <- F[, rownames(output), drop = FALSE]
      gsea <- list(F = as.matrix(F), gset = gset, output = output)
      return(gsea)
    })

    ## ================================================================================
    ## Overlap/similarity
    ## ================================================================================

    getOverlapTable <- shiny::reactive({
      ##
      ##
      ##
      if (is.null(pgx)) {
        return(NULL)
      }

      markers <- head(rownames(pgx$X), 100)
      markers <- getCurrentMarkers()
      if (is.null(markers)) {
        return(NULL)
      }

      ## fold change just for ranking of genes
      ## F <- sapply(pgx$gx.meta$meta, function(x) unclass(x$fc)[,"trend.limma"])
      F <- sapply(pgx$gx.meta$meta, function(x) x$meta.fx)
      rownames(F) <- rownames(pgx$gx.meta$meta[[1]])
      fx <- rowMeans(F**2)

      ## fisher test
      ii <- setdiff(match(toupper(markers), colnames(GSETxGENE)), NA)
      N <- cbind(
        k1 = Matrix::rowSums(GSETxGENE != 0), n1 = ncol(GSETxGENE),
        k2 = Matrix::rowSums(GSETxGENE[, ii] != 0), n2 = length(ii)
      )
      rownames(N) <- rownames(GSETxGENE)
      N <- N[which(N[, 1] > 0 | N[, 3] > 0), ]
      odds.ratio <- (N[, 3] / N[, 4]) / (N[, 1] / N[, 2])
      dim(N)

      ## WOW THIS IS FAST!!!!!!!
      pv <- corpora::fisher.pval(N[, 1], N[, 2], N[, 3], N[, 4], log.p = FALSE)
      head(pv)
      names(pv) <- rownames(N)
      pv <- pv[match(names(odds.ratio), names(pv))]
      qv <- p.adjust(pv, method = "bonferroni")
      A <- data.frame(odds.ratio = odds.ratio, p.fisher = pv, q.fisher = qv)
      dim(A)

      ## limit the list??
      A <- A[which(A$q.fisher < 0.999), ]

      ## get shared genes
      aa <- rownames(A)

      y <- 1 * (colnames(GSETxGENE) %in% toupper(markers))
      names(y) <- colnames(GSETxGENE)
      ncommon <- Matrix::colSums(Matrix::t(GSETxGENE[aa, , drop = FALSE]) * as.vector(y) != 0)
      ntotal <- Matrix::rowSums(GSETxGENE[aa, , drop = FALSE] != 0)
      A$ratio <- ncommon / ntotal
      ratio.kk <- paste0(ncommon, "/", ntotal)

      gg <- colnames(GSETxGENE)
      gset <- names(y)[which(y != 0)]
      G1 <- GSETxGENE[aa, which(y != 0)]
      commongenes <- apply(G1, 1, function(x) colnames(G1)[which(x != 0)])
      for (i in 1:length(commongenes)) {
        gg <- commongenes[[i]]
        gg <- gg[order(-abs(fx[gg]))]
        if (length(gg) > 10) {
          others <- paste0("(+", length(gg) - 10, " others)")
          gg <- c(head(gg, 10), others)
        }
        commongenes[[i]] <- paste(gg, collapse = ",")
      }
      commongenes <- unlist(commongenes)

      ## construct results dataframe
      gset.names <- substring(rownames(A), 1, 72)
      A$ratio <- round(A$ratio, digits = 3)
      A$log.OR <- round(log10(A$odds.ratio), digits = 3)
      A$odds.ratio <- round(A$odds.ratio, digits = 3)
      db <- sub(":.*", "", gset.names)
      score <- (log10(A$odds.ratio) * -log10(A$q.fisher + 1e-40))**0.5
      score <- round(score, digits = 3)
      df <- cbind(db = db, geneset = gset.names, score = score, "k/K" = ratio.kk, A, common.genes = commongenes)

      if (DEV) {
        df <- df[, c("db", "geneset", "score", "k/K", "ratio", "odds.ratio", "log.OR", "q.fisher", "common.genes")]
      } else {
        df <- df[, c("db", "geneset", "score", "k/K", "odds.ratio", "q.fisher", "common.genes")]
      }

      df <- df[order(-df$score), ]
      return(df)
    })

    ## ================================================================================
    ## Enrichment {data-height=800}
    ## ================================================================================

    getEnrichmentGeneTable <- shiny::reactive({
      shiny::req(pgx)

      gsea <- sigCalculateGSEA()
      if (is.null(gsea)) {
        return(NULL)
      }

      i <- enrichmentContrastTable$rows_selected()
      if (is.null(i) || length(i) == 0) {
        return(NULL)
      }

      meta <- playbase::pgx.getMetaFoldChangeMatrix(pgx, what = "meta")
      fc <- meta$fc
      qv <- meta$qv
      rownames(fc) <- toupper(rownames(fc))
      rownames(qv) <- toupper(rownames(qv))

      contr <- rownames(gsea$output)[i]
      fc <- fc[, contr, drop = FALSE]

      gset <- getCurrentMarkers()
      if (is.null(gset)) {
        return(NULL)
      }

      gset <- setdiff(toupper(gset), c("", NA))
      genes <- intersect(gset, rownames(fc))
      dd1 <- setdiff(genes, rownames(fc))
      dd2 <- setdiff(genes, rownames(qv))
      fc <- fc[genes, , drop = FALSE]
      qv <- qv[genes, , drop = FALSE]

      gene.tt <- substring(GENE.TITLE[toupper(rownames(fc))], 1, 40)
      names(gene.tt) <- rownames(fc)
      df <- data.frame(gene = rownames(fc), title = gene.tt, fc, check.names = FALSE)

      df
    })

    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    # Enrichment plots

    signature_plot_enplots_server(
      "enplots",
      pgx = pgx,
      sigCalculateGSEA = sigCalculateGSEA,
      enrichmentContrastTable = enrichmentContrastTable,
      watermark = WATERMARK
    )

    # Volcano plots

    signature_plot_volcano_server(
      "volcanoPlots",
      pgx = pgx,
      sigCalculateGSEA = sigCalculateGSEA,
      enrichmentContrastTable = enrichmentContrastTable,
      selected_gxmethods = selected_gxmethods,
      enrichmentGeneTable = enrichmentGeneTable,
      getEnrichmentGeneTable = getEnrichmentGeneTable,
      watermark = WATERMARK
    )

    # Signature overlap scores

    signature_plot_overlap_server(
      "overlapScorePlot",
      getOverlapTable = getOverlapTable,
      overlapTable = overlapTable,
      watermark = WATERMARK
    )

    # Overlap with other signatures

    overlapTable <- signature_table_overlap_server(
      "overlapTable",
      getOverlapTable = getOverlapTable,
      fullH = fullH,
      tabH = tabH
    )

    # Markers plot

    signature_plot_markers_server(
      "markers",
      pgx = pgx,
      getCurrentMarkers = getCurrentMarkers,
      watermark = WATERMARK
    )

    # Enrichment by contrasts

    enrichmentContrastTable <- signature_table_enrich_by_contrasts_server(
      "enrichmentContrastTable",
      sigCalculateGSEA = sigCalculateGSEA
    )

    # Genes in signature

    enrichmentGeneTable <- signature_table_genes_in_signature_server(
      "enrichmentGeneTable",
      getEnrichmentGeneTable = getEnrichmentGeneTable
    )
  })
} ## end-of-Board
