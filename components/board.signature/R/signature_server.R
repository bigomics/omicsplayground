##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

SignatureBoard <- function(id, pgx,
                           selected_gxmethods = reactive(colnames(pgx$gx.meta$meta[[1]]$fc))) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    fullH <- 800 ## full height of page
    tabH <- "70vh"

    infotext <-
      tspan("In the <strong>Signature Analysis module</strong>, users can test their gene signature by calculating an enrichment score. They can use a sample list provided on the platform or upload their own gene list. Instead of a short list, a profile can also be selected, which is a complete gene list resulted from one of the contrasts in the analysis.

<br><br>After uploading a gene list, the <strong>Markers</strong> section produces a t-SNE plot of samples for each gene, where the samples are colored with respect to the upregulation (in red) or downregulation (in blue) of that particular gene.

<br><br>The <strong>Enrichment tab</strong> performs the enrichment analysis of the gene list against all contrasts by running the GSEA algorithm and plots enrichment outputs. The enrichment statistics can be found in the corresponding table

<br><br>Under the <strong>Overlap/similarity tab</strong>, users can find the similarity of their gene list with all the gene sets and pathways in the platform, including statistics such as the total number of genes in the gene set (K), the number of intersecting genes between the list and the gene set (k), the overlapping ratio of k/K, as well as the p and q values by the Fisher’s test for the overlap test.

<br><br><br><br>
<center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=7' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>
", js = FALSE)

    ## ================================================================================
    ## ========================= INPUTS UI ============================================
    ## ================================================================================

    DEFAULT.GENES <- "MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8 HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA"

    IMMCHECK.GENES <- "ADORA2A ARHGEF5 BTLA CD160 CD244 CD27 CD274 CD276 CD47 CD80 CEACAM1 CTLA4 GEM HAVCR2 ICOS IDO1 LAG3 PDCD1 TNFSF4 VISTA VTCN1 TIGIT PVR CD28 CD40 CD40LG ICOSLG TNFRSF9 TNFSF9 CD70 TNFRSF4 TNFRSF18 TNFSF18 SIRPA LGALS9 ARG1 CD86 IDO2 PDCD1LG2 KIR2DL3"
    APOPTOSIS.GENES <- "BAD CRADD AGT FAS BCL2 PPIF S100A9 S100A8 BBC3 BCL2L11 FADD CTSH MLLT11 TRAF7 BCL2L1 HTRA2 BNIP3 BAK1 PMAIP1 LGALS9 BID"
    CELLCYCLE.GENES <- "MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8 HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA"

    DEFAULT.METABOLITES <- "16827 80470 21264 53487 29009 16843 28712 81068 2700 17111 30915 17361 16411 17858 44247 46229 17752 27849 1188  17381 16708 17052 27814 17566 27398 52222 30851 16856 17406 24813 27438 15971 16504 16082 28822 17712 53488 17442 15811 16755 65046 64557 28158"
    GLYCOLISIS.METABOLITES <- "17234 15361 17665 15946 16905 16108 29052 16001 17794 17835 18021 15422 16761 15846 24996 28602"
    CITRICACIDCYCLE.METABOLITES <- "57288 16947 32805 30887 16810 15380 15741 18012 30797 30744 15346 15846 16908 16238 17877 17552 15996 16526 15377"
    UREACYCLE.METABOLITES <- "15729 16349 16941 16467 16199 17672 17053 18012"

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
      if (DATATYPEPGX == "metabolomics") {
        shiny::updateTextAreaInput(session, "genelist", value = GLYCOLISIS.METABOLITES)
      } else {
        shiny::updateTextAreaInput(session, "genelist", value = IMMCHECK.GENES)
      }
    })
    shiny::observeEvent(input$example2, {
      if (DATATYPEPGX == "metabolomics") {
        shiny::updateTextAreaInput(session, "genelist", value = CITRICACIDCYCLE.METABOLITES)
      } else {
        shiny::updateTextAreaInput(session, "genelist", value = APOPTOSIS.GENES)
      }
    })
    shiny::observeEvent(input$example3, {
      if (DATATYPEPGX == "metabolomics") {
        shiny::updateTextAreaInput(session, "genelist", value = UREACYCLE.METABOLITES)
      } else {
        shiny::updateTextAreaInput(session, "genelist", value = CELLCYCLE.GENES)
      }
    })
    shiny::observeEvent(pgx$X, {
      if (DATATYPEPGX == "metabolomics") {
        shiny::updateTextAreaInput(session, "genelist", value = DEFAULT.METABOLITES)
        shiny::updateActionButton(session, "example1", label = "[glycolisis] ")
        shiny::updateActionButton(session, "example2", label = "[TCA cycle] ")
        shiny::updateActionButton(session, "example3", label = "[urea cycle] ")
      } else {
        shiny::updateTextAreaInput(session, "genelist", value = DEFAULT.GENES)
        shiny::updateActionButton(session, "example1", label = "[immune_chkpt] ")
        shiny::updateActionButton(session, "example2", label = "[apoptosis] ")
        shiny::updateActionButton(session, "example3", label = "[cell_cycle] ")
      }
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
        gsets <- sort(grep("HALLMARK", names(playdata::iGSETS), value = TRUE))
        names(gsets) <- tolower(gsub("_", " ", sub(".*HALLMARK_", "", gsets)))
        shiny::updateSelectInput(session, "feature", choices = gsets, selected = gsets[1])
      } else if (type == "KEGG") {
        ## collection
        gsets <- sort(grep("KEGG", names(playdata::iGSETS), value = TRUE))
        names(gsets) <- tolower(gsub("_", " ", sub(".*KEGG_", "", gsets)))
        shiny::updateSelectInput(session, "feature", choices = gsets, selected = gsets[1])
      } else if (type == "geneset") {
        ## all genesets... this is a bit too much for selectInput (DO NOT USE!!)
        gsets <- sort(names(playdata::iGSETS))
        shiny::updateSelectizeInput(session, "feature",
          choices = gsets,
          selected = gsets[1], server = TRUE
        )
      } else {
        ## custom
        shiny::updateSelectInput(session, "feature",
          choices = "<custom>",
          selected = "<custom>"
        )
      }
    })


    ## ================================================================================
    ## ======================= REACTIVE FUNCTIONS =====================================
    ## ================================================================================

    input_genelist <- shiny::eventReactive(
      {
        list(pgx$X, input$compute_button)
      },
      {
        shiny::req(input$genelist)
        gg <- input$genelist
        gg <- strsplit(as.character(gg), split = "[, \n\t]")[[1]]
        trimws(gg)
      },
      ignoreInit = FALSE,
      ignoreNULL = TRUE
    )

    getCurrentMarkers <- shiny::reactive({
      shiny::req(pgx)
      shiny::req(input$type, input$feature, input_genelist())

      ## Get current selection of markers/genes
      type <- input$type
      symbols <- toupper(pgx$genes[rownames(pgx$X), "symbol"])

      features <- NULL
      if (input$feature == "<custom>") {
        genes <- input_genelist()
        if (is.null(genes) || length(genes) == 0 || genes[1] == "") {
          return(NULL)
        }

        if (any(grepl("[*]", genes))) {
          gene.patterns <- grep("[*]", genes, value = TRUE)

          rx.genes <- c()
          for (rx in gene.patterns) {
            ## select other genes with grep-like match
            rx <- paste0("^", sub("[*]", ".*", rx))
            gg <- unique(grep(rx, symbols, value = TRUE, ignore.case = TRUE))
            rx.genes <- c(rx.genes, gg)
          }
          genes <- union(genes, rx.genes)
        }
        genes <- intersect(toupper(genes), symbols)
        ## map to probes
        features1 <- playbase::map_probes(pgx$genes, genes,
          column = "human_ortholog", ignore.case = TRUE
        )
        features2 <- playbase::map_probes(pgx$genes, genes,
          column = "symbol", ignore.case = TRUE
        )
        features <- union(features1, features2)
      } else if (input$type == "contrast" &&
        input$feature[1] %in% playbase::pgx.getContrasts(pgx)) {
        contr <- input$feature
        fx <- pgx$gx.meta$meta[[contr]]$meta.fx
        names(fx) <- rownames(pgx$gx.meta$meta[[contr]])
        ## take top 100 features
        top.genes <- fx[order(-abs(fx))]
        top.genes <- head(top.genes, min(100, length(fx) / 4))
        features <- names(top.genes)
      } else if (input$type %in% c("hallmark", "KEGG", "geneset") &&
        input$feature[1] %in% colnames(pgx$GMT)) {
        sel <- input$feature
        features <- rownames(pgx$GMT)[which(pgx$GMT[, sel] != 0)]
      } else {
        return(NULL)
      }

      ## convert to genes symbols
      symbols <- pgx$genes[features, "symbol"]

      list(
        features = unique(features),
        symbols = unique(symbols)
      )
    })

    sigCalculateGSEA <- shiny::reactive({
      ## Calculate fgsea for current marker selection and active
      ## datasets.
      if (is.null(pgx)) {
        return(NULL)
      }
      ## observe input list
      markers <- getCurrentMarkers()
      genes <- markers$symbols

      if (is.null(genes) || length(genes) == 0) {
        cat("FATAL:: sigCalculateGSEA : gset empty!\n")
        return(NULL)
      }

      ## get all logFC of this dataset
      meta <- playbase::pgx.getMetaMatrix(pgx)
      F <- meta$fc

      ## cleanup matrix
      F <- as.matrix(F)
      F <- F[, which(!duplicated(colnames(F))), drop = FALSE]

      ## convert to symbol
      fsymbol <- pgx$genes[rownames(F), "symbol"]
      F <- playbase::rowmean(F, fsymbol)
      F[is.na(F)] <- 0

      ## prioritize with quick correlation restrict to top 100
      ## comparisons (fgsea is otherwise to slow) but we do not expect
      ## many with so many comparisons
      y <- 1 * (rownames(F) %in% genes)
      ss.rank <- function(x) scale(sign(x) * rank(abs(x)), center = FALSE)[, 1]
      rho <- cor(apply(F, 2, ss.rank), y, use = "pairwise")[, 1]
      rho[is.na(rho)] <- 0
      names(rho) <- colnames(F)
      ntop <- 100
      jj <- head(order(-abs(rho)), ntop)
      F <- F[, jj, drop = FALSE]
      F <- F + 1e-4 * matrix(rnorm(length(F)), nrow(F), ncol(F))

      ## ------------- do fast GSEA
      gmt <- list("gset" = unique(genes))
      res <- NULL
      enrich_method <- "fgsea"

      if (enrich_method == "fgsea") {
        i <- 1
        shiny::withProgress(message = "Computing GSEA ...", value = 0.8, {
          res <- lapply(1:ncol(F), function(i) {
            set.seed(123)
            suppressWarnings(suppressMessages(
              res <- fgsea::fgsea(gmt, stats = F[, i], nperm = 1000, nproc = 1)
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
        message("WARNING sigCalculateGSEA:: fgsea results are corrupted?\n")
        message("WARNING sigCalculateGSEA:: got contrasts: ", res$contrast, "\n")
        message("WARNING sigCalculateGSEA:: colnames.F= ", colnames(F), "\n")
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

      gsea <- list(F = as.matrix(F), gset = genes, output = output)
      return(gsea)
    })

    ## ================================================================================
    ## Compute overlap and statistics of given list of genes and known genesets
    ## ================================================================================

    getOverlapTable <- shiny::reactive({
      shiny::req(pgx)
      shiny::req(getCurrentMarkers())

      markers <- getCurrentMarkers()
      genes <- markers$symbols

      ## fisher test
      G <- Matrix::t(pgx$GMT)
      ii <- setdiff(match(genes, colnames(G)), NA)

      N <- cbind(
        k1 = Matrix::rowSums(G != 0),
        n1 = ncol(G),
        k2 = Matrix::rowSums(G[, ii, drop = FALSE] != 0),
        n2 = length(ii)
      )
      rownames(N) <- rownames(G)
      N <- N[which(N[, 1] > 0 | N[, 3] > 0), ]
      odds.ratio <- (N[, 3] / N[, 4]) / (N[, 1] / N[, 2])

      ## WOW THIS IS FAST!!!!!!!
      ## pv <- corpora::fisher.pval(N[, 1], N[, 2], N[, 3], N[, 4], log.p = FALSE)
      pv <- try(
        corpora::fisher.pval(N[, 1], N[, 2], N[, 3], N[, 4], log.p = FALSE),
        silent = TRUE
      )
      if (class(pv) != "try-error") {
        names(pv) <- rownames(N)
        pv <- pv[match(names(odds.ratio), names(pv))]
        qv <- p.adjust(pv, method = "bonferroni")
      } else {
        message("[signature_server] corpora::fisher.pval failed. Using standard fisher.")
        pv <- c(rep(NA, nrow(N)))
        i <- 1
        for (i in 1:nrow(N)) {
          tt <- matrix(NA, nrow = 2, ncol = 2)
          tt[1, 1] <- N[i, 1]
          tt[1, 2] <- N[i, 2] - N[i, 1]
          tt[2, 1] <- N[i, 3]
          tt[2, 2] <- N[i, 4] - N[i, 3]
          ## pv[i] <- stats::fisher.test(tt, alternative = "greater")$p.value ## ??
          pv[i] <- stats::fisher.test(tt)$p.value
        }
        names(pv) <- rownames(N)
        pv <- pv[match(names(odds.ratio), names(pv))]
        qv <- p.adjust(pv, method = "bonferroni")
      }

      A <- data.frame(odds.ratio = odds.ratio, p.fisher = pv, q.fisher = qv)

      ## get shared genes
      aa <- rownames(A)
      y <- 1 * (colnames(G) %in% genes)
      names(y) <- colnames(G)
      ncommon <- Matrix::colSums(Matrix::t(G[aa, , drop = FALSE]) * as.vector(y) != 0)
      ntotal <- Matrix::rowSums(G[aa, , drop = FALSE] != 0)
      A$ratio <- ncommon / ntotal
      ratio.kk <- paste0(ncommon, "/", ntotal)

      ## determine top genes
      meta <- playbase::pgx.getMetaMatrix(pgx)
      fx <- sqrt(rowMeans(meta$fc**2, na.rm = TRUE))
      gg <- colnames(G)
      fx <- fx[match(gg, names(fx))]
      names(fx) <- gg

      gset <- names(y)[which(y != 0)]
      G1 <- G[aa, which(y != 0), drop = FALSE]

      commongenes <- apply(G1, 1, function(x) colnames(G1)[which(x != 0)])

      shiny::validate(shiny::need(
        length(commongenes) > 0,
        "No results. Perhaps you need to try on a bigger datasets with more features."
      ))

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
      df <- cbind(
        db = db, geneset = gset.names, score = score, "k/K" = ratio.kk, A,
        common.genes = commongenes
      )

      df <- df[, c("db", "geneset", "score", "k/K", "odds.ratio", "q.fisher", "common.genes")]
      df <- df[order(-df$score), ]
      return(df)
    })

    ## ================================================================================
    ## Enrichment {data-height=800}
    ## ================================================================================

    getEnrichmentGeneTable <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(sigCalculateGSEA())
      shiny::req(getCurrentMarkers())

      # Input vars
      gsea <- sigCalculateGSEA()
      markers <- getCurrentMarkers()
      gset <- markers$symbols
      features <- markers$features

      i <- enrichmentContrastTable$rows_selected()
      if (is.null(i) || length(i) == 0) {
        return(NULL)
      }
      contr <- rownames(gsea$output)[i]

      # Get metaFC mat
      meta <- playbase::pgx.getMetaMatrix(pgx)
      fc <- meta$fc
      qv <- meta$qv

      # Get contrasts, FC, and gene info
      fc <- fc[features, contr]
      qv <- qv[features, contr]
      gene.info <- pgx$genes[features, c("feature", "symbol")]

      # Return df
      df <- data.frame(gene.info, log2FC = fc, q.value = qv, check.names = FALSE)
      df <- df[order(-df$log2FC), ]

      return(df)
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
      organism = pgx$organism,
      getEnrichmentGeneTable = getEnrichmentGeneTable
    )
  })
} ## end-of-Board
