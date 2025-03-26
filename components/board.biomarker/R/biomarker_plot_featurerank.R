##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


biomarker_plot_featurerank_ui <- function(
    id,
    title,
    info.text,
    info.methods,
    info.references,
    info.extra_link,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  clust_featureRank.opts <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(ns("clust_featureRank_method"), "Method:",
        choices = c("p-value", "correlation", "meta"),
        inline = TRUE
      ),
      "Choose ranking method: p-value based or correlation-based.",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("pltmod"),
    label = label,
    plotlib = "plotly",
    title = title,
    caption = caption,
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    options = clust_featureRank.opts,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height
  )
}

biomarker_plot_featurerank_server <- function(id,
                                              pgx,
                                              ft_level,
                                              samplefilter,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    calcFeatureRanking <- shiny::reactive({
      pgx <- pgx
      ft_level <- ft_level()

      shiny::req(pgx$X, pgx$Y, pgx$gsetX, pgx$genes)      
      features <- NULL
      X <- NULL
      if (ft_level == "geneset") {
        gset_collections <- playbase::pgx.getGeneSetCollections(gsets = rownames(pgx$gsetX))
        features <- gset_collections
        X <- pgx$gsetX
        ## if(any(is.na(X))) X <- X[complete.cases(X), ]
      } else {
        features <- pgx$families
        X <- pgx$X
        if (any(is.na(X))) X <- pgx$impX
      }
      
      ## ------------ intersect features, set minimum set size
      genes <- toupper(pgx$genes$symbol)
      features <- lapply(features, toupper)
      features <- lapply(features, function(f) intersect(toupper(f), genes))
      features <- features[sapply(features, length) >= 10]

      ##------------- Supercell (scRNAseq very slow otherwise)
      ##------------- Or random downsampling by cell type?
      Y <- NULL
      c1 <- pgx$datatype %in% c("scRNAseq", "scRNA-seq")
      c2 <- ncol(pgx$counts) > 1000
      if (c1 & c2) {
        message("[biomarkers: feature-set scores] scRNAseq. >2K cells. Computing supercells.")
        group <- pgx$Y[, "celltype"]
        nb <- round(ncol(pgx$counts) / 600)
        message("[pgx.createSingleCellPGX]=======================================")
        message("[pgx.createSingleCellPGX] running SuperCell. nb = ", nb)    
        sc <- playbase::pgx.supercell(pgx$counts, pgx$Y, group = group, gamma = nb)
        message("[pgx.createSingleCellPGX] SuperCell done: ", ncol(counts), " -> ", ncol(sc$counts))
        message("[pgx.createSingleCellPGX]=======================================")
        message("[pgx.createSingleCellPGX] Normalizing supercell matrix (logCPM)")
        X <- playbase::logCPM(sc$counts, total = 1e4, prior = 1)
        X <- as.matrix(X)
        Y <- sc$meta
        remove(counts, group, sc); gc()
      }
      
      if (is.null(Y)) Y <- pgx$Y

      ##------------- rm redundant/unneeded (scRNAseq)
      if (all(c("nCount_RNA", "nCount_SCT") %in% colnames(Y))) {
        jj <- match("nCount_SCT", colnames(Y))
        Y <- Y[, -jj, drop = FALSE]
      }

      if (all(c("nFeature_RNA", "nFeature_SCT") %in% colnames(Y))) {
        jj <- match("nFeature_SCT", colnames(Y))
        Y <- Y[, -jj, drop = FALSE]
      }

      if (all(c(".cell_cycle", "Phase") %in% colnames(Y))) {
        jj <- match("Phase", colnames(Y))
        Y <- Y[, -jj, drop = FALSE]
      }

      if (all(c(".cell_cycle", "Phase") %in% colnames(Y))) {
        jj <- match("Phase", colnames(Y))
        Y <- Y[, -jj, drop = FALSE]
      }

      jj <- match("SCT_snn_res.1", colnames(Y), nomatch = FALSE)
      if (any(jj)) Y <- Y[, -jj, drop = FALSE]

      jj <- match("orig.ident", colnames(Y), nomatch = FALSE)
      if (any(jj)) Y <- Y[, -jj, drop = FALSE]
      
            
      ## ------------ Just to get current samples
      samples <- playbase::selectSamplesFromSelectedLevels(Y, samplefilter())
      X <- X[, samples, drop = FALSE]

      ## the code below overwrittes user input, and should be removed
      cvar <- playbase::pgx.getCategoricalPhenotypes(Y, max.ncat = 999)
      ss <- "sample|patient|years|days|months|gender"
      cvar <- grep(ss, cvar, invert = TRUE, value = TRUE) ## no sample IDs
      Y <- Y[colnames(X), cvar, drop = FALSE]
      kk <- which(apply(Y, 2, function(y) length(unique(y)) > 1))
      shiny::validate(shiny::need(length(kk) > 0, "Samples' filters too strict. Change 'Filter samples'."))
      Y <- Y[, kk, drop = FALSE]
      dim(Y)

      ## ------------ Note: this takes a while. Maybe better precompute off-line...
      sdx <- matrixStats::rowSds(X, na.rm = TRUE)
      names(sdx) <- rownames(X)
      S <- matrix(NA, nrow = length(features), ncol = ncol(Y))
      rownames(S) <- names(features)
      colnames(S) <- colnames(Y)

      ## ------------ Create a Progress object
      if (!interactive()) {
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        if(c1 & c2) {
          msg <- "Calculating feature-set scores on supercell data"
        } else {
          msg <- "Calculating feature-set scores"
        }
        progress$set(message = msg, value = 0)
      }

      gene.level <- TRUE
      gene.level <- (ft_level == "gene")

      if (gene.level) {
        annot <- pgx$genes
        annot <- annot[match(unique(annot$symbol), annot$symbol), ]
        rownames(annot) <- annot$symbol
      }
      
      for (i in 1:ncol(Y)) {
        if (!interactive()) progress$inc(1 / ncol(Y))

        grp <- as.character(Y[, i])

        score <- rep(NA, length(features))
        names(score) <- names(features)

        for (j in 1:length(features)) {

          pp <- features[[j]]

          if (gene.level)
            pp <- playbase::filterProbes(annot, features[[j]])

          if (c1 & c2) sdtop=300 else sdtop=1000
          pp <- head(pp[order(-sdx[pp])], sdtop)

          X1 <- playbase::rename_by(X, pgx$genes, "symbol")
          pp <- intersect(pp, rownames(X1))
          X1 <- X1[pp, , drop = FALSE]
          if (nrow(X1) == 0) { score[j]=NA; next }

          s1 <- s2 <- 1
          method <- input$clust_featureRank_method
          if (method %in% c("correlation", "meta")) {
            mx <- t(apply(X1, 1, function(x) tapply(x, grp, mean, na.rm = TRUE)))
            if (nrow(mx) == 0 || ncol(mx) == 0) next
            D <- 1 - cor(mx, use = "pairwise")
            diag(D) <- NA
            s1 <- mean(D, na.rm = TRUE)
          }

          if (method %in% c("p-value", "meta")) {
            jj <- which(!is.na(grp))
            design <- model.matrix(~ grp[jj])
            fit <- tryCatch(
              {
                suppressWarnings(limma::eBayes(limma::lmFit(X1[, jj, drop = FALSE], design)))
              },
              error = function(w) { NA }
            )
            if (all(is.na(fit))) { score[j]=NA; next }

            suppressWarnings(suppressMessages(top <- limma::topTable(fit)))
            s2 <- mean(-log10(1e-99 + top$adj.P.Val), na.rm = TRUE)
          }

          f <- 1
          f <- (1 - exp(-(length(pp) / 20)**2)) ## penalize smaller sets
          score[j] <- f * (s1 * s2) ** ifelse(method == "meta", 0.5, 1)
        }

        S[, i] <- score

      }

      S[is.na(S)] <- 0

      return(S)

    })

    render_featureRank <- function() {
      S <- calcFeatureRanking()
      c1 <- (is.null(S) || nrow(S) == 0 || ncol(S) == 0)
      if (c1) return(NULL)

      ## top scoring
      S <- tail(S[order(rowSums(S, na.rm = TRUE)), , drop = FALSE], 25)
      rownames(S) <- paste(substring(rownames(S), 1, 50), "  ")

      playbase::pgx.stackedBarplot(
        x = t(S),
        showlegend = TRUE,
        xlab = "Discriminant score",
        ylab = "",
        horiz = TRUE
      )
    }

    clust_featureRank.RENDER <- function() {
      render_featureRank() %>%
        plotly_default() %>%
        plotly::layout(legend = list(orientation = "h"))
    }

    clust_featureRank.RENDER2 <- function() {
      render_featureRank() %>% plotly_modal_default()
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = clust_featureRank.RENDER,
      func2 = clust_featureRank.RENDER2,
      csvFunc = calcFeatureRanking, ##  *** downloadable data as CSV
      res = c(72, 90), ## resolution of plots
      pdf.width = 8, pdf.height = 10,
      add.watermark = watermark
    )
  })
}
