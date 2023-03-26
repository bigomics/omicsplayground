##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


clustering_plot_featurerank_ui <- function(id,
                                           label = "",
                                           height,
                                           width) {
  ns <- shiny::NS(id)

  clust_featureRank_info <- "Ranked discriminant score for top feature sets. The plot ranks the discriminitive power of the feature set (genes) as a cumulative discriminant score for all phenotype variables. In this way, we can find which feature set (or gene family/set) can explain the variance in the data the best. <p>Correlation-based discriminative power is calculated as the average '(1-cor)' between the groups. Thus, a feature set is highly discriminative if the between-group correlation is low. P-value based scoring is computed as the average negative log p-value from the ANOVA. The 'meta' method combines the score of the former methods in a multiplicative manner."


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
    title = "Feature-set ranking",
    info.text = clust_featureRank_info,
    options = clust_featureRank.opts,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
    height = height
  )
}

clustering_plot_featurerank_server <- function(id,
                                               pgx,
                                               hm_level,
                                               hm_samplefilter,
                                               watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    calcFeatureRanking <- shiny::reactive({
      pgx <- pgx
      hm_level <- hm_level()

      shiny::req(pgx$X, pgx$Y, pgx$gsetX, pgx$genes)

      features <- X <- NULL
      if (hm_level == "geneset") {
        features <- COLLECTIONS
        X <- pgx$gsetX
      } else {
        features <- pgx$families
        X <- pgx$X
      }

      ## ------------ intersect features, set minimum set size
      rownames(X) <- toupper(rownames(X))
      genes <- toupper(rownames(X))
      features <- lapply(features, toupper)
      features <- lapply(features, function(f) intersect(toupper(f), genes))
      features <- features[sapply(features, length) >= 10]

      ## ------------ Just to get current samples
      ## samples = colnames(X)
      samples <- selectSamplesFromSelectedLevels(pgx$Y, hm_samplefilter())
      X <- X[, samples]
      cvar <- pgx.getCategoricalPhenotypes(pgx$Y, max.ncat = 999)
      cvar <- grep("sample|patient|years|days|months|gender",
        cvar,
        invert = TRUE, value = TRUE
      ) ## no sample IDs
      cvar
      Y <- pgx$Y[colnames(X), cvar, drop = FALSE]
      kk <- which(apply(Y, 2, function(y) length(unique(y)) > 1))
      Y <- Y[, kk, drop = FALSE]
      dim(Y)

      ## ------------ Note: this takes a while. Maybe better precompute off-line...
      sdx <- apply(X, 1, sd)
      names(sdx) <- rownames(X)
      S <- matrix(NA, nrow = length(features), ncol = ncol(Y))
      rownames(S) <- names(features)
      colnames(S) <- colnames(Y)

      ## ------------ Create a Progress object
      if (!interactive()) {
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Calculating feature-set scores", value = 0)
      }

      gene.level <- TRUE
      gene.level <- (hm_level == "gene")
      i <- 1
      for (i in 1:ncol(Y)) {
        if (!interactive()) progress$inc(1 / ncol(Y))

        grp <- Y[, i]
        grp <- as.character(grp)

        cat("[calcFeatureRanking] head(grp)=", head(grp), "\n")

        score <- rep(NA, length(features))
        names(score) <- names(features)
        j <- 1
        for (j in 1:length(features)) {
          pp <- features[[j]]
          if (gene.level) {
            pp <- filterProbes(pgx$genes, features[[j]])
          }
          pp <- head(pp[order(-sdx[pp])], 1000) ## how many top SD??
          pp <- intersect(pp, rownames(X))
          X1 <- X[pp, , drop = FALSE]
          dim(X1)
          ## cat("<clust_featureRank> dim(X1)=",dim(X1),"\n")
          ## if( nrow(X1)

          s1 <- s2 <- 1
          method <- input$clust_featureRank_method
          if (method %in% c("correlation", "meta")) {
            mx <- t(apply(X1, 1, function(x) tapply(x, grp, mean)))
            if (nrow(mx) == 0 || ncol(mx) == 0) next
            D <- 1 - cor(mx, use = "pairwise")
            diag(D) <- NA
            s1 <- mean(D, na.rm = TRUE)
          }

          if (method %in% c("p-value", "meta")) {
            jj <- which(!is.na(grp))
            design <- model.matrix(~ grp[jj])
            suppressWarnings(fit <- limma::eBayes(limma::lmFit(X1[, jj], design)))
            suppressWarnings(suppressMessages(top <- limma::topTable(fit)))
            ## s2 = mean(-log10(top$P.Value))  ## as score
            s2 <- mean(-log10(1e-99 + top$adj.P.Val), na.rm = TRUE) ## as score
          }

          f <- 1
          f <- (1 - exp(-(length(pp) / 20)**2)) ## penalize smaller sets
          score[j] <- f * (s1 * s2)**ifelse(method == "meta", 0.5, 1)
        }
        S[, i] <- score
      }
      S[is.na(S)] <- 0 ## missing values
      return(S)
    })

    clust_featureRank.RENDER <- shiny::reactive({

      S <- calcFeatureRanking()
      if (is.null(S) || nrow(S) == 0 || ncol(S) == 0) {
        return(NULL)
      }

      ## top scoring
      S <- tail(S[order(rowSums(S)), , drop = FALSE], 35)
      rownames(S) <- substring(rownames(S), 1, 80)
      
      pgx.stackedBarplot(
        x = t(S),
        showlegend = TRUE,
        xlab = "Discriminant score",
        ylab = "",
        horiz = TRUE
      )
      #%>%
      #  layout(legend = list(orientation = "h",   # show entries horizontally
      #    xanchor = "center",  # use center of legend as anchor
      #    x = 0.5))
    })


    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      ## plotlib2 = "plotly",
      func = clust_featureRank.RENDER,
      csvFunc = calcFeatureRanking, ##  *** downloadable data as CSV
      ## renderFunc = plotly::renderPlotly,
      ## renderFunc2 = plotly::renderPlotly,
      res = c(72, 90), ## resolution of plots
      pdf.width = 8, pdf.height = 10,
      add.watermark = watermark
    )
  })
}
