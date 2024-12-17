##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


##' Clustering board server module
##'
##' .. content for \details{} ..
##' @title
##' @param id
##' @param pgx
##' @return
##' @author kwee
ClusteringBoard <- function(id, pgx, labeltype = shiny::reactive("feature")) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    fullH <- 850 ## full height of page

    clust_infotext <-
      '<center><iframe width="1120" height="630" src="https://www.youtube.com/embed/hyDEk_MCaTk"
       title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write;
       encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></center>'

    ## ===================================================================================
    ## ======================== OBSERVERS ================================================
    ## ===================================================================================

    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Heatmap" = list(
        enable = NULL,
        disable = c("hm_clustmethod")
      ),
      "PCA/tSNE" = list(
        enable = NULL,
        disable = c("hm_features", "hm_splitby", "hm_level", "hm_filterXY", "hm_filterMitoRibo", "hm_topmode", "hm_ntop", "hm_clustk", "hm_scale", "cluster_bar", "spliby_bar")
      ),
      "Parallel" = list(
        enable = NULL,
        disable = c("selected_phenotypes", "hm_clustmethod", "pheno_bar")
      )
    )
    shiny::observeEvent(input$tabs1, {
      bigdash::update_tab_elements(input$tabs1, tab_elements)
    })

    # Board info
    shiny::observeEvent(input$board_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Clustering Board</strong>"),
        shiny::HTML(clust_infotext),
        easyClose = TRUE, size = "xl"
      ))
    })

    shiny::observeEvent(pgx$Y, {
      shiny::req(pgx$Y)
      ## input$menuitem  ## upon menuitem change
      var.types <- playbase::pgx.getCategoricalPhenotypes(pgx$samples, min.ncat = 2, max.ncat = 999)
      var.types <- var.types[grep("sample|patient", var.types, invert = TRUE)]
      vv <- c(var.types, rep("<none>", 10))
      var.types0 <- c("<none>", "<cluster>", var.types)
      var.types0 <- c("<none>", var.types)
      var.types1 <- c("<none>", var.types)
      sel <- grep("^[.]", var.types, invert = TRUE, value = TRUE)[1]
      if ("condition" %in% tolower(var.types)) {
        sel <- grep("condition", var.types,
          value = TRUE,
          ignore.case = TRUE
        )
      }
      if ("group" %in% tolower(var.types)) {
        sel <- grep("group", var.types,
          value = TRUE,
          ignore.case = TRUE
        )
      }
      shiny::updateSelectInput(session, paste0("PCAplot-", "hmpca.colvar"), choices = var.types0, selected = sel)
      shiny::updateSelectInput(session, paste0("PCAplot-", "hmpca.shapevar"), choices = var.types1, selected = "<none>")
      shiny::updateSelectInput(session, "selected_phenotypes", choices = var.types, selected = head(var.types, 6))
      choices <- c("none", "phenotype", "gene")
      choices_names <- c("none", "phenotype", tspan("gene", js = FALSE))
      names(choices) <- choices_names
      shiny::updateRadioButtons(session, "hm_splitby", choices = choices)
      shiny::updateTextAreaInput(session, "hm_customfeatures", placeholder = tspan("Paste your custom gene list", js = FALSE))
      hm_level_choices <- c("gene", "geneset")
      hm_level_choices_names <- c(tspan("gene", js = FALSE), tspan("geneset", js = FALSE))
      names(hm_level_choices) <- hm_level_choices_names
      shiny::updateSelectInput(session, "hm_level", choices = hm_level_choices)
    })

    ## update filter choices upon change of data set
    shiny::observeEvent(
      {
        list(pgx$X, pgx$Y, pgx$samples)
      },
      {
        shiny::req(pgx$X, pgx$Y, pgx$samples)
        levels <- playbase::getLevels(pgx$Y)
        shiny::updateSelectInput(session, "hm_samplefilter", choices = levels)

        ## update defaults??
        n1 <- nrow(pgx$samples) - 1
        groupings <- colnames(pgx$samples)
        groupings <- sort(groupings)
        contrasts <- playbase::pgx.getContrasts(pgx)
        shiny::updateSelectInput(session, "hm_contrast", choices = contrasts)

        ## get clusterings methods
        clustmethods <- grep("2d$", names(pgx$cluster$pos), value = TRUE)
        clustmethods <- sort(unique(sub("2d$", "", clustmethods)))
        selmethod <- ifelse("umap" %in% clustmethods, "umap", clustmethods[1])
        shiny::updateSelectInput(session, "hm_clustmethod",
          choices = clustmethods, sel = selmethod
        )
      }
    )

    ## update choices upon change of level
    shiny::observeEvent(
      {
        c(input$hm_splitvar, input$hm_level)
      },
      {
        shiny::req(pgx$families, pgx$gsetX)
        shiny::req(input$hm_level)
        choices <- names(pgx$families)
        if (input$hm_level == "geneset") {
          gset_collections <- playbase::pgx.getGeneSetCollections(gsets = rownames(pgx$gsetX))
          nk <- sapply(gset_collections, function(k) sum(k %in% rownames(pgx$gsetX)))
          choices <- names(gset_collections)[nk >= 5]
          subset_choices <- sapply(choices, function(x) any(x == input$hm_splitvar))
          choices <- names(subset_choices)
        }
        choices <- c("<custom>", "<contrast>", choices)
        choices <- sort(unique(choices))
        shiny::updateSelectInput(session, "hm_features", choices = choices)
      }
    )

    # reactive functions ##############

    shiny::observeEvent(
      {
        list(input$hm_splitby, pgx$X, pgx$samples)
      },
      {
        shiny::req(pgx$X, pgx$samples, input$hm_splitby)
        if (input$hm_splitby == "none") {
          return()
        }
        if (input$hm_splitby == "gene") {
          xgenes <- sort(rownames(getFilteredMatrix()$zx))
          shiny::updateSelectizeInput(session, "hm_splitvar", choices = xgenes, server = TRUE)
        }
        if (input$hm_splitby == "phenotype") {
          cvar <- sort(playbase::pgx.getCategoricalPhenotypes(pgx$samples, min.ncat = 2, max.ncat = 999))
          sel <- cvar[1]
          cvar0 <- grep("^[.]", cvar, value = TRUE, invert = TRUE) ## no estimated vars
          sel <- head(c(
            grep("type|family|class|stat", cvar0, ignore.case = TRUE, value = TRUE),
            cvar0, cvar
          ), 1)
          shiny::updateSelectInput(session, "hm_splitvar", choices = cvar, selected = sel)
        }
      }
    )

    ## update filter choices upon change of data set
    shiny::observeEvent(pgx$X, {
      shiny::req(pgx$X)
      shiny::updateRadioButtons(session, "hm_splitby", selected = "none")
    })

    ## ===================================================================================
    ## ============================= REACTIVES ===========================================
    ## ===================================================================================

    ##' Returns filtered matrix ready for clustering. Filtering based
    ##' on user selected geneset/features or custom list of genes.
    ##'
    ##' @title
    ##' @param id
    ##' @param pgx
    ##' @return
    getFilteredMatrix <- shiny::reactive({
      shiny::req(pgx$X, pgx$Y, pgx$gsetX, pgx$families, pgx$genes)

      ## NEED RETHINK!!!!! THIS CREATED PROBLEMS.
      if (!pgx$organism %in% c("Human", "human")) {
        genes <- pgx$genes[rownames(pgx$X), c("gene_name", "human_ortholog")]
        genes <- ifelse(genes$human_ortholog == "" | is.na(genes$human_ortholog),
          rownames(genes), genes$human_ortholog
        )
      } else {
        genes <- as.character(rownames(pgx$genes[rownames(pgx$X), ]))
      }
      genesets <- rownames(pgx$gsetX)

      ft <- input$hm_features
      shiny::req(ft)

      if (input$hm_level == "geneset") {
        ## Gene set level features #########

        gsets <- rownames(pgx$gsetX)
        gset_collections <- playbase::pgx.getGeneSetCollections(gsets = rownames(pgx$gsetX))
        gsets <- unique(gset_collections[[ft]])
        zx <- pgx$gsetX
        if (input$hm_customfeatures != "") {
          gsets1 <- genesets[grep(input$hm_customfeatures, genesets, ignore.case = TRUE)]
          if (length(gsets1) > 2) gsets <- gsets1
        }
        zx <- zx[intersect(gsets, rownames(zx)), ]
      }

      idx <- NULL
      if (input$hm_level == "gene") {
        ## Gene level features ###########
        gg <- pgx$families[[1]]

        if (ft == "<all>") {
          gg <- rownames(pgx$X)
        } else if (ft == "<contrast>") {
          ct <- input$hm_contrast
          shiny::req(ct)
          shiny::req(input$hm_ntop)
          fc <- names(sort(playbase::pgx.getMetaMatrix(pgx)$fc[, ct]))
          n1 <- floor(as.integer(input$hm_ntop) / 2)
          gg <- unique(c(head(fc, n1), tail(fc, n1)))
          if (input$hm_splitby == "gene") {
            if (!(input$hm_splitvar %in% gg)) {
              gg <- c(input$hm_splitvar, gg)
            }
          }
        } else if (ft %in% names(pgx$families)) {
          gg <- pgx$families[[ft]]
        } else if (ft == "<custom>" && ft != "") {
          message("[getFilteredMatrix] selecting for <custom> features")
          customfeatures <- "ADORA2A ARHGEF5 BTLA CD160 CD244 CD27 CD274 CD276 CD47 CD80 CEACAM1 CTLA4 GEM HAVCR2 ICOS IDO1 LAG3"
          customfeatures <- "CTLA4 GEM HAVCR2 ICOS IDO1 LAG3 PDCD1 TNFSF4 VISTA VTCN1 TIGIT PVR --- CD28 CD40 CD40LG ICOSLG TNFRSF9 TNFSF9 CD70 TNFRSF4 TNFRSF18 --- TNFSF18 SIRPA LGALS9 ARG1 CD86 IDO2 PDCD1LG2 KIR2DL3"
          customfeatures <- input$hm_customfeatures
          gg1 <- strsplit(customfeatures, split = "[, ;\n\t]")[[1]]

          is.regx <- grepl("[*.?\\[]", gg1[1])
          if (length(gg1) == 1 && is.regx) {
            gg1 <- grep(gg1, genes, ignore.case = TRUE, value = TRUE)
          }
          if (length(gg1) == 1 && !is.regx) {
            gg1 <- c(gg1, gg1) ## heatmap does not like single gene
          }

          gg1 <- gg1[toupper(gg1) %in% toupper(genes) | grepl("---", gg1)]
          idx <- NULL
          if (any(grepl("^---", gg1))) {
            message("[getFilteredMatrix] <custom> groups detected")
            idx <- rep("F1", length(gg1))
            names(idx) <- gg1
            kk <- c(1, grep("^---", gg1), length(gg1) + 1)
            for (i in 1:(length(kk) - 1)) {
              ii <- kk[i]:(kk[i + 1] - 1)
              idx[ii] <- paste0("F", i)
            }
            gg1 <- gg1[grep("---", gg1, invert = TRUE)]
            idx <- idx[gg1]
          }
          gg <- gg1
        } else {
          warning("[getFilteredMatrix] ERROR!!:: switch error : ft= ", ft)
          gg <- NULL
          return(NULL)
        }

        gg <- gg[which(toupper(gg) %in% toupper(genes))]
        if (length(gg) == 0) {
          warning("[getFilteredMatrix] warning: no genes overlap with filter")
          return(NULL)
        }

        if (input$hm_splitby == "gene") {
          if (!(input$hm_splitvar %in% gg)) {
            gg <- c(input$hm_splitvar, gg)
          }
        }
        jj <- match(toupper(gg), toupper(genes))
        pp <- rownames(pgx$X)[jj]
        zx <- pgx$X[pp, , drop = FALSE]
        if (!is.null(idx)) {
          idx <- idx[gg]
          names(idx) <- rownames(zx)
        }
      }

      if (nrow(zx) == 0) {
        return(NULL)
      }

      kk <- playbase::selectSamplesFromSelectedLevels(pgx$Y, input$hm_samplefilter)
      zx <- zx[, kk, drop = FALSE]

      if (input$hm_level == "gene" &&
        "chr" %in% names(pgx$genes) &&
        input$hm_filterXY) {
        ## Filter out X/Y chromosomes before clustering
        chr.col <- grep("^chr$|^chrom$", colnames(pgx$genes))
        chr <- pgx$genes[rownames(zx), chr.col]
        not.xy <- !(chr %in% c("X", "Y", 23, 24)) & !grepl("^X|^Y|chrX|chrY", chr)
        table(not.xy)
        zx <- zx[which(not.xy), ]
        if (!is.null(idx)) idx <- idx[rownames(zx)]
      }

      if (input$hm_level == "gene" && input$hm_filterMitoRibo) {
        ## Filter out X/Y chromosomes before clustering
        is.ribomito <- grepl("^RP[LS]|^MT-", rownames(zx), ignore.case = TRUE)
        table(is.ribomito)
        zx <- zx[which(!is.ribomito), , drop = FALSE]
        if (!is.null(idx)) idx <- idx[rownames(zx)]
      }
      shiny::validate(shiny::need(
        ncol(zx) > 0, "Filtering too restrictive. Please change 'Filter samples' settings."
      ))

      flt <- list(
        zx = zx,
        idx = idx
      )

      return(flt) ## end of getFilteredMatrix
    }) %>% bindEvent(
      input$hm_samplefilter,
      input$hm_features,
      input$hm_level,
      input$hm_customfeatures,
      input$hm_samplefilter,
      input$hm_filterXY,
      input$hm_filterMitoRibo,
      input$hm_contrast,
      pgx$X,
      ## input$hm_group,
      input$hm_ntop
    )


    ##' .. content for \description{} (no empty lines) ..
    ##'
    ##' .. content for \details{} ..
    ##' @title
    ##' @param id
    ##' @param pgx
    ##' @return
    ##' @author kwee
    getTopMatrix <- shiny::reactive({
      shiny::req(pgx$X, pgx$samples)

      flt <- getFilteredMatrix()

      zx <- flt$zx
      if (is.null(flt)) {
        return(NULL)
      }
      if (is.null(zx) || nrow(zx) == 0) {
        return(NULL)
      }

      nmax <- 4000
      nmax <- as.integer(input$hm_ntop)
      idx <- NULL
      splitvar <- "none"
      splitvar <- input$hm_splitvar
      splitby <- input$hm_splitby
      do.split <- splitby != "none"

      if (splitby == "gene" && !splitvar %in% rownames(pgx$X)) {
        return(NULL)
      }
      if (splitby == "phenotype" && !splitvar %in% colnames(pgx$samples)) {
        return(NULL)
      }

      grp <- NULL
      ## split on a phenotype variable
      if (do.split && splitvar %in% colnames(pgx$samples)) {
        grp <- pgx$samples[colnames(zx), splitvar]
        table(grp)
      }

      ## split on gene expression value: hi vs. low
      if (do.split && splitvar %in% rownames(pgx$X)) {
        gx <- pgx$X[1, ]
        gx <- pgx$X[splitvar, colnames(zx)]

        ## TODO if this code is revived again, for some datasets
        ## the number of unique values in gx can be equal or
        ## lower than k, on those cases it will crash

        ## estimate best K
        # within.ssratio <- sapply(1:4, function(k) {
        #   km <- kmeans(gx, k)
        #   km$tot.withinss / km$totss
        # })
        # within.ssratio
        # k.est <- min(which(within.ssratio < 0.10))
        # k.est <- min(which(abs(diff(within.ssratio)) < 0.10))
        # k.est
        # k.est <- pmax(pmin(k.est, 3), 2)
        k.est <- 2 ## for now...

        if (k.est == 2) {
          km <- kmeans(gx, centers = 2)
          km.rnk <- rank(km$centers, ties.method = "random")
          grp.labels <- c("low", "high")[km.rnk]
          grp <- grp.labels[km$cluster]
        } else if (k.est == 3) {
          km <- kmeans(gx, centers = 3)
          km.rnk <- rank(km$centers, ties.method = "random")
          grp.labels <- c("low", "mid", "high")[km.rnk]
          grp <- grp.labels[km$cluster]
        } else if (k.est == 4) {
          km <- kmeans(gx, centers = 4)
          km.rnk <- rank(km$centers, ties.method = "random")
          grp.labels <- c("low", "mid-low", "mid-high", "high")[km.rnk]
          grp <- grp.labels[km$cluster]
        }
        grp <- paste0(splitvar, ":", grp)
        names(grp) <- colnames(zx)
      }

      if (do.split && length(grp) == 0) {
        return(NULL)
      }

      ## Any BMC scaling?? ##########
      if (do.split && input$hm_scale == "BMC") {
        for (g in unique(grp)) {
          jj <- which(grp == g)
          zx1 <- zx[, jj, drop = FALSE]
          zx[, jj] <- zx1 - rowMeans(zx1, na.rm = TRUE)
        }
      }

      ## Create reduced matrix according to topmode #######
      topmode <- "marker"
      topmode <- "sd"
      topmode <- input$hm_topmode
      if (topmode == "marker" && length(table(grp)) <= 1) {
        topmode <- "sd"
      }
      if (!do.split && topmode == "marker") {
        topmode <- "sd"
      }

      addsplitgene <- function(gg) {
        if (do.split && splitvar %in% rownames(pgx$X)) {
          gg <- unique(c(splitvar, gg))
        }
        gg
      }

      grp.zx <- NULL
      if (topmode == "pca") {
        NPCA <- 5
        svdres <- irlba::irlba(zx - rowMeans(zx, na.rm = TRUE), nv = NPCA)
        ntop <- 12
        ntop <- as.integer(input$hm_ntop) / NPCA
        gg <- rownames(zx)
        sv.top <- lapply(1:NPCA, function(i) gg[head(order(-abs(svdres$u[, i])), ntop)])
        gg.top <- unlist(sv.top)
        for (i in 1:length(sv.top)) {
          sv.top[[i]] <- paste0("PC", i, ":", sv.top[[i]])
        }
        sv.top1 <- unlist(sv.top)
        zx <- zx[gg.top, , drop = FALSE]
        dim(zx)
        idx <- sub(":.*", "", sv.top1)
        table(idx)
      } else if (topmode == "marker" && splitby != "none") {
        ## sample cluster marker gene prioritazion
        grp.zx <- t(apply(zx, 1, function(x) tapply(x, grp, mean)))
        if (length(table(grp)) == 1) {
          grp.zx <- t(grp.zx)
          colnames(grp.zx) <- paste0(splitvar, ":", grp[1])
        }
        grp.dx <- grp.zx * 0
        nc <- ncol(grp.dx)
        for (i in 1:nc) {
          grp.dx[, i] <- grp.zx[, i] - rowMeans(grp.zx[, -i, drop = FALSE], na.rm = TRUE)
        }
        gg <- rownames(zx)
        ntop <- 12
        ntop <- ceiling(as.integer(input$hm_ntop) / ncol(grp.dx))
        grp.top <- lapply(1:nc, function(i) gg[head(order(-grp.dx[, i]), ntop)])
        idx <- unlist(mapply(rep, 1:nc, sapply(grp.top, length)))
        idx <- paste0("M", idx)
        table(idx)
        gg.top <- unlist(grp.top)
        zx <- zx[gg.top, , drop = FALSE]
        dim(zx)
      } else {
        ## Order by SD
        gg <- rownames(zx)[order(-apply(zx, 1, sd, na.rm = TRUE))]
        gg <- head(gg, nmax)
        gg <- addsplitgene(gg)
        zx <- zx[gg, , drop = FALSE] ## order
      }

      ## ------------- cluster the genes???
      if (!is.null(flt$idx)) {
        idx <- flt$idx[rownames(zx)] ## override
      }

      CLUSTK <- 4 ## number of gene groups (NEED RETHINK)
      CLUSTK <- as.integer(input$hm_clustk)
      if (is.null(idx)) {
        D <- as.dist(1 - cor(t(zx), use = "pairwise"))
        D[which(is.nan(D) | is.na(D))] <- 1
        system.time(hc <- fastcluster::hclust(D, method = "ward.D2"))
        ngrp <- min(CLUSTK, nrow(zx)) ## how many default groups???
        idx <- paste0("S", cutree(hc, ngrp))
      }

      ## ------------- matched annotation
      annot <- pgx$samples[colnames(zx), , drop = FALSE]
      kk <- grep("sample|patient", colnames(annot), invert = TRUE)
      annot <- annot[, kk, drop = FALSE] ## no group??
      samples <- colnames(zx) ## original sample list

      ## ----------------------------------------------------
      ## ------------ calculate group summarized ------------
      ## ----------------------------------------------------
      grp.zx <- NULL
      grp.var <- "group"
      grp.average <- input$hm_average_group
      is.grouped <- FALSE

      if (grp.average && input$hm_splitby != "none") {
        grp.zx <- tapply(colnames(zx), grp, function(k) {
          rowMeans(zx[, k, drop = FALSE], na.rm = TRUE)
        })
        grp.zx <- do.call(cbind, grp.zx)
        most.freq <- function(x) names(sort(-table(x)))[1]
        grp.annot <- tapply(rownames(annot), grp, function(k) {
          f <- apply(annot[k, , drop = FALSE], 2, function(x) most.freq(x))
          w.null <- sapply(f, is.null)
          if (any(w.null)) f[which(w.null)] <- NA
          unlist(f)
        })
        grp.annot <- data.frame(do.call(rbind, grp.annot))
        grp.annot <- grp.annot[colnames(grp.zx), , drop = FALSE]
        grp <- colnames(grp.zx)
        is.grouped <- TRUE
      } else {
        grp.zx <- zx
        grp.annot <- annot
      }

      ## input$top_terms
      filt <- list(
        mat = grp.zx,
        annot = grp.annot,
        grp = grp, ## sample grouping
        idx = idx, ## gene grouping
        samples = samples,
        is.grouped = is.grouped
      )

      return(filt)
    }) ## end of getTopMatrix


    ##' Same as getTopMatrix but always averaged by group
    ##'
    ##' .. content for \details{} ..
    ##' @title
    ##' @param id
    ##' @param pgx
    ##' @return
    ##' @author kwee
    getTopMatrixGrouped <- shiny::reactive({
      topmat <- getTopMatrix()
      shiny::req(topmat$mat)

      if (topmat$is.grouped || is.null(topmat$grp)) {
        return(topmat)
      }

      ## ----------------------------------------------------
      ## ------------ calculate group summarized ------------
      ## ----------------------------------------------------

      grp <- topmat$grp
      zx <- topmat$mat
      annot <- topmat$annot
      idx <- topmat$idx
      samples <- topmat$samples

      ## take most frequent term as group annotation value
      grp.zx <- tapply(colnames(zx), grp, function(k) {
        rowMeans(zx[, k, drop = FALSE], na.rm = TRUE)
      })
      grp.zx <- do.call(cbind, grp.zx)
      most.freq <- function(x) names(sort(-table(x)))[1]
      grp.annot <- tapply(rownames(annot), grp, function(k) {
        f <- apply(annot[k, , drop = FALSE], 2, function(x) most.freq(x))
        w.null <- sapply(f, is.null)
        if (any(w.null)) f[which(w.null)] <- NA
        unlist(f)
      })
      grp.annot <- data.frame(do.call(rbind, grp.annot))
      grp.annot <- grp.annot[colnames(grp.zx), , drop = FALSE]
      grp <- colnames(grp.zx)

      res <- list(
        mat = grp.zx,
        annot = grp.annot,
        grp = grp, ## sample grouping
        idx = idx, ## gene grouping
        samples = samples,
        is.grouped = TRUE
      )
      return(res)
    })


    getClustAnnotCorrelation <- shiny::reactive({
      shiny::req(pgx$X, pgx$Y, pgx$gsetX, pgx$families)
      filt <- getTopMatrix()
      shiny::req(filt)
      zx <- filt$mat
      idx <- filt$idx
      samples <- filt$samples

      if (nrow(zx) <= 1) {
        return(NULL)
      }

      ann.level <- "geneset"
      ann.refset <- "Hallmark collection"
      ann.level <- clusterannot$xann_level()
      ann.refset <- clusterannot$xann_refset()
      shiny::req(clusterannot$xann_level(), clusterannot$xann_refset())

      ref <- NULL
      ref <- pgx$gsetX[, , drop = FALSE]
      ref <- pgx$X[, , drop = FALSE]
      if (ann.level == "gene" && ann.refset %in% names(pgx$families)) {
        gg <- pgx$families[[ann.refset]]
        jj <- match(toupper(gg), toupper(pgx$genes$symbol))
        jj <- setdiff(jj, NA)
        pp <- rownames(pgx$genes)[jj]
        ref <- pgx$X[intersect(pp, rownames(pgx$X)), , drop = FALSE]
      }
      gset_collections <- playbase::pgx.getGeneSetCollections(gsets = rownames(pgx$gsetX))
      if (ann.level == "geneset" && ann.refset %in% names(gset_collections)) {
        ss <- gset_collections[[ann.refset]]
        ss <- intersect(ss, rownames(pgx$gsetX))
        length(ss)
        ref <- pgx$gsetX[ss, ]
      }
      if (ann.level == "phenotype") {
        ref <- t(playbase::expandAnnotationMatrix(pgx$Y))
      }
      if (is.null(ref)) {
        cat("<clustering:getClustAnnotCorrelation> WARNING:: ref error\n")
        return(NULL)
      }

      ## -----------  restrict to top??
      dim(ref)
      if (nrow(ref) > 1000) {
        ref <- head(ref[order(-apply(ref, 1, sd, na.rm = TRUE)), ], 1000)
      }

      ## -----------  get original data level
      X <- pgx$X
      if (input$hm_level == "geneset") X <- pgx$gsetX

      ## ----------- for each gene cluster compute average correlation
      hm_topmode <- "sd"
      hm_topmode <- input$hm_topmode
      idxx <- setdiff(idx, c(NA, " ", "   "))
      rho <- matrix(NA, nrow(ref), length(idxx))
      colnames(rho) <- idxx
      rownames(rho) <- sub(".*:", "", rownames(ref))

      if (nrow(ref) > 0) {
        for (i in 1:length(idxx)) {
          gg <- rownames(zx)[which(idx == idxx[i])]
          aa <- t(X[gg, samples, drop = FALSE])
          bb <- t(ref[, samples, drop = FALSE])
          rr <- cor(apply(aa, 2, rank), apply(bb, 2, rank), use = "pairwise")
          if (hm_topmode == "pca") rr <- abs(rr)
          rho[, i] <- colMeans(rr, na.rm = TRUE)
        }
      }

      if (input$hm_level == "gene" && ann.level == "geneset" && clusterannot$xann_odds_weighting()) {
        table(idx)
        grp <- tapply(toupper(rownames(zx)), idx, list) ## toupper for mouse!!
        gmt <- playbase::getGSETS_playbase(sub("_", ":", rownames(rho)))
        bg.genes <- toupper(rownames(X))
        P <- c()
        for (i in 1:ncol(rho)) {
          k <- colnames(rho)[i]
          res <- playbase::gset.fisher(
            grp[[k]], gmt,
            fdr = 1, min.genes = 0, max.genes = Inf,
            background = bg.genes
          )
          r <- res[, "odd.ratio"]
          odd.prob <- r / (1 + r)
          P <- cbind(P, odd.prob)
        }
        colnames(P) <- colnames(rho)
        rownames(P) <- sub(":", "_", names(gmt))
        rho <- rho[rownames(P), ]
        rho <- rho * (P / max(P))
      }

      return(rho)
    })


    selected_samples <- reactive({
      playbase::selectSamplesFromSelectedLevels(pgx$Y, input$hm_samplefilter)
    })

    # plots ##########
    clustering_plot_splitmap_server(
      id = "splitmap",
      pgx = pgx,
      getTopMatrix = getTopMatrix,
      selected_phenotypes = shiny::reactive(input$selected_phenotypes),
      hm_level = shiny::reactive(input$hm_level),
      hm_ntop = shiny::reactive(input$hm_ntop),
      hm_scale = shiny::reactive(input$hm_scale),
      hm_topmode = shiny::reactive(input$hm_topmode),
      hm_clustk = shiny::reactive(input$hm_clustk),
      watermark = WATERMARK,
      labeltype = labeltype
    )

    clustering_plot_clustpca_server("PCAplot",
      pgx = pgx,
      selected_samples = selected_samples,
      clustmethod = shiny::reactive(input$hm_clustmethod),
      watermark = WATERMARK,
      parent = ns
    )

    clustering_plot_table_parcoord_server(
      id = "parcoord",
      getTopMatrix = getTopMatrixGrouped,
      watermark = WATERMARK
    )

    clustering_plot_genemodule_server(
      id = "genemodule",
      getTopMatrix = getTopMatrixGrouped,
      watermark = WATERMARK
    )

    clustering_plot_phenoplot_server(
      id = "clust_phenoplot",
      pgx = pgx,
      selected_phenotypes = shiny::reactive(input$selected_phenotypes),
      clustmethod = shiny::reactive(input$hm_clustmethod),
      selected_samples = selected_samples,
      watermark = WATERMARK
    )

    clusterannot <- clustering_plot_clusterannot_server(
      id = "plots_clustannot",
      pgx,
      getClustAnnotCorrelation = getClustAnnotCorrelation,
      watermark = WATERMARK
    )

    # tables ##########
    clustering_table_clustannot_server(
      id = "tables_clustannot",
      getClustAnnotCorrelation = getClustAnnotCorrelation,
      xann_level = clusterannot$xann_level,
      scrollY = "calc(40vh - 236px)",
      watermark = WATERMARK
    )
  }) ## end of moduleServer
} ## end of Board
