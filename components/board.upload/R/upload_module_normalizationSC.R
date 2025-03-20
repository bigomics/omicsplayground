##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


## =========================================================================
## ===================== NORMALIZATION UI/SERVER ===========================
## =========================================================================

upload_module_normalizationSC_ui <- function(id, height = "100%") {
  ns <- shiny::NS(id)
  uiOutput(ns("normalization"), fill = TRUE)
}

upload_module_normalizationSC_server <- function(id,
                                                 r_counts,
                                                 r_samples,
                                                 r_contrasts,
                                                 upload_datatype,
                                                 is.count = FALSE,
                                                 height = 720) {

  shiny::moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
      
      ## ------------------------------------------------------------------
      ## Plot UI
      ## ------------------------------------------------------------------

      output$normalization <- shiny::renderUI({

        shiny::req(r_counts())
        counts <- r_counts()        
        nFeature_RNA <- Matrix::colSums(counts > 0, na.rm = TRUE)
        qF <- quantile(nFeature_RNA, probs = c(0.01, 0.99))
        
        shiny::req(r_samples())
        samples <- r_samples()

        metadata_vars0 <- colnames(samples)
        metadata_vars1 <- c(
          "orig.ident", "seurat_clusters",
          "nCount_RNA", "nFeature_RNA",
          "percent.mt", "percent.ribo",
          "percent.hb", "G2M.Score",
          "S.Score", "celltype.azimuth"
        )
        metadata_vars <- unique(c(metadata_vars0, metadata_vars1))
        
        dimred.infotext <- "Dimensionality reduction enables to simplify large datasets by computing representative data points capable of preserving the biological information while reducing the dimensionality of the data. Here we employ the two most widely used non-linear methods for dimensional reduction of single-cell RNA-seq data: T-distributed stochastic neighbor embedding (t-SNE), and Unifold Manifold Approximation and Projection (UMAP). https://omicsplayground.readthedocs.io/en/latest/methods/#clustering"

        cellqc.infotext <- "Data QC. Violin plots of total number of cDNA molecules (e.g., UMI) detected in each cell (nCount_RNA), number of unique genes detected in each cell (nFeature_RNA), percentage of mitochondrial gene expression in each cell (percent_mt), percentage of ribosomal gene expression in each cell (percent_ribo), percentage of globin gene expression in each cell (percent_hb), G2M cell cysle score, S cell cycle score."
        
        dimred.options <- tagList(
          shiny::radioButtons(
            ns("dimred_plottype"),
            label = "Dimensional reduction method:",
            choices = c("PCA", "tSNE", "UMAP"),
            selected = "tSNE",
            inline = FALSE
          )
        )

        cellstats.options <- tagList(
          shiny::checkboxInput(
            ns("groupby_celltype"),
            label = "Group by cell type",
            value = FALSE
          )
        )

        navmenu <- tagList(
          
          bslib::card(bslib::card_body(
            style = "padding: 0px;",
            bslib::accordion(
              multiple = FALSE,
              style = "background-color: #F7FAFD99;",

              bslib::accordion_panel(
                title = HTML("<span style='font-size: 1em;'> Single-cell transcriptomics</span>"),
                shiny::div(style = "display: flex; align-items: center; justify-content: space-between;"),
                shiny::div(
                  style = "background-color: #f0f0f0; padding: 2px; border-radius: 5px; margin-top: 0.01px;",
                  HTML(paste0("<span style='font-size: 1em;'>Number of features: </span>",
                    "<span style='font-size: 1em;'>", nrow(counts), "</span><br>",
                    "<span style='font-size: 1em;'>Number of cells: </span>",
                    "<span style='font-size: 1em;'>", ncol(counts), "</span>"))
                )
              ),
              
              bslib::accordion_panel(
                title = HTML("<span style='font-size: 1em;'> Define cell types with Azimuth</span>"),
                shiny::div(style = "display: flex; align-items: center; justify-content: space-between;"),
                shiny::selectInput(
                  ns("ref_atlas"),
                  label = "Select reference atlas",
                  choices = c(
                    "adiposeref", "bonemarrowref", "fetusref",
                    "heartref", "humancortexref", "kidneyref",
                    "lungref", "mousecortexref", "pancreasref",
                    "pbmcref", "tonsilref", "<select>"
                  ),
                  selected = "<select>"
                ),
                shiny::br()
              ),
              
              bslib::accordion_panel(
                title = HTML("<span style='font-size: 1em;'> Visualize cell clusters</span>"),
                shiny::div(style = "display: flex; align-items: center; justify-content: space-between;"),
                shiny::selectInput(
                  ns("clusterBy"),
                  label = "Visualize cell cluster by",
                  choices = metadata_vars,
                  multiple = TRUE,
                  selected = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")
                ),
                shiny::br()
              ),

              bslib::accordion_panel(
                title = HTML("<span style='font-size: 1em;'> Cell filtering</span>"),
                shiny::p("Remove cells based on QC variables"),
                shiny::checkboxInput(
                  ns("remove_cells"),
                  label = "Remove cells",
                  value = FALSE
                ),
                shiny::conditionalPanel(
                  "input.remove_cells == true",
                  ns = ns,
                  shiny::sliderInput(ns("nfeature_threshold"), "Detected genes per cell:",
                    min = 0, max = max(nFeature_RNA), value = c(qF[1], qF[2])),
                  shiny::sliderInput(ns("mt_threshold"), "MT expresssion threshold (%):", 1, 100, 10, 0),
                  shiny::sliderInput(ns("hb_threshold"), "HB expression threshold (%):", 1, 100, 10, 0)
                ),
                shiny::br()
              )
            ))
          )

        )

        ## ---------------------------- UI ----------------------------------
        ui <- div(
          bslib::as_fill_carrier(),
          style = "width: 100%; display: flex; ",
          bslib::layout_columns(
            col_widths = c(2, 10),
            style = "margin-bottom: 20px;",
            heights_equal = "row",
            ## ----------- menu ------------
            navmenu,
            ## ----------- canvas ------------
            bslib::layout_columns(
              col_widths = c(6, 6),
              row_heights = c(3, 3),
              heights_equal = "row",
              PlotModuleUI(
                ns("plot1"),
                title = "Cell statistics",
                info.text = cellqc.infotext,
                caption = cellqc.infotext,
                options = cellstats.options,
                height = c("auto", "100%"),
                show.maximize = TRUE,
                ),
              PlotModuleUI(
                ns("plot2"),
                title = "Dimensional reduction",
                info.text = dimred.infotext,
                caption = dimred.infotext,
                options = dimred.options,
                height = c("auto", "100%"),
                show.maximize = TRUE
              ),
              )
          ),
          div(shiny::checkboxInput(ns("normalizationUI"), NULL, TRUE),
            style = "visibility:hidden")
        )            
        return(ui)
      })


      ## ------------------------------------------------------------------
      ## Object reactive chain
      ## ------------------------------------------------------------------

      ## downsampled & normalized data
      ## downsample if ncells > threshold
      ds_norm_Counts <- shiny::reactive({ 

        options(future.globals.maxSize= 4*1024^100)

        ## shiny::req(input$infercelltypes)
        shiny::req(r_counts())
        shiny::req(r_samples())
        counts <- r_counts()
        samples <- r_samples()
        if (is.null(counts)) { return(NULL) }

        kk <- intersect(colnames(counts), rownames(samples))
        counts <- counts[, kk, drop = FALSE]
        samples <- samples[kk, , drop = FALSE]

        ncells <- ncol(counts)
        cells_trs <- 1000
        dbg("[normalizationSC_server:ds_norm_Counts:] N.cells in dataset:", ncells)
        if (ncells > cells_trs) {
          dbg("[normalizationSC_server:ds_norm_Counts:] Random sampling of:", cells_trs, "cells.")
          kk <- sample(colnames(counts), cells_trs)
          counts <- counts[, kk, drop = FALSE] 
          samples <- samples[kk, , drop = FALSE]
        }

        ref_tissue <- input$ref_atlas
        if(ref_tissue != "<select>") {
          dbg("[normalizationSC_server:ds_norm_Counts:] Inferring cell types with Azimuth!")
          dbg("[normalizationSC_server:ds_norm_Counts:] Reference atlas:", ref_tissue)
          counts <- as(counts, "dgCMatrix")
          msg <- "Inferring cell types with Azimuth..."
          shiny::withProgress(message = msg, value = 0.1, {
            azm <- playbase::pgx.runAzimuth(counts = counts, reference = ref_tissue)
            dbg("[normalizationSC_server:ds_norm_Counts:] Cell types inferred.")
          })
          if (class(azm) %in% c("matrix", "data.frame")) {

            kk <- grep("^predicted.*l*2$", colnames(azm))
            if (any(kk)) {
              celltype.azimuth <- azm[, kk]
              if (length(unique(celltype.azimuth)) > 15) {
                ## allows better representation
                kk <- grep("^predicted.*l*1$", colnames(azm))
                celltype.azimuth <- azm[, kk]
              }
            } else {
              kk <- grep("^predicted.*subclass*$", colnames(azm))
              celltype.azimuth <- azm[, kk]
              if (length(unique(celltype.azimuth)) > 15) {
                ## allows better representation
                kk <- grep("^predicted.class*$", colnames(azm))
                celltype.azimuth <- azm[, kk]
              }
            }
            samples <- cbind(samples, celltype.azimuth = celltype.azimuth)          

          } else if (is.vector(azm)) {
            dbg("[normalizationSC_server:ds_norm_Counts:] Your selected Azimuth reference atlas *might* be incorrect. Please double check.")
            samples <- cbind(samples, celltype.azimuth = azm)
          }

          ## Normalization & Dimensional reduction
          dbg("[normalizationSC_server] Performing logCPM normalization...")
          nX <- playbase::logCPM(as.matrix(counts), prior = 1, total = 1e4)
          jj <- head(order(-matrixStats::rowSds(nX, na.rm = TRUE)), 250)
          nX1 <- nX[jj, , drop = FALSE]
          nX1 <- nX1 - rowMeans(nX1, na.rm = TRUE)
          nb <- ceiling(min(15, dim(nX) / 8))
          pos.list <- list()
          dbg("[normalizationSC_server] Performing PCA, tSNE and UMAP...")
          msg <- "Performing PCA, tSNE, UMAP..."
          shiny::withProgress(message = msg, value = 0.5, {
            pos.list[["pca"]] <- irlba::irlba(nX1, nv = 2, nu = 0)$v
            pos.list[["tsne"]] <- Rtsne::Rtsne(t(nX1), perplexity = 2*nb, check_duplicates=F)$Y
            pos.list[["umap"]] <- uwot::umap(t(nX1), n_neighbors = max(2, nb))
            pos.list <- lapply(pos.list, function(x) {
              rownames(x) <- colnames(nX1);
              return(x)
            })
          })
          dbg("[normalizationSC_server] PCA, tSNE & UMAP completed.")

          dbg("[normalizationSC_server] Creating & preprocessing Seurat object..")
          options(Seurat.object.assay.calcn = TRUE)
          getOption("Seurat.object.assay.calcn")
          counts <- as(counts, "dgCMatrix")
          msg <- "Creating & preprocessing Seurat object..."
          shiny::withProgress(message = msg, value = 0.9, {
            SO <- playbase::pgx.createSeuratObject(counts, samples,
              batch = NULL, filter = FALSE, preprocess = FALSE)
            SO <- playbase::seurat.preprocess(SO, sct = FALSE, tsne = FALSE, umap = FALSE)
          })
          dbg("[normalizationSC_server] Seurat object created & preprocessed.")
          kk <- setdiff(colnames(samples), colnames(SO@meta.data))
          if (length(kk) > 1) {
            SO@meta.data <- cbind(SO@meta.data, samples[, kk, drop = FALSE])
          }
          LL <- list(
            SO = SO,
            samples = SO@meta.data,
            pos.pca = pos.list[["pca"]],
            pos.tsne = pos.list[["tsne"]],
            pos.umap = pos.list[["umap"]]
          )
          rm(counts, nX, nX1, samples, SO)
          return(LL)
        } else {
          rm(counts, samples)
          return(NULL)
        }
        
      })
      
      ## ------------------------------------------------------------------
      ## Plot functions
      ## ------------------------------------------------------------------
      plot1 <- function() {
        shiny::req(dim(ds_norm_Counts()$SO))
        SO <- ds_norm_Counts()$SO
        meta <- SO@meta.data
        require(scplotter)
        require(ggplot2)
        library(ggpubr)
        require(vioplot)
        vars <- input$clusterBy        
        shiny::validate(shiny::need(
          !is.null(vars),
          "For QC, please select a QC variable from the options."
        ))        
        shiny::validate(shiny::need(
          length(vars) <= 4,
          "Please select up to 4 QC variables for visualization."
        ))        
        if ("seurat_clusters" %in% colnames(meta)) {
          meta$seurat_clusters <- as.character(meta$seurat_clusters) 
        }
        i <- 1
        class.vars <- c()
        for(i in 1:length(vars)) {
          class.vars <- c(class.vars, class(meta[, vars[i]]))
        }
        names(class.vars) <- vars
        num.vars <- vars[which(class.vars %in% c("numeric","integer"))]
        char.vars <- vars[which(class.vars %in% c("character"))]

        if (!input$groupby_celltype) {
          meta$IDENT0 <- "IDENT"
          plist <- list()
          i = 1
          for(i in 1:length(vars)) {
            v <- vars[i]
            if(v %in% num.vars) {
              if (all(range(meta[, v]) %in% c(0,0))) {
                meta[, v] <- runif(nrow(meta), min = 0, max = 1e-5)
              }
              pp <- ggplot(meta, aes_string(x = "IDENT0", y = v))
              pp <- pp + geom_point(col = "white", shape = "diamond") + geom_jitter(col = "dim gray")
              ylab <- v
              if (grepl("percent", v)) {
                pp <- pp + ylim(0, 100)
                ylab <- "Percentage"
              } else {
                if (min(meta[, v]) >= 0) {
                  pp <- pp + scale_y_continuous(limits = c(0, NA))
                }
              }

              ## make yintercept from reactive values
              if (v == "percent.mt") {
                pp <- pp + geom_hline(yintercept = input$mt_threshold, col = "firebrick2")
              }
              if (v == "percent.hb") {
                pp <- pp + geom_hline(yintercept = input$hb_threshold, col = "firebrick2")
              }
              if (v == "nFeature_RNA") {
                minv <- input$nfeature_threshold[1]
                maxv <- input$nfeature_threshold[2]
                pp <- pp + geom_hline(yintercept = c(minv, maxv), col = "firebrick2")
              }
              pp <- pp + scale_x_discrete(labels = "Cells")
              plist[[v]] <- pp + labs(title = v) + ylab(ylab)
            }
            if (v %in% char.vars) {
              tt <- data.frame(table(meta[, v]))
              pp <- ggplot(tt, aes(x = Var1, y = Freq))
              pp <- pp + geom_bar(stat = "identity", fill = "dim gray")
              plist[[v]] <- pp + labs(title = v) + ylab("Number of cells")
            }
          }
          plist <- lapply(plist, function(x) {
            x <- x + theme(axis.text.x = element_text(size = 13))
            x <- x + theme(axis.text.y = element_text(size = 13))
            x <- x + RotatedAxis() + xlab("")
            x <- x + theme(panel.border = element_rect(color = "black", fill = NA,
              size = 1, linewidth = 2))
          })
          i <- 1
          if (length(plist) == 1) {
            plist[[1]]
          } else if (length(plist) == 2) {
            ggpubr::ggarrange(plist[[1]], plist[[2]], nrow = 1, ncol = 2)
          } else if (length(plist) == 3) {
            ggpubr::ggarrange(plist[[1]], plist[[2]], plist[[3]], nrow = 1, ncol = 3)
          } else if (length(plist) == 4) {
            ggpubr::ggarrange(plist[[1]], plist[[2]], plist[[3]], plist[[4]], nrow = 2, ncol = 2)
          }
        } else {
          grp <- "celltype.azimuth"
          plist <- list()
          i = 1
          for(i in 1:length(vars)) {
            v <- vars[i]
            if (v %in% num.vars) {
              pp <- ggplot(meta, aes_string(y = v, x = grp))
              pp <- pp + geom_boxplot() + RotatedAxis() + xlab("")
              ylab <- v
              pp <- pp + theme(legend.position = "none")
              if (grepl("percent", v)) {
                pp <- pp + ylim(0, 100)
                ylab <- "Percentage"
              } else if (min(meta[, v]) >= 0) {
                pp <- pp + scale_y_continuous(limits = c(0, NA))
              }
              if(v == "percent.mt") {
                pp <- pp + geom_hline(yintercept = input$mt_threshold, col = "firebrick2")
              }
              if(v == "percent.hb") {
                pp <- pp + geom_hline(yintercept = input$hb_threshold, col = "firebrick2")
              }
              if (v == "nFeature_RNA") {
                minv <- input$nfeature_threshold[1]
                maxv <- input$nfeature_threshold[2]
                pp <- pp + geom_hline(yintercept = c(minv, maxv), col = "firebrick2")
              }
              plist[[v]] <- pp + labs(title = v) + ylab(ylab)
            }
            if (v %in% char.vars) {
              if (v == "celltype.azimuth") {
                tt <- data.frame(table(meta[, v]))
                pp <- ggplot(tt, aes(x = Var1, y = Freq))
                pp <- pp + geom_bar(stat = "identity", fill = "dim gray")
                plist[[v]] <- pp + ylab("Number of cells")
              } else {
                tt <- data.frame(t(table(meta[, v], meta[, grp])))
                colnames(tt)[1:2] <- c("celltype.azimuth", v)
                pp <- ggplot(tt,  aes_string(x = v, y = "Freq", fill = "celltype.azimuth"))
                pp <- pp + geom_bar(stat = "identity")
                if (i < length(vars)) { pp <- pp + theme(legend.position = "none") }
                plist[[v]] <- pp + labs(title = v) + ylab("Number of cells")
              }
            }
          }
          plist <- lapply(plist, function(x) {
            x <- x + theme(axis.text.x = element_text(size = 13))
            x <- x + theme(axis.text.y = element_text(size = 13))
            x <- x + RotatedAxis() + xlab("")
            x <- x + theme(panel.border = element_rect(color = "black", fill = NA,
              size = 1, linewidth = 2))
          })
          i <- 1
          if (length(plist) == 1) {
            plist[[1]]
          } else if (length(plist) == 2) {
            ggpubr::ggarrange(plist[[1]], plist[[2]], nrow = 1, ncol = 2)
          } else if (length(plist) == 3) {
            ggpubr::ggarrange(plist[[1]], plist[[2]], plist[[3]], nrow = 1, ncol = 3)
          } else if (length(plist) == 4) {
            ggpubr::ggarrange(plist[[1]], plist[[2]], plist[[3]], plist[[4]], nrow = 2, ncol = 2)
          }
        }
      }

      plot2 <- function() {
        shiny::req(
          dim(ds_norm_Counts()$samples),
          dim(ds_norm_Counts()$pos.pca),
          dim(ds_norm_Counts()$pos.tsne),
          dim(ds_norm_Counts()$pos.umap)
        )
        pos.list <-  list(
          pca = ds_norm_Counts()$pos.pca,
          tsne = ds_norm_Counts()$pos.tsne,
          umap = ds_norm_Counts()$pos.umap
        )
        samples <- ds_norm_Counts()$samples
        vars <- input$clusterBy
        shiny::validate(shiny::need(
          !is.null(vars),
          "For clustering, please select a metadata variable from the menu on the left."
        ))        
        m <- tolower(input$dimred_plottype)
        if (length(vars) <= 2) {
          par(mfrow = c(1, length(vars)))
        } else if (length(vars) > 2 & length(vars) <= 4) {
          par(mfrow = c(2, 2))
        } else if (length(vars) > 4 & length(vars) <= 6) {
          par(mfrow = c(2, 3))
        } else  if (length(vars) > 6 & length(vars) <= 8) {
          par(mfrow = c(3, 3))
        }
        i <- 1
        for(i in 1:length(vars)) {
          v <- samples[, vars[i]]
          ss <- c("nFeature_RNA", "nCount_RNA")
          ## if (vars[i] %in% ss) { v <- log2(v + 1) }
          playbase::pgx.scatterPlotXY.BASE(
            pos = pos.list[[m]], var = v, title = vars[i],
            xlab = "Dim1", ylab = "Dim2"
          )
        }
      }

      ## ------------------------------------------------------------------
      ## Plot modules
      ## ------------------------------------------------------------------
      PlotModuleServer(
        "plot1",
        plotlib = "base",
        func = plot1,
        res = c(75, 120),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )
      PlotModuleServer(
        "plot2",
        plotlib = "base",
        func = plot2,
        res = c(75, 120),
        pdf.width = 12,
        pdf.height = 6,
        add.watermark = FALSE
      )

      ## counts
      counts <- shiny::reactive({
        shiny::req(r_counts())
        counts <- r_counts()
        return(counts)
      })

      ## Normalized counts
      X <- shiny::reactive({
        shiny::req(r_counts())
        counts <- r_counts()
        X <- playbase::logCPM(as.matrix(counts), total = 1e4, prior = 1)
        return(X)
      })

      ## samples
      samples <- shiny::reactive({
        shiny::req(r_samples())
        samples <- r_samples()
        return(samples)
      })
      
      azimuth_ref <- shiny::reactive({
        ref <- input$ref_atlas
        return(ref)
      })

      nfeat_thr <- shiny::reactive({
        shiny::req(input$remove_cells)
        nfeat_thr <- input$nfeature_threshold
        return(nfeat_thr)
      })

      mt_thr <- shiny::reactive({
        shiny::req(input$remove_cells)
        mt_thr <- input$mt_threshold
        return(mt_thr)
      })

      hb_thr <- shiny::reactive({
        shiny::req(input$remove_cells)
        hb_thr <- input$hb_threshold
        return(hb_thr)
      })
      
      LL <- list(
        counts = counts,
        samples = samples,
        X = X,
        impX = shiny::reactive(NULL),
        azimuth_ref = shiny::reactive(azimuth_ref()),
        nfeature_threshold = shiny::reactive(nfeat_thr()),
        mt_threshold = shiny::reactive(mt_thr()),
        hb_threshold = shiny::reactive(hb_thr()),
        norm_method = shiny::reactive("CPM")
      )

      return(LL) ## pointing to reactive

    } ## end-of-server
  )
}
