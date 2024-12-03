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

        shiny::req(r_samples())
        samples <- r_samples()
        metadata_vars <- c(
          "celltype.azimuth", "orig.ident",
          "nCount_RNA", "nFeature_RNA",
          "percent.mt", "percent.ribo",
          "percent.hb", "seurat_clusters",
          "G2M.Score", "S.Score", "Phase"
        )
        metadata_vars <- unique(c(metadata_vars, colnames(samples)))
        
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
                title = HTML("<span style='font-size: 0.9em;'> Cell type inference</span>"),
                shiny::div(
                  style = "display: flex; align-items: center; justify-content: space-between;"
                ),
                shiny::checkboxInput(
                  ns("infercelltypes"),
                  label = "Infer cell types with Azimuth",
                  value = FALSE
                ),
                shiny::conditionalPanel(
                  "input.infercelltypes == true",
                  ns = ns,
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
                  )
                ),
                shiny::br()
              ),

              bslib::accordion_panel(
                title = HTML("<span style='font-size: 0.9em;'> Visualize cell clusters</span>"),
                shiny::div(
                  style = "display: flex; align-items: center; justify-content: space-between;"
                ),
                shiny::selectInput(
                  ns("clusterBy"),
                  label = "Visualize cell cluster by",
                  choices = metadata_vars,
                  multiple = TRUE,
                  selected = c("nCount_RNA", "nFeature_RNA", "percent.mt", "celltype.azimuth")
                ),
                shiny::br()
              ),

              ## bslib::accordion_panel(
              ##   title = HTML("<span style='font-size: 0.9em;'> Phenotype of interest</span>"),
              ##   shiny::div(
              ##     style = "display: flex; align-items: center; justify-content: space-between;"
              ##   ),
              ##   shiny::selectInput(
              ##     ns("pheno"),
              ##     label = "Phenotype of interest",
              ##     choices = colnames(samples), ## reactive
              ##     selected = "<select>"
              ##   ),
              ##   shiny::br()
              ## ),

              bslib::accordion_panel(
                title = HTML("<span style='font-size: 0.9em;'> Cell filtering</span>"),
                shiny::p("Remove cells based on QC variables"),
                shiny::checkboxInput(
                  ns("remove_cells"),
                  label = "Remove cells",
                  value = FALSE
                ),
                shiny::conditionalPanel(
                  "input.remove_cells == true",
                  ns = ns,
                  shiny::sliderInput(ns("nfeature_threshold"), "Detected genes per cell:", 0, 25000, 2000, 100),
                  shiny::sliderInput(ns("mt_threshold"), "MT expresssion threshold (%):", 1, 100, 10, 0),
                  shiny::sliderInput(ns("hb_threshold"), "HB expression threshold (%):", 1, 100, 10, 0)
                ),
                shiny::br()
              )
              ##---------------------------
              ## bslib::accordion_panel(
              ##   title = HTML("<span style='font-size: 0.9em;'> Filter #1: Mitochondrial expression rate</span>"),
              ##   shiny::p("Remove cells with high MT expression"),
              ##   shiny::checkboxInput(
              ##     ns("remove_cells_1"),
              ##     label = "Remove cells",
              ##     value = FALSE
              ##   ),
              ##   shiny::conditionalPanel(
              ##     "input.remove_cells_1 == true",
              ##     ns = ns,
              ##     shiny::sliderInput(
              ##       ns("mt_threshold"), "Select threshold (%):", 1, 100, 10, 0
              ##     )
              ##   ),
              ##   shiny::br()
              ## ),
              ##---------------------------
              ## bslib::accordion_panel(
              ##   title = HTML("<span style='font-size: 0.9em;'> Filter #2: Globins' expression rate</span>"),
              ##   shiny::p("Remove cells with high globins' expression"),
              ##   shiny::checkboxInput(
              ##     ns("remove_cells_2"),
              ##     label = "Remove cells",
              ##     value = FALSE
              ##   ),
              ##   shiny::conditionalPanel(
              ##     "input.remove_cells_2 == true",
              ##     ns = ns,
              ##     shiny::sliderInput(
              ##       ns("hb_threshold"), "Select threshold (%):", 1, 100, 10, 0
              ##     )
              ##   ),
              ##   shiny::br()
              ## ),
              ##---------------------------
              ## bslib::accordion_panel(
              ##   title = HTML("<span style='font-size: 0.9em;'> Filter #3: Number of detected genes per cell</span>"),
              ##   shiny::p("Remove cells with too few or too many detected genes"),
              ##   shiny::checkboxInput(
              ##     ns("remove_cells_3"),
              ##     label = "Remove cells",
              ##     value = FALSE
              ##   ),
              ##   shiny::conditionalPanel(
              ##     "input.remove_cells_3 == true",
              ##     ns = ns,
              ##     shiny::sliderInput(
              ##       ns("nFeatureRNA_threshold"), "Select threshold:", 50, 500, 200, 100
              ##     )
              ##   ),
              ##   shiny::br()
              ## )
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

      ## downsampled and normalized data
      ## downsample or still use full dataset if ncells < threshold
      ds_norm_Counts <- shiny::reactive({ 

        options(future.globals.maxSize= 4*1024^100)

        shiny::req(input$infercelltypes)
        shiny::req(r_counts())
        shiny::req(r_samples())
        counts <- r_counts()
        samples <- r_samples()
        if (is.null(counts)) { return(NULL) }

        kk <- intersect(colnames(counts), rownames(samples))
        counts <- counts[, kk, drop = FALSE]
        samples <- samples[kk, , drop = FALSE]

        ncells <- ncol(counts)
        cells_trs <- 9000
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

          nX <- playbase::logCPM(as.matrix(counts), 1, total = 1e4)
          return(list(counts = nX, samples = samples))
          rm(counts, samples)
          return(NULL)
        } else {
          rm(counts, samples)
          return(NULL)
        }
        
      })

      ##------------PUT IN SAME REACTIVE OBJECT: IK
      ## Dim reductions
      dimred_norm_Counts <- shiny::reactive({ 

        shiny::req(dim(ds_norm_Counts()$counts))
        shiny::req(dim(ds_norm_Counts()$samples))
        shiny::req(input$infercelltypes)        
        counts <- ds_norm_Counts()$counts
        ## counts <- as.matrix(counts)
        samples <- ds_norm_Counts()$samples

        jj <- head(order(-matrixStats::rowSds(counts, na.rm = TRUE)), 250)
        counts1 <- counts[jj, , drop = FALSE]
        counts1 <- counts1 - rowMeans(counts1, na.rm = TRUE)
        nb <- ceiling(min(15, dim(counts) / 8))
        pos.list <- list()
        dbg("[normalizationSC_server:dimred_norm_Counts:] Performing PCA, tSNE and UMAP...")
        msg <- "Performing PCA, tSNE, UMAP..."
        shiny::withProgress(message = msg, value = 0.5, {
            pos.list[["pca"]] <- irlba::irlba(counts1, nv = 2, nu = 0)$v
            pos.list[["tsne"]] <- Rtsne::Rtsne(t(counts1), perplexity = 2*nb, check_duplicates=F)$Y
            pos.list[["umap"]] <- uwot::umap(t(counts1), n_neighbors = max(2, nb))
            pos.list <- lapply(pos.list, function(x) {
              rownames(x) <- colnames(counts1);
              return(x)
            })
          })
        dbg("[normalizationSC_server:dimred_norm_Counts:] PCA, tSNE & UMAP completed.")
        
        options(Seurat.object.assay.calcn = TRUE)
        getOption("Seurat.object.assay.calcn")
        dbg("[normalizationSC_server:dimred_norm_Counts:] Creating & preprocessing Seurat object..")
        options(future.globals.maxSize= 4*1024^5)
        counts <- as(counts, "dgCMatrix")
        msg <- "Creating & preprocessing Seurat object..."
        shiny::withProgress(message = msg, value = 0.9, {
          SO <- playbase::pgx.createSeuratObject(counts, samples,
            batch = NULL, filter = FALSE, preprocess = FALSE)
          SO <- playbase::seurat.preprocess(SO, sct = FALSE, tsne = FALSE, umap = FALSE)
        })
        dbg("[normalizationSC_server:dimred_norm_Counts:] Seurat object created & preprocessed.")
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
        rm(counts, samples, SO)
        return(LL)

      })

      ## ------------------------------------------------------------------
      ## Plot functions
      ## ------------------------------------------------------------------
      plot1 <- function() {
        shiny::req(dim(dimred_norm_Counts()$SO))
        SO <- dimred_norm_Counts()$SO
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
              ## yintercept make from reactive values
              if(v %in% c("percent.hb", "percent.mt")) {
                pp <- pp + geom_hline(yintercept = 10, col = "firebrick2")
              }
              if(v == "nFeature_RNA") {
                pp <- pp + geom_hline(yintercept = c(100, 2000), col = "firebrick2")
              }
              pp <- pp + scale_x_discrete(labels = "Cells")
              plist[[v]] <- pp + labs(title = v) + ylab(ylab)
            }
            if(v %in% char.vars) {
              tt <- data.frame(table(meta[, v]))
              pp <- ggplot(tt, aes(x = Var1, y = Freq))
              pp <- pp + geom_bar(stat = "identity", fill = "dim gray")
              plist[[v]] <- pp + labs(title = v) + ylab("Number of cells")
            }
          }
          plist <- lapply(plist, function(x) {
            x <- x + theme(axis.text.x = element_text(size = 13));
            x <- x + theme(axis.text.y = element_text(size = 13));
            x <- x + RotatedAxis() + xlab("")
            x <- x + theme(panel.border = element_rect(color = "black", fill = NA, size = 1));
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
            if(v %in% num.vars) {
              pp <- ggplot(meta, aes_string(y = v, x = grp, fill = grp))
              pp <- pp + geom_boxplot() + RotatedAxis() + xlab("")
              ylab <- v
              ## if (i < length(vars)) { pp <- pp + theme(legend.position = "none") }
              pp <- pp + theme(legend.position = "none")
              if (grepl("percent", v)) {
                pp <- pp + ylim(0, 100)
                ylab <- "Percentage"
              } else if (min(meta[, v]) >= 0) {
                pp <- pp + scale_y_continuous(limits = c(0, NA))
              }
              if (v %in% c("percent.hb","percent.mt")) {
                pp <- pp + geom_hline(yintercept = 10, col = "firebrick2")
              }
              if (v == "nFeature_RNA") {
                pp <- pp + geom_hline(yintercept = c(100, 2000), col = "firebrick2")
              }
              ## pp <- pp + scale_x_discrete(labels = "Cells")
              plist[[v]] <- pp + labs(title = v) + ylab(ylab)
            }
            if (v %in% char.vars) {
              if (v == "celltype.azimuth") {
                tt <- data.frame(table(meta[, v]))
                pp <- ggplot(tt, aes(x = Var1, y = Freq))
                pp <- pp + geom_bar(stat = "identity", fill = "dim gray")
                pp <- pp + ylab("Number of cells")
                plist[[v]] <- pp
              } else {
                tt <- data.frame(t(table(meta[, v], meta[, grp])))
                colnames(tt)[1:2] <- c("celltype.azimuth", v)
                pp <- ggplot(tt,  aes_string(x = v, y = "Freq", fill = "celltype.azimuth"))
                pp <- pp + geom_bar(stat = "identity", position = position_dodge())
                pp <- pp + ylab("Number of cells")
                if (i < length(vars)) { pp <- pp + theme(legend.position = "none") }
                plist[[v]] <- pp
              }
            }
          }
          plist <- lapply(plist, function(x) {
            x <- x + theme(axis.text.x = element_text(size = 13));
            x <- x + theme(axis.text.y = element_text(size = 13));
            x <- x + RotatedAxis() + xlab("")
            x <- x + theme(panel.border = element_rect(color = "black", fill = NA, size = 1));
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
          dim(dimred_norm_Counts()$samples),
          dim(dimred_norm_Counts()$pos.pca),
          dim(dimred_norm_Counts()$pos.tsne),
          dim(dimred_norm_Counts()$pos.umap)
        )
        pos.list <-  list(
          pca = dimred_norm_Counts()$pos.pca,
          tsne = dimred_norm_Counts()$pos.tsne,
          umap = dimred_norm_Counts()$pos.umap
        )
        samples <- dimred_norm_Counts()$samples
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

      ## CPM normalized counts
      X <- shiny::reactive({
        shiny::req(r_counts())
        counts <- r_counts()
        X <- playbase::logCPM(as.matrix(counts), total = 1e4, prior = 1)
        return(X)
      })

      ## smples
      samples <- shiny::reactive({
        shiny::req(r_samples())
        samples <- r_samples()
        return(samples)
      })
      
      azimuth_ref <- shiny::reactive({
        shiny::req(input$infercelltypes)
        if (input$ref_atlas != "<select>") {
          ref <- input$ref_atlas
          return(ref)
        } else {
          return(NULL)
        }
      })
      
      ## sc_pheno <- shiny::reactive({
      ##   sc_pheno <- input$pheno
      ##  return(sc_pheno)
      ## })

      LL <- list(
        counts = counts,
        samples = samples,
        X = X,
        impX = shiny::reactive(NULL),
        azimuth_ref = shiny::reactive(azimuth_ref()),
        ## sc_pheno = shiny::reactive(sc_pheno()),
        ## nfeature_threshold = shiny::reactive(nfeature_threshold()),
        ## mt_threshold = shiny::reactive(mt_threshold()),
        ## hb_threshold = shiny::reactive(hb_threshold()),
        norm_method = shiny::reactive("CPM")
      )
      
      return(LL) ## pointing to reactive

    } ## end-of-server
  )
}

## ADD WILCOXON TEST
## ADD LOW AND HIGH THRESHOLD FOR nFeature_RNA
## REMOVE ALL THOSE GSET METHODS. USE Fisher, fgsea, spearman.
## FIX REMAINING ISSUES AT UPLOAD BOARD.
## TEST ON MORE DATASETS.
## ADD COMPUTATION OPTION. WE DO METACELL IF NCOL(X)>20K.
## IDENTIFY SLOW MODULES.

