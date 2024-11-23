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
          "seurat_clusters",
          "G2M.Score", "S.Score"
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

        ## cellqc.options <- tagList(
        ##   shiny::radioButtons(
        ##     ns("sample_var"),
        ##     label = "",
        ##     choices = c(
        ##       "nFeature_RNA", "nCount_RNA",
        ##       "percent.mt", "percent.ribo", "percent.hb",
        ##       "G2M.Score", "S.Score"
        ##     ),
        ##     selected = "nFeature_RNA",
        ##     inline = FALSE
        ##   )
        ## )

        cellqc.options <- tagList(
          shiny::checkboxGroupInput(
            ns("qc_var"),
            label = "Select variable for QC:",
            choices = c(
              "nFeature_RNA", "nCount_RNA",
              "percent.mt", "percent.ribo", "percent.hb",
              "G2M.Score", "S.Score"
            ),
            selected = c("nFeature_RNA", "nCount_RNA"),
            inline = FALSE
          )
        )

        navmenu <- tagList(
          
          bslib::card(bslib::card_body(
            style = "padding: 0px;",
            bslib::accordion(
              multiple = FALSE,
              style = "background-color: #F7FAFD99;",

              bslib::accordion_panel(
                title = "Cell type inference",
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
                br()
              ),

              bslib::accordion_panel(
                title = "Visualize cell clusters",
                shiny::div(
                  style = "display: flex; align-items: center; justify-content: space-between;"
                ),
                shiny::selectInput(
                  ns("clusterBy"),
                  label = "Visualize cell cluster by",
                  choices = metadata_vars, ## reactive
                  multiple = TRUE,
                  selected = c("celltype.azimuth", "stim")
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
                title = "Samples' QC by cell type",
                info.text = cellqc.infotext,
                caption = cellqc.infotext,
                options = cellqc.options,
                height = c("auto", "100%"),
                show.maximize = FALSE,
                ),
              PlotModuleUI(
                ns("plot2"),
                title = "Dimensional reduction",
                info.text = dimred.infotext,
                caption = dimred.infotext,
                options = dimred.options,
                height = c("auto", "100%"),
                show.maximize = FALSE
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

        options(future.globals.maxSize= 4*1024^4) 

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
        cells_trs <- 1500
        if (ncells > cells_trs) {
          dbg("[normalizationSC_server:ds_norm_Counts:] N.cells in dataset:", ncells)
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
            celltype.azimuth <- azm[, grep("^predicted.*l*2$", colnames(azm))]
            samples <- cbind(samples, celltype.azimuth = celltype.azimuth)          
          } else if (is.vector(azm)) {
            dbg("[normalizationSC_server:ds_norm_Counts:] Your selected Azimuth reference atlas *might* be incorrect. Please double check.")
            samples <- cbind(samples, celltype.azimuth = azm)
          }
          nX <- playbase::logCPM(as.matrix(counts), 1, total = 1e4)
          return(list(counts = nX, samples = samples, ref_tissue = ref_tissue))
          rm(counts, samples)
          return(NULL)
        } else {
          rm(counts, samples)
          return(NULL)
        }
        
      })

      ## Dim reductions
      dimred_norm_Counts <- shiny::reactive({ 

        shiny::req(dim(ds_norm_Counts()$counts))
        shiny::req(dim(ds_norm_Counts()$samples))
        shiny::req(input$infercelltypes)
        counts <- ds_norm_Counts()$counts
        counts <- as.matrix(counts)
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
        require(scplotter)
        require(ggplot2)
        vars <- input$qc_var
        shiny::validate(shiny::need(
          !is.null(vars),
          "For QC, please select a QC variable from the options."
        ))        
        shiny::validate(shiny::need(
          length(vars)<=4,
          "Please select up to 4 QC variables for visualization."
        ))        
        pp <- scplotter::FeatureStatPlot(
          SO, features = vars, ident = "celltype.azimuth",
          facet_scales = "free_y",  drop = FALSE,
          theme_args = list(base_size = 15))
        size <- ifelse(length(vars)==1, 10, 8)
        pp + ggplot2::theme(axis.text.x = element_text(size = size))
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
          if (vars[i] %in% c("nFeature_RNA", "nCount_RNA")) {
            v <- log2(v + 1)
          }
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
        counts
      })

      ## normalized counts
      X <- shiny::reactive({
        shiny::req(r_counts())
        counts <- r_counts()
        X <- playbase::logCPM(as.matrix(counts), total = 1e4, prior = 1)
        X
      })

      ## smples
      samples <- shiny::reactive({
        shiny::req(r_samples())
        samples <- r_samples()
        samples
      })

      ## Azimuth ref
      ## azimuth_ref <- shiny::reactive({
      ##  shiny::req(ds_norm_Counts())
      ##  ref <- as.character(ds_norm_Counts()$ref_tissue())
      ##  dbg("-----MNT1:", ref)
      ##  return(ref)
      ## })
      ## observe({ if(!is.null(samples)) { dbg("MNT2----", azimuth_ref()) } })

      return(
        list(
          counts = counts,
          samples = samples, ## shiny::reactive(ds_norm_Counts()$samples),
          X = X,
          impX = shiny::reactive(NULL),
          ## azimuth_ref = azimuth_ref, ## NEW AZ # TODO
          norm_method = shiny::reactive("CPM")
        )
      ) ## pointing to reactive
    } ## end-of-server
  )
}
