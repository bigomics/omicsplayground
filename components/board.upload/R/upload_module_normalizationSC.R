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
      ## Object reactive chain
      ## ------------------------------------------------------------------

      ## downsampled and normalized data
      ## downsample or still use full dataset if ncells < threshold
      ds_norm_Counts <- shiny::reactive({ 

        shiny::req(r_counts())
        shiny::req(r_samples())
        counts <- r_counts()
        samples <- r_samples()
        if (is.null(counts)) { return(NULL) }

        kk <- intersect(colnames(counts), rownames(samples))
        counts <- counts[, kk, drop = FALSE]
        samples <- samples[kk, , drop = FALSE]
        
        ncells <- ncol(counts)
        cells_trs <- 2000
        if (ncells > cells_trs) {
          dbg("[normalizationSC_server:ds_norm_Counts:] N.cells in dataset ", ncells)
          dbg("[normalizationSC_server:ds_norm_Counts:] Performing random sampling of 1K cells.")
          kk <- sample(colnames(counts), 1000)
          counts <- counts[, kk, drop = FALSE] 
          samples <- samples[kk, , drop = FALSE]
        }

        celltype_compute <- TRUE
        ref_tissue <- input$ref_atlas
        if (!is.null(ref_tissue) && celltype_compute) {
          ## if (input$infercelltypes) {
          dbg("[normalizationSC_server:ds_norm_Counts:] Inferring cell types with Azimuth!")
          dbg("[normalizationSC_server:ds_norm_Counts:] Reference atlas:", ref_tissue)
          shiny::withProgress(
            message = "Inferring cell types with Azimuth...",
            value = 0.4,
            {
              azm <- playbase::pgx.runAzimuth(counts = counts, reference = ref_tissue)
              dbg("[normalizationSC_server:ds_norm_Counts:] Cell types inference completed.")
              celltype.azm <- azm[,grep("^predicted.*l2$", colnames(azm))]
              samples <- cbind(samples, celltype.azm = celltype.azm)
            }
          )
        } else {
          dbg("[normalizationSC_server:ds_norm_Counts:] Reference missing or cell types presents.")
        }

        nX <- playbase::logCPM(as.matrix(counts), 1, total = 1e4)
        return(list(counts = nX, samples = samples))

      })

      ## Dim reductions: top 1000 features
      dimred_norm_Counts <- shiny::reactive({ 
        
        shiny::req(dim(ds_norm_Counts()$counts), dim(ds_norm_Counts()$samples)) 
        counts <- ds_norm_Counts()$counts
        samples <- ds_norm_Counts()$samples

        counts <- as.matrix(counts)
        jj <- head(order(-matrixStats::rowSds(counts, na.rm = TRUE)), 1000)
        counts1 <- counts[jj, , drop = FALSE]
        counts1 <- counts1 - rowMeans(counts1, na.rm = TRUE)
        nb <- ceiling(min(15, dim(counts) / 8))
        pos.list <- list()
        shiny::withProgress(
          message = "Performing PCA, tSNE, UMAP...",
          value = 0.4,
          {
            pos.list[["pca"]] <- irlba::irlba(counts1, nv = 2, nu = 0)$v
            pos.list[["tsne"]] <- Rtsne::Rtsne(t(counts1), perplexity = 2*nb, check_duplicates=F)$Y
            pos.list[["umap"]] <- uwot::umap(t(counts1), n_neighbors = max(2, nb))
            pos.list <- lapply(pos.list, function(x) { rownames(x)=colnames(counts1); return(x) })
          }
        )
        LL <- list(
          samples = samples,
          pos.pca = pos.list[["pca"]],
          pos.tsne = pos.list[["tsne"]],
          pos.umap = pos.list[["umap"]]
        )

        return(LL)

      })

      ## ------------------------------------------------------------------
      ## Plot functions
      ## ------------------------------------------------------------------

      plot1 <- function() {
        shiny::req(
          dim(dimred_norm_Counts()$samples),
          dim(dimred_norm_Counts()$pos.pca),
          dim(dimred_norm_Counts()$pos.tsne),
          dim(dimred_norm_Counts()$pos.umap)
        ) 
        samples <- dimred_norm_Counts()$samples
        pos.list <-  list(
          pca = dimred_norm_Counts()$pos.pca,
          tsne = dimred_norm_Counts()$pos.tsne,
          umap = dimred_norm_Counts()$pos.umap
        )
        vars <- input$clusterBy
        m <- tolower(input$dimred_plottype)
        if (length(vars) <= 2) {
          par(mfrow = c(1, length(vars)))
        } else if (length(vars) > 2 & length(vars) <= 4) {
          par(mfrow = c(2, 2))
        } else if (length(vars) > 4) {
          par(mfrow = c(2, length(vars)))
        }
        i <- 1
        for(i in 1:length(vars)) {
          v <- samples[, vars[i]]
          playbase::pgx.scatterPlotXY.BASE(
            pos = pos.list[[m]], var = v, title = m,
            xlab = "Dim1", ylab = "Dim2"
          )
        }
      }
      
      plot2 <- function() {
        shiny::req(r_counts())
        counts <- r_counts()
        shiny::req(dim(ds_norm_Counts()$samples)) 
        samples <- ds_norm_Counts()$samples
        kk <- intersect(colnames(counts), rownames(samples))
        counts <- counts[, kk, drop = FALSE]
        samples <- samples[kk, , drop = FALSE]
        options(Seurat.object.assay.calcn = TRUE)
        getOption("Seurat.object.assay.calcn")
        shiny::withProgress(
          message = "Creating and pre-processing Seurat Object...",
          value = 0.4,
          {
            SO <- Seurat::CreateSeuratObject(counts = counts, meta.data = samples)
            SO <- Seurat::PercentageFeatureSet(SO, pattern = "^MT-|^Mt-", col.name = "percent.mt")
            SO <- Seurat::PercentageFeatureSet(SO, pattern = "^HB|^Hb", col.name = "percent.hb")
            SO <- playbase::seurat.preprocess(SO, sct = FALSE, tsne = FALSE, umap = FALSE)
          }
        )
        kk <- setdiff(colnames(samples), colnames(SO@meta.data))
        if (length(kk) > 1) {
          SO@meta.data <- cbind(SO@meta.data, samples[, kk])
        }
        require(scplotter)
        scplotter::FeatureStatPlot(SO,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"),
          ident = "celltype.azm", facet_scales = "free_y",
          theme_args = list(base_size = 18))
      }

      ## ------------------------------------------------------------------
      ## Plot UI
      ## ------------------------------------------------------------------
      output$normalization <- shiny::renderUI({

        shiny::req(dim(ds_norm_Counts()$samples))
        samples <- ds_norm_Counts()$samples
        metadata_vars <- colnames(samples)
                
        dimred.infotext <- "Dimensionality reduction enables to simplify large datasets by computing representative data points capable of preserving the biological information while reducing the dimensionality of the data. Here we employ the two most widely used non-linear methods for dimensional reduction of single-cell RNA-seq data: T-distributed stochastic neighbor embedding (t-SNE), and Unifold Manifold Approximation and Projection (UMAP). https://omicsplayground.readthedocs.io/en/latest/methods/#clustering"

        dotplot.infotext <- "Dot plot of top 10 highly variable features in the data. The dot plot allows to visualize the expression changes across different user-defined identity classes (e.g., cell clusters, cell types, phenotype classes). The size of each dot reflects the percentage of cells within a class; the color of each dot indicates the average expression level (computed using Seurat's AverageExpression function) across all cells within a class (blue is high)."
        
        dimred.options <- tagList(
          shiny::radioButtons(
            ns("dimred_plottype"),
            label = "Plot type:",
            choices = c("PCA", "tSNE", "UMAP"),
            selected = "tSNE", inline = FALSE
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
                  value = TRUE
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
                      "pbmcref", "tonsilref"
                    ),
                    selected = "pbmcref"
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
                  selected = "celltype"
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
                title = "Dimensional reduction",
                info.text = dimred.infotext,
                caption = dimred.infotext,
                options = dimred.options,
                height = c("auto", "100%"),
                show.maximize = FALSE
              ),
              PlotModuleUI(
                ns("plot2"),
                title = "Sequencing depth & samples' QC by cell type",
                ## info.text = dotplot.infotext, ## update
                ## caption = dotplot.infotext,
                ## options = dotplot.options,
                height = c("auto", "100%"),
                show.maximize = FALSE,
                )
            )
          ),
          div(shiny::checkboxInput(ns("normalizationUI"), NULL, TRUE),
            style = "visibility:hidden")
        )            
        return(ui)
      })
      
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
        X <- playbase::logCPM(as.matrix(counts), 1, total = 1e4)
        X
      })

      ## smples
      samples <- shiny::reactive({
        shiny::req(r_samples())
        samples <- r_samples()
        samples
      })
      
      return(
        list(
          counts = counts,
          samples = samples, ## shiny::reactive(ds_norm_Counts()$samples),
          X = X,
          impX = shiny::reactive(NULL),
          norm_method = shiny::reactive("CPM")
        )
      ) ## pointing to reactive
    } ## end-of-server
  )
}
