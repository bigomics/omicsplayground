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

upload_module_normalizationSC_server <- function(
                                                 id,
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

      normalizedCounts <- shiny::reactive({

        shiny::req(r_counts())
        counts <- r_counts()
        samples <- r_samples()
        if (is.null(counts)) { return(NULL) }
        kk <- intersect(colnames(counts), rownames(samples))
        counts <- counts[, kk, drop = FALSE]
        samples <- samples[kk, , drop = FALSE]
        
        cells_trs <- 1500
        ref_tissue <- input$ref_atlas
        celltype_compute <- TRUE
        if (celltype_compute) {
        ## if (input$infercelltypes) {

          ss <- c("celltype","cell_type","cell.type","CELL_TYPE")
          kk <- which(!colnames(samples) %in% ss)
          samples <- samples[, kk, drop = FALSE]

          dbg("[normalizationSC_server:normalizedCounts:] N.cells in dataset ", ncol(counts))

          if (ncol(counts) > cells_trs) {
            dbg("[normalizationSC_server:normalizedCounts:] Performing random sampling of 1000 cells.")
            kk <- sample(colnames(counts), 1000)
            counts1 <- counts[, kk, drop = FALSE] 
            samples1 <- samples[kk, , drop = FALSE]
          } else {
            counts1 <- counts
            samples1 <- samples
          }

          dbg("[normalizationSC_server:normalizedCounts:] Inferring cell types with Azimuth!")
          dbg("[normalizationSC_server:normalizedCounts:] Reference atlas:", ref_tissue)
          azm <- playbase::pgx.runAzimuth(counts = counts1, reference = ref_tissue)
          dbg("[normalizationSC_server:normalizedCounts:] Cell types inference completed.")
          celltype <- azm[,grep("^predicted.*l2$", colnames(azm))]
          samples1 <- cbind(samples1, celltype = celltype)
          nX <- playbase::logCPM(as.matrix(counts1), 1, total = 1e4)
          return(list(counts = nX, samples = samples1))
          
        } else {

          dbg("[normalizationSC_server:normalizedCounts:] Cell type already pre-defined.")
          nX <- playbase::logCPM(as.matrix(counts), 1, total = 1e4)
          return(list(counts = nX, samples = samples))

        }

      })

      ## ------------------------------------------------------------------
      ## Plot functions
      ## ------------------------------------------------------------------

      plot1 <- function() {
        shiny::req(dim(normalizedCounts()$counts), dim(normalizedCounts()$samples)) 
        counts <- normalizedCounts()$counts
        samples <- normalizedCounts()$samples
        cluster_vars <- input$clusterBy
        kk <- intersect(rownames(samples), colnames(counts))
        samples <- samples[kk, , drop = FALSE]
        counts <- counts[, kk, drop = FALSE]
        method <- tolower(input$dimred_plottype)
        par(mfrow = c(1,length(cluster_vars)))
        for(i in 1:length(cluster_vars)) {
          v <- cluster_vars[i]
          playbase::pgx.dimPlot(counts, samples[, v], method = method)
        }
      }

      plot2 <- function() {
        shiny::req(dim(normalizedCounts()$counts), dim(normalizedCounts()$samples)) 
        counts <- normalizedCounts()$counts
        samples <- normalizedCounts()$samples
        cluster_vars <- input$clusterBy
        kk <- intersect(rownames(samples), colnames(counts))
        samples <- samples[kk, , drop = FALSE]
        counts <- counts[, kk, drop = FALSE]
        SO <- playbase::pgx.justSeuratObject(counts, samples)
        ## SO <- CreateSe
        SO <- Seurat::NormalizeData(SO)
        SO <- Seurat::FindVariableFeatures(SO)
        SO <- Seurat::ScaleData(SO)
        hvf <- Seurat::VariableFeatures(SO)[1:10]
        if(!is.null(SO)) dbg("----MNT1: ", paste0(colnames(SO@meta.data), collapse=", "))
        if(!is.null(SO)) dbg("----MNT2: ", paste0(colnames(samples), collapse=", "))
        plist <- list()
        for(i in 1:length(cluster_vars)) {
          v <- as.character(cluster_vars[i])
          if (! v %in% colnames(SO@meta.data)) {
            SO@meta.data <- cbind(SO@meta.data, samples[, v])
            colnames(SO@meta.data)[ncol(SO@meta.data)] <- v
          }
          ug <- length(unique(samples[, v]))
          scaling <- ifelse(ug < 3, FALSE, TRUE)
          pl <- Seurat::DotPlot(SO, features = hvf, group.by = v, scale = scaling)
          pl <- pl + RotatedAxis() + ylab("") + xlab("")
          plist[[i]] <- pl
        }
        cowplot::plot_grid(plotlist = plist)
      }


      ## ------------------------------------------------------------------
      ## Plot UI
      ## ------------------------------------------------------------------
      output$normalization <- shiny::renderUI({

        shiny::req(dim(normalizedCounts()$samples))
        samples <- normalizedCounts()$samples
        metadata_vars <- colnames(samples)
                
        dimred.infotext <- "Dimensionality reduction enables to simplify large datasets by computing representative data points capable of preserving the biological information while reducing the dimensionality of the data. Here we employ the two most widely used non-linear methods for dimensional reduction of single-cell RNA-seq data: T-distributed stochastic neighbor embedding (t-SNE), and Unifold Manifold Approximation and Projection (UMAP). https://omicsplayground.readthedocs.io/en/latest/methods/#clustering"

        dotplot.infotext <- "Dot plot of top 10 highly variable features in the data. The dot plot allows to visualize the expression changes across different user-defined identity classes (e.g., cell clusters, cell types, phenotype classes). The size of each dot reflects the percentage of cells within a class; the color of each dot indicates the average expression level (computed using Seurat's AverageExpression function) across all cells within a class (blue is high)."
        
        dimred.options <- tagList(
          shiny::radioButtons(
            ns("dimred_plottype"),
            label = "Plot type:",
            choices = c("tSNE", "UMAP"),
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
                  value = TRUE),
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
                title = "Highly variable features",
                info.text = dotplot.infotext,
                caption = dotplot.infotext,
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

      ## log2 counts
      cX <- shiny::reactive({
        shiny::req(dim(normalizedCounts()$counts))
        log2(normalizedCounts()$counts + 1)
      })

      return(
        list(
          counts = shiny::reactive(normalizedCounts()$counts),
          samples = shiny::reactive(normalizedCounts()$samples),
          X = cX,
          impX = shiny::reactive(NULL),
          norm_method = shiny::reactive("CPM")
        )
      ) ## pointing to reactive
    } ## end-of-server
  )
}
