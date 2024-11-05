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
        if (is.null(counts)) {
          return(NULL)
        }
        
        ## ref.tissue <- input$azmRef
        ## celltype.compute <- input$celltypeCompute
        cells.trs <- 5000
        ref.tissue <- "pbmcref"
        celltype.compute <- TRUE
        ncells <- ncol(counts)
        counts.sc <- samples.sc <- NULL
        
        if (celltype.compute) {

          dbg("[normalizationSC_server:normalizedCounts:] Inferring cell types with Azimuth!")
          
          if (ncells <= cells.trs) {
            
            shiny::withProgress(
              message = "Your dataset contains ", ncells, " cells. Inferring cell types with Azimuth.",
              value = 0.3, {
                azm <- playbase::pgx.runAzimuth(counts = counts, reference = ref.tissue)
              }
            )
            if(!is.null(azm)) { dbg("----------MNT1: OK") }
            celltype <- azm[,grep("^predicted.*l2$", colnames(azm))]
            samples$azm.celltype <- celltype
            nX <- as.matrix(playbase::logCPM(counts, 1, total = 1e5))
            return(list(counts = nX, samples = samples))

          } else {

            shiny::withProgress(
              message = "Your dataset contains > 5K cells. Computing metacells.",
              value = 0.3, {
                SC <- playbase::pgx.supercell(counts = counts, meta = samples)
                counts.sc <- SC$counts
                samples.sc <- SC$meta
                ## samples$sc.membership <- SC$membership
              }
            )
            shiny::withProgress(
              message = "Inferring cell types with Azimuth on metacells.",
              value = 0.3, {
                azm <- playbase::pgx.runAzimuth(counts = counts.sc, reference = ref.tissue)
                celltype <- azm[,grep("^predicted.*l2$", colnames(azm))]
                samples.sc$azm.celltype <- celltype
              }
            )
            nX <- as.matrix(playbase::logCPM(counts.sc, 1, total = 1e5))
            return(list(counts = nX, samples = samples.sc))

          }

        } else {

          dbg("[normalizationSC_server:normalizedCounts:] Cell type already pre-defined.")
          if (ncells > cells.trs) {
            shiny::withProgress(
              message = "Your dataset contains > 5K cells. Computing metacells.",
              value = 0.3, {
                SC <- playbase::pgx.supercell(counts = counts, meta = samples)
                counts.sc <- SC$counts
                samples.sc <- SC$meta
                ## samples$sc.membership <- SC$membership
              }
            )
          }
          nX <- as.matrix(playbase::logCPM(counts.sc, 1, total = 1e5))
          return(list(counts = nX, samples = samples.sc))

        }

      })

      ## ------------------------------------------------------------------
      ## Plot functions
      ## ------------------------------------------------------------------

      plot1 <- function() {
        shiny::req(dim(normalizedCounts()$counts), dim(normalizedCounts()$samples)) 
        counts <- normalizedCounts()$counts
        samples <- normalizedCounts()$samples
        kk <- intersect(rownames(samples), colnames(counts))
        samples <- samples[kk, , drop = FALSE]
        counts <- counts[, kk, drop = FALSE]
        Vars <- c("celltype","stim")
        for(i in 1:length(Vars)) {
          y <- samples[, Vars[i]]
          playbase::pgx.dimPlot(counts, y, method="umap")
        }
      }

      ## ------------------------------------------------------------------
      ## Plot UI
      ## ------------------------------------------------------------------
      output$normalization <- shiny::renderUI({

        dimred.infotext <- "T-distributed stochastic neighbor embedding (t-SNE)"
        
        dropout.options <- tagList(
          shiny::radioButtons(
            ns("dropout_plottype"),
            label = "Plot type:",
            choices = c("boxplot", "histogram"),
            selected = "histogram", inline = FALSE
          )
        )

        navmenu <- tagList(
          bslib::card(bslib::card_body(
            style = "padding: 0px;",
            bslib::accordion(
              multiple = FALSE,
              style = "background-color: #F7FAFD99;",
              bslib::accordion_panel(
                title = "1. Step 1",
                shiny::div(
                  style = "display: flex; align-items: center; justify-content: space-between;",
                  shiny::p("wip")
                ),
                shiny::checkboxInput(ns("zero_as_na"), label = "Treat zero as NA", value=FALSE),
                br()
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
                ## info.text = dimred.infotext,
                ## caption = dropout.infotext,
                ## options = dimred.options,
                height = c("auto", "100%"),
                show.maximize = FALSE),
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
