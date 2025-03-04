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

        SC <- playbase::pgx.supercell(
          counts = counts,
          meta = samples,
          gamma = 10,
          group = NULL
        )

        nX <- playbase::logCPM(SC[["counts"]], 1, total = 1e5, log = FALSE)
        list(counts = nX, samples = SC[["meta"]])

      })

      ## ------------------------------------------------------------------
      ## Plot functions
      ## ------------------------------------------------------------------

      plot1 <- function() {
        ## counts <- r_counts()
        ## samples <- r_samples()
        counts <- normalizedCounts()$counts
        samples <- normalizedCounts()$samples
        kk <- intersect(rownames(samples), colnames(counts))
        samples <- samples[kk, , drop = FALSE]
        counts <- counts[, kk, drop = FALSE]
        ## plot(1:100, log(1:100), cex=1)
        pgx.dimPlot(log2(counts+1), samples[,"stim"])
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
