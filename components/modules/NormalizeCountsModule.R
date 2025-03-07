##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


## =====================================================================================
## ===================== NORMALIZE GADGET UI ===========================================
## =====================================================================================

NormalizeCountsGadget <- function(X, pheno, height = 720) {
  gadgetize(NormalizeCountsUI, NormalizeCountsServer,
    title = "NormalizeCounts",
    X = X, pheno = pheno, height = height
  )
}

NormalizeCountsUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("UI"))
}

NormalizeCountsServerRT <- function(id, counts, height = 720) {
  ##
  ## counts:  reactive input object
  ##
  shiny::moduleServer(
    id,
    function(input, output, session) {
      all.methods <- c("none", "scale", "quantile", "CPM", "TMM", "RLE")
      nc_info <- "normalization module"
      nc_info <- ""

      output$UI <- shiny::renderUI({
        ns <- session$ns
        shiny::sidebarLayout(
          shiny::sidebarPanel(
            shiny::helpText(nc_info),
            withTooltip(
              shiny::selectInput(ns("selectmethod"), "Select normalization:",
                choices = all.methods, selected = "CPM"
              ),
              "Select initial normalization method.",
              placement = "top", options = list(container = "body")
            ),
            withTooltip(
              shiny::checkboxInput(ns("postqn"), "Post quantile normalization"),
              "Apply additional quantile normalization after scaling method.",
              placement = "top", options = list(container = "body")
            ),
            withTooltip(
              shiny::checkboxInput(ns("addnoise"), "Simulate unnormalized"),
              "Simulated unnormalized data by adding random scaling to raw data",
              placement = "top", options = list(container = "body")
            ),
            width = 2
          ),
          shiny::mainPanel(
            shiny::plotOutput(ns("canvas"), width = "100%", height = height) %>% shinycssloaders::withSpinner(),
            width = 10
          )
        )
      })
      shiny::outputOptions(output, "UI", suspendWhenHidden = FALSE) ## important!!!

      pgx <- shiny::reactive({
        if (is.null(input$addnoise)) {
          return(NULL)
        }
        pgx <- list(counts = counts())

        ## Set missing values to zero !!!!!!!!!!!!!!!
        is.missing <- is.na(pgx$counts) | is.infinite(pgx$counts)
        pgx$counts[is.missing] <- 0

        ## Simulate raw signals
        if (input$addnoise) {
          mult <- 1e4 * runif(ncol(pgx$counts))
          pgx$counts <- t(t(pgx$counts) * mult)
        }
        pgx
      })

      output$canvas <- shiny::renderPlot({
        shiny::req(counts())
        ## Show all methods
        postqn <- input$postqn
        viz.NormalizeCounts(
          pgx(),
          methods = all.methods,
          post.qn = postqn,
          title = NULL, subtitle = NULL, caption = NULL
        )
      })

      normalized_counts <- shiny::reactive({
        method <- input$selectmethod
        message("[normalized_counts] length(method) = ", length(method))
        message("[normalized_counts] method = ", method)
        if (length(method) == 0) {
          message("[normalized_counts] setting method to 'none'")
          method <- "none"
        }
        norm.counts <- NULL
        pgx <- pgx()
        if (method == "none") {
          message("[normalized_counts] >>> no normalization")
          norm.counts <- pgx$counts ## no-normalization
        } else {
          message("[normalized_counts] >>> normalizing counts with ", method)
          norm.counts <- playbase::pgx.countNormalization(pgx$counts, method)
        }
        norm.counts
      })

      return(normalized_counts) ## pointing to reactive
    } ## end-of-server
  )
}
