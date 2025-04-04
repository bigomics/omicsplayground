##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

PcsfBoard <- function(id, pgx, board_observers=NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    fullH <- 800
    tabH <- "70vh"

    pcsf_info <- div(
      "This PCSF analysis module..."
    )

    ## ========================================================================
    ## ============================ OBSERVERS =================================
    ## ========================================================================

    my_observers <- list()
    
    my_observers[[1]] <- observeEvent(input$pcsf_info, {
      showModal(
        modalDialog(
          title = tags$strong("PCSF Network Analysis"),
          pcsf_info,
          easyClose = TRUE,
          size = "xl"
        )
      )
    })

    my_observers[[2]] <- observe({
      if (is.null(pgx)) {
        return(NULL)
      }
      comparisons <- playbase::pgx.getContrasts(pgx)
      comparisons <- sort(comparisons[!grepl("^IA:", comparisons)])
      updateSelectInput(
        session,
        "contrast",
        choices = comparisons,
        selected = head(comparisons, 1)
      )
    })

    ## add to list global of observers. suspend by default.
    my_observers <- my_observers[!sapply(my_observers,is.null)]
    # lapply( my_observers, function(b) b$suspend() )
    if(!is.null(board_observers)) board_observers[[id]] <- my_observers
    
    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    ## PCSF  analysis
    pcsf_compute <- shiny::eventReactive(
      {
        list(pgx$X, input$contrast, input$pcsf_beta, input$pcsf_ntop)
      },
      {
        shiny::req(pgx$X, input$contrast)
        comparisons <- colnames(pgx$model.parameters$contr.matrix)
        shiny::req(input$contrast %in% comparisons)

        beta <- as.numeric(input$pcsf_beta)
        ntop <- as.integer(input$pcsf_ntop)
        contrast <- input$contrast
        
        pcsf <- playbase::pgx.computePCSF(
          pgx,
          contrast,
          level = "gene",
          ntop = ntop,
          ncomp = 2,
          beta = 10^beta,
          use.corweight = TRUE,
          dir = "both",
          rm.negedge = TRUE,
          as.name = c("mx")
        )
        if (is.null(pcsf)) {
          validate()
          shiny::validate(
            !is.null(pcsf),
            "No PCSF solution found. Beta value is probably too small. Please adjust beta or increase network size."
          )
          return(NULL)
        }

        pcsf
      }
    )

    pcsf_plot_network_server(
      "pcsf_network",
      pgx,
      pcsf_compute = pcsf_compute,
      r_layout = shiny::reactive(input$layout),
      watermark = WATERMARK
    )

    pcsf_table_centrality_server(
      "centrality_table",
      pgx,
      r_contrast = shiny::reactive(input$contrast),
      r_pcsf = pcsf_compute
    )

    ## pcsf_plot_heatmap_server(
    ##   "pcsf_heatmap",
    ##   pgx,
    ##   pcsf_compute = pcsf_compute,
    ##   watermark = WATERMARK
    ## )
  })
}
