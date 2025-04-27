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

    tab_elements <- list(
      "Gene PCSF" = list(disable = c("gset_accordion")),      
      "Geneset PCSF" = list(disable = c("pcsf_accordion"))
    )

    my_observers[[1]] <- shiny::observeEvent(input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })
    
    
    my_observers[[2]] <- observeEvent(input$pcsf_info, {
      showModal(
        modalDialog(
          title = tags$strong("PCSF Network Analysis"),
          pcsf_info,
          easyClose = TRUE,
          size = "xl"
        )
      )
    })

    my_observers[[3]] <- observe({
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
    
    ## =========================================================================
    ## =========================== FUNCTIONS ===================================
    ## =========================================================================

    ## PCSF  analysis
    pcsf_compute <- shiny::eventReactive(
      {
        list(pgx$X, input$contrast, input$pcsf_beta, input$pcsf_ntop)
      },
      {
        shiny::req(pgx$X, input$contrast, input$pcsf_ntop)
        comparisons <- colnames(pgx$model.parameters$contr.matrix)

        shiny::req(input$contrast %in% comparisons)

        beta <- as.numeric(input$pcsf_beta)
        ntop <- as.integer(input$pcsf_ntop)
        contrast <- input$contrast

        if(pgx$datatype == "multi-omics") {
          ## Multi-omics PCSF
          info("[PcsfBoard:pcsf_compute] computing multi-omics PCSF...")
          datatypes <- unique(playbase::mofa.get_prefix(rownames(pgx$X)))
          if(all(c("gx","px") %in% datatypes)) {
            datatypes <- setdiff(datatypes, c("gx"))
          }
          info("[PcsfBoard:pcsf_compute] datatypes =", datatypes)
          pcsf <- playbase::pgx.computePCSF_multiomics(
            pgx,
            contrast,
            datatypes = datatypes,
            ntop = ntop,
            ncomp = 3,
            beta = 10^beta,
            rm.negedge = TRUE,
            highcor = 0.8, 
            dir = "both",
            ppi = c("STRING", "GRAPHITE"),
            as.name = c("mx")
          ) 
        } else {
          ## Single-omics PCSF
          info("[PcsfBoard:pcsf_compute] computing single-omics PCSF...")
          pcsf <- playbase::pgx.computePCSF(
            pgx,
            contrast,
            level = "gene",
            ntop = ntop,
            ncomp = 2,
            beta = 10^beta,
            dir = "both",
            rm.negedge = TRUE,
            as.name = c("mx")
          )
        }
        
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

    ## =========================================================================
    ## =========================== PANELS ======================================
    ## =========================================================================
    
    pcsf_genepanel_server(
      "genepanel",
      pgx,
      r_contrast = shiny::reactive(input$contrast),
      r_ntop = shiny::reactive(input$pcsf_ntop),
      r_beta = shiny::reactive(input$pcsf_beta),
      r_cut = shiny::reactive(input$pcsf_cut),
      r_nclust = shiny::reactive(input$pcsf_nclust),      
      watermark = WATERMARK
    )

    pcsf_gset_server(
      "gsetpanel",
      pgx,
      r_contrast = shiny::reactive(input$contrast),
      r_ntop = shiny::reactive(input$gset_ntop),
      r_beta = shiny::reactive(input$gset_beta),
      r_cut = shiny::reactive(input$gset_cut),
      r_nclust = shiny::reactive(input$gset_nclust),      
      watermark = WATERMARK
    )

  })
}
