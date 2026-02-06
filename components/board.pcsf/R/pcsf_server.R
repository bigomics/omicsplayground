##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

PcsfBoard <- function(id, pgx) {
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
      "Geneset PCSF" = list(disable = c("pcsf_accordion")),
      "AI Summary" = list(disable = c("pcsf_accordion", "gset_accordion"))
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

    ## =========================================================================
    ## =========================== FUNCTIONS ===================================
    ## =========================================================================


    ## =========================================================================
    ## =========================== PANELS ======================================
    ## =========================================================================

    genepanel_out <- pcsf_genepanel_server(
      "genepanel",
      pgx,
      r_contrast = shiny::reactive(input$contrast),
      watermark = WATERMARK
    )

    pcsf_gsetpanel_server(
      "gsetpanel",
      pgx,
      r_contrast = shiny::reactive(input$contrast),
      watermark = WATERMARK
    )

    # AI PCSF summary
    pcsf_ai_summary_server(
      "pcsfAISummary",
      pgx = pgx,
      getCentralityTable = genepanel_out$table_data,
      getPcsfGraph = shiny::reactive({
        res <- genepanel_out$gene_pcsf()
        shiny::req(res)
        res$graph
      }),
      r_contrast = shiny::reactive(input$contrast),
      session = session,
      watermark = WATERMARK
    )
  })
}
