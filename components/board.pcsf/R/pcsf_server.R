##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

PcsfBoard <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    fullH <- 800
    tabH <- "70vh"

    pcsf_info <- div(
      "This PCSF analysis module..."
    )

    ## Visibility gating: pause this board's own observers while its tab is
    ## off screen, resume (re-running if invalidated while suspended) when
    ## shown again. See components/ui/ui-board-visibility.R.
    is_visible <- board_is_visible(input, label = "PcsfBoard")
    observers <- board_observer_registry()

    ## ========================================================================
    ## ============================ OBSERVERS =================================
    ## ========================================================================

    tab_elements <- list(
      "Gene PCSF" = list(disable = c("gset_accordion")),
      "Geneset PCSF" = list(disable = c("pcsf_accordion"))
    )

    observers$add(shiny::observeEvent(input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    }))

    observers$add(observeEvent(input$pcsf_info, {
      showModal(
        modalDialog(
          title = tags$strong("PCSF Network Analysis"),
          pcsf_info,
          easyClose = TRUE,
          size = "xl"
        )
      )
    }))

    observers$add(observe({
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
    }))

    ## =========================================================================
    ## =========================== FUNCTIONS ===================================
    ## =========================================================================


    ## =========================================================================
    ## =========================== PANELS ======================================
    ## =========================================================================

    pcsf_genepanel_server(
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

    board_pause_resume_observers(is_visible, observers, label = "PcsfBoard")
  })
}
