##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## DEAN ATTALI code recommendation:
## example of how a board module should be written (TcgaBoard)
##
## https://github.com/bigomics/omicsplayground/pull/20/commits/bd943d84d316d76dca9140f2fd3610b3d1dfc950

TcgaBoard <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    fullH <- 800
    tabH <- "70vh"

    tcga_tcgasurv_info <- div(
      "This", tags$strong("TCGA analysis module"),
      "computes the survival probability in (more than 10000) cancer patients of 32 TCGA cancer types, for your selected contrast.",
      "Each cohort is dichotomized into positively and negatively correlated with your signature.",
      "The survival probabilities are computed and tested using the Kaplan-Meier method."
    )

    ## ========================================================================
    ## ============================ OBSERVERS =================================
    ## ========================================================================

    observeEvent(input$tcga_info, {
      showModal(
        modalDialog(
          title = tags$strong("TCGA Analysis Board"),
          tcga_tcgasurv_info,
          easyClose = TRUE,
          size = "l"
        )
      )
    })

    observe({
      if (is.null(pgx)) {
        return(NULL)
      }
      comparisons <- colnames(pgx$model.parameters$contr.matrix)
      comparisons <- sort(comparisons)
      updateSelectInput(session, "contrast", choices = comparisons, selected = head(comparisons, 1))
    })

    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    # TCGA survival analysis

    tcga_plot_survival_server(
      "tcga_tcgasurv",
      pgx,
      contrast = shiny::reactive(input$contrast),
      sigtype = shiny::reactive(input$sigtype),
      genelist = shiny::reactive(input$genelist),
      watermark = WATERMARK
    )

    # AI TCGA summary

    tcga_ai_summary_server(
      "tcgaAISummary",
      pgx = pgx,
      contrast = shiny::reactive(input$contrast),
      sigtype = shiny::reactive(input$sigtype),
      genelist = shiny::reactive(input$genelist),
      session = session,
      watermark = WATERMARK
    )
  })
}
