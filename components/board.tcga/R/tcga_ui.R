##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## DEAN ATTALI code recommendation:
## example of how a board module should be written (TcgaBoard)
##
## https://github.com/bigomics/omicsplayground/pull/20/commits/bd943d84d316d76dca9140f2fd3610b3d1dfc950


TcgaInputs <- function(id) {
  ns <- NS(id)

  bigdash::tabSettings(
    withTooltip(
      radioButtons(
        ns("sigtype"),
        "Signature type:",
        choices = c("contrast", "genelist"),
        selected = "contrast",
        inline = TRUE
      ),
      "Choose the type of signature as input",
      placement = "right",
      options = list(container = "body")
    ),
    conditionalPanel(
      "input.sigtype == 'contrast'",
      ns = ns,
      withTooltip(
        selectInput(ns("contrast"), NULL, choices = NULL, multiple = FALSE),
        "Select the contrast that you want to correlate with survival.",
        placement = "right",
        options = list(container = "body")
      ),
    ),
    conditionalPanel(
      "input.sigtype == 'genelist'",
      ns = ns,
      withTooltip(
        textAreaInput(
          ns("genelist"),
          NULL,
          value = NULL,
          height = "100px",
          width = "100%",
          rows = 4,
          placeholder = "Paste your custom gene list"
        ),
        "Paste a custom list of genes to be used as features.",
        placement = "bottom"
      )
    ),
    br(),
    withTooltip(
      actionLink(ns("tcga_options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "top",
      options = list(container = "body")
    )
  )
}

tcga_info <- "This analysis module computes the survival probability in (more than 10000) cancer patients of 32 TCGA cancer types, for your selected contrast. Each cohort is dichotomized into positively and negatively correlated with your signature. The survival probabilities are computed and tested using the Kaplan-Meier method."

TcgaUI <- function(id) {
  ns <- NS(id)

  div(
    boardHeader(
      title = "TCGA",
      info_link = ns("tcga_info")
    ),
    tabsetPanel(
      id = ns("tabs1"),
      tabPanel(
        "TCGA survival",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          tcga_plot_survival_ui(
            ns("tcga_tcgasurv"),
            caption = paste(
              "TCGA survival analysis",
              "Survival probability of cancer patients in 32 TCGA cancer types.",
              "Each cohort is dichotomized into positively and negatively correlated with your signature.",
              "The survival probabilities are computed and tested using the Kaplan-Meier method."
            ),
            info.text = tcga_info,
            height = c("100%", "70vh"),
            width = c("auto", "100%")
          )
        )
      ),
      tabPanel(
        "AI Summary",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          AiTextCardUI(
            ns("tcgaAISummary"),
            title = "AI TCGA Summary",
            info.text = "AI-generated summary of the TCGA survival analysis results for the selected contrast.",
            caption = "AI-generated TCGA survival summary.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      )
    )
  )
}
