##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## EWAS settings panel.
##
## bigdash::tabSettings is per bigTabItem, not per nav_panel, so all three
## EWAS sub-tabs share one settings block. Rather than showing every control
## everywhere - which invites re-running the whole model from a screen that
## only consumes its output - each group is shown only on the sub-tabs that
## actually use it, keyed on the navset id.
##
##   hits    : the model itself, plus the threshold
##   regions : threshold (it defines the hit list gometh tests), DMR and
##             gene-set controls
##   qq      : threshold, and how many CpGs to chart
##
## The fitted-model summary stays visible everywhere: on the regions and QQ
## screens it is the answer to "what are these results based on?".

methylome_ewas_inputs <- function(id) {
  ns <- shiny::NS(id)
  on_tab <- function(..., x) {
    shiny::conditionalPanel(
      condition = paste(sprintf("input.ewas_subtab == '%s'", c(...)), collapse = " || "),
      ns = ns, x
    )
  }

  bigdash::tabSettings(

    ## ---- the model: only where it is run ----
    on_tab("hits", x = shiny::tagList(
      withTooltip(
        shiny::selectInput(ns("ewas_contrast"), "Outcome:", choices = NULL),
        "What to test against methylation: a two-group contrast, or a continuous variable from the sample sheet such as age, BMI or pack-years. Changing it refits the model.",
        placement = "top"
      ),
      withTooltip(
        shiny::selectizeInput(ns("ewas_covars"), "Adjust for:", choices = NULL,
                              multiple = TRUE,
                              options = list(placeholder = "none (unadjusted)")),
        "Covariates added to the model. In bulk tissue an unadjusted comparison is usually a test of cell composition rather than of the phenotype, so adjusting is the norm rather than the exception.",
        placement = "top"
      ),
      withTooltip(
        shiny::checkboxInput(ns("ewas_adjust_cells"), "Adjust for cell composition", FALSE),
        "Add the estimated cell proportions as covariates. Estimate them first on the Cell composition screen.",
        placement = "top"
      ),
      withTooltip(
        shiny::checkboxGroupInput(
          ns("ewas_mask"), "Mask probes:",
          choices = c("Common SNP in probe" = "snp", "Cross-reactive" = "xreactive"),
          selected = c("snp", "xreactive")
        ),
        "Probes carrying a common SNP measure genotype rather than methylation; cross-reactive probes map to more than one locus. Both are standard QC exclusions.",
        placement = "top"
      ),
      shiny::actionButton(ns("run_ewas"), "Run EWAS",
                          class = "btn btn-primary btn-sm", width = "100%")
    )),

    ## ---- what was fitted: context everywhere ----
    shiny::uiOutput(ns("ewas_model")),
    shiny::tags$hr(style = "margin:10px 0;"),

    ## ---- threshold: defines the hit list on every sub-tab ----
    withTooltip(
      shiny::selectInput(ns("thresh_type"), "Threshold on:",
                         choices = c("FDR (q-value)" = "q", "Nominal p-value" = "p"),
                         selected = "q"),
      "Call hits on the multiple-testing corrected q-value, or on the raw p-value.",
      placement = "top"
    ),
    withTooltip(
      shiny::numericInput(ns("thresh_value"), "Cut-off:", value = 0.05,
                          min = 0, max = 1, step = 0.01),
      "Probes at or below this value are called hits.",
      placement = "top"
    ),
    withTooltip(
      shiny::numericInput(ns("min_dbeta"), "Minimum |delta-beta|:", value = 0,
                          min = 0, max = 1, step = 0.01),
      "Effect-size filter on the beta difference between groups. 0 disables it. A hard cut deletes the real signal of a population EWAS, where true differences are typically under 5%.",
      placement = "top"
    ),

    ## ---- charting: only where the stripcharts are ----
    on_tab("qq", x = withTooltip(
      shiny::numericInput(ns("top_n"), "Top CpGs to chart:", value = 6,
                          min = 2, max = 12, step = 1),
      "How many of the most significant CpGs to draw as per-sample beta stripcharts.",
      placement = "top"
    )),

    ## ---- regions and gene sets: only on their own sub-tab ----
    on_tab("regions", x = shiny::tagList(
      shiny::tags$hr(style = "margin:10px 0;"),
      withTooltip(
        shiny::numericInput(ns("dmr_maxgap"), "Max gap between CpGs (bp):",
                            value = 500, min = 100, max = 5000, step = 100),
        "CpGs further apart than this are not joined into the same region.",
        placement = "top"
      ),
      shiny::actionButton(ns("run_dmr"), "Call regions",
                          class = "btn btn-primary btn-sm", width = "100%"),
      withTooltip(
        shiny::selectInput(ns("dmr_pick"), "Region to draw:", choices = NULL),
        "Which called region the detail plot below the table shows. Populated once regions have been called.",
        placement = "top"
      ),
      withTooltip(
        shiny::selectInput(ns("gs_collection"), "Gene-set collection:",
                           choices = c("GO", "KEGG"), selected = "GO"),
        "Which annotation collection to test the hit list against.",
        placement = "top"
      ),
      shiny::actionButton(ns("run_gometh"), "Test gene sets",
                          class = "btn btn-primary btn-sm", width = "100%")
    ))
  )
}
