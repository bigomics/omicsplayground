##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Settings for the Methylome landscape screen - the panels moved over from the
## old Epigenomics board, whose plot modules are unchanged and still carry their
## epigenomics_ names.
##
## The original showed the chromosome selector only on the ideogram tab, but the
## chromosomal boxplot reads the same input, so on the landscape tab there was no
## way to change the chromosomes it drew. All three controls are shown on both
## sub-tabs instead, which also means no conditionalPanel is needed.

methylome_landscape_inputs <- function(id) {
  ns <- shiny::NS(id)
  bigdash::tabSettings(
    withTooltip(
      shiny::selectizeInput(ns("select_chromosome"), "Chromosome:",
                            choices = NULL, multiple = TRUE),
      "Which chromosomes to draw. The ideogram takes at most six at a time; the chromosomal boxplot has no such limit.",
      placement = "top"
    ),
    withTooltip(
      shiny::selectInput(ns("data_samplefilter"), "Filter out samples:",
                         choices = NULL, multiple = TRUE),
      "Samples to exclude from both panels.",
      placement = "top"
    ),
    withTooltip(
      shiny::selectInput(ns("select_pheno"), "Select phenotype:", choices = NULL),
      "Split the profiles by a sample variable. With <ungrouped> everything is pooled; with a two-level phenotype the ideogram gains a per-group delta-beta track.",
      placement = "top"
    )
  )
}
