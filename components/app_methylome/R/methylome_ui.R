##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Methylome Profiler shell. Same construction as qsee_ui(): a bigdash page
## with one sidebarItem per screen and one bigTabItem holding it.

methylome_ui <- function(id = "methylome", height = "100%") {
  ns <- shiny::NS(id)

  bigdash::bigPage(
    id = id,
    navbar = shiny::div(style = "visibility: hidden; display: none;"),
    sidebar = bigdash::sidebar(
      "Methylome Profiler",
      id = id,
      bigdash::sidebarItem("Sample ledger", ns("ledger-tab")),
      bigdash::sidebarItem("Epigenetic age", ns("age-tab")),
      bigdash::sidebarItem("Methylome character", ns("character-tab")),
      bigdash::sidebarItem("EWAS", ns("ewas-tab"))
    ),
    bigdash::bigTabs(
      bigdash::bigTabItem("ledger-tab", methylome_ledger_ui(ns("ledger"))),
      bigdash::bigTabItem("age-tab", methylome_age_ui(ns("age"))),
      bigdash::bigTabItem("character-tab", methylome_character_ui(ns("character"))),
      bigdash::bigTabItem("ewas-tab", methylome_ewas_ui(ns("ewas")))
    )
  )
}
