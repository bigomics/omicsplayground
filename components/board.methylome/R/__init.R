##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Methylome board registration. Same shape as MODULE.wgcna: module_ui declares
## the tab shells with their settings panels, module_ui2 supplies the content
## that the platform inserts on demand, and module_server wires one server for
## the whole menu group.
##
## Tab ids are global (the platform inserts them into its own bigTabs), but the
## module namespace stays "methylome" so every child module id resolves exactly
## as it did when this was a standalone app.

MODULE.methylome <- list(
  module_menu = function() {
    c(
      mledger = "Sample ledger",
      mage = "Epigenetic age",
      mchar = "Methylome character",
      mland = "Methylome landscape",
      mdeconv = "Cell composition",
      mewas = "EWAS"
    )
  },
  module_ui = function() {
    list(
      bigdash::bigTabItem(
        "mledger-tab",
        create_loader("mledger-loader")
      ),
      bigdash::bigTabItem(
        "mage-tab",
        methylome_age_inputs("methylome"),
        create_loader("mage-loader")
      ),
      bigdash::bigTabItem(
        "mchar-tab",
        create_loader("mchar-loader")
      ),
      bigdash::bigTabItem(
        "mland-tab",
        methylome_landscape_inputs("methylome"),
        create_loader("mland-loader")
      ),
      bigdash::bigTabItem(
        "mdeconv-tab",
        methylome_deconv_inputs("methylome"),
        create_loader("mdeconv-loader")
      ),
      bigdash::bigTabItem(
        "mewas-tab",
        methylome_ewas_inputs("methylome"),
        create_loader("mewas-loader")
      )
    )
  },
  module_ui2 = function() {
    list(
      list("mledger-tab", methylome_ui_ledger("methylome")),
      list("mage-tab", methylome_ui_age("methylome")),
      list("mchar-tab", methylome_ui_character("methylome")),
      list("mland-tab", methylome_ui_landscape("methylome")),
      list("mdeconv-tab", methylome_ui_deconv("methylome")),
      list("mewas-tab", methylome_ui_ewas("methylome"))
    )
  },
  module_server = function(PGX) {
    methylome_server("methylome", pgx = PGX)
  },
  module_help = function() {
    list(
      bigdash::sidebarTabHelp(
        "mledger-tab", "Sample ledger",
        tspan("Per-sample quality, identity and epigenetic age. Needs no contrast.")
      ),
      bigdash::sidebarTabHelp(
        "mage-tab", "Epigenetic age",
        tspan("Epigenetic clocks, age acceleration and per-clock coverage.")
      ),
      bigdash::sidebarTabHelp(
        "mchar-tab", "Methylome character",
        tspan("Where methylation sits: CpG island context and position in gene.")
      ),
      bigdash::sidebarTabHelp(
        "mland-tab", "Methylome landscape",
        tspan("Genome-wide overview across chromosomes and samples, and ideograms.")
      ),
      bigdash::sidebarTabHelp(
        "mdeconv-tab", "Cell composition",
        tspan("Estimated cell proportions - the dominant confounder in bulk tissue.")
      ),
      bigdash::sidebarTabHelp(
        "mewas-tab", "EWAS",
        tspan("Covariate-adjusted epigenome-wide association, regions and gene sets.")
      )
    )
  }
)
