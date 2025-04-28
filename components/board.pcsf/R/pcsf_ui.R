##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

PcsfInputs <- function(id) {
  ns <- NS(id)

  bigdash::tabSettings(
    withTooltip(
      selectInput(ns("contrast"), "Select contrast:",
        choices = NULL, multiple = FALSE
      ),
      "Select contrast.",
      placement = "right"
    ),
    hr(),
    br(),
    bslib::accordion(
      id = ns("pcsf_accordion"),
      open = FALSE,
      bslib::accordion_panel(
        "Options",
        icon = icon("cog", lib = "glyphicon"),
        pcsf_genepanel_settings_ui(ns("genepanel"))
        ## withTooltip(
        ##   shiny::radioButtons(
        ##     ns("pcsf_ntop"), "Network size:",
        ##     choices = c("S" = 250, "M" = 500, "L" = 1000, "XL" = 2000),
        ##     selected = 500, inline = TRUE
        ##   ),
        ##   "Select initial network size (number of top genes) for ."
        ## ),
        ## hr(),
        ## withTooltip(
        ##   shiny::checkboxInput(ns("pcsf_cut"), "Cut clusters", FALSE),
        ##   "Cut network into smaller clusters"
        ## ),
        ## shiny::conditionalPanel(
        ##   "input.pcsf_cut == true",
        ##   ns = ns,
        ##   withTooltip(
        ##     shiny::selectInput(ns("pcsf_nclust"), "Number of clusters",
        ##                        choices = c(1,4,9,16,25,99), selected=9),
        ##     "Maximum number of components"
        ##   ),
        ##   shiny::selectInput(ns("pcsf_resolution"), "Resolution",
        ##                        choices = c(0.01,0.05,0.1,0.2,0.5,1), selected=0.1)
        ## ),
        ## hr()        
      )
    ),
    bslib::accordion(
      id = ns("gset_accordion"),
      open = FALSE,
      bslib::accordion_panel(
        "Options",
        icon = icon("cog", lib = "glyphicon"),
        pcsf_gsetpanel_settings_ui(ns("gsetpanel"))
        ## withTooltip(
        ##   shiny::radioButtons(ns("gset_ntop"), "Network size:",
        ##     choices = c("S" = 250, "M" = 500, "L" = 1000, "XL" = 2000),
        ##     selected = 1000, inline = TRUE
        ##   ),
        ##   "Select initial network size (number of top genes) for ."
        ## ),
        ## hr(),
        ## withTooltip(
        ##   shiny::checkboxInput(ns("gset_cut"), "Cut clusters", TRUE),
        ##   "Cut network into smaller clusters"
        ## ),
        ## shiny::conditionalPanel(
        ##   "input.gset_cut == true",
        ##   ns = ns,
        ##   withTooltip(
        ##     shiny::selectInput(ns("gset_nclust"), "Number of clusters",
        ##       choices = c(1,4,9,16,25,99), selected=9),
        ##     "Maximum number of components"
        ##   ),
        ##   shiny::selectInput(ns("gset_resolution"), "Resolution",
        ##     choices = c(0.01,0.05,0.1,0.2,0.5,1), selected=0.1)
        ## )
      )
    )
  )
}

pcsf_module_info <- "The PCSF network analysis uses the Prize-collection Steiner Forest algorithm to determine high-correlated subnetworks of the most differentially expressed genes. Interactions from STRING and pathway databases are used as template. The PCSF solution may be used to identify 'driver' genes that appear as hubs in the computed network."

pcsf_graph_info <- "Prize-collection Steiner Forest solution for the top differential genes using the STRING database as backbone. 'Driver' genes appear as hubs in the network computed using a page-rank centrality measure."

PcsfUI <- function(id) {
  ns <- NS(id)
  tagList(
    boardHeader(
      title = "Prize-Collecting Steiner Forest", info_link = ns("pcsf_info")
    ),
    shiny::tabsetPanel(
      id = ns("tabs"),
      shiny::tabPanel(
        "Gene PCSF",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          ##--------- begin tab content ------------
          bs_alert(pcsf_module_info),
          bslib::layout_columns(
            col_widths = c(6, 6),
            height = "calc(100vh - 181px)",
            bslib::layout_columns(
              col_widths = c(12),
              pcsf_genepanel_seriesplot_ui(
                ns("genepanel"),
                caption = "",
                info.text = pcsf_graph_info,
                height = c("100%", "75vh"),
                width = c("auto", "100%")
              ),
              pcsf_genepanel_table_ui(
                ns("genepanel"),
                title = "Centrality score",
                info.text = "",
                caption = "Table showing the centrality score of genes.",
                width = c("100%", "100%"),
                height = c("100%", TABLE_HEIGHT_MODAL)
              )
            ),
            pcsf_genepanel_networkplot_ui(
              ns("genepanel"),
              caption = "PCSF network analysis. Functional analysis of biological networks using Prize-collection Steiner Forest algorithm that determines high-confidence subnetworks.",
              info.text = pcsf_graph_info,
              height = c("100%", "75vh"),
              width = c("auto", "100%")
            )            
          )
          ##--------- end tab content ------------            
        )
      ),  ## end tabpanel2

      shiny::tabPanel(
        "Geneset PCSF",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          ##--------- begin tab content ------------
          bs_alert(pcsf_module_info),
          bslib::layout_columns(
            col_widths = c(6, 6),
            height = "calc(100vh - 181px)",
            bslib::layout_columns(
              col_widths = c(12),
              pcsf_gsetpanel_seriesplot_ui(
                ns("gsetpanel"),
                title = "Centrality score",
                caption = "",
                info.text = pcsf_graph_info,
                height = c("100%", "75vh"),
                width = c("auto", "100%")
              ),
              pcsf_gsetpanel_table_ui(
                ns("gsetpanel"),
                title = "Centrality score",
                info.text = "",
                caption = "Table showing the centrality score of genes.",
                width = c("100%", "100%"),
                height = c("100%", TABLE_HEIGHT_MODAL)
              )
            ),
            pcsf_gsetpanel_networkplot_ui(
              ns("gsetpanel"),
              caption = "PCSF network analysis. Functional analysis of biological networks using Prize-collection Steiner Forest algorithm that determines high-confidence subnetworks.",
              info.text = pcsf_graph_info,
              height = c("100%", "75vh"),
              width = c("auto", "100%")
            )            
          )
          ##--------- end tab content ------------            
        )
      ) ## end tabpanel2

    ) ## end tabsetpanel
  )
}
