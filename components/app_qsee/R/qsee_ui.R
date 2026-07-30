##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


qsee_ui.save <- function(id, height = "100%") {
  ns <- shiny::NS(id)


  title <- div("Qsee/Bsee", style = "font-size: 18px;")
  
  ui <- bslib::page_fillable(
    padding = 0,
    div(class = "navbar navbar-static-top", div(title, class = "container-fluid"),
      style = "margin-top: 24px; height: 55px;"),
  
    bslib::navset_pill_list(
      widths = c(2, 10),
      bslib::nav_panel(
        title = "Normalization",
        qsee_normalization_ui(ns("normalize"))
      ),
      bslib::nav_panel(
        title = "Missing values",
        qsee_imputation_ui(ns("impute"))
      ),
      bslib::nav_panel(
        title = "Batch-effects",
        qsee_bsee_ui(ns("bsee"))
      ),      
      bslib::nav_panel(
        title = "Filtering",
        qsee_filtering_ui(ns("filtering"))
      ),
      bslib::nav_spacer(),
      selectInput(ns("testinput"),"test", choices=NULL)
    )        
  )
  
  return(ui)
}

qsee_ui <- function(id, height = "100%") {
  ns <- shiny::NS(id)

  title <- div("Qsee/Bsee", style = "font-size: 18px;")
  
  ui <- bigdash::bigPage(
    id = id,
    navbar = bigdash::navbar("Qsee/Bsee", style="margin-top: 20px;"),
    sidebar = bigdash::sidebar("Menu",
      id = id,
      style = "max-height: 100vh !important;",      
      bigdash::sidebarItem("Normalization", ns("normalize-tab")),
      bigdash::sidebarItem("Missing values", ns("impute-tab")),
      bigdash::sidebarItem("Outlier analysis", ns("outlier-tab")),
      bigdash::sidebarItem("Batch-effects", ns("bsee-tab")),
      bigdash::sidebarItem("Filtering", ns("filtering-tab")),
      br(),
      shiny::actionButton(ns("upload"),"upload CSV")      
    ),
    settings = bigdash::settings("Settings", p("Settings will appear here."), id = id),
    bigdash::sidebarHelp(
      id = id
      #sidebarTabHelp("home", "Welcome!","This is the homepage, welcome!")
    ),
    bigdash::bigTabs(
      id = id,
      bigdash::bigTabItem(ns("normalize-tab"), qsee_normalization_ui(ns("normalize"))),
      bigdash::bigTabItem(ns("impute-tab"), qsee_imputation_ui(ns("impute"))),
      bigdash::bigTabItem(ns("outlier-tab"), qsee_outlier_ui(ns("outlier"))),
      bigdash::bigTabItem(ns("bsee-tab"), qsee_bsee_ui(ns("bsee"))),
      bigdash::bigTabItem(ns("filtering-tab"), qsee_filtering_ui(ns("filtering")))      
    )
  )

  return(ui)
}
