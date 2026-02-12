##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mox_ui <- function(id) {
  ns <- shiny::NS(id) ## namespace

  ui <- bigdash::bigPage(
    ##----------------------------------------------------------------------
    navbar = bigdash::navbar(
      tags$img(
        src = "assets/img/bigomics.png",
        width = "110"
      ),
      tags$div(
        "title in navbar",
        style = 'text-align:center;width: 100%;'
      ),
      bigdash::navbarDropdown(
        "Support",
        bigdash::navbarDropdownItem(
          "Documentation"
        ),
        bigdash::navbarDropdownItem(
          "Contact"
        )
      ),
      bigdash::navbarDropdown(
        "Tutorials",
        bigdash::navbarDropdownItem(
          "Get started"
        ),
        bigdash::navbarDropdownItem(
          "Advanced"
        )
      )
    ),
    ##----------------------------------------------------------------------
    sidebar = bigdash::sidebar(
      "Menu",
      bigdash::sidebarItem(
        "Home",
        "xhome"
      ),
      bigdash::sidebarMenu(
        "Upload",
        bigdash::sidebarMenuItem(
          "Data",
          "xtab1"
        ),
        bigdash::sidebarMenuItem(
          "Upload",
          "xtab2"
        )
      ),
      bigdash::sidebarItem(
        "TabSetPanel",
        "xtab3"
      )
    ),
    ##----------------------------------------------------------------------
    settings = bigdash::settings(
      "Settings",
      p(
        "Settings will appear here."
      )
    ),
    ##----------------------------------------------------------------------
    bigdash::sidebarHelp(
      bigdash::sidebarTabHelp(
        "home",
        "Welcome!",
        "This is the homepage, welcome!"
      ),
      bigdash::sidebarTabHelp(
        "tab1",
        "Upload",
        "This is the first tab!"
      )
    ),
    bigdash::bigTabs(
      bigdash::bigTabItem(
        "xhome",
        bigdash::fullPage(
          .class = "bg-secondary text-center",
          tags$img(
            src = "assets/img/mascotte-sc.png",
            class = "img-fluid",
            style = "max-height: 20rem;"
          )
        )
      ),
      bigdash::bigTabItem(
        "xtab1",
        div(
          class = "p-4",
          h3("TabSetPanel"),
          shiny::tabsetPanel(
            shiny::tabPanel(
              "First subtab",
              tags$img(
                src = "www/tab1.png",
                width = "100%"
              )
            ),
            shiny::tabPanel(
              "Second subtab",
              tags$img(
                src = "www/tab2.png",
                width = "100%"
              )
            )
          )
        )
      ),
      bigdash::bigTabItem(
        "xtab2",
        div(
          class = "p-4",
          h2("Tab2")
        )
      ),
      bigdash::bigTabItem(
        "xtab3",
        div(
          class = "p-4",
          h2("Tab3")
        )
      )
    )
  )
  
  return(ui)
}
