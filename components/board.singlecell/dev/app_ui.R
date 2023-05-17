#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {

  # header
  header <- shiny::tagList(
    golem_add_external_resources('board.singlecell')
  )

  # sidebar
  sidebar <- bigdash::sidebar(
    "Menu",
    div(class="sidebar-item",
        bigdash::sidebarItem(
          "Cell profiling",
          "cell-tab"
        )
    )
  )

  big_theme <- bslib::bs_add_variables(
    bigdash::big_theme(),
    "grid-breakpoints" = "map-merge($grid-breakpoints, ('xxxl': 2400px))",
    .where = "declarations"
  )

  bigdash::bigPage(
    header,
    title = "Omics Playground v3",
    theme = big_theme,
    sidebar = sidebar,
    navbar = bigdash::navbar(
      tags$img(
        id = "logo-bigomics",
        src = "assets/img/bigomics.png",
        width = "110",
      )
    ),
    settings = bigdash::settings(
      "Settings"
    ),
    bigdash::bigTabs(
      bigdash::bigTabItem(
        "cell-tab",
        SingleCellInputs("cell"),
        SingleCellUI("cell")
      )
    )
  )

}


