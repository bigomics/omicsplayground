##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


admin_table_users_ui <- function(
  id,
  title,
  height,
  width = c("auto", "100%"),
  caption,
  info.text
) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("tbl"),
    width = width,
    height = height,
    title = title,
    info.text = info.text,
    caption = caption
  )
}

admin_table_users_server <- function(id, auth,
                                     scrollY = "calc(100vh - (240px + 140px))") {
  moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      shiny::req(isTRUE(auth$ADMIN))
      dbg("[AdminBoard::user_stats] table_data")

      user_dirs <- list.dirs(PGX.DIR, full.names = FALSE, recursive = FALSE)
      user_dirs <- grep("@", user_dirs, value = TRUE)

      dt <- data.frame(
        Metric = c("Total User Directories", "Data Directory"),
        Value = c(length(user_dirs), PGX.DIR),
        check.names = FALSE
      )
      dt
    })

    table.RENDER <- function() {
      dt <- table_data()
      req(dt)
      DT::datatable(dt,
        class = "compact hover",
        rownames = FALSE,
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = 1),
        options = list(
          dom = "lfrtip",
          scroller = TRUE,
          scrollX = TRUE,
          scrollY = scrollY,
          scrollResize = TRUE,
          deferRender = TRUE
        )
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    modal_table.RENDER <- function() {
      dt <- table_data()
      req(dt)
      DT::datatable(dt,
        class = "compact hover",
        rownames = FALSE,
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        options = list(
          dom = "lfrtip",
          scroller = TRUE,
          scrollX = TRUE,
          scrollY = SCROLLY_MODAL,
          deferRender = TRUE
        )
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "20px", lineHeight = "70%")
    }

    TableModuleServer(
      "tbl",
      func = table.RENDER,
      func2 = modal_table.RENDER,
      selector = "none"
    )
  }) ## end of moduleServer
} ## end of server
