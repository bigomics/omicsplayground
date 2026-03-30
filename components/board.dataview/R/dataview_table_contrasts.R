##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


dataview_table_contrasts_ui <- function(
  id,
  title,
  info.text,
  caption,
  width,
  height
) {
  ns <- shiny::NS(id)

  opts <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(
        ns("ctbygroup"),
        "Show by:",
        choices = c("sample", "group")
      ),
      "Show contrasts by group or by samples.",
      placement = "right", options = list(container = "body")
    )
  )

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    width = width,
    height = height,
    title = title,
    caption = caption,
    options = opts
  )
}


dataview_table_contrasts_server <- function(id,
                                            pgx,
                                            r.samples = reactive(""),
                                            scrollY) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    contrasts_data <- shiny::reactive({
      shiny::req(pgx$Y, pgx$model.parameters)
      shiny::req(r.samples(), !is.null(input$ctbygroup))

      ## dereference reactives
      samples <- r.samples()
      ctbygroup <- input$ctbygroup

      dt <- NULL
      if (ctbygroup == "group") {
        ct <- pgx$model.parameters$contr.matrix
        kk <- which(rownames(ct) %in% pgx$model.parameters$group[samples])
        dt <- ct[kk, , drop = FALSE]
      } else {
        dt <- pgx$model.parameters$exp.matrix[samples, , drop = FALSE]
      }
      dt <- sign(dt)

      # We have pgx$contrasts, however this works also for old pgx files +
      # for ctbygroup == "group" we need to do this conversion anyway,
      # so we do it for all cases.
      cm <- dt
      result <- matrix("",
        nrow = nrow(cm), ncol = ncol(cm),
        dimnames = dimnames(cm)
      )
      for (col in colnames(cm)) {
        treatment <- sub("_vs_.*", "", col)
        control <- sub(".*_vs_", "", col)
        result[cm[, col] == -1, col] <- control
        result[cm[, col] == 1, col] <- treatment
      }
      dt <- result
      dt
    }) ## %>% bindCache(pgx$Y, r.samples(), ctbygroup)


    table.RENDER <- function() {
      dt <- contrasts_data()
      req(dt)

      colnames(dt) <- sub("[_. ]vs[_. ]", "\nvs ", colnames(dt))

      DT::datatable(
        data = dt,
        class = "compact hover",
        rownames = TRUE,
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = list(
          mode = "single",
          target = "row",
          selected = 1
        ),
        options = list(
          dom = "lfrtip",
          scroller = TRUE,
          scrollX = TRUE,
          scrollY = scrollY,
          scrollResize = TRUE,
          deferRender = TRUE,
          autoWidth = TRUE
        )
      ) %>%
        DT::formatStyle(
          columns = 0,
          target = "row",
          fontSize = "14px",
          lineHeight = "70%"
        ) %>%
        DT::formatStyle(
          columns = colnames(dt),
          backgroundColor = DT::styleEqual(
            levels = {
              all_labels <- unique(as.vector(dt))
              all_labels[all_labels != ""]
            },
            values = {
              all_labels <- unique(as.vector(dt))
              all_labels <- all_labels[all_labels != ""]
              controls <- unique(sub(".*\nvs ", "", colnames(dt)))
              ifelse(all_labels %in% controls,
                omics_colors("orange"),
                omics_colors("brand_blue")
              )
            }
          ),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    }

    table.RENDER_modal <- shiny::reactive({
      dt <- table.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL ## nice!
      dt
    })

    table_csv <- function() {
      dt <- table.RENDER()
      dt <- dt$x$data
      colnames(dt) <- make.names(colnames(dt))
      return(dt)
    }

    TableModuleServer(
      "datasets",
      func = table.RENDER,
      csvFunc = table_csv,
      func2 = table.RENDER_modal,
      selector = "none"
    )
  }) ## end of moduleServer
} ## end of server
