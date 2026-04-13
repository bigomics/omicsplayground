#' Evidence Panel Sub-module
#'
#' Right panel: dataset context + plot + table slots for agent output.
#' The main plot viewer uses PlotModuleUI/PlotModuleServer for consistent
#' download and enlarge behavior.  Plot history is kept in-memory for the
#' current session and exposed as a compact carousel strip.

copilot_panel_evidence_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(
    style = "padding: 8px; height: 100%; overflow-y: auto;",

    ## Dataset context card (Task 4.2)
    shiny::conditionalPanel(
      condition = paste0("output['", ns("has_dataset"), "']"),
      bslib::card(
        class = "mb-2",
        bslib::card_header("Dataset Info", class = "py-1 px-2 fw-bold"),
        bslib::card_body(
          class = "p-2",
          shiny::uiOutput(ns("dataset_info"))
        )
      )
    ),

    ## Plot card — three PlotModuleUI instances, one per plotlib kind.
    ## Only the instance matching the active plot kind is visible.
    shiny::conditionalPanel(
      condition = paste0("output['", ns("has_plot"), "']"),
      shiny::conditionalPanel(
        condition = paste0("output['", ns("is_ggplot_active"), "']"),
        PlotModuleUI(
          ns("evidence_ggplot"),
          plotlib = "ggplot",
          height = c(350, 700),
          caption = "Evidence Plot",
          info.text = "Copilot-generated evidence plot",
          download.fmt = c("png", "pdf"),
          translate = FALSE
        )
      ),
      shiny::conditionalPanel(
        condition = paste0("output['", ns("is_plotly_active"), "']"),
        PlotModuleUI(
          ns("evidence_plotly"),
          plotlib = "plotly",
          height = c(350, 700),
          caption = "Evidence Plot",
          info.text = "Copilot-generated evidence plot",
          download.fmt = c("png", "pdf"),
          translate = FALSE
        )
      ),
      shiny::conditionalPanel(
        condition = paste0("output['", ns("is_iheatmapr_active"), "']"),
        PlotModuleUI(
          ns("evidence_iheatmapr"),
          plotlib = "iheatmapr",
          height = c(350, 700),
          caption = "Evidence Plot",
          info.text = "Copilot-generated evidence plot",
          download.fmt = c("png"),
          translate = FALSE
        )
      )
    ),

    ## Plot history carousel — visible only when 2+ plots in history
    shiny::conditionalPanel(
      condition = paste0("output['", ns("has_history"), "']"),
      shiny::div(
        class = "mb-2",
        shiny::div(
          class = "d-flex align-items-center mb-1",
          shiny::tags$small(class = "text-muted fw-bold", "Plot History")
        ),
        shiny::div(
          style = "overflow-x: auto; white-space: nowrap; padding: 4px 0;",
          shiny::uiOutput(ns("plot_carousel"))
        )
      )
    ),

    ## Table card
    shiny::conditionalPanel(
      condition = paste0("output['", ns("has_table"), "']"),
      bslib::card(
        class = "mb-2",
        bslib::card_header("Results", class = "py-1 px-2"),
        bslib::card_body(
          class = "p-1",
          DT::dataTableOutput(ns("evidence_table"))
        )
      )
    ),

    ## Empty state
    shiny::conditionalPanel(
      condition = paste0("!output['", ns("has_dataset"), "']"),
      shiny::div(
        style = "text-align: center; padding: 40px; color: #999;",
        shiny::icon("chart-line", class = "fa-3x", style = "color: #ddd;"),
        shiny::br(), shiny::br(),
        shiny::p("Evidence and visualizations will appear here when you interact with the Copilot.")
      )
    )
  )
}

copilot_detect_plot_kind <- function(plot_obj) {
  if (inherits(plot_obj, "plotly")) return("plotly")
  if (inherits(plot_obj, "iheatmapr")) return("iheatmapr")
  if (inherits(plot_obj, "ggplot")) return("ggplot")
  NULL
}

copilot_panel_evidence_server <- function(id, local_pgx = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## --- Plot history state ---
    plot_history <- shiny::reactiveVal(list())
    active_plot_idx <- shiny::reactiveVal(NULL)

    active_record <- shiny::reactive({
      idx <- active_plot_idx()
      hist <- plot_history()
      if (is.null(idx) || idx < 1L || idx > length(hist)) return(NULL)
      hist[[idx]]
    })

    active_kind <- shiny::reactive({
      rec <- active_record()
      if (is.null(rec)) return(NULL)
      rec$kind
    })

    table_data <- shiny::reactiveVal(NULL)

    ## --- Conditional panel flags ---
    output$has_plot <- shiny::reactive({ !is.null(active_record()) })
    shiny::outputOptions(output, "has_plot", suspendWhenHidden = FALSE)

    output$is_ggplot_active <- shiny::reactive({ identical(active_kind(), "ggplot") })
    shiny::outputOptions(output, "is_ggplot_active", suspendWhenHidden = FALSE)

    output$is_plotly_active <- shiny::reactive({ identical(active_kind(), "plotly") })
    shiny::outputOptions(output, "is_plotly_active", suspendWhenHidden = FALSE)

    output$is_iheatmapr_active <- shiny::reactive({ identical(active_kind(), "iheatmapr") })
    shiny::outputOptions(output, "is_iheatmapr_active", suspendWhenHidden = FALSE)

    output$has_history <- shiny::reactive({ length(plot_history()) > 1L })
    shiny::outputOptions(output, "has_history", suspendWhenHidden = FALSE)

    output$has_table <- shiny::reactive({ !is.null(table_data()) })
    shiny::outputOptions(output, "has_table", suspendWhenHidden = FALSE)

    output$has_dataset <- shiny::reactive({
      !is.null(local_pgx) && !is.null(local_pgx())
    })
    shiny::outputOptions(output, "has_dataset", suspendWhenHidden = FALSE)

    ## --- PlotModuleServer instances (one per plotlib) ---
    PlotModuleServer(
      "evidence_ggplot",
      plotlib = "ggplot",
      func = function() {
        rec <- active_record()
        shiny::validate(shiny::need(
          !is.null(rec) && identical(rec$kind, "ggplot"), ""
        ))
        rec$plot
      },
      pdf.width = 8,
      pdf.height = 6
    )

    PlotModuleServer(
      "evidence_plotly",
      plotlib = "plotly",
      func = function() {
        rec <- active_record()
        shiny::validate(shiny::need(
          !is.null(rec) && identical(rec$kind, "plotly"), ""
        ))
        rec$plot
      },
      pdf.width = 8,
      pdf.height = 6
    )

    PlotModuleServer(
      "evidence_iheatmapr",
      plotlib = "iheatmapr",
      func = function() {
        rec <- active_record()
        shiny::validate(shiny::need(
          !is.null(rec) && identical(rec$kind, "iheatmapr"), ""
        ))
        rec$plot
      }
    )

    ## --- Carousel rendering ---
    output$plot_carousel <- shiny::renderUI({
      hist <- plot_history()
      if (length(hist) <= 1L) return(NULL)
      idx <- active_plot_idx()

      items <- lapply(seq_along(hist), function(i) {
        rec <- hist[[i]]
        is_active <- identical(i, idx)
        border_style <- if (is_active) "border-color: #0d6efd;" else "border-color: #dee2e6;"
        bg_style <- if (is_active) "background: rgba(13,110,253,0.08);" else "background: #f8f9fa;"
        time_label <- format(rec$timestamp, "%H:%M")

        shiny::tags$div(
          style = paste0(
            "display: inline-block; cursor: pointer; padding: 4px 10px;",
            "margin-right: 4px; border-radius: 6px; border: 2px solid;",
            "min-width: 70px; text-align: center; vertical-align: top;",
            border_style, bg_style
          ),
          onclick = sprintf(
            "Shiny.setInputValue('%s', %d, {priority: 'event'})",
            ns("select_plot"), i
          ),
          shiny::tags$div(
            style = "font-size: 0.75em; font-weight: 600;",
            rec$label
          ),
          shiny::tags$div(
            style = "font-size: 0.65em; color: #888;",
            time_label
          )
        )
      })

      shiny::tagList(items)
    })

    ## --- Carousel click handler ---
    shiny::observeEvent(input$select_plot, {
      idx <- input$select_plot
      hist <- plot_history()
      if (!is.null(idx) && idx >= 1L && idx <= length(hist)) {
        active_plot_idx(idx)
      }
    })

    ## --- Table output ---
    output$evidence_table <- DT::renderDataTable({
      df <- table_data()
      shiny::validate(shiny::need(!is.null(df), "No results to display yet."))
      DT::datatable(df,
        rownames = FALSE,
        options = list(pageLength = 10, dom = "ftip", scrollX = TRUE)
      )
    })

    ## Task 4.2 — Dataset context display
    output$dataset_info <- shiny::renderUI({
      shiny::validate(shiny::need(
        !is.null(local_pgx) && !is.null(local_pgx()),
        "No dataset loaded yet."
      ))
      pgx <- local_pgx()

      n_samples <- if (!is.null(pgx$X)) ncol(pgx$X) else "?"
      n_genes <- if (!is.null(pgx$X)) nrow(pgx$X) else "?"
      organism <- if (!is.null(pgx$organism)) pgx$organism else "unknown"

      contrasts <- "none"
      if (!is.null(pgx$contrasts)) {
        ct <- colnames(pgx$contrasts)
        if (length(ct) > 3) {
          contrasts <- paste(c(ct[1:3], paste0("+ ", length(ct) - 3, " more")), collapse = ", ")
        } else {
          contrasts <- paste(ct, collapse = ", ")
        }
      }

      shiny::tags$div(
        style = "font-size: 0.85em; line-height: 1.6;",
        shiny::tags$div(shiny::icon("dna"), shiny::strong(" Organism: "), organism),
        shiny::tags$div(shiny::icon("vials"), shiny::strong(" Samples: "), n_samples),
        shiny::tags$div(shiny::icon("bars"), shiny::strong(" Genes: "), n_genes),
        shiny::tags$div(shiny::icon("arrows-left-right"), shiny::strong(" Contrasts: "), contrasts)
      )
    })

    ## --- Return API for copilot_server.R ---
    list(
      append_plot = function(record) {
        hist <- plot_history()
        hist[[length(hist) + 1L]] <- record
        plot_history(hist)
        active_plot_idx(length(hist))
      },
      clear_plots = function() {
        plot_history(list())
        active_plot_idx(NULL)
      },
      update_table = function(df) { table_data(df) },
      clear_table = function() { table_data(NULL) }
    )
  })
}
