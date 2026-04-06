#' Evidence Panel Sub-module
#'
#' Right panel: dataset context + plot + table slots for agent output

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

    ## Plot card
    shiny::conditionalPanel(
      condition = paste0("output['", ns("has_plot"), "']"),
      bslib::card(
        class = "mb-2",
        bslib::card_header("Plot", class = "py-1 px-2"),
        bslib::card_body(
          class = "p-1",
          shiny::uiOutput(ns("plot_container"))
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
  if (inherits(plot_obj, "plotly")) {
    return("plotly")
  }

  if (inherits(plot_obj, "iheatmapr")) {
    return("iheatmapr")
  }

  if (inherits(plot_obj, "ggplot")) {
    return("ggplot")
  }

  NULL
}

copilot_panel_evidence_server <- function(id, local_pgx = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    plot_obj <- shiny::reactiveVal(NULL)
    plot_kind <- shiny::reactiveVal(NULL)
    table_data <- shiny::reactiveVal(NULL)

    ## Conditional panel flags
    output$has_plot <- shiny::reactive({ !is.null(plot_obj()) })
    shiny::outputOptions(output, "has_plot", suspendWhenHidden = FALSE)

    output$has_table <- shiny::reactive({ !is.null(table_data()) })
    shiny::outputOptions(output, "has_table", suspendWhenHidden = FALSE)

    output$has_dataset <- shiny::reactive({
      !is.null(local_pgx) && !is.null(local_pgx())
    })
    shiny::outputOptions(output, "has_dataset", suspendWhenHidden = FALSE)

    output$plot_container <- shiny::renderUI({
      kind <- plot_kind()
      shiny::validate(shiny::need(!is.null(kind), "No plot available yet."))

      if (identical(kind, "plotly")) {
        return(plotly::plotlyOutput(session$ns("evidence_plotly"), height = "400px"))
      }

      if (identical(kind, "iheatmapr")) {
        return(iheatmapr::iheatmaprOutput(session$ns("evidence_iheatmap"), height = "400px"))
      }

      shiny::plotOutput(session$ns("evidence_plot"), height = "400px")
    })

    output$evidence_plot <- shiny::renderPlot({
      shiny::validate(shiny::need(
        identical(plot_kind(), "ggplot"),
        "This plot cannot be shown in the ggplot renderer."
      ))
      plot_value <- plot_obj()
      shiny::validate(shiny::need(!is.null(plot_value), "No plot available yet."))
      print(plot_value)
    })

    output$evidence_plotly <- plotly::renderPlotly({
      shiny::validate(shiny::need(
        identical(plot_kind(), "plotly"),
        "This plot cannot be shown in the plotly renderer."
      ))
      plot_value <- plot_obj()
      shiny::validate(shiny::need(!is.null(plot_value), "No plot available yet."))
      plot_value
    })

    output$evidence_iheatmap <- iheatmapr::renderIheatmap({
      shiny::validate(shiny::need(
        identical(plot_kind(), "iheatmapr"),
        "This plot cannot be shown in the iheatmap renderer."
      ))
      plot_value <- plot_obj()
      shiny::validate(shiny::need(!is.null(plot_value), "No plot available yet."))
      plot_value
    })

    ## Table output
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

      ## Get contrasts
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

    ## Return update functions for the server to call
    list(
      update_plot = function(value) {
        kind <- copilot_detect_plot_kind(value)
        if (is.null(kind)) {
          stop("Unsupported plot object for evidence panel.", call. = FALSE)
        }
        plot_obj(value)
        plot_kind(kind)
      },
      clear_plot = function() {
        plot_obj(NULL)
        plot_kind(NULL)
      },
      update_table = function(df) { table_data(df) },
      clear_table = function() { table_data(NULL) }
    )
  })
}
