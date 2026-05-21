# CopilotEvidenceServer.R — Evidence module server
#
# Displays agent-generated plots (ggplot, plotly, iheatmapr), a results table,
# and a clickable plot-history carousel. Single `active_artifact` reactiveVal
# replaces the three per-kind reactiveVals from the old implementation.
#
# Entry point for new evidence records is `$append_artifact(record)` — the ONLY
# write path into plot_history. All other callers (run controller, restore
# controller) call `$clear_plots()` or `$update_table()`, never write to
# plot_history directly.

# ---- Null-coalescing operator (local) ----
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

#' Copilot Evidence Module Server
#'
#' Manages the artifact display panel: plot renderers, carousel, results table,
#' and dataset-context card. Provides a public API consumed by the orchestrator
#' and the run-bindings plot_callback.
#'
#' @param id Character scalar — Shiny module namespace.
#' @param local_pgx Reactive returning the current PGX list, or NULL (optional,
#'   used only for the dataset-context card).
#' @return Named list: append_artifact, clear_plots, update_table, clear_table,
#'   clear (alias for clear_plots), active (reactive returning active_artifact()).
#' @export
CopilotEvidenceServer <- function(id, local_pgx = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- Internal state ----
    active_artifact <- shiny::reactiveVal(NULL)
    # shape: list(
    #   kind             = "ggplot" | "plotly" | "iheatmapr",
    #   plot             = <R plot object>,
    #   prerendered_path = character(1) | NULL,
    #   plot_type        = character(1),
    #   args             = list,
    #   artifact         = <ArtifactRecord | NULL>,
    #   label            = character(1),
    #   timestamp        = POSIXct(1)
    # )

    plot_history    <- shiny::reactiveVal(list())
    active_idx      <- shiny::reactiveVal(NULL)
    table_data      <- shiny::reactiveVal(NULL)

    # ---- append_artifact ----
    # Called by the plot_callback closure in copilot_bindings.R. Receives a
    # fully-rendered record list built by the bindings closure (see B1).
    .append_artifact <- function(record) {
      history <- c(plot_history(), list(record))
      plot_history(history)
      idx <- length(history)
      active_idx(idx)
      active_artifact(record)
      invisible(NULL)
    }

    # ---- clear_plots ----
    .clear_plots <- function() {
      # Clean up temp PNG files before discarding history
      paths <- vapply(
        plot_history(),
        function(r) r$prerendered_path %||% NA_character_,
        character(1L)
      )
      paths <- paths[!is.na(paths)]
      if (length(paths)) copilot_prerender_cleanup(paths)
      plot_history(list())
      active_artifact(NULL)
      active_idx(NULL)
      invisible(NULL)
    }

    # ---- Conditional panel outputs (suspendWhenHidden = FALSE required) ----
    output$has_plot <- shiny::reactive({
      !is.null(active_artifact())
    })
    shiny::outputOptions(output, "has_plot", suspendWhenHidden = FALSE)

    output$is_ggplot_active <- shiny::reactive({
      isTRUE(active_artifact()$kind == "ggplot")
    })
    shiny::outputOptions(output, "is_ggplot_active", suspendWhenHidden = FALSE)

    output$is_plotly_active <- shiny::reactive({
      isTRUE(active_artifact()$kind == "plotly")
    })
    shiny::outputOptions(output, "is_plotly_active", suspendWhenHidden = FALSE)

    output$is_iheatmapr_active <- shiny::reactive({
      isTRUE(active_artifact()$kind == "iheatmapr")
    })
    shiny::outputOptions(output, "is_iheatmapr_active", suspendWhenHidden = FALSE)

    output$has_history <- shiny::reactive({
      length(plot_history()) > 1L
    })
    shiny::outputOptions(output, "has_history", suspendWhenHidden = FALSE)

    output$has_table <- shiny::reactive({
      !is.null(table_data())
    })
    shiny::outputOptions(output, "has_table", suspendWhenHidden = FALSE)

    output$has_dataset <- shiny::reactive({
      !is.null(local_pgx) && !is.null(local_pgx())
    })
    shiny::outputOptions(output, "has_dataset", suspendWhenHidden = FALSE)

    # ---- Render outputs ----

    output$evidence_ggplot <- shiny::renderImage({
      rec <- active_artifact()
      shiny::req(
        !is.null(rec),
        identical(rec$kind, "ggplot"),
        !is.null(rec$prerendered_path),
        file.exists(rec$prerendered_path)
      )
      list(
        src         = rec$prerendered_path,
        contentType = "image/png",
        width       = "100%",
        alt         = "Evidence Plot"
      )
    }, deleteFile = FALSE)

    output$evidence_plotly <- plotly::renderPlotly({
      rec <- active_artifact()
      shiny::req(!is.null(rec), identical(rec$kind, "plotly"))
      rec$plot
    })
    shiny::outputOptions(output, "evidence_plotly", suspendWhenHidden = FALSE)

    ## renderIheatmap internally calls to_widget() which only has an S4
    ## method for "Iheatmap". omicsplots::pgx.plot_heatmap already converts
    ## to an "iheatmapr" htmlwidget via to_widget(), so renderIheatmap would
    ## try to convert *again* and fail with:
    ##   unable to find an inherited method for function 'to_widget'
    ##     for signature '"iheatmapr"'
    ## Use shinyRenderWidget directly so we can handle both cases:
    ## - an "Iheatmap" S4 object  -> call to_widget() once here
    ## - an "iheatmapr" htmlwidget already converted -> pass through as-is
    output$evidence_iheatmapr <- htmlwidgets::shinyRenderWidget(
      quote({
        rec <- active_artifact()
        shiny::req(!is.null(rec), identical(rec$kind, "iheatmapr"))
        p <- rec$plot
        if (methods::is(p, "Iheatmap")) iheatmapr::to_widget(p) else p
      }),
      iheatmapr::iheatmaprOutput,
      environment(),
      quoted = TRUE
    )

    # ---- Download handler ----
    # Dispatches on artifact kind:
    #   ggplot    -> PNG copy of the pre-rendered file
    #   plotly    -> standalone HTML via htmlwidgets::saveWidget
    #   iheatmapr -> standalone HTML via htmlwidgets::saveWidget
    .download_ext <- function(kind) {
      switch(kind, ggplot = "png", plotly = "html", iheatmapr = "html", "bin")
    }
    output$evidence_download <- shiny::downloadHandler(
      filename = function() {
        rec <- active_artifact()
        kind <- if (!is.null(rec)) rec$kind else "plot"
        ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
        sprintf("copilot-%s-%s.%s", kind, ts, .download_ext(kind))
      },
      content = function(file) {
        rec <- active_artifact()
        shiny::req(rec)
        if (identical(rec$kind, "ggplot")) {
          shiny::req(!is.null(rec$prerendered_path), file.exists(rec$prerendered_path))
          file.copy(rec$prerendered_path, file, overwrite = TRUE)
        } else if (identical(rec$kind, "plotly")) {
          htmlwidgets::saveWidget(rec$plot, file, selfcontained = TRUE)
        } else if (identical(rec$kind, "iheatmapr")) {
          p <- rec$plot
          if (methods::is(p, "Iheatmap")) p <- iheatmapr::to_widget(p)
          htmlwidgets::saveWidget(p, file, selfcontained = TRUE)
        }
      }
    )

    # ---- Maximize modal slot ----
    # The modal itself is declared statically in CopilotEvidenceUI (Bootstrap
    # data-bs-toggle, no R round-trip). Here we emit ONLY the active kind's
    # output element so plotly/iheatmapr don't mount inside a display:none
    # conditionalPanel (which would freeze them at zero width and produce the
    # squeezed-column render in the modal).
    output$evidence_modal_slot <- shiny::renderUI({
      rec <- active_artifact()
      if (is.null(rec)) return(NULL)
      mh <- "calc(80vh - 100px)"
      switch(rec$kind,
        ggplot    = shiny::imageOutput(ns("evidence_ggplot_modal"), width = "100%", height = mh),
        plotly    = plotly::plotlyOutput(ns("evidence_plotly_modal"), width = "100%", height = mh),
        iheatmapr = iheatmapr::iheatmaprOutput(ns("evidence_iheatmapr_modal"), width = "100%", height = mh),
        NULL
      )
    })

    output$evidence_ggplot_modal <- shiny::renderImage({
      rec <- active_artifact()
      shiny::req(
        !is.null(rec),
        identical(rec$kind, "ggplot"),
        !is.null(rec$prerendered_path),
        file.exists(rec$prerendered_path)
      )
      list(
        src         = rec$prerendered_path,
        contentType = "image/png",
        width       = "100%",
        alt         = "Evidence Plot (maximized)"
      )
    }, deleteFile = FALSE)

    output$evidence_plotly_modal <- plotly::renderPlotly({
      rec <- active_artifact()
      shiny::req(!is.null(rec), identical(rec$kind, "plotly"))
      # Clear any width/height baked in by plotly_build() at inline-card size.
      # Belt-and-braces: the modal_slot uiOutput already avoids the zero-width
      # mount, but cached dimensions in $x$layout would still override autosize.
      p <- rec$plot
      if (!is.null(p$x) && !is.null(p$x$layout)) {
        p$x$layout$width  <- NULL
        p$x$layout$height <- NULL
      }
      p$width  <- NULL
      p$height <- NULL
      p
    })

    output$evidence_iheatmapr_modal <- htmlwidgets::shinyRenderWidget(
      quote({
        rec <- active_artifact()
        shiny::req(!is.null(rec), identical(rec$kind, "iheatmapr"))
        p <- rec$plot
        if (methods::is(p, "Iheatmap")) iheatmapr::to_widget(p) else p
      }),
      iheatmapr::iheatmaprOutput,
      environment(),
      quoted = TRUE
    )

    # ---- Carousel rendering ----
    output$plot_carousel <- shiny::renderUI({
      hist <- plot_history()
      if (length(hist) <= 1L) return(NULL)
      idx <- active_idx()

      items <- lapply(seq_along(hist), function(i) {
        rec          <- hist[[i]]
        is_active    <- identical(i, idx)
        border_style <- if (is_active) "border-color: #0d6efd;" else "border-color: #dee2e6;"
        bg_style     <- if (is_active) "background: rgba(13,110,253,0.08);" else "background: #f8f9fa;"
        time_label   <- format(rec$timestamp, "%H:%M")

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
            rec$label %||% rec$kind
          ),
          shiny::tags$div(
            style = "font-size: 0.65em; color: #888;",
            time_label
          )
        )
      })

      shiny::tagList(items)
    })

    # ---- Carousel click handler ----
    shiny::observeEvent(input$select_plot, {
      idx  <- input$select_plot
      hist <- plot_history()
      if (is.null(idx) || idx < 1L || idx > length(hist)) return()
      active_idx(idx)
      active_artifact(hist[[idx]])
    })

    # ---- Table output ----
    output$evidence_table <- DT::renderDataTable({
      df <- table_data()
      shiny::validate(shiny::need(!is.null(df), "No results to display yet."))
      DT::datatable(df,
        rownames = FALSE,
        options  = list(pageLength = 10, dom = "ftip", scrollX = TRUE)
      )
    })

    # ---- Dataset context card ----
    output$dataset_info <- shiny::renderUI({
      shiny::validate(shiny::need(
        !is.null(local_pgx) && !is.null(local_pgx()),
        "No dataset loaded yet."
      ))
      pgx <- local_pgx()

      n_samples <- if (!is.null(pgx$X)) ncol(pgx$X) else "?"
      n_genes   <- if (!is.null(pgx$X)) nrow(pgx$X) else "?"
      organism  <- if (!is.null(pgx$organism)) pgx$organism else "unknown"

      contrasts <- "none"
      if (!is.null(pgx$contrasts)) {
        ct <- colnames(pgx$contrasts)
        if (length(ct) > 3L) {
          contrasts <- paste(c(ct[1:3], paste0("+ ", length(ct) - 3L, " more")), collapse = ", ")
        } else {
          contrasts <- paste(ct, collapse = ", ")
        }
      }

      shiny::tags$div(
        style = "font-size: 0.85em; line-height: 1.6;",
        shiny::tags$div(shiny::icon("dna"),               shiny::strong(" Organism: "),  organism),
        shiny::tags$div(shiny::icon("vials"),             shiny::strong(" Samples: "),   n_samples),
        shiny::tags$div(shiny::icon("bars"),              shiny::strong(" Genes: "),     n_genes),
        shiny::tags$div(shiny::icon("arrows-left-right"), shiny::strong(" Contrasts: "), contrasts)
      )
    })

    # ---- Public API ----
    list(
      # Primary write path — called by plot_callback in copilot_bindings.R
      append_artifact = .append_artifact,
      # Called by run/restore controller on new chat or restore
      clear_plots     = .clear_plots,
      clear           = .clear_plots,   # alias used by run controller (.do_reset)
      # Table management — called by chat server when agent emits a table
      update_table    = function(df) { table_data(df) },
      clear_table     = function()   { table_data(NULL) },
      # Reactive for orchestrator / save controller consumption
      active          = shiny::reactive(active_artifact())
    )
  })
}
