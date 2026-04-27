#' Evidence Panel — Server
#'
#' Displays agent-generated plots (ggplot, plotly, iheatmapr), a results
#' table, and a clickable plot-history carousel.  Uses direct Shiny
#' renderers (renderImage / renderPlotly / renderIheatmap) for minimal
#' flush-cycle overhead.

#' @param id Character scalar, Shiny module namespace.
#' @param local_pgx Reactive returning the current PGX object, or `NULL`.
#' @return A named list with the following callables:
#'   \describe{
#'     \item{`append_plot(record)`}{Add a plot record to the viewer and carousel.}
#'     \item{`clear_plots()`}{Reset the viewer, carousel, and plot cache.}
#'     \item{`update_table(df)`}{Display a data frame in the results table.}
#'     \item{`clear_table()`}{Remove the results table.}
#'     \item{`plot_rendered`}{Reactive wrapping `input$plot_rendered` (client
#'       confirmation that a plot frame was painted).}
#'   }
copilot_panel_evidence_server <- function(id, local_pgx = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## --- Plot history state ---
    plot_history    <- shiny::reactiveVal(list())
    plot_cache      <- new.env(parent = emptyenv())
    active_plot_idx <- shiny::reactiveVal(NULL)

    ## Per-kind reactiveVals: only the matching renderer fires on update
    active_ggplot_record    <- shiny::reactiveVal(NULL)
    active_plotly_record    <- shiny::reactiveVal(NULL)
    active_iheatmapr_record <- shiny::reactiveVal(NULL)

    log_render_handoff <- function(rec, stage) {
      if (is.null(rec)) return(invisible(NULL))
      now <- Sys.time()
      queued_secs <- if (!is.null(rec$queued_at)) {
        round(as.numeric(difftime(now, rec$queued_at, units = "secs")), 3)
      } else {
        NA_real_
      }
      built_secs <- if (!is.null(rec$built_at)) {
        round(as.numeric(difftime(now, rec$built_at, units = "secs")), 3)
      } else {
        NA_real_
      }
      appended_secs <- if (!is.null(rec$appended_at)) {
        round(as.numeric(difftime(now, rec$appended_at, units = "secs")), 3)
      } else {
        NA_real_
      }
      if (getOption("copilot.trace", FALSE)) {
        dbg(
          "[CopilotEvidence] ", stage,
          " plot_id=", rec$id %||% "<none>",
          " kind=", rec$kind %||% "<none>",
          " label=", rec$label %||% "<none>",
          " queue_to_now=", queued_secs, "s",
          " build_to_now=", built_secs, "s",
          " append_to_now=", appended_secs, "s"
        )
      }
      invisible(NULL)
    }

    table_data <- shiny::reactiveVal(NULL)

    shiny::observe({
      rec <- active_ggplot_record() %||% active_plotly_record() %||% active_iheatmapr_record()
      if (is.null(rec)) return()
      if (getOption("copilot.trace", FALSE)) {
        dbg(
          "[CopilotEvidence] active_record",
          " plot_id=", rec$id %||% "<none>",
          " kind=", rec$kind %||% "<none>",
          " label=", rec$label %||% "<none>"
        )
      }
    })

    ## --- Conditional panel flags ---
    output$has_plot <- shiny::reactive({
      !is.null(active_ggplot_record()) || !is.null(active_plotly_record()) || !is.null(active_iheatmapr_record())
    })
    shiny::outputOptions(output, "has_plot", suspendWhenHidden = FALSE)

    output$is_ggplot_active <- shiny::reactive({ !is.null(active_ggplot_record()) })
    shiny::outputOptions(output, "is_ggplot_active", suspendWhenHidden = FALSE)

    output$is_plotly_active <- shiny::reactive({ !is.null(active_plotly_record()) })
    shiny::outputOptions(output, "is_plotly_active", suspendWhenHidden = FALSE)

    output$is_iheatmapr_active <- shiny::reactive({ !is.null(active_iheatmapr_record()) })
    shiny::outputOptions(output, "is_iheatmapr_active", suspendWhenHidden = FALSE)

    output$has_history <- shiny::reactive({ length(plot_history()) > 1L })
    shiny::outputOptions(output, "has_history", suspendWhenHidden = FALSE)

    output$has_table <- shiny::reactive({ !is.null(table_data()) })
    shiny::outputOptions(output, "has_table", suspendWhenHidden = FALSE)

    output$has_dataset <- shiny::reactive({
      !is.null(local_pgx) && !is.null(local_pgx())
    })
    shiny::outputOptions(output, "has_dataset", suspendWhenHidden = FALSE)

    ## --- Direct render outputs ---
    output$evidence_ggplot <- shiny::renderImage({
      rec <- active_ggplot_record()
      shiny::req(rec, rec$prerendered_path, file.exists(rec$prerendered_path))
      log_render_handoff(rec, "render_handoff")
      list(src = rec$prerendered_path, contentType = "image/png",
           width = "100%", alt = "Evidence Plot")
    }, deleteFile = FALSE)

    ## Eliminates hidden-output resume latency when the evidence panel first appears.
    output$evidence_plotly <- plotly::renderPlotly({
      rec <- active_plotly_record()
      shiny::req(rec)
      log_render_handoff(rec, "render_handoff")
      rec$plot
    })
    shiny::outputOptions(output, "evidence_plotly", suspendWhenHidden = FALSE)

    ## renderIheatmap internally calls to_widget() which only has an S4

    ## method for "Iheatmap".  omicsplots::pgx.plot_heatmap already converts
    ## to an "iheatmapr" htmlwidget via to_widget(), so renderIheatmap would
    ## try to convert *again* and fail with:
    ##   unable to find an inherited method for function 'to_widget'
    ##     for signature '"iheatmapr"'
    ## Use shinyRenderWidget directly so we can handle both cases.
    output$evidence_iheatmapr <- htmlwidgets::shinyRenderWidget(
      quote({
        rec <- active_iheatmapr_record()
        shiny::req(rec)
        log_render_handoff(rec, "render_handoff")
        p <- rec$plot
        if (methods::is(p, "Iheatmap")) iheatmapr::to_widget(p) else p
      }),
      iheatmapr::iheatmaprOutput,
      environment(),
      quoted = TRUE
    )

    ## --- plot_rendered event handler ---
    shiny::observeEvent(input$plot_rendered, {
      rec <- active_ggplot_record() %||% active_plotly_record() %||% active_iheatmapr_record()
      if (is.null(rec)) return()
      evt <- input$plot_rendered
      now <- Sys.time()
      appended_secs <- if (!is.null(rec$appended_at)) {
        round(as.numeric(difftime(now, rec$appended_at, units = "secs")), 3)
      } else {
        NA_real_
      }
      built_secs <- if (!is.null(rec$built_at)) {
        round(as.numeric(difftime(now, rec$built_at, units = "secs")), 3)
      } else {
        NA_real_
      }
      dbg(
        "[CopilotEvidence] client_rendered",
        " output_id=", evt$output_id %||% "<unknown>",
        " plot_id=", rec$id %||% "<none>",
        " kind=", rec$kind %||% "<none>",
        " build_to_client=", built_secs, "s",
        " append_to_client=", appended_secs, "s",
        " client_ts=", evt$ts %||% "<none>"
      )
    }, ignoreNULL = TRUE)

    ## --- Carousel rendering ---
    output$plot_carousel <- shiny::renderUI({
      hist <- plot_history()
      if (length(hist) <= 1L) return(NULL)
      idx <- active_plot_idx()

      items <- lapply(seq_along(hist), function(i) {
        rec <- hist[[i]]
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
      idx  <- input$select_plot
      hist <- plot_history()
      if (is.null(idx) || idx < 1L || idx > length(hist)) return()
      active_plot_idx(idx)
      rec <- hist[[idx]]
      plot_id <- rec$id %||% ""
      if (nzchar(plot_id) && exists(plot_id, envir = plot_cache, inherits = FALSE)) {
        rec$plot <- get(plot_id, envir = plot_cache, inherits = FALSE)
      }
      kind <- rec$kind
      if (identical(kind, "ggplot")) {
        active_ggplot_record(rec)
        active_plotly_record(NULL)
        active_iheatmapr_record(NULL)
      } else if (identical(kind, "plotly")) {
        active_plotly_record(rec)
        active_ggplot_record(NULL)
        active_iheatmapr_record(NULL)
      } else if (identical(kind, "iheatmapr")) {
        active_iheatmapr_record(rec)
        active_ggplot_record(NULL)
        active_plotly_record(NULL)
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

    ## --- Dataset context display ---
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
        if (length(ct) > 3) {
          contrasts <- paste(c(ct[1:3], paste0("+ ", length(ct) - 3, " more")), collapse = ", ")
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

    ## --- Return API for copilot_server.R ---
    list(
      append_plot = function(record) {
        appended_at <- Sys.time()
        plot_id <- record$id %||% paste0(
          "copilot-plot-",
          format(appended_at, "%Y%m%d%H%M%S"),
          "-",
          sprintf("%06d", sample.int(999999L, 1L))
        )
        ## Cache plot object for carousel navigation; keep it alive
        ## in the active record so renderers can read rec$plot.
        assign(plot_id, record$plot, envir = plot_cache)
        record$id         <- plot_id
        record$appended_at <- appended_at

        kind <- record$kind
        if (identical(kind, "ggplot")) {
          active_ggplot_record(record)
          active_plotly_record(NULL)
          active_iheatmapr_record(NULL)
        } else if (identical(kind, "plotly")) {
          active_plotly_record(record)
          active_ggplot_record(NULL)
          active_iheatmapr_record(NULL)
        } else if (identical(kind, "iheatmapr")) {
          active_iheatmapr_record(record)
          active_ggplot_record(NULL)
          active_plotly_record(NULL)
        }

        dbg(
          "[CopilotEvidence] append_plot",
          " plot_id=", plot_id,
          " kind=", record$kind %||% "<none>",
          " label=", record$label %||% "<none>"
        )

        queued_secs <- if (!is.null(record$queued_at)) {
          round(as.numeric(difftime(appended_at, record$queued_at, units = "secs")), 3)
        } else {
          NA_real_
        }
        built_secs <- if (!is.null(record$built_at)) {
          round(as.numeric(difftime(appended_at, record$built_at, units = "secs")), 3)
        } else {
          NA_real_
        }

        ## Defer history update past the plot flush — carousel renderUI would compete
        ## for the same flush cycle.
        ## Drop the plot ref from the history copy to avoid duplication with cache.
        history_record <- record
        history_record$plot <- NULL
        session$onFlushed(function() {
          hist <- shiny::isolate(plot_history())
          hist[[length(hist) + 1L]] <- history_record
          plot_history(hist)
          active_plot_idx(length(hist))
          if (getOption("copilot.trace", FALSE)) {
            dbg(
              "[CopilotEvidence] post_flush",
              " plot_id=", plot_id,
              " kind=", history_record$kind %||% "<none>",
              " label=", history_record$label %||% "<none>",
              " history_len=", length(hist),
              " queue_to_flush=", queued_secs, "s",
              " build_to_flush=", built_secs, "s"
            )
          }
        }, once = TRUE)
      },
      clear_plots = function() {
        ## Clean up pre-rendered temp PNGs before discarding history
        hist <- shiny::isolate(plot_history())
        png_paths <- vapply(hist, function(r) r$prerendered_path %||% NA_character_, character(1))
        png_paths <- png_paths[!is.na(png_paths)]
        if (length(png_paths)) copilot_prerender_cleanup(png_paths)
        cache_keys <- ls(envir = plot_cache, all.names = TRUE)
        if (length(cache_keys)) rm(list = cache_keys, envir = plot_cache)
        plot_history(list())
        active_plot_idx(NULL)
        active_ggplot_record(NULL)
        active_plotly_record(NULL)
        active_iheatmapr_record(NULL)
      },
      update_table = function(df) { table_data(df) },
      clear_table  = function()   { table_data(NULL) },
      ## Expose plot_rendered input to parent module for save deferral
      plot_rendered = shiny::reactive({ input$plot_rendered })
    )
  })
}
