##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' WGCNA AI Report Tab UI
#'
#' Entry point for the AI Report tab body.
#' @param id Module namespace ID
#' @return Shiny UI element
wgcna_ai_report_ui <- function(id) {
  ns <- shiny::NS(id)
  controls_ns <- shiny::NS(ns("controls"))

  ## "Show prompt" checkbox: namespaced to controls so the controls server
  ## can read it via input$show_prompt, but rendered in the card hamburger menu.
  show_prompt_input <- shiny::checkboxInput(
    controls_ns("show_prompt"),
    "Show prompt",
    FALSE
  )

  ## Infographic style/layout controls: namespaced to controls so the controls
  ## server can read them, but rendered in the image card's hamburger menu.
  ## The include_infographic checkbox lives in the sidebar (ai_report_controls.R).
  infographic_options <- shiny::tagList(
    shiny::tags$hr(),
    shiny::selectInput(
      controls_ns("image_style"),
      "Infographic Style:",
      choices = NULL,
      width = "100%"
    ),
    shiny::radioButtons(
      controls_ns("image_blocks"),
      "Layout:",
      choices = c("1 Panel" = "1", "2 Panels" = "2", "3 Panels" = "3"),
      selected = "1",
      inline = TRUE
    )
  )

  ai_report_layout_ui(
    ns("layout"),
    text_title = "AI Report",
    diagram_title = "Module Diagram",
    infographic_title = "Graphical Abstract",
    text_options = show_prompt_input,
    infographic_options = infographic_options
  )
}

#' WGCNA AI Report Sidebar UI
#'
#' Entry point for the AI Report sidebar controls.
#' @param id Module namespace ID
#' @return Shiny UI element
wgcna_ai_report_inputs_ui <- function(id) {
  ns <- shiny::NS(id)
  ai_report_controls_ui(ns("controls"))
}

#' AI Report Layout - Always 3 Cards
#'
#' Creates a fixed 3-card layout: text (left 50%) + diagram and infographic
#' (right 50% stacked). The text card is shared between summary and report
#' modes — the parent server decides what content feeds into it.
#'
#' @param id Module namespace ID
#' @param text_title Title for text panel (default "AI Report")
#' @param diagram_title Title for diagram panel (default "Board Diagram")
#' @param infographic_title Title for infographic panel (default "Graphical Abstract")
#'
#' @return Shiny UI tagList
ai_report_layout_ui <- function(id,
                                text_title = "AI Report",
                                diagram_title = "Board Diagram",
                                infographic_title = "Graphical Abstract",
                                text_options = NULL,
                                infographic_options = NULL) {
  ns <- shiny::NS(id)

  bslib::layout_columns(
    col_widths = 12,
    height = "calc(100vh - 180px)",
    row_heights = c("auto", 1, "auto"),

    # Disclaimer banner
    bs_alert(
      HTML(paste0(
        "\u26a0\ufe0f <b>Disclaimer:</b> This page contains AI-generated content. ",
        "Please verify important information independently."
      )),
      translate = FALSE
    ),

    # Always 3 cards: text (left) | diagram + infographic (right stacked)
    bslib::layout_columns(
      col_widths = c(6, 6),
      height = "calc(100vh - 220px)",

      # LEFT: Text card (shared between summary and report modes)
      PlotModuleUI(
        ns("report_text"),
        outputFunc = htmlOutput,
        title = text_title,
        options = text_options,
        caption = "AI-generated analysis content",
        height = c("100%", TABLE_HEIGHT_MODAL),
        width = c("auto", "100%"),
        download.fmt = c("pdf", "docx", "md")
      ),

      # RIGHT: Diagram + Infographic stacked
      bslib::layout_columns(
        col_widths = 12,
        row_heights = c(1, 1),

        AiDiagramCardUI(
          ns("diagram"),
          title = diagram_title,
          caption = "AI-generated module diagram",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        ),

        AiImageCardUI(
          ns("infographic"),
          title = infographic_title,
          caption = "AI-generated graphical abstract",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%"),
          extra_options = infographic_options
        )
      )
    ),

    ## ── Deep Report evidence trail (only rendered when the last fired
    ## generation was a Deep Report and the agent emitted tool calls).
    ## Re-uses the collapsible <details>/<summary> pattern from
    ## board.copilot (.format_tool_request in copilot_agent.R) without
    ## pulling the full plot-history evidence panel — Deep Report
    ## produces no plots, only tool-call summaries plus the agent's
    ## final bibliography (which is rendered inline in the main text
    ## card under the ## Bibliography heading).
    shiny::conditionalPanel(
      condition = paste0("output['", ns("has_deep_trace"), "']"),
      bslib::card(
        class = "mt-2",
        bslib::card_header(
          shiny::tags$details(
            shiny::tags$summary(
              shiny::tags$strong("Evidence trail (Deep Report)"),
              shiny::tags$small(class = "text-muted ms-2",
                                "Tool calls used to ground citations")
            ),
            shiny::div(
              class = "mt-2",
              shiny::uiOutput(ns("deep_trace"))
            )
          )
        )
      )
    )
  )
}

#' AI Report Layout - Server (mode-agnostic)
#'
#' Renders the 3 cards with whatever content the parent provides.
#' No mode awareness — the parent server handles mode-switching logic.
#'
#' @param id Module namespace ID
#' @param text_reactive Reactive returning markdown text string.
#'   Fed by parent: summary text in summary mode, report text in report mode.
#' @param diagram_result_reactive Deprecated; unused
#' @param infographic_reactive Deprecated; unused
#' @param watermark Logical; add watermark to outputs
#'
#' @return NULL (no return values needed)
ai_report_layout_server <- function(id,
                                    text_reactive,
                                    diagram_result_reactive = NULL,
                                    infographic_reactive = NULL,
                                    deep_turns_reactive = NULL,
                                    last_deep_reactive = NULL,
                                    watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    # ============================================================
    # CARD 1: Text (left panel — shared between summary and report)
    # ============================================================
    text.RENDER <- function() {
      text <- text_reactive()
      shiny::validate(shiny::need(!is.null(text), "Click 'Generate!' to create AI content"))
      shiny::div(
        class = "gene-info",
        shiny::HTML(opg_markdown_to_html(text))
      )
    }

    text.RENDER2 <- function() {
      text <- text_reactive()
      shiny::validate(shiny::need(!is.null(text), "Content not available"))
      shiny::div(
        style = "font-size: 18px;",
        shiny::HTML(opg_markdown_to_html(text))
      )
    }

    dl_text <- shiny::reactive({
      txt <- text_reactive()
      shiny::req(!is.null(txt))
      .aicards_text_markdown(txt)
    })

    PlotModuleServer(
      "report_text",
      plotlib = "generic",
      plotlib2 = "generic",
      func = text.RENDER,
      func2 = text.RENDER2,
      renderFunc = shiny::renderUI,
      renderFunc2 = shiny::renderUI,
      download.fmt = c("pdf", "docx", "md"),
      download.pdf  = .aicards_download_handler(dl_text, "wgcna-ai-report", "WGCNA AI Report", "pdf"),
      download.docx = .aicards_download_handler(dl_text, "wgcna-ai-report", "WGCNA AI Report", "docx"),
      download.md   = .aicards_download_handler(dl_text, "wgcna-ai-report", "WGCNA AI Report", "md"),
      pdf.width = 8, pdf.height = 11,
      res = c(75, 100),
      add.watermark = watermark
    )

    # ============================================================
    # Deep Report evidence trail (tool-call summaries)
    # ============================================================
    output$has_deep_trace <- shiny::reactive({
      if (is.null(last_deep_reactive) || is.null(deep_turns_reactive)) return(FALSE)
      if (!isTRUE(last_deep_reactive())) return(FALSE)
      turns <- deep_turns_reactive()
      !is.null(turns) && length(turns) > 0
    })
    shiny::outputOptions(output, "has_deep_trace", suspendWhenHidden = FALSE)

    output$deep_trace <- shiny::renderUI({
      if (is.null(deep_turns_reactive)) return(NULL)
      turns <- deep_turns_reactive()
      shiny::req(turns)
      .render_deep_tool_trace(turns)
    })

    invisible(NULL)
  })
}

#' Render the Deep Report tool-call trace as a list of collapsible blocks.
#'
#' One <details> per ContentToolRequest plus a small block per matching
#' ContentToolResult. Mirrors `.format_tool_request` in
#' `board.copilot/R/copilot_agent.R` but renders results inline instead
#' of dropping them from the chat stream.
.render_deep_tool_trace <- function(turns) {
  trunc <- function(s, n = 600L) {
    s <- as.character(s)
    if (length(s) == 0L || is.na(s)) return("")
    if (nchar(s) > n) paste0(substr(s, 1L, n), "\n… [truncated]") else s
  }
  scalar <- function(x) {
    if (is.null(x)) return("")
    paste(vapply(x, function(v) as.character(v)[1L], character(1L)), collapse = ", ")
  }

  blocks <- list()
  for (turn in turns) {
    contents <- tryCatch(turn@contents, error = function(e) list())
    for (item in contents) {
      if (S7::S7_inherits(item, ellmer::ContentToolRequest)) {
        tool_name <- tryCatch(item@name, error = function(e) "<unknown>")
        args      <- tryCatch(item@arguments, error = function(e) list())
        args_txt <- if (length(args) > 0L) {
          paste(names(args),
                vapply(args, function(v) trunc(scalar(v), 200L), character(1L)),
                sep = " = ", collapse = "\n")
        } else {
          "(no arguments)"
        }
        blocks[[length(blocks) + 1L]] <- shiny::tags$details(
          class = "mb-1",
          shiny::tags$summary(
            shiny::icon("wrench"),
            shiny::tags$code(tool_name)
          ),
          shiny::tags$pre(
            style = "font-size: 0.85em; white-space: pre-wrap;",
            args_txt
          )
        )
      } else if (S7::S7_inherits(item, ellmer::ContentToolResult)) {
        result_txt <- tryCatch(
          {
            v <- item@value
            if (is.character(v)) v else paste(utils::capture.output(print(v)), collapse = "\n")
          },
          error = function(e) ""
        )
        blocks[[length(blocks) + 1L]] <- shiny::tags$div(
          class = "ms-3 mb-2 small text-muted",
          shiny::tags$pre(
            style = "font-size: 0.8em; white-space: pre-wrap; background: #f7f7f7; padding: 4px;",
            trunc(result_txt, 800L)
          )
        )
      }
    }
  }
  if (length(blocks) == 0L) {
    return(shiny::tags$em("No tool calls were made during this run."))
  }
  shiny::tagList(blocks)
}
