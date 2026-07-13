##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

if(0) {
  source("~/Playground/playbase/dev/include.R", chdir=TRUE)
  devtools::load_all("~/Playground/playbase")
  pgx <- playbase::pgx.load("~/Playground/omicsplayground/data/example-data-reports.pgx")
  pgx <- playbase::pgx.load("~/Playground/omicsplayground/data/mox-daf2.pgx")
  pgx <- playbase::pgx.load("~/Playground/omicsplayground/data/klidel-multiomics2.pgx")
  names(pgx)
  names(pgx$wgcna)
  names(pgx$mofa)  
  names(pgx$wgcna$report)
  names(pgx$wgcna_mox$report)  
  names(pgx$mofa)
}

AiReportSettings <- function(id) {
  ns <- shiny::NS(id)

  info ="AI reports for the current dataset. These reports are AI-generated and can be inaccurate; please always double-check its responses."
  
  shiny::div(
    style = "padding: 10px 15px;",
    h3("AI reports"),
    shiny::br(),
    shiny::div(info),
    shiny::br(),
    br(),
    br(),
    shiny::actionButton(ns("generate"), "Generate reports",
      class="btn btn-primary", style="margin-bottom: 6px;"),
    shiny::actionButton(ns("save_edits"), "Save edits",
      class="btn btn-secondary", style="margin-bottom: 6px;"),
    shiny::actionButton(ns("clear"), "Clear",
      class="btn btn-warning-outline", style="margin-bottom: 6px;"),
    br(),
    shiny::checkboxInput(ns("force"), "Force regenerate", FALSE),
    br(),
    shiny::fileInput(ns("logofile"), "Replace logo:", accept=c(".png",".jpg"))
  )
}

  
AiReportUI <- function(id) {
  ns <- shiny::NS(id)

  preview_ui <- bslib::navset_tab(
    id = ns("navset_preview"),
    selected = "Summary",
    bslib::nav_panel(title = "Summary",
      shiny::div( shiny::htmlOutput(ns("summary")) %>% bigLoaders::useSpinner(),
        style="align-items: center;")
    ),
    bslib::nav_panel(title = "WGCNA",
      shiny::div( shiny::htmlOutput(ns("wgcna"))%>% bigLoaders::useSpinner(),
        style="align-items: center;")      
    ),
    bslib::nav_panel(title = "moxWGCNA",
      shiny::div( shiny::htmlOutput(ns("wgcna2")) %>% bigLoaders::useSpinner(),
        style="align-items: center;")
    ),
    bslib::nav_panel(title = "MOFA",
      shiny::div( shiny::htmlOutput(ns("mofa")) %>% bigLoaders::useSpinner(),
        style="align-items: center;")
    ),
    bslib::nav_panel(title = "DE",
      shiny::div( shiny::htmlOutput(ns("de")) %>% bigLoaders::useSpinner(),
        style="align-items: center;")
    ),    
    bslib::nav_panel(title = "Enrichment",
      shiny::div( shiny::htmlOutput(ns("enrichment")) %>% bigLoaders::useSpinner(),
        style="align-items: center;")
    )
  )

  ## Local UI factory for the edit tabs.
  ## Kept inside AiReportUI because it only depends on this module namespace
  ## and avoids repeating the same textAreaInput boilerplate seven times.
  text_area <- function(id) {
    div( shiny::textAreaInput(ns(id), NULL, value = "",
      height = "calc(100vh - 160px)", width = "100%"),
      id="reportTA", style='font-size: 1em;')
  }
  
  edit_ui <- bslib::navset_tab(
    id = ns("navset_edit"),
    selected = "Summary",
    bslib::nav_panel(title = "Summary", text_area("edit_summary")),
    bslib::nav_panel(title = "WGCNA", text_area("edit_wgcna")),
    bslib::nav_panel(title = "moxWGCNA", text_area("edit_wgcna2")),
    bslib::nav_panel(title = "MOFA", text_area("edit_mofa")),
    bslib::nav_panel(title = "DE", text_area("edit_de")),
    bslib::nav_panel(title = "Enrichment", text_area("edit_enrichment"))
  )

  ui <- bslib::navset_hidden(
    id = ns("navset"),
    selected = "Preview",
    bslib::nav_panel(title = "Preview", preview_ui),    
    bslib::nav_panel(title = "Edit", edit_ui),
    header = tagList(
      bslib::nav_spacer(),
      div( bslib::input_switch(ns("edit_mode"), "edit", FALSE),
        style = "position: absolute; width: 100px; right: 0px; font-size: 1.1em;")
      )
  )
  return(ui)
}


AiReportServer <- function(id, pgx, save_pgx = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    
    pdf_tempdir <- tempdir()
    opg.logo = file.path(OPG,"components/app/R/www/bigomics-logo-small.png")
    logopath = file.path(pdf_tempdir,"logo.png")
    file.copy(opg.logo, logopath)
    
    de_pdf <- file.path(pdf_tempdir,"report-de.pdf")
    enrichment_pdf <- file.path(pdf_tempdir,"report-enrichment.pdf")      
    mofa_pdf <- file.path(pdf_tempdir,"report-mofa.pdf")
    wgcna_pdf <- file.path(pdf_tempdir,"report-wgcna.pdf")
    wgcna2_pdf <- file.path(pdf_tempdir,"report-wgcna2.pdf")      
    summary_pdf <- file.path(pdf_tempdir,"report-summary.pdf")      

    shiny::addResourcePath('pdf', pdf_tempdir)
    on.exit({
      unlink(de_pdf)
      unlink(enrichment_pdf)
      unlink(mofa_pdf)
      unlink(wgcna_pdf)
      unlink(wgcna2_pdf)
      unlink(summary_pdf)      
      unlink(Sys.glob(file.path(pdf_tempdir, "report-drugs_*.pdf")))
      removeResourcePath('pdf')
    })

    dynamic_drug_tabs <- character(0)

    ## Local edit-control factory for drug tabs inserted at runtime.
    ## It mirrors the static edit text areas but must run inside the server
    ## module so the dynamically generated ids receive the session namespace.
    drug_text_area <- function(id) {
      div(shiny::textAreaInput(ns(id), NULL, value = "",
        height = "calc(100vh - 160px)", width = "100%"),
        id = "reportTA", style = "font-size: 1em;")
    }
    
    ## Capture the current reactive PGX state as a plain list.
    ## Async report generation cannot safely carry reactiveValues into futures,
    ## so all downstream helpers operate on this immutable snapshot.
    pgx_snapshot <- function() {
      shiny::reactiveValuesToList(pgx)
    }

    ## Read one AI report body from the normalized pgx$ai module schema.
    ## Missing reports are represented as empty strings so Shiny text areas
    ## can be refreshed without null checks at every call site.
    report_text <- function(module, pgx_list = pgx_snapshot()) {
      entry <- ai_report_get_module(pgx_list, module)
      if (is.null(entry)) "" else entry$report
    }

    ## Insert one Preview/Edit tab pair per concrete pgx$ai$drugs_* slot.
    ## Drug reports are generated as multiple independent reports, so the UI
    ## must not use a synthetic concatenated "drugs" report.
    sync_drug_tabs <- function(pgx_list = pgx_snapshot()) {
      for (slot in dynamic_drug_tabs) {
        try(bslib::nav_remove("navset_preview", target = slot), silent = TRUE)
        try(bslib::nav_remove("navset_edit", target = slot), silent = TRUE)
      }
      dynamic_drug_tabs <<- character(0)

      drug_slots <- ai_report_drug_slots(pgx_list)
      target <- "moxWGCNA"
      for (slot in drug_slots) {
        label <- ai_report_drug_label(pgx_list, slot)
        output_id <- paste0("preview_", slot)
        input_id <- paste0("edit_", slot)
        pdfname <- paste0("report-", slot, ".pdf")

        bslib::nav_insert("navset_preview",
          target = target,
          position = "after",
          nav = bslib::nav_panel(
            title = label,
            value = slot,
            shiny::div(shiny::htmlOutput(ns(output_id)) %>%
              bigLoaders::useSpinner(), style = "align-items: center;")
          )
        )
        bslib::nav_insert("navset_edit",
          target = target,
          position = "after",
          nav = bslib::nav_panel(
            title = label,
            value = slot,
            drug_text_area(input_id)
          )
        )

        local({
          slot_i <- slot
          input_i <- input_id
          output_i <- output_id
          pdf_i <- pdfname
          label_i <- label
          output[[output_i]] <- renderUI({
            updatePDF(input[[input_i]], slot_i, pdf_i)
            shiny::validate(need(file.exists(file.path(pdf_tempdir, pdf_i)),
              paste("missing", label_i, "report")))
            return(pdf.iframe(pdf_i))
          })
        })

        dynamic_drug_tabs <<- c(dynamic_drug_tabs, slot)
        target <- slot
      }
      invisible(drug_slots)
    }

    ## Push current report text into the edit-mode text areas.
    ## This centralizes module name mapping, because static UI tab names differ
    ## from pgx$ai module keys for summary, enrichment, and moxWGCNA.
    update_report_inputs <- function(pgx_list = pgx_snapshot()) {
      updateTextAreaInput(session, "edit_wgcna",
        value = report_text("wgcna", pgx_list))
      updateTextAreaInput(session, "edit_wgcna2",
        value = report_text("wgcna_mox", pgx_list))
      updateTextAreaInput(session, "edit_mofa",
        value = report_text("mofa", pgx_list))
      updateTextAreaInput(session, "edit_summary",
        value = report_text("combined", pgx_list))
      updateTextAreaInput(session, "edit_de",
        value = report_text("de", pgx_list))
      updateTextAreaInput(session, "edit_enrichment",
        value = report_text("pathways", pgx_list))
      for (slot in ai_report_drug_slots(pgx_list)) {
        updateTextAreaInput(session, paste0("edit_", slot),
          value = ai_report_get(pgx_list, slot)$report)
      }
    }

    ## Pick the Generate/Reset button label from report availability.
    ## Keeping this as one helper prevents load, clear, success, and error
    ## observers from drifting in their button state logic.
    report_button_label <- function(pgx_list = pgx_snapshot()) {
      if (ai_report_has(pgx_list)) "Reset reports" else "Generate reports"
    }

    edited_reports <- function(pgx_list = pgx_snapshot()) {
      reports <- c(
        combined = input$edit_summary,
        wgcna = input$edit_wgcna,
        wgcna_mox = input$edit_wgcna2,
        mofa = input$edit_mofa,
        de = input$edit_de,
        pathways = input$edit_enrichment
      )
      for (slot in ai_report_drug_slots(pgx_list)) {
        reports[[slot]] <- input[[paste0("edit_", slot)]]
      }
      reports[ai_report_slots(pgx_list)]
    }

    report_module_labels <- function(modules) {
      labels <- c(
        combined = "Summary",
        wgcna = "WGCNA",
        wgcna_mox = "moxWGCNA",
        mofa = "MOFA",
        drugs = "Drug reports",
        de = "Differential Expression",
        pathways = "Enrichment"
      )
      out <- labels[modules]
      out[is.na(out)] <- modules[is.na(out)]
      unname(out)
    }

    selected_existing_modules <- function(pgx_list, modules) {
      selected <- modules[vapply(modules, function(module) {
        ai_report_has(pgx_list, module)
      }, logical(1))]
      if (length(selected)) selected else modules
    }

    ## Show only tabs backed by reports present in pgx$ai.
    ## This replaced the old slot-specific logic so the UI follows the new
    ## report schema instead of legacy pgx$report and module$report paths.
    update_nav <- function(pgx_list = pgx_snapshot()) {
      drug_slots <- sync_drug_tabs(pgx_list)
      drug_reports <- lapply(drug_slots, function(slot) ai_report_get(pgx_list, slot))
      names(drug_reports) <- drug_slots

      reports <- c(
        list(
          Summary = ai_report_get_module(pgx_list, "combined"),
          WGCNA = ai_report_get_module(pgx_list, "wgcna"),
          moxWGCNA = ai_report_get_module(pgx_list, "wgcna_mox")
        ),
        drug_reports,
        list(
          MOFA = ai_report_get_module(pgx_list, "mofa"),
          DE = ai_report_get_module(pgx_list, "de"),
          Enrichment = ai_report_get_module(pgx_list, "pathways")
        )
      )
      available <- names(reports)[!vapply(reports, is.null, logical(1))]
      selected <- if (length(available)) available[[1L]] else NULL

      for (target in names(reports)) {
        if (is.null(reports[[target]])) {
          bslib::nav_hide("navset_preview", target)
          bslib::nav_hide("navset_edit", target)
        } else {
          bslib::nav_show("navset_preview", target,
            select = identical(target, selected))
          bslib::nav_show("navset_edit", target,
            select = identical(target, selected))
        }
      }
      if (!is.null(selected)) {
        bslib::nav_select("navset",
          if (isTRUE(input$edit_mode)) "Edit" else "Preview")
      }
      invisible(selected)
    }

    ## Refresh all visible report controls from one PGX snapshot.
    ## Used after dataset load, pgx$ai updates, and async completion so those
    ## paths cannot accidentally diverge in UI state handling.
    refresh_report_ui <- function(pgx_list = pgx_snapshot()) {
      update_nav(pgx_list)
      clear_files()
      update_report_inputs(pgx_list)
      updateCheckboxInput(session, "force", value = FALSE)
      updateActionButton(session, "generate",
        label = report_button_label(pgx_list))
    }

    has.warned=0
    
    shiny::observeEvent( input$edit_mode, {
      if(input$edit_mode) {
        bslib::nav_select("navset_edit", input$navset_preview)                  
        bslib::nav_select("navset", "Edit")
        ## shinyalert::shinyalert(
        ##   text = HTML("Warning. Edits are currently <b>not saved</b> after logout!"),
        ##   html = TRUE,
        ##   closeOnClickOutside = TRUE,
        ##   timer = ifelse(has.warned, 1000, 0)
        ## )
        has.warned <<- has.warned + 1
      } else {
        bslib::nav_select("navset_preview", input$navset_edit)          
        bslib::nav_select("navset", "Preview")
      }
    })
    
    ## ----------------------- PDF outputs -------------------------
    ## Track letter counts to detect text changes before re-rendering PDFs.
    ## This is local because the cache only serves the current Shiny session
    ## and depends on transient edit text, temp files, and logo state.
    countLetters <- function(text) {
      az <- c(letters, LETTERS)
      cc <- sapply(az, function(s) sum(gregexpr(s,text)[[1]]>0))
      cc[az]
    }

    char_counts <- c( "wgcna" = 0, "wgcna2" = 0, "mofa" = 0,
      "summary" = 0, "de" = 0, "enrichment" = 0)    
    letter_counts = rep(list(countLetters("")), length(names(char_counts)))
    names(letter_counts) <- names(char_counts)
    
    ## Decide whether a report needs a fresh PDF render.
    ## The helper updates the local cache counters as a side effect, keeping
    ## the renderUI blocks short and consistent across all report tabs.
    has.changed <- function(what, now) {
      if(is.null(now)) now <- ""
      if (is.na(char_counts[what])) char_counts[what] <<- 0
      if (is.null(letter_counts[[what]])) {
        letter_counts[[what]] <<- countLetters("")
      }
      ct <- nchar(now)
      changed1 <- (!is.null(now) && ct != char_counts[what])
      char_counts[what] <<- ct
      lt <- countLetters(now)
      changed2 <- !all((lt - letter_counts[[what]]) == 0)
      letter_counts[[what]] <<- lt
      return(changed1 || changed2)
    }

    ## Render markdown to the temp PDF only when report text changed.
    ## Kept local because it closes over the session tempdir, logo path, and
    ## per-tab cache counters used by has.changed().
    updatePDF <- function(rpt, rptname, pdfname) {
      if(has.changed(rptname, rpt)) {
        shiny::withProgress(message = "Rendering PDF...", value = 0.7, {
          file <- file.path(pdf_tempdir, pdfname)
          playbase::markdownToPDF(rpt, file=file, logo=logopath, quiet=TRUE)
        })
      }
    }

    ## Build the PDF preview iframe for a generated temp PDF.
    ## This isolates the resource-path URL convention so each renderUI block
    ## only names the PDF it wants to display.
    pdf.iframe <- function(pdfname) {
      HTML(paste0('<embed style="height: calc(100vh - 200px); width: calc(100% - 20px)" src="pdf/',pdfname,'" type="application/pdf" />'))      
    }
    
    output$summary <- renderUI({
      updatePDF(input$edit_summary, "summary", "report-summary.pdf")
      shiny::validate(need(file.exists(summary_pdf), "missing summary report"))
      return(pdf.iframe("report-summary.pdf"))
    })

    output$wgcna <- renderUI({
      updatePDF(input$edit_wgcna, "wgcna", "report-wgcna.pdf")      
      shiny::validate(need(file.exists(wgcna_pdf),"missing WGCNA report"))      
      return(pdf.iframe("report-wgcna.pdf"))
    })

    output$wgcna2 <- renderUI({
      updatePDF(input$edit_wgcna2, "wgcna2", "report-wgcna2.pdf")            
      shiny::validate(need(file.exists(wgcna2_pdf),"missing moxWGCNA report"))      
      return(pdf.iframe("report-wgcna2.pdf"))
    })

    output$mofa <- renderUI({
      updatePDF(input$edit_mofa, "mofa", "report-mofa.pdf")                        
      shiny::validate(need(file.exists(mofa_pdf),"missing mofa report"))            
      return(pdf.iframe("report-mofa.pdf"))
    })

    output$de <- renderUI({
      updatePDF(input$edit_de, "de", "report-de.pdf")
      shiny::validate(need(file.exists(de_pdf),"missing DE report"))            
      return(pdf.iframe("report-de.pdf"))
    })

    output$enrichment <- renderUI({
      updatePDF(input$edit_enrichment, "enrichment", "report-enrichment.pdf")
      shiny::validate(need(file.exists(enrichment_pdf),"missing Enrichment report"))
      return(pdf.iframe("report-enrichment.pdf"))
    })

    report_jobs <- shiny::reactiveValues(
      running = FALSE,
      run_id = 0L,
      token = NULL,
      total = 0L,
      done = 0L,
      failed = 0L,
      queue = character(0),
      progress = NULL,
      changed = FALSE
    )

    report_progress <- function(detail = NULL) {
      if (is.null(report_jobs$progress)) return(invisible(NULL))
      value <- if (report_jobs$total > 0L) report_jobs$done / report_jobs$total else 0
      if (is.null(detail)) {
        detail <- paste0(report_jobs$done, " / ",
          report_jobs$total, " reports complete")
      }
      report_jobs$progress$set(
        message = "Regenerating AI reports",
        detail = detail,
        value = value
      )
      invisible(NULL)
    }

    finish_report_jobs <- function(message = NULL, type = "message") {
      info("[AiReportServer] report jobs finish: done=", report_jobs$done,
        " total=", report_jobs$total,
        " failed=", report_jobs$failed,
        " changed=", report_jobs$changed,
        " message=", message)
      if (!is.null(report_jobs$progress)) {
        report_jobs$progress$close()
        report_jobs$progress <- NULL
      }
      report_jobs$running <- FALSE
      report_jobs$queue <- character(0)
      updateActionButton(session, "generate",
        label = report_button_label(pgx_snapshot()))
      if (isTRUE(report_jobs$changed) && !is.null(save_pgx)) save_pgx(pgx)
      refresh_report_ui(pgx_snapshot())
      if (!is.null(message)) {
        shiny::showNotification(message, type = type, session = session)
      }
    }

    launch_next_report_job <- function(llm_model, run_id, cred_fn = NULL) {
      if (!isTRUE(report_jobs$running) ||
          !identical(run_id, report_jobs$run_id)) {
        return(NULL)
      }
      if (!length(report_jobs$queue)) return(NULL)

      module <- report_jobs$queue[[1L]]
      report_jobs$queue <- report_jobs$queue[-1L]
      launch_report_job(module, llm_model, run_id, cred_fn)
    }

    launch_report_job <- function(module, llm_model, run_id, cred_fn = NULL) {
      pgx_list <- pgx_snapshot()
      token <- report_jobs$token
      label <- report_module_labels(module)
      report_progress(paste("Starting", label))

      promise <- promises::future_promise({
        list(
          module = module,
          token = token,
          pgx = ai_report_generate(
            pgx_list,
            llm_model = llm_model,
            force = TRUE,
            select = module,
            img_model = NULL,
            report_type = "normal",
            on_error = "warn",
            credentials = cred_fn
          )
        )
      }, seed = TRUE)

      promises::then(
        promise,
        onFulfilled = function(result) {
          if (!isTRUE(report_jobs$running) ||
              !identical(run_id, report_jobs$run_id)) {
            info("[AiReportServer] report job ignored: module=", result$module,
              " stale run_id=", run_id,
              " current_run_id=", report_jobs$run_id)
            return(NULL)
          }

          current_pgx <- pgx_snapshot()
          if (!identical(result$token, ai_report_dataset_token(current_pgx))) {
            info("[AiReportServer] report job stale dataset: module=", result$module)
            finish_report_jobs(
              "AI reports finished for a previous dataset; results were ignored.",
              type = "warning"
            )
            return(NULL)
          }

          updated <- shiny::isolate(
            ai_report_merge_into_reactive(pgx, result$pgx)
          )
          if (isTRUE(updated)) {
            report_jobs$changed <- TRUE
          }

          report_jobs$done <- report_jobs$done + 1L

          # Telemetry: log one event per report generation on the main thread.
          # TODO(telemetry): tokens unavailable until playbase returns usage from
          # pgx.update_reports() (computed in the future worker, not returned here).
          # TODO(telemetry): thread auth$email for attribution (no auth in StudioServer scope).
          tryCatch(
            ai_telemetry_record(
              source     = "report",
              session_id = paste0(session$token, ":", result$module),
              user_email = NA_character_,
              model      = llm_model,
              usage      = NULL
            ),
            error = function(e) NULL
          )

          report_progress(paste("Finished", report_module_labels(result$module)))
          info("[AiReportServer] report job complete: module=", result$module,
            " done=", report_jobs$done, "/", report_jobs$total,
            " updated=", updated)

          if (report_jobs$done >= report_jobs$total) {
            if (report_jobs$failed > 0L) {
              finish_report_jobs("AI reports completed with warnings.",
                type = "warning")
            } else {
              finish_report_jobs("AI reports ready.")
              shinyalert::shinyalert(
                title = "AI reports ready",
                text = "Your AI reports are ready.",
                type = "success",
                confirmButtonText = "OK"
              )
            }
          } else {
            launch_next_report_job(llm_model, run_id, cred_fn)
          }
          NULL
        },
        onRejected = function(err) {
          if (!isTRUE(report_jobs$running) ||
              !identical(run_id, report_jobs$run_id)) {
            info("[AiReportServer] report job error ignored: module=", module,
              " stale run_id=", run_id,
              " current_run_id=", report_jobs$run_id)
            return(NULL)
          }
          info("[AiReportServer] report job failed: module=", module,
            " error=", conditionMessage(err))
          report_jobs$failed <- report_jobs$failed + 1L
          report_jobs$done <- report_jobs$done + 1L
          report_progress(paste("Failed", label))
          shiny::showNotification(conditionMessage(err),
            type = "error", session = session)

          if (report_jobs$done >= report_jobs$total) {
            finish_report_jobs("AI report generation completed with errors.",
              type = "error")
          } else {
            launch_next_report_job(llm_model, run_id, cred_fn)
          }
          NULL
        }
      )
      invisible(promise)
    }

    start_report_jobs <- function(modules, llm_model, cred_fn = NULL) {
      modules <- unique(as.character(modules))
      modules <- modules[!is.na(modules) & nzchar(modules)]
      if (!length(modules)) {
        shiny::showNotification("Select at least one report to regenerate.",
          type = "warning", session = session)
        return(NULL)
      }

      report_jobs$run_id <- report_jobs$run_id + 1L
      run_id <- report_jobs$run_id
      report_jobs$running <- TRUE
      report_jobs$token <- ai_report_dataset_token(pgx_snapshot())
      report_jobs$total <- length(modules)
      report_jobs$done <- 0L
      report_jobs$failed <- 0L
      report_jobs$changed <- FALSE

      first_modules <- setdiff(modules, "combined")
      report_jobs$queue <- c(first_modules, intersect(modules, "combined"))
      report_jobs$progress <- shiny::Progress$new(session, min = 0, max = 1)
      report_progress("Starting selected reports")
      updateActionButton(session, "generate", label = "Generating...")

      launch_next_report_job(llm_model, run_id, cred_fn)
      invisible(NULL)
    }

    ##---------------------------------------------------------------
    ## --------- refresh observer: pgx load -------------------------
    ##---------------------------------------------------------------
    ## On dataset change, only refresh the Studio UI from whatever
    ## reports the pgx already carries. Auto-generation on load is
    ## owned by LoadingBoard's yes/no prompt
    ## (board.loading/R/loading_server.R::maybe_offer_ai_reports).

    shiny::observeEvent(list(pgx$X, pgx$name), {
      shiny::req(pgx$X, pgx$name)
      dbg("[AiReportServer] pgx load refresh (no auto-generation)")
      refresh_report_ui(pgx_snapshot())
    })

    shiny::observeEvent(pgx$ai, {
      shiny::req(pgx$X, pgx$name)
      dbg("[AiReportServer] pgx$ai refresh")
      refresh_report_ui(pgx_snapshot())
    }, ignoreInit = TRUE, ignoreNULL = FALSE)

    ##---------------------------------------------------------------
    ## --------- main observer: generate reports --------------------
    ##---------------------------------------------------------------
    ## Only the Generate button triggers report computation. PGX
    ## load is handled by the refresh observer above and must not
    ## re-enter this path, otherwise we re-introduce the silent
    ## auto-generation that LoadingBoard's prompt is meant to gate.

    shiny::observeEvent(input$generate, {
      shiny::req(pgx$X, pgx$name)

      llm_model <- getUserOption(session, "llm_model")
      if (is.null(llm_model) || llm_model == "") {
        shiny::showNotification("Please enable AI/LLM in user settings.",
          type = "warning", session = session)
        return(NULL)
      }

      if (isTRUE(report_jobs$running)) {
        shiny::showNotification("AI report generation is already running.",
          type = "message", session = session)
        return(NULL)
      }

      pgx_list <- pgx_snapshot()
      report_modules <- ai_report_modules_for_pgx(pgx_list)
      if (!length(report_modules)) {
        shiny::showNotification("No reportable modules found in this dataset.",
          type = "warning", session = session)
        return(NULL)
      }

      selected <- selected_existing_modules(pgx_list, report_modules)
      choices <- stats::setNames(report_modules,
        report_module_labels(report_modules))
      shiny::showModal(shiny::modalDialog(
        title = "Regenerate AI reports",
        shiny::checkboxGroupInput(ns("report_select"), NULL,
          choices = choices, selected = selected),
        footer = shiny::tagList(
          shiny::modalButton("Cancel"),
          shiny::actionButton(ns("confirm_generate"),
            "Regenerate selected", class = "btn btn-primary")
        ),
        easyClose = TRUE
      ))
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$confirm_generate, {
      modules <- input$report_select
      if (!length(modules)) {
        shiny::showNotification("Select at least one report to regenerate.",
          type = "warning", session = session)
        return(NULL)
      }
      llm_model <- getUserOption(session, "llm_model")
      if (is.null(llm_model) || llm_model == "") {
        shiny::showNotification("Please enable AI/LLM in user settings.",
          type = "warning", session = session)
        return(NULL)
      }
      if (isTRUE(report_jobs$running)) {
        shiny::showNotification("AI report generation is already running.",
          type = "message", session = session)
        return(NULL)
      }
      # Capture credential closure here — in the reactive context — before it
      # enters the future_promise worker (where session is unavailable).
      cred_fn <- get_ai_credentials(session)

      shiny::removeModal()
      info("[AiReportServer] report modal confirmed: modules=",
        paste(modules, collapse = ","),
        " llm_model=", llm_model)
      start_report_jobs(modules, llm_model, cred_fn)
    }, ignoreInit = TRUE)

    shiny::observeEvent(pgx$name, {
      if (isTRUE(report_jobs$running)) {
        info("[AiReportServer] report jobs cancelled: dataset changed")
        finish_report_jobs(
          "AI report generation was cancelled because the dataset changed.",
          type = "warning"
        )
      }
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$save_edits, {
      pgx_list <- pgx_snapshot()
      if (!ai_report_has(pgx_list)) {
        shiny::showNotification("No AI reports to save.",
          type = "warning", session = session)
        return(NULL)
      }

      updated_pgx <- ai_report_update_text(pgx_list, edited_reports(pgx_list))
      updated <- shiny::isolate(ai_report_copy_into_reactive(pgx, updated_pgx))
      if (isTRUE(updated) && !is.null(save_pgx)) save_pgx(pgx)

      refresh_report_ui(pgx_snapshot())
      shiny::showNotification("AI report edits saved.",
        type = "message", session = session)
    }, ignoreInit = TRUE)
    
    ##-------------------------------------------------
    ##------------ other observers --------------------
    ##-------------------------------------------------

    observe({
      file <- input$logofile
      req(file)
      ext <- tools::file_ext(file$datapath)
      validate(need(ext %in% c('jpg','png'), "Please upload a png or jpg file"))
      file.copy(file$datapath, logopath, overwrite=TRUE)
    })

    shiny::observeEvent(input$clear, {
      clear_files()
      update_nav()
      updateCheckboxInput(session, "force", value = FALSE)      
      updateActionButton(session, "generate", label = "Generate reports")
    }, ignoreInit = TRUE)

    ## Clear edit text, cached PDFs, and local change counters.
    ## This is session-local cleanup, not PGX mutation, so the helper stays
    ## inside the server module where temp paths and UI ids are available.
    clear_files <- function() {
      dbg("[clear_files] clearing files!")
      updateTextAreaInput(session, "edit_wgcna", value = "")
      updateTextAreaInput(session, "edit_wgcna2", value = "")
      updateTextAreaInput(session, "edit_mofa", value = "")
      updateTextAreaInput(session, "edit_summary", value = "")
      updateTextAreaInput(session, "edit_de", value = "")
      updateTextAreaInput(session, "edit_enrichment", value = "")
      for (slot in dynamic_drug_tabs) {
        updateTextAreaInput(session, paste0("edit_", slot), value = "")
      }
      unlink(file.path(pdf_tempdir,"report-wgcna.pdf"))
      unlink(file.path(pdf_tempdir,"report-wgcna2.pdf"))
      unlink(file.path(pdf_tempdir,"report-mofa.pdf"))
      unlink(file.path(pdf_tempdir,"report-summary.pdf"))
      unlink(file.path(pdf_tempdir,"report-de.pdf"))
      unlink(file.path(pdf_tempdir,"report-enrichment.pdf"))
      unlink(file.path(pdf_tempdir,
        paste0("report-", dynamic_drug_tabs, ".pdf")))
      char_counts <<- char_counts * 0
      letter_counts <<- lapply(letter_counts, function(x) x*0)
    }

    
  }) ## end of moduleServer
}
