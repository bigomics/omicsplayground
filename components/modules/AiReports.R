##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# PGX AI report access helpers
# =============================================================================

.ai_report_valid_entry <- function(x) {
  is.list(x) &&
    is.character(x$report) &&
    length(x$report) > 0L &&
    !is.na(x$report[[1L]]) &&
    nzchar(x$report[[1L]])
}

ai_report_slots <- function(pgx) {
  ai <- playbase::ai_report_ai_slot(pgx)
  if (is.null(ai)) return(character(0))

  slots <- setdiff(names(ai), "meta")
  slots[vapply(slots, function(slot) {
    .ai_report_valid_entry(ai[[slot]])
  }, logical(1))]
}

ai_report_get <- function(pgx, slot) {
  if (is.null(slot) || length(slot) != 1L) return(NULL)
  slot <- tryCatch(as.character(slot)[[1L]], error = function(e) NA_character_)
  if (is.na(slot) || !nzchar(slot)) {
    return(NULL)
  }

  ai <- playbase::ai_report_ai_slot(pgx)
  if (is.null(ai)) return(NULL)

  entry <- ai[[slot]]
  if (!.ai_report_valid_entry(entry)) return(NULL)

  list(
    slot = slot,
    report = as.character(entry$report)[[1L]],
    prompt = if (is.character(entry$prompt) && length(entry$prompt) > 0L) {
      as.character(entry$prompt)[[1L]]
    } else {
      NULL
    },
    meta = if (is.null(ai$meta)) NULL else ai$meta
  )
}

ai_report_has <- function(pgx, select = NULL) {
  slots <- ai_report_slots(pgx)
  if (is.null(select)) return(length(slots) > 0L)
  select <- tryCatch(as.character(select), error = function(e) character(0))
  select <- select[!is.na(select) & nzchar(select)]
  if (!length(select)) return(FALSE)
  all(vapply(select, function(slot) {
    slot %in% slots ||
      (identical(slot, "drugs") && any(startsWith(slots, "drugs_")))
  }, logical(1)))
}

ai_report_needs_generation <- function(pgx) {
  length(ai_report_modules_for_pgx(pgx)) > 0L && !ai_report_has(pgx)
}

ai_report_drug_slots <- function(pgx) {
  slots <- ai_report_slots(pgx)
  slots[startsWith(slots, "drugs_")]
}

ai_report_drug_label <- function(pgx, slot) {
  slot <- tryCatch(as.character(slot)[[1L]], error = function(e) NA_character_)
  if (is.na(slot) || !nzchar(slot)) return("")

  label <- gsub("_", " ", sub("^drugs_", "", slot))
  dbs <- tryCatch(names(pgx$drugs), error = function(e) character(0))
  if (length(dbs)) {
    safe <- paste0("drugs_", gsub("[^A-Za-z0-9]+", "_", dbs))
    match_idx <- match(slot, safe)
    if (!is.na(match_idx)) label <- dbs[[match_idx]]
  }

  if (grepl("^L1000([_/ -]|$)", label, ignore.case = TRUE)) {
    if (grepl("activ", label, ignore.case = TRUE)) return("L1000 Activity")
    if (grepl("gene", label, ignore.case = TRUE)) return("L1000 Gene")
  }
  label
}

ai_report_modules_for_pgx <- function(pgx) {
  if (is.null(pgx) || !is.list(pgx)) return(character(0))

  modules <- c(character(0),
    if (!is.null(pgx$wgcna)) "wgcna",
    if (!is.null(pgx$wgcna_mox)) "wgcna_mox",
    if (!is.null(pgx$mofa)) "mofa",
    if (!is.null(pgx$drugs) && length(pgx$drugs) > 0L) "drugs",
    if (!is.null(pgx$gx.meta)) "de",
    if (!is.null(pgx$gset.meta)) "pathways"
  )
  if (length(modules) > 0L) modules <- c(modules, "combined")
  unique(modules)
}

ai_report_get_module <- function(pgx, module) {
  if (is.null(module) || length(module) != 1L) return(NULL)
  module <- tryCatch(
    as.character(module)[[1L]],
    error = function(e) NA_character_
  )
  if (is.na(module) || !nzchar(module)) return(NULL)

  module <- switch(module,
    summary = "combined",
    enrichment = "pathways",
    module
  )
  ai_report_get(pgx, module)
}

ai_report_dataset_token <- function(pgx) {
  if (is.null(pgx) || !is.list(pgx)) return("")
  x <- tryCatch(pgx$X, error = function(e) NULL)
  xdim <- if (is.null(dim(x))) "" else paste(dim(x), collapse = "x")
  name <- if (is.null(pgx$name)) "" else as.character(pgx$name)[[1L]]
  date <- if (is.null(pgx$date)) "" else as.character(pgx$date)[[1L]]
  paste(name, xdim, date, sep = "|")
}

ai_report_generate <- function(pgx,
                               llm_model,
                               force = FALSE,
                               select = NULL,
                               img_model = NULL,
                               report_type = "normal",
                               on_error = "warn",
                               credentials = NULL) {
  if (is.null(llm_model) || !nzchar(llm_model)) return(pgx)
  if (is.null(select)) select <- ai_report_modules_for_pgx(pgx)
  select <- tryCatch(as.character(select), error = function(e) character(0))
  select <- select[!is.na(select) & nzchar(select)]
  if (!length(select)) return(pgx)

  playbase::pgx.update_reports(
    pgx,
    ai = list(
      llm_model = llm_model,
      img_model = img_model,
      select = select,
      report_type = report_type,
      force = isTRUE(force),
      on_error = on_error,
      credentials = credentials
    )
  )
}

# =============================================================================
# Asynchronous, parallel AI report generation
# =============================================================================
#
# ai_report_generate() above is the synchronous path: it blocks the caller for
# the whole run, and pgx.update_reports() issues one provider call at a time.
# Measured on a 29-contrast multi-omics pgx that is ~14 minutes for four of up
# to eight modules, essentially all of it spent waiting on the provider (prompt
# assembly for the same four: 5s).
#
# The async path below splits that in two. Prompts are assembled in the Shiny
# process - cheap, and it has to happen here because it needs the pgx - and
# only the resulting strings (<1MB per job, against a 269MB pgx) are handed to
# mirai workers, which do nothing but make the provider call. That is the same
# shape InfographicServer already uses for image jobs.
#
# Modules run concurrently, with one exception: `combined` summarises the other
# reports, so it is built and run in a second phase once phase one has been
# folded back into the pgx.

#' Number of mirai daemons to use for AI work
#'
#' These jobs are network-bound, not CPU-bound, so the useful count is not tied
#' to core count; it is capped low because each daemon is a full R process.
ai_report_daemon_count <- function() {
  n <- suppressWarnings(as.integer(Sys.getenv("OPG_AI_DAEMONS", "4")))
  if (is.na(n) || n < 1L) n <- 4L
  min(n, 8L)
}

#' Ensure mirai daemons exist, without disturbing an existing pool
#'
#' Called lazily on first use rather than at app start so a session that never
#' touches AI never pays for the worker processes.
ai_report_ensure_daemons <- function(n = ai_report_daemon_count()) {
  status <- tryCatch(mirai::status(), error = function(e) NULL)
  ## `status$daemons` is the socket URL once a pool exists, not a count - the
  ## count lives in `connections`. Reading the wrong field made this always
  ## fall through to daemons(), which tears the pool down and kills every
  ## in-flight mirai.
  current <- tryCatch(as.integer(status$connections), error = function(e) 0L)
  if (length(current) != 1L || is.na(current)) current <- 0L
  ## Never re-call daemons() on a live pool: it resets it. Every caller asks
  ## for the same default, so there is nothing to grow towards anyway.
  if (current > 0L) return(invisible(current))
  tryCatch(mirai::daemons(n), error = function(e) NULL)
  invisible(n)
}

#' Run one report job in a mirai worker
#'
#' The worker gets plain strings and a credential closure - never the pgx, and
#' never the job's `finalize`, which stays here and is applied on the way back.
#'
#' @param job An `ai_report_job` from `playbase::pgx.build_report_jobs()`.
#' @param llm_model Model id.
#' @param credentials Nullary credential closure, or NULL.
#' @param timeout_s Per-request deadline handed to omicsai.
#' @param retries Maximum attempts per job.
#' @return A promise resolving to a result list with `ok`, `slot` and either
#'   `text`/`usage` or `error`.
ai_report_job_promise <- function(job, llm_model, credentials = NULL,
                                  timeout_s = 240L, retries = 2L,
                                  reasoning_effort = "low") {
  promises::as.promise(mirai::mirai(
    {
      suppressPackageStartupMessages(library(omicsai))
      started <- Sys.time()
      out <- tryCatch({
        cfg_args <- list(
          model           = llm_model,
          system_prompt   = system_prompt,
          credentials     = credentials,
          timeout_seconds = timeout_s,
          retries         = retries
        )
        # Strict `extra` validation rejects the key outright on models that do
        # not declare it, so ask the registry first.
        if (!is.null(reasoning_effort) && nzchar(reasoning_effort) &&
            omicsai::omicsai_model_accepts_extra(llm_model, "reasoning_effort")) {
          cfg_args$reasoning_effort <- reasoning_effort
        }
        cfg <- do.call(omicsai::omicsai_config, cfg_args)
        res <- omicsai::omicsai_gen_text(board, config = cfg)
        usage <- res$metadata$usage
        if (!is.null(usage)) usage$model <- llm_model
        list(ok = TRUE, text = res$text, usage = usage)
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
      out$slot <- slot
      out$secs <- as.numeric(difftime(Sys.time(), started, units = "secs"))
      out
    },
    slot             = job$slot,
    system_prompt    = job$system,
    board            = job$board,
    llm_model        = llm_model,
    credentials      = credentials,
    timeout_s        = timeout_s,
    retries          = retries,
    reasoning_effort = reasoning_effort
  ))
}

#' Human-readable label for a report slot
#'
#' Slots are not module names: drugs fan out into `drugs_<db>`. Used for
#' progress text, so it must never fail on an unexpected slot.
#'
#' @param slot Slot name, or NA for a generic "overall summary" step.
#' @param pgx_list Optional pgx, used to recover the original drug DB name.
ai_report_slot_label <- function(slot, pgx_list = NULL) {
  if (is.null(slot) || length(slot) != 1L || is.na(slot)) return("overall summary")
  labels <- c(
    combined  = "Summary",
    wgcna     = "WGCNA",
    wgcna_mox = "moxWGCNA",
    mofa      = "MOFA",
    de        = "Differential Expression",
    pathways  = "Enrichment"
  )
  if (slot %in% names(labels)) return(unname(labels[[slot]]))
  if (startsWith(slot, "drugs_")) {
    label <- if (!is.null(pgx_list)) ai_report_drug_label(pgx_list, slot) else ""
    if (!nzchar(label)) label <- gsub("_", " ", sub("^drugs_", "", slot))
    return(paste("Drugs -", label))
  }
  slot
}

#' Rebuild the stored prompt record for a job
#'
#' Kept identical to what the synchronous playbase path writes, so reports look
#' the same in the Studio prompt view regardless of which path produced them.
ai_report_job_prompt_text <- function(job) {
  paste0("# SYSTEM\n\n", job$system, "\n\n---\n\n# BOARD\n\n", job$board)
}

#' Build the report options list shared by both generation paths
ai_report_options <- function(llm_model, select, force = FALSE,
                              credentials = NULL, report_type = "normal",
                              timeout_s = 240L, retries = 2L,
                              reasoning_effort = "low") {
  list(
    llm_model        = llm_model,
    select           = select,
    report_type      = report_type,
    force            = isTRUE(force),
    on_error         = "warn",
    credentials      = credentials,
    timeout_seconds  = timeout_s,
    retries          = retries,
    reasoning_effort = reasoning_effort
  )
}

#' Generate AI reports asynchronously, running modules in parallel
#'
#' @param pgx_list Plain (non-reactive) pgx list.
#' @param llm_model Model id.
#' @param select Modules to generate; defaults to everything the pgx supports.
#' @param force Regenerate modules that already have a report.
#' @param credentials Nullary credential closure, or NULL.
#' @param on_progress `function(done, total, slot, ok)`, called on the main
#'   thread as jobs complete; `slot` is NA for lifecycle notices such as the
#'   start of the summary phase. Cheap work only - it runs between fold-ins.
#' @param on_phase `function(phase, ai)` called once after each phase has been
#'   folded in - the right place to persist, since saving a large pgx costs
#'   seconds and doing it per module would hand back the stall we just removed.
#' @param still_valid `function()` returning FALSE to abandon the run (dataset
#'   switched, newer run started). Checked before each fold-in and before
#'   phase two.
#' @return A promise resolving to `list(ai=, done=, failed=, failures=)`, where
#'   `ai` is the updated `pgx$ai` slot. Individual module failures resolve
#'   normally and are reported in `failed`/`failures`; the promise only rejects
#'   if the run could not be set up at all.
ai_report_generate_async <- function(pgx_list, ...) {
  ## Prompt assembly happens synchronously, before any promise exists, so an
  ## error there would escape the caller's promise chain entirely and strand
  ## whatever progress UI it had opened. Funnel it into a rejection instead.
  tryCatch(
    .ai_report_generate_async(pgx_list, ...),
    error = function(e) promises::promise_reject(e)
  )
}

.ai_report_generate_async <- function(pgx_list,
                                     llm_model,
                                     select = NULL,
                                     force = FALSE,
                                     credentials = NULL,
                                     on_progress = NULL,
                                     on_phase = NULL,
                                     still_valid = NULL,
                                     timeout_s = 240L,
                                     retries = 2L,
                                     reasoning_effort = "low") {
  if (is.null(llm_model) || !nzchar(llm_model)) {
    return(promises::promise_resolve(NULL))
  }
  if (is.null(select)) select <- ai_report_modules_for_pgx(pgx_list)
  select <- tryCatch(as.character(select), error = function(e) character(0))
  select <- select[!is.na(select) & nzchar(select)]
  if (!length(select)) return(promises::promise_resolve(NULL))

  ai_report_ensure_daemons()
  opts <- ai_report_options(llm_model, select, force, credentials,
    timeout_s = timeout_s, retries = retries,
    reasoning_effort = reasoning_effort)

  ## `force` in playbase wipes pgx$ai wholesale before regenerating. The async
  ## path regenerates only what was asked for and merges, so clear exactly the
  ## requested slots here and let the builders see the rest.
  work <- pgx_list
  if (isTRUE(force) && is.list(work$ai)) {
    for (slot in names(work$ai)) {
      ## drugs fan out into drugs_<db> slots; everything else is 1:1.
      base <- sub("^(drugs)_.*$", "\\1", slot)
      if (base %in% select) work$ai[[slot]] <- NULL
    }
  }

  state <- new.env(parent = emptyenv())
  state$pgx      <- work
  state$done     <- 0L
  state$failures <- character(0)

  alive <- function() is.null(still_valid) || isTRUE(still_valid())

  jobs <- playbase::pgx.build_report_jobs(work, opts)
  total <- length(jobs) + as.integer("combined" %in% select)

  # Fold one worker result back into the pgx being assembled. Runs on the main
  # thread via promises, so no locking is needed around `state`.
  absorb <- function(job, result) {
    state$done <- state$done + 1L
    if (isTRUE(result$ok)) {
      state$pgx <- playbase::pgx.apply_report_result(
        state$pgx, job,
        list(report = result$text,
             prompt = ai_report_job_prompt_text(job),
             usage  = result$usage)
      )
    } else {
      state$failures <- c(state$failures, job$slot)
      warning("[ai_report_generate_async] ", job$slot, ": ", result$error,
        call. = FALSE)
    }
    if (is.function(on_progress)) {
      tryCatch(on_progress(state$done, total, job$slot, isTRUE(result$ok)),
        error = function(e) NULL)
    }
    invisible(NULL)
  }

  launch <- function(job) {
    promises::then(
      ai_report_job_promise(job, llm_model, credentials, timeout_s, retries,
        reasoning_effort),
      onFulfilled = function(result) {
        if (!alive()) return(NULL)
        absorb(job, result)
        NULL
      },
      onRejected = function(err) {
        if (!alive()) return(NULL)
        absorb(job, list(ok = FALSE, error = conditionMessage(err)))
        NULL
      }
    )
  }

  finish <- function() {
    list(
      ai       = state$pgx$ai,
      done     = state$done,
      failed   = length(state$failures),
      failures = state$failures
    )
  }

  phase_one <- if (!length(jobs)) {
    promises::promise_resolve(NULL)
  } else {
    promises::promise_all(.list = lapply(jobs, launch))
  }

  promises::then(phase_one, function(...) {
    if (!alive()) return(finish())
    if (is.function(on_phase)) {
      tryCatch(on_phase("modules", state$pgx$ai), error = function(e) NULL)
    }
    if (!"combined" %in% select) return(finish())

    ## Built now, not earlier: the combined prompt is assembled from the
    ## reports phase one just produced.
    cjobs <- playbase::pgx.build_combined_report_job(state$pgx, opts)
    if (!length(cjobs)) return(finish())
    if (is.function(on_progress)) {
      tryCatch(on_progress(state$done, total, NA_character_, TRUE),
        error = function(e) NULL)
    }
    promises::then(launch(cjobs[[1L]]), function(...) {
      if (is.function(on_phase)) {
        tryCatch(on_phase("combined", state$pgx$ai), error = function(e) NULL)
      }
      finish()
    })
  })
}

ai_report_update_text <- function(pgx, reports) {
  if (is.null(pgx) || !is.list(pgx)) return(pgx)
  ai <- playbase::ai_report_ai_slot(pgx)
  if (is.null(ai)) return(pgx)

  reports <- reports[!is.na(names(reports)) & nzchar(names(reports))]
  if (!length(reports)) return(pgx)

  for (slot in names(reports)) {
    if (!is.list(ai[[slot]])) next
    value <- tryCatch(as.character(reports[[slot]])[[1L]],
      error = function(e) NA_character_)
    if (is.na(value)) next
    ai[[slot]]$report    <- value
    ai[[slot]]$edited    <- TRUE
    ai[[slot]]$edited_at <- as.numeric(Sys.time())
  }

  pgx$ai <- ai
  pgx
}

ai_report_merge_into_reactive <- function(pgx_rv, ai) {
  if (is.null(pgx_rv)) return(invisible(FALSE))

  ai_slot <- if (is.list(ai) && !is.null(ai$ai) && is.list(ai$ai)) ai$ai else ai
  probe <- list(ai = ai_slot)
  slots <- ai_report_slots(probe)
  if (is.null(ai_slot) || !is.list(ai_slot) || !length(slots)) {
    return(invisible(FALSE))
  }

  current <- pgx_rv$ai
  if (is.null(current) || !is.list(current)) current <- list()
  for (slot in slots) {
    current[[slot]] <- ai_slot[[slot]]
  }
  if (!is.null(ai_slot$meta)) current$meta <- ai_slot$meta
  pgx_rv$ai <- current
  invisible(TRUE)
}

ai_report_copy_into_reactive <- function(pgx_rv, ai) {
  if (is.null(pgx_rv)) return(invisible(FALSE))

  ai_slot <- if (is.list(ai) && !is.null(ai$ai) && is.list(ai$ai)) ai$ai else ai
  probe <- list(ai = ai_slot)
  if (is.null(ai_slot) || !is.list(ai_slot) || !ai_report_has(probe)) {
    return(invisible(FALSE))
  }

  pgx_rv$ai <- ai_slot
  invisible(TRUE)
}

#' Convert a stored AI infographic entry into a Shiny renderImage value
#'
#' @param img Infographic entry returned by \code{playbase::ai_infographic_get()}.
#' @param tmpdir Directory used to materialize stored raw bytes.
#' @param name Base filename for materialized image bytes.
#' @return List suitable for \code{shiny::renderImage()}, or \code{NULL}.
ai_infographic_render_value <- function(img, tmpdir, name = "infographic") {
  if (is.null(img) || !is.list(img)) return(NULL)
  if (identical(img$status, "error")) {
    # Keep legacy/raw stored errors friendly before they reach the UI.
    msg <- playbase::ai_infographic_friendly_error(img$error)
    shiny::validate(shiny::need(FALSE, msg))
  }

  content_type <- if (is.null(img$content_type)) "image/png" else img$content_type
  ext <- if (identical(content_type, "image/jpeg")) ".jpg" else ".png"
  src <- NULL
  if (is.raw(img$bytes) && length(img$bytes) > 0L) {
    dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
    src <- file.path(tmpdir, paste0(name, ext))
    writeBin(img$bytes, src)
  } else if (is.character(img$path) && length(img$path) > 0L &&
      file.exists(img$path[[1L]])) {
    src <- img$path[[1L]]
  }

  if (is.null(src)) return(NULL)
  list(src = src, height = "auto", width = "100%", contentType = content_type)
}
