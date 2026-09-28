##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# App-scope AI report runs
# =============================================================================
#
# ai_report_generate_async() is the engine: it fans report modules out across
# mirai workers and folds the results back in. Driven directly from a module
# server, though, the run is owned by one Shiny session - so closing the tab or
# loading another dataset abandons work that has already been paid for, and the
# only record of a finished report is the pgx the session happens to hold.
#
# This manager owns runs at app scope instead. A run is keyed by dataset, user
# and model, not by session; sessions merely subscribe to it for progress. Two
# things follow:
#
#   * a run outlives the session that started it, and a second session (or the
#     same user in another tab) joins the run in progress rather than starting
#     a duplicate one - but only when the second request is equivalent, since
#     the joiner would otherwise get reports written by the starter's model and
#     BYOK credentials;
#   * the pgx is committed at every phase boundary, so a crash costs at most
#     the phase still in flight rather than the whole run. Each report is also
#     written to a sidecar file the moment it lands, which is finer-grained
#     still and is what a future resume pass would read.
#
# Committing per phase is not free: ai_report_commit_to_disk() re-reads and
# re-writes the pgx on the Shiny thread and that stalls every session in the
# process for seconds on a large dataset. Per phase (two of them) is the price
# of durability; per module would be several times that, which is why the
# engine offers a phase hook at all.
#
# What this deliberately does NOT do is survive a restart of the R process:
# mirai daemons die with their host, so in-flight calls are lost.

#' Directory holding the sidecar reports for one run
#'
#' Sits next to the pgx it belongs to - for permissions and locality only.
#' Nothing copies or deletes it with the dataset (dataset management in
#' `loading_table_datasets.R` touches the `.pgx` alone), so the manager clears
#' it itself once the reports are safely in the pgx.
#'
#' The dataset token is part of the directory name: two runs on the same path
#' but different versions of the dataset must not overwrite each other's
#' reports, exactly as they do not share a run key.
#'
#' @param save_path Absolute path of the target `.pgx` file.
#' @param token Dataset token from `ai_report_dataset_token()`.
#' @return Directory path (not created).
ai_report_sidecar_dir <- function(save_path, token) {
  stem <- sub("[.]pgx$", "", basename(save_path))
  ## Hashed rather than sanitised: the token is built from the dataset name,
  ## which is free text and cannot be trusted as a path component.
  file.path(dirname(save_path),
            ".ai_reports",
            paste0(stem, "-", substr(digest::digest(token), 1L, 8L)))
}

#' Write one finished report to its sidecar
#'
#' Best-effort: a failure here must never take down the run, since the report
#' is still in memory and will be folded into the pgx at the next phase.
ai_report_write_sidecar <- function(save_path, token, slot, text, usage = NULL) {
  dir <- ai_report_sidecar_dir(save_path, token)
  tryCatch({
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(text, file.path(dir, paste0(slot, ".md")))
    meta <- list(slot = slot, written = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
                 chars = nchar(text), usage = usage)
    writeLines(jsonlite::toJSON(meta, auto_unbox = TRUE, null = "null"),
               file.path(dir, paste0(slot, ".json")))
    TRUE
  }, error = function(e) {
    info("[AiReportManager] sidecar write failed: slot=", slot,
         " error=", conditionMessage(e))
    FALSE
  })
}

## Drop the sidecars once the pgx holds the same reports: they are a crash
## recovery buffer, not an archive, and nothing else ever cleans them up.
.ai_report_clear_sidecars <- function(save_path, token) {
  unlink(ai_report_sidecar_dir(save_path, token), recursive = TRUE)
  invisible(NULL)
}

## Registry of live runs, one entry per key. App scope on purpose: this file is
## sourced once, so the environment is shared by every session in the process.
.ai_report_runs <- new.env(parent = emptyenv())

#' Identity of a run
#'
#' The save target rather than the dataset name, because the path already
#' encodes whose copy of a dataset is being written. User and model are part of
#' the key as well: a run carries the starter's model and BYOK credentials, so
#' another user joining it would silently have their reports generated - and
#' billed - against someone else's account.
#'
#' @param save_path Absolute path of the target `.pgx` file.
#' @param token Dataset token from `ai_report_dataset_token()`.
#' @param user_key Stable per-user id (email, or whatever auth offers).
#' @param llm_model Model id.
#' @return Run key.
ai_report_run_key <- function(save_path, token, user_key, llm_model) {
  if (is.null(user_key)) user_key <- ""
  if (is.null(llm_model)) llm_model <- ""
  paste(normalizePath(save_path, mustWork = FALSE), token,
        as.character(user_key)[[1L]], as.character(llm_model)[[1L]],
        sep = "::")
}

#' Is a run for this key already in flight?
ai_report_run_active <- function(key) {
  run <- .ai_report_runs[[key]]
  !is.null(run) && isTRUE(run$running)
}

#' Snapshot of a run, or NULL if no run is registered under this key
#'
#' A plain list so a caller can read it outside any reactive context - used to
#' tell "already generating" apart from "never started" without touching the
#' run itself.
#'
#' @return `list(running, done, total, failed, slot, started)`, or NULL.
ai_report_run_status <- function(key) {
  run <- .ai_report_runs[[key]]
  if (is.null(run)) return(NULL)
  list(running = run$running, done = run$done, total = run$total,
       failed = run$failed, slot = run$slot, started = run$started)
}

#' Subscribe a session to a run's progress
#'
#' @param key Run key.
#' @param on_progress `function(done, total, slot, ok)`, or NULL.
#' @param on_done `function(result)` called once when the run finishes, where
#'   `result` carries `ai`, `done`, `failed` and `failures`.
#' @return A subscriber id to pass to `ai_report_run_unsubscribe()`, or NULL if
#'   there was no run to subscribe to.
ai_report_run_subscribe <- function(key, on_progress = NULL, on_done = NULL) {
  run <- .ai_report_runs[[key]]
  if (is.null(run)) return(NULL)
  ## Ids are minted here rather than supplied by the caller: two tabs of one
  ## Shiny session share a session token, and the second subscription would
  ## then silently replace the first one's callbacks.
  run$next_sub <- run$next_sub + 1L
  id <- paste0("sub", run$next_sub)
  run$subscribers[[id]] <- list(on_progress = on_progress, on_done = on_done)
  id
}

#' Drop a subscriber
#'
#' Called from `session$onSessionEnded`. The run itself is untouched - that is
#' the whole point of this module.
ai_report_run_unsubscribe <- function(key, sub_id) {
  run <- .ai_report_runs[[key]]
  if (is.null(run) || is.null(sub_id)) return(invisible(FALSE))
  run$subscribers[[sub_id]] <- NULL
  invisible(TRUE)
}

## Fan a callback out to every subscriber, guarding each one separately: a
## session that died mid-run will throw from its progress handle, and that must
## not stop the others from being told.
.ai_report_notify <- function(run, what, ...) {
  for (id in names(run$subscribers)) {
    fn <- run$subscribers[[id]][[what]]
    if (is.function(fn)) tryCatch(fn(...), error = function(e) NULL)
  }
  invisible(NULL)
}

## Commit what the run has produced so far, remembering which slots reached
## disk. The final commit is then free when the last phase already wrote
## exactly those slots - one pgx.load + pgx.save less on the Shiny thread.
.ai_report_commit_run <- function(run, save_path, ai) {
  if (identical(run$committed, names(ai))) return(TRUE)
  ok <- ai_report_commit_to_disk(save_path, ai)
  if (ok) run$committed <- names(ai)
  ok
}

## Telemetry is recorded here, not in the starting session, because the run may
## well outlive it. record_reports() only reads `name` and `ai`, so the run's
## own results stand in for the pgx - re-reading it from disk would cost
## another multi-second load for two fields. Event ids are derived from each
## slot's created_at, so a session that also records is deduplicated.
.ai_report_record_telemetry <- function(run, ai) {
  if (is.null(run$user_email)) return(invisible(NULL))
  tryCatch(
    ai_telemetry_record_reports(list(name = run$dataset, ai = ai),
                                user_email = run$user_email),
    error = function(e) NULL
  )
  invisible(NULL)
}

## The tail of a run: commit, tell the subscribers, unregister. Kept out of
## ai_report_run_start() so the lifetime of `run$running` is readable in one
## screen - every path through here settles it exactly once, before notifying.
.ai_report_run_finish <- function(run, key, save_path, token, promise) {
  settle <- function(result) {
    ## Not in finally(): promises schedules catch() and finally() on separate
    ## later() turns, and a join landing in that window would attach to a run
    ## that has already notified everyone and then wait forever.
    run$running <- FALSE
    .ai_report_notify(run, "on_done", result)
    NULL
  }

  promises::finally(
    promises::catch(
      promises::then(promise, function(result) {
        ## The engine resolves NULL when there was nothing to generate. That is
        ## still the end of the run, and subscribers are holding a progress bar
        ## open until they hear so.
        if (is.null(result)) {
          info("[AiReportManager] run produced nothing: key=", key)
          return(settle(list(ai = NULL, done = run$done, total = run$total,
                             failed = run$failed, failures = character(0),
                             error = "no reports were generated")))
        }
        ok <- .ai_report_commit_run(run, save_path, result$ai)
        if (ok) {
          .ai_report_record_telemetry(run, result$ai)
          .ai_report_clear_sidecars(save_path, token)
        }
        info("[AiReportManager] run finished: key=", key,
             " done=", result$done, " failed=", result$failed,
             " saved=", ok)
        settle(result)
      }),
      function(err) {
        info("[AiReportManager] run failed: key=", key,
             " error=", conditionMessage(err))
        settle(list(ai = NULL, done = run$done, total = run$total,
                    failed = run$failed, failures = character(0),
                    error = conditionMessage(err)))
      }
    ),
    function() {
      run$subscribers <- list()
      ## Only unregister ourselves: settling happens a turn earlier, so a
      ## session may already have started a fresh run under the same key.
      if (identical(.ai_report_runs[[key]], run)) {
        rm(list = key, envir = .ai_report_runs)
      }
    }
  )
}

#' Start an app-scope report run, or join one already in progress
#'
#' @param pgx_list Plain (non-reactive) pgx.
#' @param save_path Absolute path of the pgx to write. Resolved by the caller
#'   while it still has session context; the run never touches the session.
#' @param llm_model Model id.
#' @param select Modules to generate; NULL means everything the pgx supports.
#' @param credentials Nullary credential closure, or NULL.
#' @param user_key Stable per-user id; part of the run key.
#' @param user_email Captured at dispatch for telemetry, since the session that
#'   knows it may be gone by the time the run finishes.
#' @param on_progress,on_done Callbacks for the starting session; identical to
#'   the ones `ai_report_run_subscribe()` takes.
#' @return `list(key=, sub_id=)`, or NULL if there was nothing to generate or
#'   the request could not join the run already in flight.
ai_report_run_start <- function(pgx_list,
                                save_path,
                                llm_model,
                                select = NULL,
                                credentials = NULL,
                                user_key = NULL,
                                user_email = NULL,
                                force = FALSE,
                                on_progress = NULL,
                                on_done = NULL,
                                timeout_s = 600L,
                                retries = 3L) {
  if (is.null(save_path) || !nzchar(save_path)) return(NULL)
  if (is.null(llm_model) || !nzchar(llm_model)) return(NULL)

  ## Resolved here rather than left to the engine so that a joiner can be
  ## compared against what the running job will actually produce.
  if (is.null(select)) select <- ai_report_modules_for_pgx(pgx_list)
  select <- tryCatch(as.character(select), error = function(e) character(0))
  select <- select[!is.na(select) & nzchar(select)]
  if (!length(select)) return(NULL)

  token <- ai_report_dataset_token(pgx_list)
  key <- ai_report_run_key(save_path, token, user_key, llm_model)

  ## Join rather than duplicate. Two tabs on the same dataset would otherwise
  ## both pay for a full run and then race to save the pgx, with the loser's
  ## reports silently lost. Only a request the running job already covers may
  ## join, though: anything else would be told a run finished that never
  ## generated what it asked for.
  if (ai_report_run_active(key)) {
    run <- .ai_report_runs[[key]]
    if (!identical(isTRUE(force), isTRUE(run$force)) ||
        !all(select %in% run$select)) {
      info("[AiReportManager] declined join, request differs: key=", key)
      return(NULL)
    }
    sub_id <- ai_report_run_subscribe(key, on_progress, on_done)
    info("[AiReportManager] joined run in progress: key=", key)
    return(list(key = key, sub_id = sub_id))
  }

  run <- new.env(parent = emptyenv())
  run$running     <- TRUE
  run$done        <- 0L
  run$total       <- NA_integer_
  run$failed      <- 0L
  run$slot        <- NA_character_
  run$started     <- Sys.time()
  run$select      <- select
  run$force       <- isTRUE(force)
  run$dataset     <- pgx_list$name
  run$user_email  <- user_email
  run$subscribers <- list()
  run$next_sub    <- 0L
  run$committed   <- NULL
  .ai_report_runs[[key]] <- run
  sub_id <- ai_report_run_subscribe(key, on_progress, on_done)

  ## The run is not bound to the dataset still being on screen: results are
  ## written to the pgx they were generated from, addressed by path. Switching
  ## dataset mid-run therefore no longer throws the work away.
  promise <- ai_report_generate_async(
    pgx_list,
    llm_model   = llm_model,
    select      = select,
    force       = force,
    credentials = credentials,
    timeout_s   = timeout_s,
    retries     = retries,
    on_result = function(slot, text, usage) {
      ai_report_write_sidecar(save_path, token, slot, text, usage)
    },
    on_phase = function(phase, ai) {
      ## Durability boundary: a crash after this point costs the next phase,
      ## not the whole run. See the note on cost at the top of this file.
      ok <- .ai_report_commit_run(run, save_path, ai)
      info("[AiReportManager] phase committed: key=", key,
           " phase=", phase, " saved=", ok)
    },
    on_progress = function(done, total, slot, ok) {
      run$done <- done
      run$total <- total
      run$slot <- slot
      if (!isTRUE(ok) && !is.na(slot)) run$failed <- run$failed + 1L
      .ai_report_notify(run, "on_progress", done, total, slot, ok)
    }
  )

  .ai_report_run_finish(run, key, save_path, token, promise)

  info("[AiReportManager] run started: key=", key, " model=", llm_model)
  list(key = key, sub_id = sub_id)
}

#' Merge a finished `ai` slot into the pgx on disk
#'
#' Re-reads the pgx rather than writing the caller's in-memory copy: the run may
#' have outlived the session it started in, and anything else written to that
#' dataset in the meantime must not be clobbered by a stale snapshot.
#'
#' This runs on the Shiny main thread and a pgx.load + pgx.save pair costs
#' seconds on a large dataset, blocking every session in the process meanwhile.
#' That stall is the price of the run being durable; callers should invoke it
#' only at phase boundaries.
#'
#' @return TRUE if something was written, FALSE otherwise - including when
#'   there is nothing to merge.
ai_report_commit_to_disk <- function(save_path, ai) {
  if (!is.list(ai)) return(FALSE)
  slots <- names(ai)
  slots <- slots[!is.na(slots) & nzchar(slots)]
  ## A NULL element would delete the slot it names on assignment, quietly
  ## dropping a report that is already on disk.
  slots <- slots[!vapply(ai[slots], is.null, logical(1))]
  if (!length(slots) || !file.exists(save_path)) return(FALSE)
  tryCatch({
    pgx <- playbase::pgx.load(save_path)
    if (is.null(pgx)) return(FALSE)
    existing <- pgx$ai
    if (!is.list(existing)) existing <- list()
    for (slot in slots) existing[[slot]] <- ai[[slot]]
    pgx$ai <- existing
    playbase::pgx.save(pgx, file = save_path)
    TRUE
  }, error = function(e) {
    info("[AiReportManager] commit failed: path=", save_path,
         " error=", conditionMessage(e))
    FALSE
  })
}
