#' Copilot History Tools
#'
#' Helper functions for session transcript projection, JSON export,
#' and replay policy decisions. All functions are pure (no reactives)
#' except \code{copilot_replay_turns()} which calls \code{shinychat}.

#' Build a transcript record from an ellmer Turn object
#'
#' @param turn One ellmer Turn object.
#' @param idx Integer turn index.
#' @return List record with text/tool/thinking channels.
copilot_record_from_turn <- function(turn, idx) {
  role <- tryCatch(turn@role, error = function(e) NA_character_)
  contents <- tryCatch(turn@contents, error = function(e) list())

  text_parts <- character(0)
  tool_requests <- list()
  tool_results <- list()
  thinking <- character(0)

  for (ct in contents) {
    if (S7::S7_inherits(ct, ellmer::ContentText)) {
      text_parts <- c(text_parts, ct@text %||% "")
    } else if (S7::S7_inherits(ct, ellmer::ContentToolRequest)) {
      tool_requests[[length(tool_requests) + 1L]] <- list(
        id = ct@id %||% "",
        name = ct@name %||% "",
        arguments = ct@arguments %||% list()
      )
    } else if (S7::S7_inherits(ct, ellmer::ContentToolResult)) {
      tool_results[[length(tool_results) + 1L]] <- list(
        request_id = tryCatch(ct@request@id, error = function(e) NA_character_),
        value = ct@value %||% ""
      )
    } else if (S7::S7_inherits(ct, ellmer::ContentThinking)) {
      thinking <- c(thinking, ct@thinking %||% "")
    }
  }

  text <- paste(text_parts[nzchar(text_parts)], collapse = "\n")
  content_types <- vapply(contents, function(x) class(x)[1], character(1))

  list(
    idx = as.integer(idx),
    role = role,
    text = text,
    has_text = nzchar(text),
    has_tool_request = length(tool_requests) > 0L,
    has_tool_result = length(tool_results) > 0L,
    has_tool = (length(tool_requests) + length(tool_results)) > 0L,
    has_thinking = any(nzchar(thinking)),
    content_types = content_types,
    tool_requests = tool_requests,
    tool_results = tool_results,
    thinking = thinking[nzchar(thinking)]
  )
}

#' Build a transcript record from one flattened SQLite turn row
#'
#' @param row Named list with turns table columns.
#' @return List record with text/tool/thinking channels.
copilot_record_from_row <- function(row) {
  ## %||% guards NULL (missing list element); has_val/is.na guard NA from SQL.
  has_val <- function(x) !is.null(x) && !is.na(x) && nzchar(x)

  role <- row$role %||% NA_character_
  is_tool_request <- has_val(row$tool_name)
  is_tool_result <- !is_tool_request && identical(role, "user") && has_val(row$tool_call_id)

  tool_requests <- list()
  if (is_tool_request) {
    args <- if (has_val(row$tool_args)) {
      tryCatch(jsonlite::fromJSON(row$tool_args, simplifyVector = FALSE), error = function(e) list())
    } else {
      list()
    }
    tool_requests[[1L]] <- list(
      id = row$tool_call_id %||% "",
      name = row$tool_name %||% "",
      arguments = args
    )
  }

  tool_results <- list()
  if (is_tool_result) {
    tool_results[[1L]] <- list(
      request_id = row$tool_call_id %||% "",
      value = row$content %||% ""
    )
  }

  text <- if (is_tool_result) "" else {
    val <- row$content %||% ""
    if (is.na(val)) "" else as.character(val)
  }
  think <- row$thinking %||% ""
  if (is.na(think)) think <- ""

  list(
    idx = as.integer(row$idx %||% NA_integer_),
    role = role,
    text = text,
    has_text = nzchar(text),
    has_tool_request = is_tool_request,
    has_tool_result = is_tool_result,
    has_tool = is_tool_request || is_tool_result,
    has_thinking = nzchar(think),
    content_types = if (is_tool_request) {
      if (nzchar(text)) c("ellmer::ContentText", "ellmer::ContentToolRequest") else "ellmer::ContentToolRequest"
    } else if (is_tool_result) {
      "ellmer::ContentToolResult"
    } else {
      if (nzchar(text)) "ellmer::ContentText" else character(0)
    },
    tool_requests = tool_requests,
    tool_results = tool_results,
    thinking = if (nzchar(think)) think else character(0)
  )
}

#' Decide whether a turn record should be shown in the chat replay.
#'
#' "current" policy hides any turn containing tool content (request or result).
#' "expected" policy only hides tool-result turns, so assistant messages that
#' mixed text + a tool call still show their text portion.
#'
#' @param rec Turn record from \code{copilot_record_from_turn/row()}.
#' @param policy \code{"current"} or \code{"expected"}.
#' @return List with \code{show} (logical) and \code{hidden_reasons} (character).
.copilot_check_turn_visibility <- function(rec, policy = c("current", "expected")) {
  policy <- match.arg(policy)
  role_ok <- !is.null(rec$role) && !is.na(rec$role) && rec$role %in% c("user", "assistant")

  hidden_reasons <- character(0)
  if (!role_ok) hidden_reasons <- c(hidden_reasons, "non_chat_role")
  if (!rec$has_text) hidden_reasons <- c(hidden_reasons, "no_text_content")
  if (rec$has_thinking) hidden_reasons <- c(hidden_reasons, "contains_thinking")

  if (identical(policy, "current")) {
    if (rec$has_tool) hidden_reasons <- c(hidden_reasons, "contains_tool_content")
  } else {
    if (rec$has_tool_result) hidden_reasons <- c(hidden_reasons, "internal_tool_result")
  }

  ## show is derived solely from hidden_reasons — single source of truth
  list(show = length(hidden_reasons) == 0L, hidden_reasons = unique(hidden_reasons))
}

#' Build a chat transcript projection from turn records
#'
#' @param records List of records from \code{copilot_record_from_turn()} or
#'   \code{copilot_record_from_row()}.
#' @param policy Replay policy: \code{"current"} reproduces existing behavior,
#'   \code{"expected"} keeps user-visible text while hiding tool-result payloads.
#' @return List containing counts, replay messages, and foldable transcript entries.
copilot_build_transcript_from_records <- function(records, policy = c("current", "expected")) {
  policy <- match.arg(policy)
  entries <- lapply(records, function(rec) {
    vis <- .copilot_check_turn_visibility(rec, policy = policy)
    list(
      idx = rec$idx %||% NA_integer_,
      role = rec$role %||% NA_character_,
      render = list(
        chat_message = if (rec$has_text) rec$text else NULL,
        show_by_default = vis$show,
        expandable = list(
          has_tool_request = rec$has_tool_request,
          has_tool_result = rec$has_tool_result,
          has_thinking = rec$has_thinking
        )
      ),
      hidden = list(
        reasons = vis$hidden_reasons,
        thinking = rec$thinking %||% character(0),
        tool_requests = rec$tool_requests %||% list(),
        tool_results = rec$tool_results %||% list()
      ),
      raw = list(
        has_text = rec$has_text,
        has_tool = rec$has_tool,
        content_types = rec$content_types %||% character(0)
      )
    )
  })

  replay_messages <- lapply(Filter(function(x) x$render$show_by_default, entries), function(x) {
    list(
      idx = x$idx,
      role = x$role,
      content = x$render$chat_message %||% ""
    )
  })

  list(
    policy = policy,
    counts = list(
      total_turns = length(entries),
      visible_messages = length(replay_messages),
      hidden_messages = length(entries) - length(replay_messages),
      turns_with_tool_request = sum(vapply(records, function(x) x$has_tool_request, logical(1))),
      turns_with_tool_result = sum(vapply(records, function(x) x$has_tool_result, logical(1))),
      turns_with_thinking = sum(vapply(records, function(x) x$has_thinking, logical(1)))
    ),
    replay_messages = replay_messages,
    messages = entries
  )
}

#' Build a transcript projection directly from ellmer Turn objects
#'
#' @param turns List of ellmer Turn objects.
#' @param policy Replay policy (\code{"current"} or \code{"expected"}).
#' @return Transcript projection list.
copilot_build_transcript_from_turns <- function(turns, policy = c("current", "expected")) {
  records <- lapply(seq_along(turns), function(i) copilot_record_from_turn(turns[[i]], i))
  copilot_build_transcript_from_records(records, policy = policy)
}

#' Build a transcript projection directly from a persisted session in SQLite
#'
#' Provisional — no callers yet. Intended for future session-export and
#' debugging workflows. Opens its own DBI connection; the column set must
#' stay in sync with the omicsagentovi turns table schema.
#'
#' @param session_id Session identifier.
#' @param session_dir Directory containing \code{sessions.sqlite}.
#' @param policy Replay policy (\code{"current"} or \code{"expected"}).
#' @return Transcript projection list.
copilot_build_transcript_from_session <- function(session_id, session_dir, policy = c("current", "expected")) {
  policy <- match.arg(policy)
  db_path <- file.path(session_dir, "sessions.sqlite")
  if (!file.exists(db_path)) {
    stop("sessions.sqlite not found in session_dir: ", session_dir, call. = FALSE)
  }

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  rows <- DBI::dbGetQuery(
    con,
    "SELECT idx, role, content, tool_name, tool_call_id, tool_args, thinking, json_payload, token_count
     FROM turns WHERE session_id = ? ORDER BY idx",
    params = list(session_id)
  )

  records <- lapply(seq_len(nrow(rows)), function(i) copilot_record_from_row(as.list(rows[i, , drop = FALSE])))
  transcript <- copilot_build_transcript_from_records(records, policy = policy)
  transcript$session_id <- session_id
  transcript
}

#' Serialize a transcript projection to JSON
#'
#' Provisional — no callers yet. Intended for future session-export UI.
#'
#' @param transcript Transcript list from \code{copilot_build_transcript_*()}.
#' @param pretty Logical; pretty-print JSON.
#' @return JSON string.
copilot_transcript_to_json <- function(transcript, pretty = TRUE) {
  jsonlite::toJSON(
    transcript,
    auto_unbox = TRUE,
    pretty = pretty,
    null = "null"
  )
}

#' Compute replayable chat messages without touching Shiny UI
#'
#' @param turns List of ellmer Turn objects.
#' @param policy Replay policy (\code{"current"} or \code{"expected"}).
#' @return List of message objects with \code{idx}, \code{role}, \code{content}.
copilot_collect_replay_messages <- function(turns, policy = c("current", "expected")) {
  policy <- match.arg(policy)
  transcript <- copilot_build_transcript_from_turns(turns, policy = policy)
  transcript$replay_messages
}

#' Replay persisted user/assistant text turns into a shinychat instance
#'
#' Appends visible messages to the shinychat widget via
#' \code{shinychat::chat_append_message()}. This is the only function in
#' this file that touches Shiny UI.
#'
#' Policy controls which turns are shown:
#' \itemize{
#'   \item \code{"current"} — hides any turn containing tool content (request
#'     or result). This was the original behaviour and drops assistant messages
#'     that mixed text with a tool call.
#'   \item \code{"expected"} — only hides tool-result turns. Assistant messages
#'     that contain both text and a tool request still show their text, so the
#'     full user-visible conversation is restored.
#' }
#'
#' @param chat_id The shinychat widget id to append messages to.
#' @param turns List of ellmer Turn objects from a restored session.
#' @param policy Replay policy (\code{"current"} or \code{"expected"}).
#' @param replay_messages Optional precomputed messages from
#'   \code{copilot_collect_replay_messages()} to avoid recomputation.
#' @param mode Replay mode: \code{"single"} appends one message at a time via
#'   \code{chat_append_message()}, \code{"batch"} sends all messages in one
#'   custom payload and lets the browser append them in chunks.
#' @param session Shiny session used for \code{mode="batch"}.
#' @param batch_size Number of messages appended per browser frame in batch mode.
#' @param done_input_id Optional input id set by the browser when batch replay
#'   finishes. Use a fully-qualified id (e.g. \code{session$ns("chat_replay_done")}).
#' @return Integer count of user turns replayed.
copilot_replay_turns <- function(chat_id,
                                 turns,
                                 policy = c("current", "expected"),
                                 replay_messages = NULL,
                                 mode = c("single", "batch"),
                                 session = shiny::getDefaultReactiveDomain(),
                                 batch_size = 32L,
                                 done_input_id = NULL) {
  policy <- match.arg(policy)
  mode <- match.arg(mode)
  if (is.null(replay_messages)) {
    replay_messages <- copilot_collect_replay_messages(turns, policy = policy)
  }
  n_user <- 0L

  if (identical(mode, "batch") && !is.null(session) && is.function(session$sendCustomMessage)) {
    resolved_id <- if (is.function(session$ns)) session$ns(chat_id) else chat_id
    payload_messages <- lapply(replay_messages, function(msg) {
      list(
        role = msg$role,
        content = as.character(msg$content %||% "")
      )
    })
    session$sendCustomMessage("copilot-chat-batch-append", list(
      id = resolved_id,
      messages = payload_messages,
      batch_size = as.integer(batch_size),
      done_input_id = done_input_id
    ))
    n_user <- sum(vapply(replay_messages, function(msg) identical(msg$role, "user"), logical(1L)))
    return(as.integer(n_user))
  }

  for (msg in replay_messages) {
    shinychat::chat_append_message(chat_id, list(role = msg$role, content = msg$content), chunk = FALSE)
    if (identical(msg$role, "user")) {
      n_user <- n_user + 1L
    }
  }
  n_user
}
