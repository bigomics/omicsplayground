## test-copilot-followups-inputs.R
##
## Unit tests for copilot_followups_inputs.R — the pure deterministic
## gatherers that build the follow-up helper's payload from an Agent.

.board_dir <- if (dir.exists("components/app_copilot/R")) {
  "components/app_copilot/R"
} else {
  "../../components/app_copilot/R"
}
.module_dir <- if (dir.exists("components/modules")) {
  "components/modules"
} else {
  "../../components/modules"
}

# AiReports + Reports server provide copilot_report_label(); context_blocks
# provides .copilot_dataset_context_text() used by build_followup_payload().
source(file.path(.module_dir, "AiReports.R"),                local = TRUE)
source(file.path(.board_dir,  "CopilotReportsServer.R"),    local = TRUE)
source(file.path(.board_dir,  "copilot_context_blocks.R"),  local = TRUE)
source(file.path(.board_dir,  "copilot_followups_inputs.R"), local = TRUE)

# ===========================================================================
# .first_sentence — leading-md strip, terminator detect, char cap
# ===========================================================================

test_that(".first_sentence strips leading markdown then takes first sentence", {
  expect_equal(.first_sentence("## Differential expression. Hub genes."), "Differential expression.")
  expect_equal(.first_sentence("- The dataset has 50 samples. Plenty of contrasts."), "The dataset has 50 samples.")
  expect_equal(.first_sentence("> quoted note. Tail."), "quoted note.")
})

test_that(".first_sentence caps at max chars when no terminator", {
  long <- paste(rep("x", 200L), collapse = "")
  expect_equal(nchar(.first_sentence(long, max = 50L)), 50L)
})

test_that(".first_sentence returns empty on NULL/NA/empty", {
  expect_equal(.first_sentence(NULL), "")
  expect_equal(.first_sentence(NA_character_), "")
  expect_equal(.first_sentence(""), "")
})

# ===========================================================================
# .trim_tool_output_body — drop Evidence/Primer, keep signal sections
# ===========================================================================

.sample_tool_output <- paste(
  "<!-- tool: query_de | case: ranked-up | dataset: example-data -->",
  "<!-- args: what=contrast:act48h_vs_notact -->",
  "",
  "Header:",
  "example-data | query_de | what=contrast:act48h_vs_notact",
  "",
  "Primer:",
  "Primer omitted by request.",
  "",
  "Headline:",
  "3169 proteins are significant at q < 0.05.",
  "",
  "Evidence:",
  "| name | description | fc | fdr |",
  "| --- | --- | --- | --- |",
  "| KCNN4 | potassium channel | 9.6 | 1.6e-11 |",
  "| ANLN | anillin | 9.3 | 1.9e-05 |",
  "",
  "Caveat:",
  "Output is bounded to top_n = 20 rows.",
  "",
  "== SUGGESTED PLOTS ==",
  "- `omicspgxmcp.show_omics_plot({\"plot_type\":\"volcano\",...})`",
  "  Why: visualize significance vs effect size",
  "",
  "Next:",
  "- `omicspgxmcp.query_genesets({\"what\":\"table:AGING\",...})`",
  "  Why: explain dominant biology",
  sep = "\n"
)

test_that(".trim_tool_output_body keeps Header/Headline/Caveat/SUGGESTED/Next", {
  out <- .trim_tool_output_body(.sample_tool_output)
  expect_match(out, "Header:",                fixed = TRUE)
  expect_match(out, "example-data | query_de", fixed = TRUE)
  expect_match(out, "Headline:",              fixed = TRUE)
  expect_match(out, "3169 proteins",          fixed = TRUE)
  expect_match(out, "Caveat:",                fixed = TRUE)
  expect_match(out, "== SUGGESTED PLOTS ==",  fixed = TRUE)
  expect_match(out, "show_omics_plot",        fixed = TRUE)
  expect_match(out, "Next:",                  fixed = TRUE)
  expect_match(out, "query_genesets",         fixed = TRUE)
})

test_that(".trim_tool_output_body drops the Evidence table and Primer", {
  out <- .trim_tool_output_body(.sample_tool_output)
  expect_false(grepl("Evidence:",            out, fixed = TRUE))
  expect_false(grepl("| KCNN4 |",            out, fixed = TRUE))
  expect_false(grepl("| ANLN |",             out, fixed = TRUE))
  expect_false(grepl("Primer omitted",       out, fixed = TRUE))
})

test_that(".trim_tool_output_body keeps the top-of-file <!-- comments", {
  out <- .trim_tool_output_body(.sample_tool_output)
  expect_match(out, "<!-- tool: query_de", fixed = TRUE)
  expect_match(out, "<!-- args:",          fixed = TRUE)
})

test_that(".trim_tool_output_body returns '' on empty/NULL", {
  expect_equal(.trim_tool_output_body(""), "")
  expect_equal(.trim_tool_output_body(NULL), "")
})

# ===========================================================================
# .extract_recent_tool_outputs — agent walking, ≤2 cap, last-user stop
# ===========================================================================

# Build a minimal agent stub: agent@chat$get_turns() returning a list of
# turn-like S7 objects. The walker only reads $role and $contents via .field,
# and treats any class containing "ContentToolResult" as a tool result.

.make_turn <- function(role, contents = list()) {
  list(role = role, contents = contents)
}
.make_tool_result <- function(name, value) {
  request <- structure(list(name = name, id = "req1"),
                       class = c("ellmer::ContentToolRequest"))
  structure(list(value = value, request = request),
            class = c("ellmer::ContentToolResult"))
}
.make_text <- function(text) {
  structure(list(text = text),
            class = c("ellmer::ContentText"))
}
# Use a real S7 class so `agent@chat` succeeds. The walker calls
# `chat$get_turns()` which works because we expose the function on an env.
.make_agent_stub <- function(turns) {
  env <- new.env(parent = emptyenv())
  env$get_turns <- function() turns
  AgC <- S7::new_class("FollowupTurnAgentStub",
                       properties = list(chat = S7::class_any))
  AgC(chat = env)
}

# .fu_field reads via S7::prop first; for our stubs it falls into the
# is.list() branch — verified by spot-checks below.

test_that(".extract_recent_tool_outputs returns list() for NULL agent", {
  expect_identical(.extract_recent_tool_outputs(NULL), list())
})

test_that(".extract_recent_tool_outputs returns list() when no tool results", {
  turns <- list(
    .make_turn("user",      list(.make_text("hi"))),
    .make_turn("assistant", list(.make_text("hello")))
  )
  agent <- .make_agent_stub(turns)
  expect_identical(.extract_recent_tool_outputs(agent), list())
})

test_that(".extract_recent_tool_outputs collects assistant tool results", {
  turns <- list(
    .make_turn("user",      list(.make_text("show de"))),
    .make_turn("assistant", list(
      .make_tool_result("query_de",
                       "Headline:\nstrong DE\n\nEvidence:\n| g | fc |\n\nNext:\n- foo")
    ))
  )
  agent <- .make_agent_stub(turns)
  out <- .extract_recent_tool_outputs(agent)
  expect_length(out, 1L)
  expect_equal(out[[1L]]$name, "query_de")
  expect_match(out[[1L]]$text, "Headline:", fixed = TRUE)
  expect_match(out[[1L]]$text, "Next:",     fixed = TRUE)
  expect_false(grepl("Evidence:", out[[1L]]$text, fixed = TRUE))
})

test_that(".extract_recent_tool_outputs caps at max_outputs", {
  contents <- list(
    .make_tool_result("query_de",       "Headline:\nA"),
    .make_tool_result("query_pathways", "Headline:\nB"),
    .make_tool_result("query_wgcna",    "Headline:\nC")
  )
  turns <- list(
    .make_turn("user",      list(.make_text("explore"))),
    .make_turn("assistant", contents)
  )
  agent <- .make_agent_stub(turns)
  out <- .extract_recent_tool_outputs(agent, max_outputs = 2L)
  expect_length(out, 2L)
  # Should keep the most recent two — pathways and wgcna.
  expect_equal(vapply(out, `[[`, character(1L), "name"),
               c("query_pathways", "query_wgcna"))
})

test_that(".extract_recent_tool_outputs stops at the last user turn", {
  # Earlier assistant turn-cluster should NOT bleed into the result.
  turns <- list(
    .make_turn("user",      list(.make_text("first"))),
    .make_turn("assistant", list(.make_tool_result("query_de", "Headline:\nold"))),
    .make_turn("user",      list(.make_text("now"))),
    .make_turn("assistant", list(.make_tool_result("query_pathways", "Headline:\nnew")))
  )
  agent <- .make_agent_stub(turns)
  out <- .extract_recent_tool_outputs(agent)
  expect_length(out, 1L)
  expect_equal(out[[1L]]$name, "query_pathways")
  expect_match(out[[1L]]$text, "new", fixed = TRUE)
})

# ===========================================================================
# .build_tool_catalog — degrades gracefully without omicsagentovi
# ===========================================================================

test_that(".build_tool_catalog returns '' when omicsagentovi is absent", {
  local_mocked_bindings(
    requireNamespace = function(pkg, ...) {
      if (identical(pkg, "omicsagentovi")) FALSE else TRUE
    },
    .package = "base"
  )
  expect_identical(.build_tool_catalog(), "")
})

test_that(".build_tool_catalog renders one line per tool", {
  skip_if_not_installed("omicsagentovi")
  out <- .build_tool_catalog()
  # Either we got an empty string (no registered tools in this test env) or a
  # well-formed multi-line bullet list.
  if (!nzchar(out)) succeed()
  else {
    expect_true(grepl("^-\\s", out))
    expect_true(grepl("\n-\\s|\n$|^[^\n]*$", out))  # ≥1 line, '- ' prefixed
  }
})

# ===========================================================================
# .build_report_inventory — walk pgx$ai, filter empty $report
# ===========================================================================

.fake_agent_with_ai <- function(ai) {
  pgx_obj <- structure(list(ai = ai), class = "pgx")
  ctx <- structure(list(pgx = pgx_obj), class = c("StubContext", "S7_object"))
  structure(list(context = ctx), class = c("StubAgent", "S7_object"))
}

# Our test agents are stub lists, but .build_report_inventory reads
# agent@context@pgx — emulate via a class with @ overload? Simpler: build a
# real-ish S7-like list and read via .fu_field's is.list() fallback.
# However .build_report_inventory uses `agent@context@pgx` directly, which
# requires S7 prop access. Wrap that inline using a tiny S7 class so the
# probe `agent@context@pgx` succeeds.

.make_real_agent_ai <- function(ai) {
  CtxC <- S7::new_class("FollowupCtxStub",
                        properties = list(pgx = S7::class_any))
  AgC  <- S7::new_class("FollowupAgentStub",
                        properties = list(context = CtxC, chat = S7::class_any))
  pgx  <- list(ai = ai)
  class(pgx) <- "pgx"
  AgC(context = CtxC(pgx = pgx), chat = NULL)
}

test_that(".build_report_inventory returns '' when no pgx$ai entries", {
  agent <- .make_real_agent_ai(NULL)
  expect_identical(.build_report_inventory(agent), "")
})

test_that(".build_report_inventory keeps entries with non-empty $report", {
  ai <- list(
    de       = list(report = "Differential expression. Hub genes detected."),
    pathways = list(report = ""),                # skipped
    wgcna    = list(report = "WGCNA modules. ME1 is biggest.")
  )
  agent <- .make_real_agent_ai(ai)
  out <- .build_report_inventory(agent)
  expect_match(out, "de",       fixed = TRUE)
  expect_match(out, "wgcna",    fixed = TRUE)
  expect_false(grepl("pathways", out, fixed = TRUE))
  # First-sentence extracts must be present.
  expect_match(out, "Differential expression.", fixed = TRUE)
  expect_match(out, "WGCNA modules.",           fixed = TRUE)
})

test_that(".build_report_inventory caps each entry's one-line summary", {
  long <- paste(rep("x", 500L), collapse = "")
  ai   <- list(combined = list(report = long))
  agent <- .make_real_agent_ai(ai)
  out   <- .build_report_inventory(agent)
  # The summary segment after the colon must be ≤120 chars.
  parts <- strsplit(out, ": ", fixed = TRUE)[[1L]]
  expect_lte(nchar(parts[[length(parts)]]), 120L)
})

test_that(".build_report_inventory filters to only_slots when non-empty", {
  ai <- list(
    de       = list(report = "Differential expression. Hub genes."),
    pathways = list(report = "Pathway enrichment. Top hallmarks."),
    wgcna    = list(report = "WGCNA modules. ME1 is biggest.")
  )
  agent <- .make_real_agent_ai(ai)

  # Only `pathways` ticked — `de` and `wgcna` must be filtered out.
  out <- .build_report_inventory(agent, only_slots = "pathways")
  expect_match(out, "pathways", fixed = TRUE)
  expect_false(grepl("^- de \\(",    out))
  expect_false(grepl("^- wgcna \\(", out))
})

test_that(".build_report_inventory with empty only_slots vector returns all", {
  ai <- list(
    de       = list(report = "DE report."),
    pathways = list(report = "Pathway report.")
  )
  agent <- .make_real_agent_ai(ai)
  # Empty filter -> no filter (back-compat for callers that don't gate)
  out <- .build_report_inventory(agent, only_slots = character(0))
  expect_match(out, "de",       fixed = TRUE)
  expect_match(out, "pathways", fixed = TRUE)
})

# ===========================================================================
# build_followup_payload — assembly + degraded inputs
# ===========================================================================

test_that("build_followup_payload returns 5 named blocks even with NULL agent", {
  out <- build_followup_payload(agent = NULL, last_text = "hello")
  expect_named(out,
               c("last_text", "dataset_context", "recent_tool_outputs",
                 "tool_catalog", "available_reports"),
               ignore.order = TRUE)
  expect_equal(out$last_text, "hello")
  expect_equal(out$dataset_context, "")
  expect_identical(out$recent_tool_outputs, list())
  expect_equal(out$available_reports, "")
})

test_that("build_followup_payload coerces NULL last_text to ''", {
  out <- build_followup_payload(agent = NULL, last_text = NULL)
  expect_equal(out$last_text, "")
})

test_that("build_followup_payload gates available_reports on attached slots", {
  ai <- list(
    de       = list(report = "Differential expression. Hub genes."),
    pathways = list(report = "Pathway enrichment. Top hallmarks."),
    combined = list(report = "Summary across modules.")
  )
  agent <- .make_real_agent_ai(ai)

  # Default (no slots ticked) → empty inventory.
  out0 <- build_followup_payload(agent = agent, last_text = "x")
  expect_equal(out0$available_reports, "")

  # `pathways` ticked → only `pathways` line, not `de` or `combined`.
  out1 <- build_followup_payload(agent = agent, last_text = "x",
                                 attached_report_slots = "pathways")
  expect_match(out1$available_reports, "pathways", fixed = TRUE)
  expect_false(grepl("^- de \\(",        out1$available_reports))
  expect_false(grepl("^- combined \\(",  out1$available_reports))

  # `de` + `combined` ticked → both lines, no `pathways`.
  out2 <- build_followup_payload(agent = agent, last_text = "x",
                                 attached_report_slots = c("de", "combined"))
  expect_match(out2$available_reports, "de",       fixed = TRUE)
  expect_match(out2$available_reports, "combined", fixed = TRUE)
  expect_false(grepl("^- pathways \\(",  out2$available_reports))
})
