## test-copilot-followups.R
##
## Tests for copilot_followups.R — parse_followup_list, format_followup_bubble,
## make_followup_generator.

.board_dir <- if (dir.exists("components/app_copilot/R")) {
  "components/app_copilot/R"
} else {
  "../../components/app_copilot/R"
}

source(file.path(.board_dir, "copilot_followups.R"), local = TRUE)

# ===========================================================================
# parse_followup_list — 9 cases
# ===========================================================================

test_that("numbered list '1. q1\\n2. q2' parses to character(2)", {
  expect_equal(
    parse_followup_list("1. First question?\n2. Second question?", n = 2L),
    c("First question?", "Second question?")
  )
})

test_that("parenthesized '1) q1\\n2) q2' parses to character(2)", {
  expect_equal(
    parse_followup_list("1) First?\n2) Second?", n = 2L),
    c("First?", "Second?")
  )
})

test_that("dash bullets '- q1\\n- q2' parse to character(2)", {
  expect_equal(
    parse_followup_list("- First?\n- Second?", n = 2L),
    c("First?", "Second?")
  )
})

test_that("markdown-wrapped items strip ** and * formatting", {
  expect_equal(
    parse_followup_list("1. **First?**\n2. *Second?*", n = 2L),
    c("First?", "Second?")
  )
})

test_that("single-line unenumerated text returns character(0)", {
  expect_equal(parse_followup_list("Just one sentence."), character(0))
})

test_that("extra items beyond n are truncated", {
  out <- parse_followup_list("1. a\n2. b\n3. c\n4. d", n = 2L)
  expect_equal(out, c("a", "b"))
})

test_that("fewer than n returns what was found, no padding", {
  expect_equal(parse_followup_list("1. only one", n = 2L), "only one")
})

test_that("empty / NULL / NA input returns character(0)", {
  expect_equal(parse_followup_list(""), character(0))
  expect_equal(parse_followup_list(NULL), character(0))
  expect_equal(parse_followup_list(character(0)), character(0))
})

test_that("bullet chars (*, •) are recognised", {
  expect_equal(
    parse_followup_list("* First?\n• Second?", n = 2L),
    c("First?", "Second?")
  )
})

# ===========================================================================
# format_followup_bubble — HTML escape (XSS smoke), shape
# ===========================================================================

test_that("format_followup_bubble produces <ul><li class='suggestion submit'>", {
  html <- format_followup_bubble(c("q1?", "q2?"))
  expect_match(html, "^<ul>")
  expect_match(html, "</ul>$")
  expect_match(html, "<li class='suggestion submit'>q1\\?</li>")
  expect_match(html, "<li class='suggestion submit'>q2\\?</li>")
})

test_that("format_followup_bubble HTML-escapes untrusted LLM output (XSS smoke)", {
  html <- format_followup_bubble("<script>alert(1)</script>")
  # Live <script> tag must not survive — escaped to &lt;script&gt;
  expect_false(grepl("<script>", html, fixed = TRUE))
  expect_match(html, "&lt;script&gt;")
})

test_that("format_followup_bubble returns empty string on length-0 input", {
  expect_identical(format_followup_bubble(character(0)), "")
})

# ===========================================================================
# make_followup_generator — runtime-optional dependency
# ===========================================================================

test_that("make_followup_generator returns NULL when omicsai is absent", {
  # Mock requireNamespace via local trace
  local_mocked_bindings(
    requireNamespace = function(pkg, ...) {
      if (identical(pkg, "omicsai")) FALSE else TRUE
    },
    .package = "base"
  )
  expect_null(make_followup_generator())
})

test_that("make_followup_generator returns a list with generate() when omicsai present", {
  skip_if_not_installed("omicsai")
  skip_if_not_installed("promises")
  gen <- make_followup_generator()
  expect_true(is.list(gen))
  expect_true(is.function(gen$generate))
})

test_that("generate(empty) resolves to character(0) without calling the LLM", {
  skip_if_not_installed("omicsai")
  skip_if_not_installed("promises")
  gen <- make_followup_generator()
  if (is.null(gen)) skip("omicsai not loadable in this env")

  result <- NULL
  p <- gen$generate("")
  promises::then(p, onFulfilled = function(v) result <<- v)
  later::run_now(2)
  expect_equal(result, character(0))
})

test_that("generate accepts a payload list with empty last_text", {
  skip_if_not_installed("omicsai")
  skip_if_not_installed("promises")
  gen <- make_followup_generator()
  if (is.null(gen)) skip("omicsai not loadable in this env")

  result <- NULL
  p <- gen$generate(list(last_text = ""))
  promises::then(p, onFulfilled = function(v) result <<- v)
  later::run_now(2)
  expect_equal(result, character(0))
})

# ===========================================================================
# render_followup_prompt — sectioned template substitution
# ===========================================================================

test_that("render_followup_prompt substitutes every named block", {
  payload <- list(
    last_text           = "The most upregulated gene is KCNN4.",
    dataset_context     = "name: example\norganism: Homo sapiens",
    tool_catalog        = "- query_de: differential expression",
    available_reports   = "- de (Differential Expression): top genes",
    recent_tool_outputs = list(list(name = "query_de", text = "Headline:\nstrong DE"))
  )
  out <- render_followup_prompt(payload, n = 2L)

  expect_match(out, "Exactly 2 questions", fixed = TRUE)
  expect_match(out, "KCNN4",                fixed = TRUE)
  expect_match(out, "Homo sapiens",         fixed = TRUE)
  expect_match(out, "query_de: differential", fixed = TRUE)
  expect_match(out, "Differential Expression", fixed = TRUE)
  expect_match(out, "## query_de",          fixed = TRUE)
  expect_match(out, "strong DE",            fixed = TRUE)
})

test_that("render_followup_prompt fills missing blocks with '(none available)'", {
  payload <- list(
    last_text = "hello",
    dataset_context = "",
    tool_catalog = "",
    available_reports = "",
    recent_tool_outputs = list()
  )
  out <- render_followup_prompt(payload, n = 2L)

  expect_match(out, "(none available)", fixed = TRUE)
  expect_match(out, "(none)",           fixed = TRUE)   # tool outputs block
  expect_match(out, "hello",            fixed = TRUE)
})

test_that("render_followup_prompt encodes the anti-hallucination rules", {
  payload <- list(last_text = "x", dataset_context = "", tool_catalog = "",
                  available_reports = "", recent_tool_outputs = list())
  out <- render_followup_prompt(payload, n = 2L)
  expect_match(out, "NEVER output JSON",            fixed = TRUE)
  expect_match(out, "Never invent",                 fixed = TRUE)
  expect_match(out, "wet-lab metadata",             fixed = TRUE)
  expect_match(out, "show_omics_plot",              fixed = TRUE)
})

test_that("render_followup_prompt encodes the plot-suggestion hard rule", {
  payload <- list(last_text = "x", dataset_context = "", tool_catalog = "",
                  available_reports = "", recent_tool_outputs = list())
  out <- render_followup_prompt(payload, n = 2L)
  # Headline of the rule
  expect_match(out, "PLOT-SUGGESTION HARD RULE", fixed = TRUE)
  # The empty-case prohibition must be present verbatim
  expect_match(out, "DO NOT SUGGEST ANY PLOTTING FOLLOW-UP", fixed = TRUE)
  # Forbidden phrases must be enumerated
  expect_match(out, "volcano plot",   fixed = TRUE)
  expect_match(out, "heatmap",        fixed = TRUE)
  expect_match(out, "show_omics_plot", fixed = TRUE)
  # Discovery-question escape hatch is also there
  expect_match(out, "list the available contrasts", fixed = TRUE)
})

test_that("render_followup_prompt encodes the reports hard rule", {
  payload <- list(last_text = "x", dataset_context = "", tool_catalog = "",
                  available_reports = "", recent_tool_outputs = list())
  out <- render_followup_prompt(payload, n = 2L)
  # Headline + gate semantics
  expect_match(out, "REPORTS HARD RULE",                 fixed = TRUE)
  expect_match(out, "Attached precomputed reports",      fixed = TRUE)
  expect_match(out, "TICKED by the user",                fixed = TRUE)
  # Both-empty case spelled out
  expect_match(out, "BOTH the `Recent tool outputs`",    fixed = TRUE)
  expect_match(out, "pure",                              fixed = TRUE)
  expect_match(out, "tool-discovery",                    fixed = TRUE)
})

test_that("render_followup_prompt enforces user-voice phrasing", {
  payload <- list(last_text = "x", dataset_context = "", tool_catalog = "",
                  available_reports = "", recent_tool_outputs = list())
  out <- render_followup_prompt(payload, n = 2L)
  # Positive voice signals
  expect_match(out, "user's own voice",        fixed = TRUE)
  expect_match(out, "Can you",                 fixed = TRUE)
  expect_match(out, "Show me",                 fixed = TRUE)
  # Anti-patterns explicitly called out
  expect_match(out, "Would you like to see",   fixed = TRUE)
  expect_match(out, "Are you interested in",   fixed = TRUE)
  expect_match(out, "user IS the asker",       fixed = TRUE)
})
