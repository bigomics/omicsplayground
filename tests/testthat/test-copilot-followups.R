## test-copilot-followups.R
##
## Tests for copilot_followups.R — parse_followup_list, format_followup_bubble,
## make_followup_generator.

.board_dir <- if (dir.exists("components/board.copilot/R")) {
  "components/board.copilot/R"
} else {
  "../../components/board.copilot/R"
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
