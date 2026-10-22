## test-copilot-plot-render.R
##
## Pure-function tests for copilot_plot_render.R helpers.
## No Shiny server loop required.
##
## copilot_build_plot is integration-only (depends on omicspgxmcp::: internals
## and a real PGX object) and is NOT tested here. It is exercised by the
## board-level smoke test when a real session is available.

.board_dir <- if (dir.exists("components/board.copilot/R")) {
  "components/board.copilot/R"
} else {
  "../../components/board.copilot/R"
}

source(file.path(.board_dir, "copilot_options.R"),  local = TRUE)
source(file.path(.board_dir, "copilot_messages.R"), local = TRUE)
source(file.path(.board_dir, "copilot_logger.R"),   local = TRUE)
source(file.path(.board_dir, "copilot_plot_render.R"), local = TRUE)

# ===========================================================================
# copilot_detect_plot_kind
# ===========================================================================

test_that("copilot_detect_plot_kind returns 'ggplot' for a ggplot object", {
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
         ggplot2::geom_point()
  expect_equal(copilot_detect_plot_kind(p), "ggplot")
})

test_that("copilot_detect_plot_kind returns 'plotly' for a plotly object", {
  p <- plotly::plot_ly(x = 1:3, y = 1:3, type = "scatter", mode = "lines")
  expect_equal(copilot_detect_plot_kind(p), "plotly")
})

test_that("copilot_detect_plot_kind returns NULL for an unrecognised object", {
  expect_null(copilot_detect_plot_kind(list(foo = "bar")))
  expect_null(copilot_detect_plot_kind(42L))
})

# ===========================================================================
# copilot_parse_features
# ===========================================================================

test_that("copilot_parse_features splits a comma-separated string", {
  out <- copilot_parse_features("a,b,c")
  expect_equal(out, c("a", "b", "c"))
})

test_that("copilot_parse_features trims whitespace around names", {
  out <- copilot_parse_features(" TP53 , BRCA1 , MYC ")
  expect_equal(out, c("TP53", "BRCA1", "MYC"))
})

test_that("copilot_parse_features returns NULL for NULL input", {
  expect_null(copilot_parse_features(NULL))
})

test_that("copilot_parse_features returns NULL for empty/whitespace string", {
  expect_null(copilot_parse_features(""))
  expect_null(copilot_parse_features("   "))
})

test_that("copilot_parse_features drops empty tokens from trailing commas", {
  out <- copilot_parse_features("TP53,,MYC,")
  expect_equal(out, c("TP53", "MYC"))
})

# ===========================================================================
# copilot_prerender_ggplot
# ===========================================================================

test_that("copilot_prerender_ggplot writes a non-empty PNG and returns path", {
  p    <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
            ggplot2::geom_point()
  path <- copilot_prerender_ggplot(p)
  on.exit(copilot_prerender_cleanup(path))
  expect_true(file.exists(path))
  expect_gt(file.size(path), 0L)
  expect_match(path, "\\.png$", ignore.case = TRUE)
})

# ===========================================================================
# copilot_prerender_cleanup
# ===========================================================================

test_that("copilot_prerender_cleanup removes existing PNG files", {
  # Create a real file to be removed
  p    <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
            ggplot2::geom_point()
  path <- copilot_prerender_ggplot(p)
  expect_true(file.exists(path))
  copilot_prerender_cleanup(path)
  expect_false(file.exists(path))
})

test_that("copilot_prerender_cleanup is silent for NA paths", {
  expect_invisible(copilot_prerender_cleanup(NA_character_))
})

test_that("copilot_prerender_cleanup is silent for missing paths", {
  expect_invisible(copilot_prerender_cleanup("/nonexistent/path/plot.png"))
})

test_that("copilot_prerender_cleanup is silent for NULL paths", {
  expect_invisible(copilot_prerender_cleanup(NULL))
})

test_that("copilot_prerender_cleanup handles a mixed vector with NAs and real paths", {
  p    <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
            ggplot2::geom_point()
  path <- copilot_prerender_ggplot(p)
  expect_true(file.exists(path))
  copilot_prerender_cleanup(c(NA_character_, path, "/nonexistent/foo.png"))
  expect_false(file.exists(path))
})
