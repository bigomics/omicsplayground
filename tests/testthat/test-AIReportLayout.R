##
## Test suite for AI Report Layout module
##

context("AI Report Layout")

test_that("ai_report_layout_ui creates proper UI structure", {
  # Source the module from new location
  source("../../components/board.wgcna/R/ai.report/ai_report_ui.R")

  # Call UI function
  ui <- ai_report_layout_ui("test_id")

  # Test that UI is returned
  expect_true(!is.null(ui))

  # Test that it's a shiny tag
  expect_true(inherits(ui, "shiny.tag") || inherits(ui, "shiny.tag.list"))
})

test_that("ai_report_layout_ui accepts custom titles", {
  source("../../components/board.wgcna/R/ai.report/ai_report_ui.R")

  # Call with custom titles
  ui <- ai_report_layout_ui(
    "test_id",
    text_title = "Custom Report",
    diagram_title = "Custom Diagram",
    infographic_title = "Custom Infographic"
  )

  # Should not error
  expect_true(!is.null(ui))
})

test_that("ai_report_layout_server function exists and has correct signature", {
  # Source all ai.report files
  for (f in list.files("../../components/board.wgcna/R/ai.report",
                      pattern = "\\.R$", full.names = TRUE)) {
    source(f)
  }

  # Test function exists
  expect_true(exists("ai_report_layout_server"))

  # Test it's a function
  expect_true(is.function(ai_report_layout_server))

  # Test it has expected parameters (updated signature)
  args <- formals(ai_report_layout_server)
  expect_true("id" %in% names(args))
  expect_true("text_reactive" %in% names(args))
  expect_true("diagram_result_reactive" %in% names(args))
  expect_true("infographic_reactive" %in% names(args))
  expect_true("diagram_colors" %in% names(args))
  expect_true("watermark" %in% names(args))

  # Test default values
  expect_equal(args$watermark, FALSE)
  expect_null(args$diagram_colors)
})

test_that("module follows Shiny module pattern", {
  source("../../components/board.wgcna/R/ai.report/ai_report_ui.R")

  # UI function should use NS
  ui_code <- deparse(body(ai_report_layout_ui))
  expect_true(any(grepl("NS\\(id\\)", ui_code)))
  expect_true(any(grepl("ns\\(", ui_code)))

  # Server function should use moduleServer
  for (f in list.files("../../components/board.wgcna/R/ai.report",
                      pattern = "\\.R$", full.names = TRUE)) {
    source(f)
  }
  server_code <- deparse(body(ai_report_layout_server))
  expect_true(any(grepl("moduleServer", server_code)))
})

test_that("module integrates with PlotModule correctly", {
  source("../../components/board.wgcna/R/ai.report/ai_report_ui.R")

  ui_code <- deparse(body(ai_report_layout_ui))

  for (f in list.files("../../components/board.wgcna/R/ai.report",
                      pattern = "\\.R$", full.names = TRUE)) {
    source(f)
  }
  server_code <- deparse(body(ai_report_layout_server))

  # UI should call PlotModuleUI
  expect_true(any(grepl("PlotModuleUI", ui_code)))

  # Server should call PlotModuleServer
  expect_true(any(grepl("PlotModuleServer", server_code)))
})

test_that("module implements three panels with correct plotlib types", {
  for (f in list.files("../../components/board.wgcna/R/ai.report",
                      pattern = "\\.R$", full.names = TRUE)) {
    source(f)
  }

  server_code <- deparse(body(ai_report_layout_server))

  # Should have PlotModuleServer calls for each panel type
  expect_true(any(grepl("plotlib.*=.*\"generic\"", server_code)))
  expect_true(any(grepl("plotlib.*=.*\"visnetwork\"", server_code)))
  expect_true(any(grepl("plotlib.*=.*\"image\"", server_code)))
})

test_that("diagram panel uses visNetwork rendering via omicsai", {
  for (f in list.files("../../components/board.wgcna/R/ai.report",
                      pattern = "\\.R$", full.names = TRUE)) {
    source(f)
  }

  server_code <- deparse(body(ai_report_layout_server))

  # Should call omicsai_diagram_render
  expect_true(any(grepl("omicsai_diagram_render", server_code)))

  # Should NOT reference DiagrammeR or svgPanZoom
  expect_false(any(grepl("DiagrammeR", server_code)))
  expect_false(any(grepl("svgPanZoom", server_code)))
})

test_that("module includes validation messages", {
  for (f in list.files("../../components/board.wgcna/R/ai.report",
                      pattern = "\\.R$", full.names = TRUE)) {
    source(f)
  }

  server_code <- deparse(body(ai_report_layout_server))

  # Should use shiny::validate and shiny::need
  expect_true(any(grepl("shiny::validate", server_code)))
  expect_true(any(grepl("shiny::need", server_code)))

  # Should have user-friendly messages
  expect_true(any(grepl("Click.*Generate", server_code)))
  expect_true(any(grepl("will appear", server_code)))
})
