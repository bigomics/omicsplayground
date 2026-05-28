## test-copilot-evidence-server.R
##
## Module-server tests for CopilotEvidenceServer using shiny::testServer.
## Covers the public API: append_artifact, clear/clear_plots, update_table,
## clear_table, carousel navigation, and conditional-panel state.
##
## Conditional-panel output flags (is_ggplot_active, etc.) are validated by
## inspecting active_artifact() directly via session$returned$active — the
## public reactive exposed by the module. output$evidence_ggplot (renderImage)
## IS accessible via testServer since it uses a proper render function.

.board_dir <- if (dir.exists("components/board.copilot/R")) {
  "components/board.copilot/R"
} else {
  "../../components/board.copilot/R"
}

source(file.path(.board_dir, "copilot_options.R"),  local = TRUE)
source(file.path(.board_dir, "copilot_messages.R"), local = TRUE)
source(file.path(.board_dir, "copilot_logger.R"),   local = TRUE)
source(file.path(.board_dir, "copilot_plot_render.R"),    local = TRUE)
source(file.path(.board_dir, "CopilotEvidenceServer.R"),  local = TRUE)

library(shiny)

# ---- Helpers ----

make_ggplot_record <- function(label = "pca plot") {
  p    <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
            ggplot2::geom_point()
  path <- copilot_prerender_ggplot(p)
  list(
    kind             = "ggplot",
    plot             = p,
    prerendered_path = path,
    plot_type        = "pca",
    args             = list(),
    artifact         = NULL,
    label            = label,
    timestamp        = Sys.time()
  )
}

make_plotly_record <- function(label = "volcano plot") {
  p <- plotly::plot_ly(x = 1:3, y = 1:3, type = "scatter", mode = "lines")
  list(
    kind             = "plotly",
    plot             = p,
    prerendered_path = NULL,
    plot_type        = "volcano",
    args             = list(),
    artifact         = NULL,
    label            = label,
    timestamp        = Sys.time()
  )
}

# ===========================================================================
# append_artifact
# ===========================================================================

test_that("append_artifact adds record to history and sets active_artifact", {
  shiny::testServer(CopilotEvidenceServer, {
    rec <- make_ggplot_record()
    session$returned$append_artifact(rec)
    on.exit(copilot_prerender_cleanup(rec$prerendered_path))

    expect_equal(length(shiny::isolate(plot_history())), 1L)
    art <- shiny::isolate(active_artifact())
    expect_equal(art$kind, "ggplot")
    expect_equal(art$label, "pca plot")
  })
})

test_that("append_artifact increments history on each call", {
  shiny::testServer(CopilotEvidenceServer, {
    rec1 <- make_ggplot_record("plot 1")
    rec2 <- make_plotly_record("plot 2")
    on.exit(copilot_prerender_cleanup(rec1$prerendered_path))

    session$returned$append_artifact(rec1)
    session$returned$append_artifact(rec2)

    expect_equal(length(shiny::isolate(plot_history())), 2L)
  })
})

test_that("active_artifact reflects the last appended record", {
  shiny::testServer(CopilotEvidenceServer, {
    rec1 <- make_ggplot_record("first")
    rec2 <- make_plotly_record("second")
    on.exit(copilot_prerender_cleanup(rec1$prerendered_path))

    session$returned$append_artifact(rec1)
    session$returned$append_artifact(rec2)

    art <- shiny::isolate(active_artifact())
    expect_equal(art$kind, "plotly")
    expect_equal(art$label, "second")
  })
})

# ===========================================================================
# Kind detection via active_artifact
# ===========================================================================

test_that("active_artifact$kind is 'ggplot' after appending a ggplot record", {
  shiny::testServer(CopilotEvidenceServer, {
    rec <- make_ggplot_record()
    on.exit(copilot_prerender_cleanup(rec$prerendered_path))
    session$returned$append_artifact(rec)

    art <- shiny::isolate(active_artifact())
    expect_equal(art$kind, "ggplot")
  })
})

test_that("active_artifact$kind is 'plotly' after appending a plotly record", {
  shiny::testServer(CopilotEvidenceServer, {
    rec <- make_plotly_record()
    session$returned$append_artifact(rec)

    art <- shiny::isolate(active_artifact())
    expect_equal(art$kind, "plotly")
  })
})

test_that("active_artifact is NULL before any append and set after", {
  shiny::testServer(CopilotEvidenceServer, {
    expect_null(shiny::isolate(active_artifact()))

    rec <- make_ggplot_record()
    on.exit(copilot_prerender_cleanup(rec$prerendered_path))
    session$returned$append_artifact(rec)

    expect_false(is.null(shiny::isolate(active_artifact())))
  })
})

test_that("has_history flag: FALSE with 0-1 plots, TRUE with 2+", {
  shiny::testServer(CopilotEvidenceServer, {
    expect_equal(length(shiny::isolate(plot_history())), 0L)

    rec1 <- make_ggplot_record()
    on.exit(copilot_prerender_cleanup(rec1$prerendered_path))
    session$returned$append_artifact(rec1)
    expect_equal(length(shiny::isolate(plot_history())), 1L)
    # Only 1 entry — carousel should not show
    expect_false(length(shiny::isolate(plot_history())) > 1L)

    rec2 <- make_plotly_record()
    session$returned$append_artifact(rec2)
    expect_true(length(shiny::isolate(plot_history())) > 1L)
  })
})

# ===========================================================================
# clear / clear_plots
# ===========================================================================

test_that("clear resets history and active_artifact", {
  shiny::testServer(CopilotEvidenceServer, {
    rec <- make_ggplot_record()
    # clear() will handle the PNG cleanup
    session$returned$append_artifact(rec)
    expect_false(is.null(shiny::isolate(active_artifact())))

    session$returned$clear()

    expect_null(shiny::isolate(active_artifact()))
    expect_equal(length(shiny::isolate(plot_history())), 0L)
  })
})

test_that("clear_plots removes prerendered PNG files", {
  shiny::testServer(CopilotEvidenceServer, {
    rec  <- make_ggplot_record()
    path <- rec$prerendered_path
    expect_true(file.exists(path))

    session$returned$append_artifact(rec)
    session$returned$clear_plots()

    expect_false(file.exists(path))
  })
})

# ===========================================================================
# Table management
# ===========================================================================

test_that("update_table sets table_data; clear_table resets it", {
  shiny::testServer(CopilotEvidenceServer, {
    expect_null(shiny::isolate(table_data()))

    session$returned$update_table(data.frame(a = 1:3, b = letters[1:3]))
    expect_false(is.null(shiny::isolate(table_data())))
    expect_equal(nrow(shiny::isolate(table_data())), 3L)

    session$returned$clear_table()
    expect_null(shiny::isolate(table_data()))
  })
})

# ===========================================================================
# renderImage — verify via active_artifact state (testServer encodes as base64)
# ===========================================================================

test_that("active ggplot record has a valid prerendered_path that exists on disk", {
  shiny::testServer(CopilotEvidenceServer, {
    rec <- make_ggplot_record()
    on.exit(copilot_prerender_cleanup(rec$prerendered_path))
    session$returned$append_artifact(rec)

    art  <- shiny::isolate(active_artifact())
    path <- art$prerendered_path
    expect_true(!is.null(path))
    expect_true(file.exists(path))
    expect_gt(file.size(path), 0L)
    expect_match(path, "\\.png$", ignore.case = TRUE)
  })
})

# ===========================================================================
# Carousel navigation
# ===========================================================================

test_that("carousel click (select_plot) updates active_artifact to the chosen entry", {
  shiny::testServer(CopilotEvidenceServer, {
    rec1 <- make_ggplot_record("plot 1")
    rec2 <- make_plotly_record("plot 2")
    on.exit(copilot_prerender_cleanup(rec1$prerendered_path))

    session$returned$append_artifact(rec1)
    session$returned$append_artifact(rec2)

    # Currently active is rec2 (plotly). Navigate back to rec1 via carousel.
    session$setInputs(select_plot = 1L)

    art <- shiny::isolate(active_artifact())
    expect_equal(art$kind, "ggplot")
    expect_equal(art$label, "plot 1")
  })
})

test_that("public active() reactive returns the same value as active_artifact()", {
  shiny::testServer(CopilotEvidenceServer, {
    rec <- make_plotly_record("test")
    session$returned$append_artifact(rec)

    from_api    <- shiny::isolate(session$returned$active())
    from_direct <- shiny::isolate(active_artifact())
    expect_identical(from_api, from_direct)
  })
})
