## test-copilot-reports-server.R

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

source(file.path(.board_dir, "copilot_logger.R"), local = TRUE)
source(file.path(.board_dir, "copilot_pgx_normalize.R"), local = TRUE)
source(file.path(.module_dir, "AiReports.R"), local = TRUE)
source(file.path(.board_dir, "CopilotReportsServer.R"), local = TRUE)

.reports_pgx <- function() {
  structure(
    list(
      name = "demo",
      date = "2026-06-11",
      X = matrix(1, nrow = 2, ncol = 2),
      drugs = list("L1000_ACTIVITYS_N20D1011" = list()),
      ai = list(
        combined = list(report = "Summary"),
        de = list(report = "DE"),
        drugs_L1000_ACTIVITYS_N20D1011 = list(report = "Drug")
      )
    ),
    class = "pgx"
  )
}

test_that("report labels use helper labels including drugs", {
  pgx <- .reports_pgx()
  expect_identical(copilot_report_label(pgx, "combined"), "Summary")
  expect_identical(copilot_report_label(pgx, "de"), "Differential Expression")
  expect_identical(
    copilot_report_label(pgx, "drugs_L1000_ACTIVITYS_N20D1011"),
    "L1000 Activity"
  )
})

test_that("reports server defaults to combined and removes consumed slots", {
  pgx_rv <- shiny::reactiveVal(.reports_pgx())

  shiny::testServer(CopilotReportsServer, args = list(pgx = pgx_rv), {
    expect_identical(selected_reports(), "combined")

    mark_consumed("combined")
    expect_false("combined" %in% selected_reports())
  })
})

test_that("reports server resets consumed reports when report text changes", {
  pgx_rv <- shiny::reactiveVal(.reports_pgx())

  shiny::testServer(CopilotReportsServer, args = list(pgx = pgx_rv), {
    mark_consumed("combined")
    expect_false("combined" %in% selected_reports())

    updated <- .reports_pgx()
    updated$ai$combined$report <- "Updated summary"
    pgx_rv(updated)
    session$flushReact()

    expect_identical(selected_reports(), "combined")
  })
})

test_that("report pgx helper reads populated reactiveValues", {
  rv <- shiny::reactiveValues()
  pgx <- .reports_pgx()
  for (nm in names(pgx)) {
    rv[[nm]] <- pgx[[nm]]
  }

  out <- shiny::isolate(.copilot_report_pgx(rv))

  expect_identical(ai_report_slots(out), ai_report_slots(pgx))
})

test_that("reports server exposes tools_enabled toggle", {
  pgx_rv <- shiny::reactiveVal(.reports_pgx())

  shiny::testServer(CopilotReportsServer, args = list(pgx = pgx_rv), {
    expect_true(tools_enabled())
    session$setInputs(tools_enabled = FALSE)
    expect_false(tools_enabled())
    session$setInputs(tools_enabled = TRUE)
    expect_true(tools_enabled())
  })
})
