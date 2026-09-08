##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

.module_dir <- if (dir.exists("components/modules")) {
  "components/modules"
} else {
  "../../components/modules"
}

source(file.path(.module_dir, "AiReports.R"), local = TRUE)

test_that("ai_report_slots discovers valid pgx$ai report slots", {
  pgx <- list(ai = list(
    combined = list(report = "Combined report", prompt = "combined prompt"),
    de = list(report = "DE report"),
    meta = list(llm_model = "test"),
    malformed = list(prompt = "missing report"),
    empty = list(report = "")
  ))

  expect_equal(ai_report_slots(pgx), c("combined", "de"))
  expect_true(ai_report_has(pgx))
  expect_true(ai_report_has(pgx, "combined"))
  expect_true(ai_report_has(pgx, c("combined", "de")))
  expect_false(ai_report_has(pgx, c("combined", "mofa")))
})

test_that("ai_report_has treats drug multi-slots as the drugs module", {
  pgx <- list(ai = list(
    drugs_L1000_activity = list(report = "Drug activity report"),
    combined = list(report = "Combined report")
  ))

  expect_true(ai_report_has(pgx, "drugs"))
  expect_true(ai_report_has(pgx, c("drugs", "combined")))
  expect_false(ai_report_has(pgx, c("drugs", "mofa")))
})

test_that("ai_report_needs_generation only requires one usable report", {
  no_reports <- list(wgcna = list(), gx.meta = data.frame())
  partial_reports <- list(
    wgcna = list(),
    gx.meta = data.frame(),
    ai = list(de = list(report = "DE report"))
  )
  no_reportable_modules <- list()

  expect_true(ai_report_needs_generation(no_reports))
  expect_false(ai_report_needs_generation(partial_reports))
  expect_false(ai_report_needs_generation(no_reportable_modules))
})

test_that("ai_report_modules_for_pgx returns available generation modules", {
  pgx <- list(
    wgcna = list(),
    drugs = list(L1000 = list()),
    gx.meta = data.frame(),
    gset.meta = data.frame()
  )

  expect_equal(
    ai_report_modules_for_pgx(pgx),
    c("wgcna", "drugs", "de", "pathways", "combined")
  )
  expect_equal(ai_report_modules_for_pgx(list()), character(0))
})

test_that("ai_report_get reads dynamic report and prompt slots", {
  pgx <- list(ai = list(
    combined = list(report = "Combined report", prompt = "combined prompt"),
    de = list(report = "DE report"),
    meta = list(llm_model = "test")
  ))

  combined <- ai_report_get(pgx, "combined")
  expect_equal(combined$slot, "combined")
  expect_equal(combined$report, "Combined report")
  expect_equal(combined$prompt, "combined prompt")
  expect_equal(combined$meta$llm_model, "test")

  de <- ai_report_get(pgx, "de")
  expect_equal(de$report, "DE report")
  expect_null(de$prompt)
})

test_that("ai_report_update_text edits pgx$ai reports only", {
  pgx <- list(
    report = list(report = "old summary"),
    wgcna = list(report = list(report = "old wgcna")),
    ai = list(
      combined = list(report = "Combined report", prompt = "combined prompt"),
      wgcna = list(report = "WGCNA report"),
      meta = list(llm_model = "test")
    )
  )

  updated <- ai_report_update_text(pgx, c(
    combined = "Edited summary",
    wgcna = "Edited WGCNA",
    missing = "Ignored report"
  ))

  expect_equal(updated$ai$combined$report, "Edited summary")
  expect_equal(updated$ai$combined$prompt, "combined prompt")
  expect_equal(updated$ai$wgcna$report, "Edited WGCNA")
  expect_null(updated$ai$missing)
  expect_equal(updated$ai$meta$llm_model, "test")
  expect_equal(updated$report$report, "old summary")
  expect_equal(updated$wgcna$report$report, "old wgcna")
})

test_that("ai_report_update_text sets edited=TRUE without touching created_at", {
  created_epoch <- 1750000000

  pgx <- list(ai = list(
    combined = list(
      report     = "original",
      prompt     = "p",
      usage      = list(total = 100L, model = "m"),
      created_at = created_epoch,
      edited     = FALSE,
      edited_at  = ""
    )
  ))

  before <- Sys.time()
  updated <- ai_report_update_text(pgx, c(combined = "edited text"))
  after  <- Sys.time()

  slot <- updated$ai$combined
  expect_equal(slot$report, "edited text")
  expect_true(isTRUE(slot$edited))
  expect_gte(slot$edited_at, as.numeric(before))
  expect_lte(slot$edited_at, as.numeric(after))
  # created_at must be untouched
  expect_equal(slot$created_at, created_epoch)
  # usage must be untouched
  expect_equal(slot$usage$total, 100L)
  expect_equal(slot$prompt, "p")
})

test_that("ai_report_merge_into_reactive preserves existing pgx$ai slots", {
  pgx_rv <- new.env(parent = emptyenv())
  pgx_rv$ai <- list(
    de = list(report = "DE report"),
    combined = list(report = "Combined report"),
    meta = list(llm_model = "old")
  )

  result <- list(ai = list(
    pathways = list(report = "Pathway report"),
    drugs_L1000 = list(report = "Drug report"),
    meta = list(llm_model = "new")
  ))

  expect_true(ai_report_merge_into_reactive(pgx_rv, result))
  expect_equal(pgx_rv$ai$de$report, "DE report")
  expect_equal(pgx_rv$ai$combined$report, "Combined report")
  expect_equal(pgx_rv$ai$pathways$report, "Pathway report")
  expect_equal(pgx_rv$ai$drugs_L1000$report, "Drug report")
  expect_equal(pgx_rv$ai$meta$llm_model, "new")

  expect_false(ai_report_merge_into_reactive(pgx_rv, list(ai = list())))
})

test_that("ai_report_get_module maps static UI modules only", {
  pgx <- list(ai = list(
    combined = list(report = "Combined report"),
    pathways = list(report = "Pathway report"),
    drugs_L1000_activity = list(report = "L1000 report"),
    drugs_GDSC_sensitivity = list(report = "GDSC report")
  ))

  expect_equal(ai_report_get_module(pgx, "summary")$report, "Combined report")
  expect_equal(ai_report_get_module(pgx, "enrichment")$report, "Pathway report")
  expect_null(ai_report_get_module(pgx, "drugs"))
  expect_null(ai_report_get_module(pgx, "cmap"))
})

test_that("drug report helpers expose independent drug slots", {
  pgx <- list(
    drugs = list(
      "L1000_ACTIVITYS_N20D1011" = list(),
      "GDSC/sensitivity" = list()
    ),
    ai = list(
      drugs_L1000_ACTIVITYS_N20D1011 = list(report = "L1000 report"),
      drugs_GDSC_sensitivity = list(report = "GDSC report")
    )
  )

  expect_equal(
    ai_report_drug_slots(pgx),
    c("drugs_L1000_ACTIVITYS_N20D1011", "drugs_GDSC_sensitivity")
  )
  expect_equal(
    ai_report_drug_label(pgx, "drugs_L1000_ACTIVITYS_N20D1011"),
    "L1000 Activity"
  )
  expect_equal(
    ai_report_drug_label(pgx, "drugs_GDSC_sensitivity"),
    "GDSC/sensitivity"
  )
  expect_equal(
    ai_report_get(pgx, "drugs_GDSC_sensitivity")$report,
    "GDSC report"
  )
})

test_that("ai report helpers tolerate invalid and old-shape inputs", {
  old_shape <- list(
    report = list(report = "old summary"),
    wgcna = list(report = list(report = "old wgcna")),
    mofa = list(report = "old mofa")
  )
  malformed <- list(ai = list(
    combined = "not a list",
    de = list(report = NA_character_),
    meta = list()
  ))

  for (pgx in list(NULL, list(), old_shape, malformed)) {
    expect_equal(ai_report_slots(pgx), character(0))
    expect_false(ai_report_has(pgx))
    expect_null(ai_report_get(pgx, "combined"))
  }
})

test_that("ai_report_copy_into_reactive copies only valid pgx$ai shapes", {
  pgx_rv <- new.env(parent = emptyenv())
  ai <- list(
    combined = list(report = "Combined report"),
    meta = list(date = "now")
  )

  expect_true(ai_report_copy_into_reactive(pgx_rv, ai))
  expect_equal(pgx_rv$ai$combined$report, "Combined report")

  pgx_rv$ai <- NULL
  expect_true(ai_report_copy_into_reactive(pgx_rv, list(ai = ai)))
  expect_equal(pgx_rv$ai$combined$report, "Combined report")

  pgx_rv$ai <- ai
  expect_false(ai_report_copy_into_reactive(pgx_rv, list(meta = list())))
  expect_equal(pgx_rv$ai$combined$report, "Combined report")
})

test_that("ai_report_copy_into_reactive does not copy old report slots", {
  pgx_rv <- new.env(parent = emptyenv())
  pgx_rv$wgcna <- list(report = list(report = "old wgcna"))

  updated <- list(
    ai = list(wgcna = list(report = "New WGCNA report")),
    wgcna = list(report = list(report = "generated old-shape report"))
  )

  expect_true(ai_report_copy_into_reactive(pgx_rv, updated))
  expect_equal(pgx_rv$ai$wgcna$report, "New WGCNA report")
  expect_equal(pgx_rv$wgcna$report$report, "old wgcna")
})

test_that("ai_infographic_set stores image bytes and metadata under pgx$ai", {
  img <- tempfile(fileext = ".png")
  writeBin(as.raw(c(137, 80, 78, 71, 13, 10, 26, 10)), img)
  result <- structure(list(
    path = img,
    prompt = "draw report",
    metadata = list(model = "gemini-test", cached = FALSE)
  ), class = "omicsai_image_result")

  pgx <- list(ai = list(combined = list(report = "Combined report")))
  updated <- playbase::ai_infographic_set(pgx, "combined", result,
    style = "bigomics", n_blocks = 2)

  stored <- updated$ai$combined$infographic
  expect_equal(stored$status, "done")
  expect_type(stored$bytes, "raw")
  expect_gt(length(stored$bytes), 0)
  expect_equal(stored$content_type, "image/png")
  expect_equal(stored$prompt, "draw report")
  expect_equal(stored$model, "gemini-test")
  expect_equal(stored$style, "bigomics")
  expect_equal(stored$n_blocks, 2)
  expect_equal(playbase::ai_infographic_slots(updated), "combined")
})

test_that("ai_infographic_get and render tolerate missing and error states", {
  pgx <- list(ai = list(
    combined = list(report = "Combined report"),
    wgcna = list(
      report = "WGCNA report",
      infographic = list(status = "error", error = "provider failed")
    )
  ))

  expect_null(playbase::ai_infographic_get(pgx, "combined"))
  expect_equal(playbase::ai_infographic_get(pgx, "wgcna")$error, "provider failed")
  expect_equal(playbase::ai_infographic_slots(pgx), "wgcna")
})

test_that("ai_infographic_friendly_error hides provider internals", {
  raw_error <- paste(
    "Error in .image_api_call(prompt = full_prompt):",
    "No image data in Gemini response"
  )

  expect_match(playbase::ai_infographic_friendly_error(raw_error),
    "server seems saturated", fixed = TRUE)
  expect_false(grepl(
    "\\.image_api_call|Gemini",
    playbase::ai_infographic_friendly_error(raw_error)
  ))
  expect_match(playbase::ai_infographic_friendly_error("HTTP 429 rate limit"),
    "rate limit", fixed = TRUE)
  expect_match(playbase::ai_infographic_friendly_error(""),
    "Please try again", fixed = TRUE)
})

test_that("ai_infographic_set stores friendly errors", {
  pgx <- list(ai = list(pathways = list(report = "Enrichment report")))
  updated <- playbase::ai_infographic_set(
    pgx,
    "pathways",
    NULL,
    status = "error",
    error = "No image data in Gemini response"
  )

  stored <- updated$ai$pathways$infographic
  expect_equal(stored$status, "error")
  expect_match(stored$error, "server seems saturated", fixed = TRUE)
  expect_false(grepl("Gemini", stored$error))
})

test_that("ai_infographic_render_value writes stored bytes to a temp file", {
  tmpdir <- tempfile()
  img <- list(
    status = "done",
    bytes = as.raw(c(1, 2, 3, 4)),
    content_type = "image/png"
  )

  rendered <- ai_infographic_render_value(img, tmpdir, "combined")
  expect_true(file.exists(rendered$src))
  expect_equal(readBin(rendered$src, "raw", n = 4), img$bytes)
  expect_equal(rendered$contentType, "image/png")
  expect_equal(rendered$width, "100%")
})
