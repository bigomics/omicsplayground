## test-custom-labels.R
##
## Regression guard: on datasets keyed by non-symbol IDs (UniProt, Ensembl),
## typing a gene symbol into the editor's "Custom labels" box must resolve to
## its feature instead of silently falling back to the default labels.

.ui_dir <- if (dir.exists("components/ui")) "components/ui" else "../../components/ui"
source(file.path(.ui_dir, "ui-EditorHelpers.R"), encoding = "UTF-8", local = TRUE)

test_that("custom labels resolve gene symbols to non-symbol feature IDs", {
  features <- c("P05067", "P02649", "P10636", "Q8WWW0")
  symbols <- c("APP", "APOE", "MAPT", "RASSF5")
  input <- list(custom_labels = TRUE, label_features = "APP\nMAPT")

  expect_equal(get_custom_labels(input, features, defaults = "top"), "top")
  expect_setequal(
    get_custom_labels(input, features, defaults = "top", alt_names = list(features, symbols)),
    c("P05067", "P10636")
  )
})
