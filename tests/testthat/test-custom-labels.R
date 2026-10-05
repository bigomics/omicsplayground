## test-custom-labels.R
##
## Regression guard: on datasets keyed by non-symbol IDs (UniProt, Ensembl),
## typing a gene symbol into the editor's "Custom labels" box must resolve to
## its feature instead of silently falling back to the default labels.

.ui_dir <- if (dir.exists("components/ui")) "components/ui" else "../../components/ui"
source(file.path(.ui_dir, "ui-EditorHelpers.R"), encoding = "UTF-8", local = TRUE)
source(file.path(.ui_dir, "ui-ColorDefaults.R"), encoding = "UTF-8", local = TRUE)

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

test_that("Volcano (methods) draws a custom label typed as a gene symbol", {
  source(file.path(.ui_dir, "../board.expression/R/expression_plot_volcanoMethods.R"), local = TRUE)
  PlotModuleServer <- function(...) NULL # ponytail: only base.plots() is under test

  features <- c("P05067", "P02649", "P10636")
  fc <- matrix(c(2, -1, 0.5, 1.8, -1.2, 0.4), 3, dimnames = list(features, c("ttest", "limma")))
  mx <- data.frame(row.names = features)
  mx$fc <- fc
  mx$q <- fc * 0 + 0.01
  pgx <- list(
    X = matrix(0, 3, 2, dimnames = list(features, NULL)),
    gx.meta = list(meta = list(A_vs_B = mx)),
    genes = data.frame(
      feature = features, symbol = c("APP", "APOE", "MAPT"),
      gene_title = c("Amyloid precursor protein", "Apolipoprotein E", "Tau"),
      row.names = features
    )
  )

  shiny::testServer(expression_plot_volcanoMethods_server,
    args = list(
      pgx = pgx, comp = shiny::reactive("A_vs_B"), fdr = shiny::reactive(0.05),
      lfc = shiny::reactive(1), show_pv = shiny::reactive(FALSE),
      genes_selected = shiny::reactive(list(sel.genes = features, lab.genes = "P02649")),
      pval_cap = shiny::reactive(1e-12)
    ),
    {
      session$setInputs(custom_labels = TRUE, label_features = "APP", color_selection = FALSE)
      expect_true("APP" %in% base.plots()$data$label)
    }
  )
})
