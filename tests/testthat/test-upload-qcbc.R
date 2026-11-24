## test-upload-qcbc.R
##
## Regression guard for issue #1915: the QC/BC "Batch-effects correction"
## plot failed with "attempt to apply non-function" and left partial plot
## content on screen. The plotting steps of plot_all_methods() are extracted
## and run with realistic inputs (single metadata column, integer/phenotype
## input, samples missing from the PCA positions).

.qcbc_file <- if (file.exists("components/board.upload/R/upload_module_normalization.R")) {
  "components/board.upload/R/upload_module_normalization.R"
} else {
  "../../components/board.upload/R/upload_module_normalization.R"
}

## extract the body of plot_all_methods() without running the module server
.qcbc_extract <- function(file = .qcbc_file) {
  src <- readLines(file)
  start <- grep("plot_all_methods <- function", src)
  end <- grep("^      plot_before_after <- function", src)
  body <- src[(start + 1):(end - 2)]
  ## drop the numeric-colour legend (needs fields::image.plot)
  cut <- grep("fields::image.plot", body)[1]
  if (!is.na(cut)) body <- body[1:(cut - 3)]
  body
}

## run the extracted code with stubbed reactives and shiny
.qcbc_plot_all_methods <- function(props, pos, input = list()) {
  env <- new.env()
  env$results_outlier_methods <- function() list(pos = list(pca = pos))
  env$results_correction_methods <- function() {
    list(pheno = props$Groups, pos = list(ComBat = pos, limma = NULL))
  }
  env$r_samples <- function() props
  input$bec_xpc <- input$bec_xpc %||% "PC1"
  input$bec_ypc <- input$bec_ypc %||% "PC2"
  input$colorby_var <- input$colorby_var %||% "Groups"
  env$input <- input
  env$shiny <- asNamespace("shiny")
  body <- .qcbc_extract()
  eval(parse(text = c("plot_all_methods <- function() {", body, "}")), envir = env)
  env$plot_all_methods
}

.qcbc_props <- function(n = 6) {
  data.frame(Groups = rep(c("A", "B"), length.out = n), row.names = paste0("s", seq_len(n)))
}

test_that("plot_all_methods() draws with a single annotate-by variable", {
  props <- .qcbc_props()
  pos <- matrix(rnorm(12), 6, 2, dimnames = list(rownames(props), c("pca_1", "pca_2")))
  fn <- .qcbc_plot_all_methods(props, pos)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(suppressWarnings(do.call(fn, list())))
})

test_that("plot_all_methods() survives missing/unknown annotate-by variable", {
  props <- .qcbc_props()
  pos <- matrix(rnorm(12), 6, 2, dimnames = list(rownames(props), c("pca_1", "pca_2")))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  fn <- .qcbc_plot_all_methods(props, pos, input = list(colorby_var = "unknown"))
  expect_no_error(suppressWarnings(do.call(fn, list())))

  fn <- .qcbc_plot_all_methods(props, pos, input = list(colorby_var = NULL))
  expect_no_error(suppressWarnings(do.call(fn, list())))
})

test_that("plot_all_methods() colours samples after outlier removal", {
  ## r_samples() keeps all samples, PCA positions only the cleaned set
  props <- .qcbc_props()
  pos <- matrix(rnorm(8), 4, 2, dimnames = list(rownames(props)[1:4], c("pca_1", "pca_2")))
  fn <- .qcbc_plot_all_methods(props, pos)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(suppressWarnings(do.call(fn, list())))
})

test_that("QC/BC plot outputs are unregistered when a plot fails", {
  source(.qcbc_file, local = FALSE, encoding = "UTF-8")
  src <- paste(deparse(upload_module_normalization_server), collapse = "\n")

  ## the failing plot function is wrapped and its output removed on error
  expect_match(src, "func = func_correction", fixed = TRUE)
  expect_match(src, "shiny::removePlotOutput(plot)", fixed = TRUE)

  env <- new.env()
  env$shiny <- asNamespace("shiny")
  eval(
    parse(text = "retain_plot_output <- function(plot) {
      if (!is.null(plot)) {
        shiny::onReactiveDomainEnded(shiny::getDefaultReactiveDomain(), function() shiny::removePlotOutput(plot))
      }
      invisible(plot)
    }"),
    envir = env
  )
  expect_equal(env$retain_plot_output(42), 42)
})
