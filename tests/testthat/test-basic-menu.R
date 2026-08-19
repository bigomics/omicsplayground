## test-basic-menu.R

.app_dir <- if (dir.exists("components/app/R")) "components/app/R" else "../../components/app/R"
source(file.path(.app_dir, "opg_ui.R"), local = TRUE)

BASIC_MENU_DEFAULT <- c("dataview", "clustersamples", "diffexpr")

.tree <- list(
  "DataView" = c(dataview = "DataView"),
  "Clustering" = c(clustersamples = "Samples", clusterfeatures = "Features"),
  "Expression" = c(diffexpr = "Differential expression", corr = "Correlation analysis")
)

test_that("basic menu keeps the admin selection, in full-menu order", {
  m <- opg_basic_menu_tree(.tree, c("corr", "dataview"))

  expect_equal(names(m), c("dataview", "corr"))
  ## key must equal the tab name, else createMenu renders a nested group
  expect_equal(names(m[["corr"]]), "corr")
  expect_equal(unname(m[["corr"]]), "Correlation analysis")
})

test_that("flattened clustering entries get a standalone title", {
  m <- opg_basic_menu_tree(.tree, "clustersamples")
  expect_equal(unname(m[["clustersamples"]]), "Cluster Samples")
})

test_that("unknown boards are dropped and an empty result falls back to the default", {
  expect_equal(names(opg_basic_menu_tree(.tree, c("dataview", "nosuchboard"))), "dataview")
  expect_equal(names(opg_basic_menu_tree(.tree, "nosuchboard")), BASIC_MENU_DEFAULT)
  expect_equal(names(opg_basic_menu_tree(.tree, character(0))), BASIC_MENU_DEFAULT)
})

## --- advanced_option(): which boards get their settings greyed out ---------

.ui_dir <- if (dir.exists("components/ui")) "components/ui" else "../../components/ui"
source(file.path(.ui_dir, "ui-utils.R"), local = TRUE)

test_that("lock_advanced greys every accordion of a locked board's settings", {
  opt <<- list(BASIC_LOCKED = c("dataview"))
  on.exit(rm(opt, envir = globalenv()), add = TRUE)

  tab <- function(name) {
    shiny::div(
      class = "big-tab", `data-name` = name,
      shiny::div(
        class = "tab-settings",
        shiny::div(class = "accordion", id = "opts"),
        shiny::div(class = "accordion keep-live", id = "primary"),
        shiny::selectInput("gene", "Gene:", "A") ## a primary control, untouched
      )
    )
  }
  marked <- function(x) {
    html <- as.character(x)
    lengths(regmatches(html, gregexpr("advanced-option", html)))
  }

  expect_equal(marked(lock_advanced(tab("dataview-tab"))), 1) ## keep-live opted out
  expect_equal(marked(lock_advanced(tab("corr-tab"))), 0)     ## board not locked

  html <- paste(as.character(lock_advanced(tab("dataview-tab"))), collapse = "")
  expect_true(grepl('class="accordion advanced-option"', html, fixed = TRUE))
  expect_true(grepl('class="accordion keep-live"', html, fixed = TRUE))
})

test_that("real bslib accordions: keep-live and action buttons stay usable", {
  ## built with bslib, not hand-rolled divs: bslib::accordion() emits several
  ## `class` attributes, which is what broke the keep-live check once already
  opt <<- list(BASIC_LOCKED = "wgcna")
  on.exit(rm(opt, envir = globalenv()), add = TRUE)

  settings <- function(...) {
    shiny::div(
      class = "big-tab", `data-name` = "wgcna-tab",
      shiny::div(class = "d-none tab-settings", ...)
    )
  }
  marked <- function(x) {
    html <- paste(as.character(x), collapse = "")
    lengths(regmatches(html, gregexpr("advanced-option", html)))
  }

  plain <- bslib::accordion(id = "a", bslib::accordion_panel("Options", "fdr"))
  expect_equal(marked(lock_advanced(settings(plain))), 1)

  ## opts out explicitly
  kept <- bslib::accordion(id = "b", class = "keep-live", bslib::accordion_panel("Options", "x"))
  expect_equal(marked(lock_advanced(settings(kept))), 0)

  ## holds a control, not just settings (WGCNA's "Recompute!" button)
  compute <- bslib::accordion(
    id = "c",
    bslib::accordion_panel(
      "Recompute",
      shiny::selectInput("ngenes", "Number genes:", c(1000, 2000)),
      shiny::actionButton("compute", "Recompute!")
    )
  )
  expect_equal(marked(lock_advanced(settings(compute))), 0)
})

test_that("an empty basic menu is survivable, not a crash", {
  ## createMenu() does `for (i in 1:length(tree))`, so an empty tree used to
  ## take down the whole UI build
  tree <- list("GeneSets" = c(enrich = "Geneset Enrichment"))
  expect_length(opg_basic_menu_tree(tree, "nosuchboard"), 0)
})

test_that("lock_advanced passes through anything that is not a tab tag", {
  opt <<- list(BASIC_LOCKED = "dataview")
  on.exit(rm(opt, envir = globalenv()), add = TRUE)
  expect_equal(lock_advanced("not a tag"), "not a tag")
})

test_that("advanced_option marks a block only for boards the admin locked", {
  opt <<- list(BASIC_LOCKED = c("dataview", "diffexpr"))
  on.exit(rm(opt, envir = globalenv()), add = TRUE)

  expect_equal(advanced_option("dataview"), "advanced-option")
  expect_null(advanced_option("enrich"))

  ## unticking every board leaves all settings usable
  opt <<- list(BASIC_LOCKED = character(0))
  expect_null(advanced_option("dataview"))
})

## --- FORCE_BASIC: pinning one user to the basic menu -----------------------
## read_user_options() drops anything not in its ALLOWED_USER_OPTS whitelist,
## so a missing entry there silently ignores the setting.

.mod_dir <- if (dir.exists("components/modules")) "components/modules" else "../../components/modules"
source(file.path(.mod_dir, "AuthenticationModule_functions.R"), local = TRUE)
dbg <- function(...) invisible(NULL) ## app global, not sourced here

test_that("a per-user OPTIONS file can pin the user to the basic menu", {
  opt <<- list(FORCE_BASIC = FALSE, WATERMARK = TRUE)
  on.exit(rm(opt, envir = globalenv()), add = TRUE)

  user_dir <- withr::local_tempdir()
  writeLines(c("FORCE_BASIC = TRUE", "MAX_SAMPLES = 10"), file.path(user_dir, "OPTIONS"))

  user_opt <- read_user_options(user_dir)
  expect_true(isTRUE(as.logical(user_opt$FORCE_BASIC)))
  expect_equal(user_opt$MAX_SAMPLES, 10)
  ## untouched globals still come through
  expect_true(user_opt$WATERMARK)
})

test_that("the deployment-wide default applies when the user sets nothing", {
  opt <<- list(FORCE_BASIC = TRUE)
  on.exit(rm(opt, envir = globalenv()), add = TRUE)

  user_dir <- withr::local_tempdir() ## no OPTIONS file
  expect_true(isTRUE(as.logical(read_user_options(user_dir)$FORCE_BASIC)))
})

## --- per-deploy file names in the shared etc/ mount -------------------------
## etc/ is an S3 bucket shared by every deploy, so two servers must never write
## the same file.

.utils_dir <- if (dir.exists("components/app/R/utils")) "components/app/R/utils" else "../../components/app/R/utils"
source(file.path(.utils_dir, "utils.R"), local = TRUE)

test_that("each deploy writes its own file, named from opt$HOSTNAME", {
  ETC <<- "/etc-mount"
  opt <<- list(HOSTNAME = "opg-prod-01")
  on.exit(rm(ETC, opt, envir = globalenv()), add = TRUE)

  expect_equal(etc_host_file("BASIC_MENU"), "/etc-mount/BASIC_MENU-opg-prod-01")

  ## a hostname with path-unsafe characters cannot escape the folder
  opt$HOSTNAME <<- "opg/../prod 02"
  expect_equal(etc_host_file("BASIC_MENU"), "/etc-mount/BASIC_MENU-opg_.._prod_02")
})

test_that("the deploy's own file wins over the shared default", {
  ETC <<- withr::local_tempdir()
  opt <<- list(HOSTNAME = "SERVER-A")
  on.exit(rm(ETC, opt, envir = globalenv()), add = TRUE)

  writeLines(c("dataview", "corr"), file.path(ETC, "BASIC_MENU"))          ## fleet default
  expect_equal(read_board_list("BASIC_MENU", "fallback"), c("dataview", "corr"))

  writeLines("diffexpr", file.path(ETC, "BASIC_MENU-SERVER-A"))            ## this deploy
  expect_equal(read_board_list("BASIC_MENU", "fallback"), "diffexpr")

  ## another server on the same bucket is unaffected by SERVER-A's choice
  opt$HOSTNAME <<- "SERVER-B"
  expect_equal(read_board_list("BASIC_MENU", "fallback"), c("dataview", "corr"))
})

test_that("no file at all falls back, an empty file means none", {
  ETC <<- withr::local_tempdir()
  opt <<- list(HOSTNAME = "SERVER-A")
  on.exit(rm(ETC, opt, envir = globalenv()), add = TRUE)

  expect_equal(read_board_list("BASIC_LOCKED", c("a", "b")), c("a", "b"))

  file.create(file.path(ETC, "BASIC_LOCKED-SERVER-A"))
  expect_equal(read_board_list("BASIC_LOCKED", c("a", "b")), character(0))
})
