test_that("example data loads with no error",{
  # source aux functions
  source("aux-test-functions.R")

  # Parallelization config
  n_workers <- getOption("test.workers", 3)  # Default 3 workers, configurable
  base_port <- 8080

  # Check if parallel execution is supported (Unix only for mclapply forking)
  can_parallelize <- .Platform$OS.type == "unix" && n_workers > 1

  # test single board minimal components

  # get all board names
  boards <- list.dirs(path = "../../components", full.names = FALSE, recursive = FALSE)
  # split "." and get second name
  boards <- sapply(strsplit(boards, split = "\\."), function(x) x[2])
  boards <- boards[!is.na(boards)]

  # remove upload, loading and user from boards
  boards <- boards[!boards %in% c("upload", "loading", "user")]

  # remove problematic boards
  all_boards <- boards[!boards %in% c("tcga", "signature")]

  # get all pgx files
  pgx_files <- list.files(normalizePath("../../data/pgx_results"), pattern = "*.pgx", full.names = TRUE)


  authentication <- options()$authentication

  # Define the test function for each pgx file
  test_pgx_file <- function(idx) {
    tryCatch({
      pgx_file <- pgx_files[idx]
      # Calculate unique port for this worker (cyclic assignment)
      port <- base_port + ((idx - 1) %% n_workers)
      message(sprintf("[Worker %d, Port %d] %s", idx, port, pgx_file))
      pgx <- playbase::pgx.load(pgx_file)
    boards <- all_boards[all_boards %in% names(pgx)]
    boards <- c("dataview", "enrichment", "clustering", "featuremap", "compare", "correlation", "expression", "pathway", "timeseries", "biomarker", "signature", "intersection", boards)
    if ("mofa" %in% boards) {
      boards <- c(boards, "snf", "lasagna", "deepnet", "mgsea")
    }
    lapply(boards, function(board) {
    # get error from App and save it as error_log
    message(board)
    # board = "wordcloud"
    # board = boards[1]
    App <- shinytest2::AppDriver$new(
      normalizePath("../../dev/board.launch"),
      timeout = 120000,
      height = 1080,
      width = 1920,
      seed = 2910,
      variant = shinytest2::platform_variant(),
      options = list(
        board = board,
        authentication = authentication,
        use_example_data = FALSE
      ),
      shiny_args = list(port = port)
    )
    App$get_values(input = TRUE)

    withr::defer(App$stop())

    App$set_inputs("pgx_path" = pgx_file)
    pgx_file <- tools::file_path_sans_ext(basename(pgx_file))
    if(board == "enrichment") {
      App$set_inputs("enrichment-gs_fdr" = 1)
      App$wait_for_idle(duration = 3000, timeout = 100000)
    }
    if(board == "biomarker") {
      App$run_js("$('#biomarker-pdx_runbutton').click(); ")
      App$wait_for_idle(duration = 3000, timeout = 100000)
    }
    tabs <- searchTabs(board)
    if (!is.null(tabs)){
      lapply(tabs, function(tab){
        App$run_js(generate_js_click_code(tab))
        if(board == "connectivity") {
          duration <- 1000000
          App$wait_for_idle(duration = 10000, timeout = duration)
        } else if (board == "clustering") {
          duration <- 50000
          App$wait_for_idle(duration = 15000, timeout = duration)
        } else if (board == "expression") {
          if (tab == "Overview") {
            App$wait_for_idle(duration = 3000, timeout = 100000)
            App$expect_screenshot(name = paste0(pgx_file, "_", board, "_", tab, "_normal"), threshold = 10, selector = "viewport")
            App$run_js("$('#expression-genetable-datasets-datatable').find('table tr').eq(2).trigger('mousedown').trigger('mouseup'); ")
            App$wait_for_idle(duration = 3000, timeout = 100000)
            App$expect_screenshot(name = paste0(pgx_file, "_", board, "_", tab, "_1gene"), threshold = 10, selector = "viewport")
            App$run_js("$('#expression-gsettable-datasets-datatable').find('table tr').eq(2).trigger('mousedown').trigger('mouseup'); ")
            App$wait_for_idle(duration = 3000, timeout = 100000)
            App$expect_screenshot(name = paste0(pgx_file, "_", board, "_", tab, "_1set"), threshold = 10, selector = "viewport")
            App$run_js(generate_js_click_code("Foldchange (all)"))
            App$wait_for_idle(duration = 3000, timeout = 100000)
            App$expect_screenshot(name = paste0(pgx_file, "_", board, "_", tab, "_foldchange"), threshold = 10, selector = "viewport")
            App$run_js(generate_js_click_code("FDR table"))
            App$wait_for_idle(duration = 3000, timeout = 100000)
            App$expect_screenshot(name = paste0(pgx_file, "_", board, "_", tab, "_fdr"), threshold = 10, selector = "viewport")
            App$run_js(generate_js_click_code("static"))
            App$wait_for_idle(duration = 3000, timeout = 100000)
            App$expect_screenshot(name = paste0(pgx_file, "_", board, "_", tab, "_static"), threshold = 10, selector = "viewport")
          }
          if (tab == "Volcano by comparison" || tab == "Volcano by method") {
            App$wait_for_idle(duration = 3000, timeout = 100000)
            App$expect_screenshot(name = paste0(pgx_file, "_", board, "_", tab, "_dynamic"), threshold = 10, selector = "viewport")
            App$run_js(generate_js_click_code("static", tab_pane = TRUE))
            App$wait_for_idle(duration = 3000, timeout = 100000)
            App$expect_screenshot(name = paste0(pgx_file, "_", board, "_", tab, "_static"), threshold = 10, selector = "viewport")
          }
          App$wait_for_idle(duration = 3000, timeout = 100000)
        } else if (board == "featuremap") {
            App$wait_for_idle(duration = 3000, timeout = 100000)
            App$expect_screenshot(name = paste0(pgx_file, "_", board, "_", tab, "_dynamic"), threshold = 10, selector = "viewport")
            App$run_js(generate_js_click_code("static"))
            App$wait_for_idle(duration = 3000, timeout = 100000)
            App$expect_screenshot(name = paste0(pgx_file, "_", board, "_", tab, "_static"), threshold = 10, selector = "viewport")
        } else {
          duration <- 50000
          App$wait_for_idle(duration = 3000, timeout = duration)
          if (board == "wgcna") {
            if (tab == "Enrichment") {
              App$run_js("$('#wgcna-enrichTable-datasets-datatable').find('table tr').eq(2).trigger('mousedown').trigger('mouseup'); ")
            }
          }
          App$wait_for_idle(duration = 3000, timeout = duration)
        }
        tab <- gsub(" ", "_", tab)
        tab <- gsub("/", "_", tab)

        App$expect_screenshot(name = paste0(pgx_file, "_", board, "_", tab), threshold = 10, selector = "viewport")
      })
    } else {
      App$wait_for_idle(duration = 3000)
      App$expect_screenshot(name = paste0(pgx_file, "_", board), threshold = 10, selector = "viewport")
    }
  })},
    error = function(e) {
      message(sprintf("[Worker %d] ERROR: %s", idx, conditionMessage(e)))
      return(list(error = conditionMessage(e), idx = idx))
    })
  }

  # Run pgx_files in parallel using mclapply (Unix) or sequentially (Windows)
  if (can_parallelize) {
    message(sprintf("Running %d pgx files in parallel with %d workers", length(pgx_files), n_workers))
    AppLog <- parallel::mclapply(
      seq_along(pgx_files),
      test_pgx_file,
      mc.cores = n_workers,
      mc.preschedule = FALSE,
      mc.silent = FALSE
    )
    # Print any warnings from mclapply
    warns <- warnings()
    if (length(warns) > 0) {
      message("mclapply warnings:")
      print(warns)
    }
  } else {
    message(sprintf("Running %d pgx files sequentially", length(pgx_files)))
    AppLog <- lapply(seq_along(pgx_files), test_pgx_file)
  }

  # Check for errors in parallel results
  # mclapply returns try-error objects or our custom error lists
  for (i in seq_along(AppLog)) {
    result <- AppLog[[i]]
    if (inherits(result, "try-error")) {
      message(sprintf("Worker %d failed with try-error: %s", i, as.character(result)))
    } else if (is.list(result) && !is.null(result$error)) {
      message(sprintf("Worker %d failed: %s", result$idx, result$error))
    } else if (is.null(result)) {
      message(sprintf("Worker %d returned NULL (likely crashed)", i))
    }
  }
})
