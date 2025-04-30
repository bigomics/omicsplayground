test_that("example data loads with no error",{
  # source aux functions
  source("aux-test-functions.R")
  
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

  AppLog <- lapply(pgx_files, function(pgx_file) {
    message(pgx_file)
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
      shiny_args = list(port = 8080)
    )
    App$get_values(input = TRUE)

    withr::defer(App$stop())

    App$set_inputs("pgx_path" = pgx_file)
    pgx_file <- tools::file_path_sans_ext(basename(pgx_file))
    if(board == "enrichment") {
      App$set_inputs("enrichment-gs_fdr" = 1)
      App$wait_for_idle(duration = 10000, timeout = 50000)
    }
    if(board == "expression") {
      App$run_js("$('#expression-genetable-datasets-datatable').find('table tr').eq(2).trigger('mousedown').trigger('mouseup'); ")
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
        } else {
          duration <- 50000
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
  })})
})
