test_that("normal acc: auth-code works",{

  source("aux-test-functions.R")

  App <- shinytest2::AppDriver$new(
    normalizePath("../../dev/board.launch"),
    timeout = 120000,
    height = 1080,
    width = 1920,
    seed = 2910,
    variant = shinytest2::platform_variant(),
    options = list(
      board = "dataview",
      authentication = "login-code",
      use_example_data = FALSE
    ),
    shiny_args = list(port = 8080)
  )

  App$set_inputs(`auth-login_email` = "joe@doe.com")
  App$click("auth-login_submit_btn")
  logs <- App$get_logs()
  login_code_log <- grep("sending login code", logs, value = TRUE)
  login_code <- sub(".*sending login code ([A-Z-]+) .*", "\\1", login_code_log)

  App$run_js('document.dispatchEvent(new KeyboardEvent("keydown", {
                key: "Escape",
                code: "Escape",
                keyCode: 27,
                which: 27,
                bubbles: true
              }));')

  App$set_inputs(`auth-login2_password` = login_code[[1]])

  App$click("auth-login2_submit_btn")

  testthat::expect_equal(App$get_value(output = "current_user"), "joe@doe.com")

  App$stop()
})

test_that("universal acc: auth-code works",{
  source("aux-test-functions.R")

  App <- shinytest2::AppDriver$new(
    normalizePath("../../dev/board.launch"),
    timeout = 120000,
    height = 1080,
    width = 1920,
    seed = 2910,
    variant = shinytest2::platform_variant(),
    options = list(
      board = "dataview",
      authentication = "login-code",
      use_example_data = FALSE
    ),
    shiny_args = list(port = 8080)
  )

  App$set_inputs(`auth-login_email` = "jane@doe.com")
  App$click("auth-login_submit_btn")
  logs <- App$get_logs()
  login_code_log <- grep("sending login code", logs, value = TRUE)
  login_code <- sub(".*sending login code ([A-Z-]+) .*", "\\1", login_code_log)

  App$run_js('document.dispatchEvent(new KeyboardEvent("keydown", {
                key: "Escape",
                code: "Escape",
                keyCode: 27,
                which: 27,
                bubbles: true
              }));')

  App$set_inputs(`auth-login2_password` = login_code[[1]])

  App$click("auth-login2_submit_btn")

  testthat::expect_equal(App$get_value(output = "current_user"), "jane@doe.com")

  App$stop()
})

test_that("universal acc: auth-code universal works",{
  source("aux-test-functions.R")

  App <- shinytest2::AppDriver$new(
    normalizePath("../../dev/board.launch"),
    timeout = 120000,
    height = 1080,
    width = 1920,
    seed = 2910,
    variant = shinytest2::platform_variant(),
    options = list(
      board = "dataview",
      authentication = "login-code",
      use_example_data = FALSE
    ),
    shiny_args = list(port = 8080)
  )

  App$set_inputs(`auth-login_email` = "jane@doe.com")
  App$click("auth-login_submit_btn")

  App$run_js('document.dispatchEvent(new KeyboardEvent("keydown", {
                key: "Escape",
                code: "Escape",
                keyCode: 27,
                which: 27,
                bubbles: true
              }));')

  App$set_inputs(`auth-login2_password` = "123456")

  App$click("auth-login2_submit_btn")

  testthat::expect_equal(App$get_value(output = "current_user"), "jane@doe.com")

  App$stop()
})
