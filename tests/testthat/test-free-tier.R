## test-free-tier.R
##
## Free accounts expire FREE_TRIAL_DAYS after first login (dated by the user
## folder's files); the upgrade link depends on whether the email is academic.

.mod_dir <- if (dir.exists("components/modules")) "components/modules" else "../../components/modules"
source(file.path(.mod_dir, "AuthenticationModule_functions.R"), encoding = "UTF-8", local = TRUE)
opt <- list(FREE_TRIAL_DAYS = 14, ACADEMIC_DOMAINS = "")
get_cro_emails <- function() "lab@cro.com"
SHARE.DIR <- tempfile()

test_that("free account expiry is dated from example-data.pgx, not received datasets", {
  user_dir <- file.path(tempfile(), "joe@doe.com")
  dir.create(user_dir, recursive = TRUE)
  pgx <- file.path(user_dir, "example-data.pgx")
  file.create(pgx)

  Sys.setFileTime(pgx, Sys.time() - 13 * 86400)
  expect_false(as.Date(free_account_expiry(user_dir)) < Sys.Date())
  received <- file.path(user_dir, "shared-by-cro.pgx")
  file.create(received)
  Sys.setFileTime(received, Sys.time() - 60 * 86400)
  expect_false(as.Date(free_account_expiry(user_dir)) < Sys.Date())

  Sys.setFileTime(pgx, Sys.time() - 15 * 86400)
  expect_true(as.Date(free_account_expiry(user_dir)) < Sys.Date())
  file.rename(pgx, paste0(pgx, "_")) ## deleted example still dates the account
  expect_true(as.Date(free_account_expiry(user_dir)) < Sys.Date())
})

test_that("a free account holding a CRO dataset never expires", {
  user_dir <- file.path(tempfile(), "joe@doe.com")
  dir.create(user_dir, recursive = TRUE)
  file.create(file.path(user_dir, "example-data.pgx"))
  Sys.setFileTime(file.path(user_dir, "example-data.pgx"), Sys.time() - 60 * 86400)
  expect_true(as.Date(free_account_expiry(user_dir)) < Sys.Date())

  ## waiting in the received-shares queue
  dir.create(SHARE.DIR)
  queued <- file.path(SHARE.DIR, "study.pgx__to__joe@doe.com__from__lab@cro.com__")
  file.create(queued)
  expect_equal(free_account_expiry(user_dir), "unlimited")
  file.remove(queued)

  ## accepted into the library
  file.create(file.path(user_dir, "study.pgx"))
  write.csv(data.frame(dataset = c("example-data", "study"), creator = c("", "lab@cro.com")),
    file.path(user_dir, "datasets-info.csv"), row.names = FALSE)
  expect_equal(free_account_expiry(user_dir), "unlimited")
})

test_that("academic emails are recognised", {
  expect_true(is_academic_email("a@mit.edu"))
  expect_true(is_academic_email("a@ox.ac.uk"))
  expect_true(is_academic_email("a@unsw.edu.au"))
  expect_false(is_academic_email("a@pfizer.com"))
  expect_false(is_academic_email("john.edu@pfizer.com"))
  expect_false(is_academic_email("a@ethz.ch"))
  expect_true(is_academic_email("a@ethz.ch", domains = "epfl.ch|ethz.ch"))
  expect_false(is_academic_email(NA))
})
