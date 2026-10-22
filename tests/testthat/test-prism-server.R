## test-prism-server.R

.board_dir <- if (dir.exists("components/app_prism/R")) {
  "components/app_prism/R"
} else {
  "../../components/app_prism/R"
}

source(file.path(.board_dir, "prism_server.R"), local = TRUE)

test_that("df_to_csv round-trips a data.frame including row names", {
  df <- data.frame(x = 1:2, y = c("a", "b"), row.names = c("r1", "r2"))
  csv <- df_to_csv(df)

  expect_type(csv, "character")
  expect_length(csv, 1)

  back <- read.csv(text = csv, row.names = 1)
  expect_equal(back$x, 1:2)
  expect_equal(as.character(back$y), c("a", "b"))
  expect_equal(rownames(back), c("r1", "r2"))
})

test_that("strip_code_fences removes r/R/plain fence markers", {
  expect_equal(strip_code_fences("```r\nggplot(mtcars)\n```"), "\nggplot(mtcars)\n")
  expect_equal(strip_code_fences("```R\nggplot(mtcars)\n```"), "\nggplot(mtcars)\n")
  expect_equal(strip_code_fences("```\nggplot(mtcars)\n```"), "\nggplot(mtcars)\n")
  expect_equal(strip_code_fences("no fences here"), "no fences here")
})

test_that("build_prism_prompt embeds the request, data context and preloaded packages", {
  p <- build_prism_prompt(
    msg = "plot mpg vs wt",
    vars = "mpg, wt",
    rows = "Mazda RX4, Datsun 710",
    last_plotcode = "",
    dataset = "mtcars",
    pointsize = 3,
    fontsize = 12,
    theme = "classic",
    packages = PRISM_WEBR_PACKAGES
  )

  expect_type(p, "character")
  expect_match(p, "plot mpg vs wt", fixed = TRUE)
  expect_match(p, "mpg, wt", fixed = TRUE)
  expect_match(p, "Mazda RX4, Datsun 710", fixed = TRUE)
  expect_match(p, PRISM_WEBR_PACKAGES, fixed = TRUE)
  expect_match(p, "never call library\\(\\)")
})

test_that("update_plotcode returns NULL when there is no last plotcode", {
  expect_null(update_plotcode("", 3, 12, "classic", "mtcars"))
})

test_that("update_plotcode patches theme, point size, font size and title onto the last code", {
  last <- "ggplot(data, aes(x=wt,y=mpg)) + geom_point()"
  out <- update_plotcode(last, 5, 20, "dark", "mtcars")

  expect_match(out, "theme_dark()", fixed = TRUE)
  expect_match(out, "geom_point(size=5)", fixed = TRUE)
  expect_match(out, "element_text(size=20)", fixed = TRUE)
  expect_match(out, "labs(title='mtcars')", fixed = TRUE)
})
