data_error_modal <- function(path, data_type) {
  con <- file(path, "rb")
  text_raw <- readBin(con, "raw", n = 1000)
  close(con)
  is_bin <- tryCatch(
    {
      paste(rawToChar(text_raw), collapse = "")
    },
    error = function(w) {
      NULL
    }
  )
  if (!is.null(is_bin)) {
    shinyalert::shinyalert(
      title = "Upload error",
      text = shiny::tagList(
        shiny::HTML(glue::glue("There seems to be an issue with your data file. It could be the data has been wrongly exported, <u><b><a href='https://omicsplayground.readthedocs.io/en/latest/dataprep/{data_type}/' target='_blank'>here</a></b></u> you can find the file specification documentation, and <u><b><a href='https://raw.githubusercontent.com/bigomics/playbase/refs/heads/main/inst/extdata/{data_type}.csv' target='_blank'>here</a></b></u> you can find a correct data file. <br><br> The 10 first lines of the input file are the following:")),
        shiny::HTML("<br><br>"),
        shiny::tags$pre(
          paste(readLines(path, n = 10), collapse = "\n")
        ),
        shiny::HTML("We suggest visualizing on your computer the data using a plain text editor (rather than Excel) to check on your end the file is correct.")
      ),
      html = TRUE,
      size = "l",
      type = "error",
      className = "left-align"
    )
  } else {
    shinyalert::shinyalert(
      title = "Upload error",
      text = shiny::tagList(
        shiny::HTML(glue::glue("There seems to be an issue with your data file. We detected it is a binary file rather than a plain text csv file. It could be the data has been wrongly exported. You can find the file specification documentation <u><b><a href='https://omicsplayground.readthedocs.io/en/latest/dataprep/{data_type}/' target='_blank'>here</a></b></u>.<br><br>")),
        shiny::HTML("We suggest visualizing on your computer the data using a plain text editor (rather than Excel) to check on your end the file is correct.")
      ),
      html = TRUE,
      size = "l",
      type = "error",
      className = "left-align"
    )
  }
}
