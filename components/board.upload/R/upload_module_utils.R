data_error_modal <- function(path, data_type) {
  shinyalert::shinyalert(
    title = "Upload error",
    text = shiny::tagList(
      shiny::HTML(glue::glue("There seems to be an issue with your data file. It could be the data has been wrongly exported, <b><a href='https://raw.githubusercontent.com/bigomics/playbase/refs/heads/main/inst/extdata/{data_type}.csv' target='_blank'>here</a></b> you can find a correct data file, and <b><a href='https://omicsplayground.readthedocs.io/en/latest/dataprep/{data_type}/' target='_blank'>here</a></b> you can find the file specification documentation. Here are the 10 first lines of your input file, we suggest visualizing it using a plain text editor (rather than Excel) to check on your end the file is correct.")),
      shiny::HTML("<br><br>"),
      shiny::tags$pre(
        paste(readLines(path, n = 10), collapse = "\n")
      )
    ),
    html = TRUE,
    size = "l",
    type = "error",
    className = "left-align"
  )
}
