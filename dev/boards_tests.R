# test single board minimal components

### board specific files ###

driver <- shinytest::ShinyDriver$new(
  path = "components/board.tcga/dev_MMM/"
  )

driver$getUrl()
driver$getDebugLog()
driver$listWidgets()$input
driver$listWidgets()$output

# changing pgx file

pgx_file <- normalizePath("data/example-data.pgx")

driver$getValue("pgx_path")

driver$setValue("pgx_path", pgx_file)

driver$getValue("pgx_path")

driver$finalize()
