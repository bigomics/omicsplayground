## Extra packages needed

install.packages("ggfittext")
pak::pak("alastairrushworth/inspectdf")

BiocManager::install("pdftools")
BiocManager::install("waiter")
BiocManager::install("sever")
BiocManager::install("shiny.i18n")

## create python venv and preinstall plotly+kaleido
BiocManager::install("reticulate")
if(!reticulate::virtualenv_exists()) {
  reticulate::virtualenv_create()
}
reticulate::use_virtualenv()
reticulate::py_install('plotly')
reticulate::py_install('kaleido')

