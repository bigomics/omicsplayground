# Sass code compilation
sass::sass(input = sass::sass_file("inst/app/www/custom.sass"), output = "inst/app/www/custom.css", cache = NULL)

# Sass code compilation
#sass::sass(input = sass::sass_file("inst/app/www/custom.sass"), output = "inst/app/www/custom.css", cache = NULL)

## set to package root!!
setwd("~/R/golem/testpackage/")  

# Set options here
options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode

# Comment this if you don't want the app to be served on a random port
##options(shiny.port = httpuv::randomPort())

# Detach all loaded packages and clean your environment
golem::detach_all_attached()
rm(list=ls(all.names = TRUE))
remove.packages('testpackage')

# Document and reload your package
golem::document_and_reload()

# Run the application
run_app2()
shiny::shinyApp(ui=app_ui, server=app_server)

pkgload::load_all()
ls()
ls("package:testpackage")
run_app2()
shiny::shinyApp(ui=app_ui, server=app_server)


library(golem)
source("R/00Headers.Rx")
ls()
run_app2()
shiny::shinyApp(ui=app_ui, server=app_server)



