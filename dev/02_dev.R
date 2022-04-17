# Building a Prod-Ready, Robust Shiny Application.
#
# README: each step of the dev files is optional, and you don't have to
# fill every dev scripts before getting started.
# 01_start.R should be filled at start.
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
#
#
###################################
#### CURRENT FILE: DEV SCRIPT #####
###################################

setwd("~/Playground/omicsplayground3/")  ## set to package root!!
usethis::proj_set()
here::here()

# Hack for subfolders
source("dev/create_headers.R")
create_headers(c('R','shiny'), add.header=TRUE, excl=c("pgx-include.R","global.R"))

## Dependencies ----
## Amend DESCRIPTION with dependencies read from package code parsing
##attachment::att_amend_desc()  ## last error 17.04.2022 (missing or mispelled)

## Update NAMESPACE
roxygen2::roxygenize()               ## writing NAMESPACE
file.show('NAMESPACE')

## Add modules ----
## Create a module infrastructure in R/
#golem::add_module(name = "module1", with_test = TRUE) # Name of the module
#golem::add_module(name = "module2", export=TRUE, with_test = TRUE) # Name of the module

## Add helper functions ----
## Creates fct_* and utils_*
#golem::add_fct("helpers", with_test = TRUE)
#golem::add_utils("helpers", with_test = TRUE)

## External resources
## Creates .js and .css files at inst/app/www
golem::add_js_file("script")
golem::add_js_handler("handlers")
golem::add_css_file("custom")
golem::add_sass_file("custom")

## Add internal datasets ----
## If you have data in your package
usethis::use_data_raw(name = "my_dataset", open = FALSE)

## Tests ----
## Add one line by test you want to create
usethis::use_test("app")

# Documentation

## Vignette ----
usethis::use_vignette("omicsplayground-vignette","OmicsPlayground Vignette")
##usethis::use_vignette("quickstart")

## ERROR!!! seems running out of disk space...
devtools::build_vignettes()          ## builds vignettes, copy to doc
devtools::document()                 ## builds Rd in man folder
devtools::build_manual(path='doc')   ## builds PDF reference manual


