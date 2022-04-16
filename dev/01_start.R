# Building a Prod-Ready, Robust Shiny Application.
#
# README: each step of the dev files is optional, and you don't have to
# fill every dev scripts before getting started.
# 01_start.R should be filled at start.
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
#
#
########################################
#### CURRENT FILE: ON START SCRIPT #####
########################################

setwd("~/Playground/omicsplayground/")  ## set to package root!!

## Fill the DESCRIPTION ----
## Add meta data about your application
##
## /!\ Note: if you want to change the name of your app during development,
## either re-run this function, call golem::set_golem_name(), or don't forget
## to change the name in the app_sys() function in app_config.R /!\
##
golem::fill_desc(
  pkg_name = "omicsplayground", # The Name of the package containing the App
  pkg_title = "Visual self-service analytics platform for big omics data", 
  pkg_description = "Omics Playground is a comprehensive self-service analytics platform for the visualization, analytics and exploration of Big Omics Data. It allows biologists to apply a multitude of state-of-the-art analysis tools to their own data to explore and discover underlying biology without coding.",
  author_first_name = "BigOmics Analytics SA", # Your First Name
  author_last_name = "(Switzerland)", # Your Last Name
  author_email = "hello@bigomics.ch", # Your Email
  repo_url = "https://github.com/bigomics/omicsplayground" 
)

## Set {golem} options ----
golem::set_golem_options()

## Create Common Files ----
## See ?usethis for more information
#usethis::use_gpl3_license() # You can set another license here
#usethis::use_readme_rmd(open = FALSE)
# Note that `contact` is required since usethis version 2.1.5
# If your {usethis} version is older, you can remove that param
#usethis::use_code_of_conduct(contact = "Golem User")
#usethis::use_lifecycle_badge("Experimental")
#usethis::use_news_md(open = FALSE)

## Use git ----
usethis::use_git()

## Init Testing Infrastructure ----
## Create a template for tests
golem::use_recommended_tests()

## Favicon ----
# If you want to change the favicon (default is golem's one)
#golem::use_favicon() # path = "path/to/ico". Can be an online file.
golem::remove_favicon() # Uncomment to remove the default favicon

## Add helper functions ----
golem::use_utils_ui(with_test = TRUE)
golem::use_utils_server(with_test = TRUE)

# You're now set! ----

# go to dev/02_dev.R
rstudioapi::navigateToFile("dev/02_dev.R")
