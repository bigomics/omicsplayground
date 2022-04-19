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
### CURRENT FILE: CDI/CD SCRIPT ###
###################################

setwd("~/Playground/omicsplayground3/")  ## set to package root!!


## Code Coverage----
## Set the code coverage service ("codecov" or "coveralls")
#usethis::use_coverage()

# Create a summary readme for the testthat subdirectory
##remotes::install_github('yonicd/covrpage')
#covrpage::covrpage()

## CI ----
## Use this part of the script if you need to set up a CI
## service for your application
##
## (You'll need GitHub there)
#usethis::use_github()

# GitHub Actions
#usethis::use_github_action()
# Chose one of the three
# See https://usethis.r-lib.org/reference/use_github_action.html
#usethis::use_github_action_check_release()
#usethis::use_github_action_check_standard()
#usethis::use_github_action_check_full()
# Add action for PR
#usethis::use_github_action_pr_commands()

# Travis CI
#usethis::use_travis()
#usethis::use_travis_badge()

# AppVeyor
#usethis::use_appveyor()
#usethis::use_appveyor_badge()

# Circle CI
#usethis::use_circleci()
#usethis::use_circleci_badge()

# Jenkins
#usethis::use_jenkins()

# GitLab CI
#usethis::use_gitlab_ci()

# You're now set! ----
# go to dev/03_deploy.R
#rstudioapi::navigateToFile("dev/03_deploy.R")
