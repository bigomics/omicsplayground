# R options for docker building
cat("sourcing .Rprofile (docker)\n")

options(Ncpus = 8L)
options(timeout = 99999)  ## download time.out

## Speed up installation using binary packages from RStudio. 
options(repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/__linux__/noble/latest"))

