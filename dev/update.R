
options(Ncpus = 8L)
options(timeout = 99999)  ## download time.out
options(repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/__linux__/jammy/latest"))
update.packages(ask=FALSE)
