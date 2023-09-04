
options(Ncpus = 8L)
options(timeout = 99999)  ## download time.out

options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])))
source("https://docs.rstudio.com/rspm/admin/check-user-agent.R")
options(repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/__linux__/jammy/latest"))
update.packages(ask=FALSE)
