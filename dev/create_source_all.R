
## Create global "source_all" file
source("dev-utils.R")
appdir <- get_appdir() 
setwd(appdir)
appdir

setwd(appdir)
create_SourceAll('components',add.comments=FALSE)

if(0) {
    ##create_allcode('R') 
    source("components/00SourceAll.R")
    source("components/00SourceAll.R",chdir=TRUE)
}
