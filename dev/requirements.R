# This file is supposed to run from the root Playground folder
if (basename(getwd()) != "omicsplayground") {
    stop("Please run from the OmicsPlayground root folder")
}

require <- function(pkg) (pkg %in% installed.packages()[,'Package'])

if(!require("renv")) install.packages("renv")
if(!require("BiocManager")) install.packages("BiocManager")
if(!require("remotes")) install.packages("remotes")
##if(!require("devtools")) install.packages("devtools")
if(!require("reticulate")) install.packages("reticulate")

## ---------------------------------------------------------------------
## Automatically scan all used packages and install
## ---------------------------------------------------------------------
## We use renv to detect dependencies. Renv is looking for library and
## require statements in the r/R source files. The script also sets:
## - pgx.imports
## - pgx.remotes
source("dev/functions.R")
source("dev/write_description.R")
pkg <- scan_description(path=NULL) 

## some are perhaps already installed by playbase or other packages
installed.pkgs <- sort(as.character(installed.packages()[,"Package"]))
pkgs <- c(pkg$imports, names(pkg$remotes))
message("Total packages (not in playbase): ", length(pkgs))
message("Packages not yet installed: ", sum(!pkgs %in% installed.pkgs))

## determine missing packagers
P <- available.packages()
missing.imports <- setdiff( pkg$imports, installed.pkgs )
missing.remotes <- pkg$remotes[!names(pkg$remotes) %in% installed.pkgs]

## suppress ultra-verbose packages that link to RccpEigen
pkg.eigen <- rownames(P)[grep("RcppEigen",as.character(P[,"LinkingTo"]))]
pkg.eigen <- c("RccpEigen",pkg.eigen)
pkg.eigen <- intersect(missing.imports,pkg.eigen)
if(length(pkg.eigen) > 0) {
  message("> Pre-installing ultra-verbose packages: ", paste(pkg.eigen,collapse=" "))
  if(!dir.exists("~/.R")) dir.create("~/.R")
  if(!file.exists("~/.R/Makevars")) file.copy("dev/Makevars","~/.R/Makevars.save")
  file.copy("dev/Makevars.no-error","~/.R/Makevars")
  BiocManager::install(pkg.eigen,ask=FALSE)
  file.remove("~/.R/Makevars")
  if(!file.exists("~/.R/Makevars.save")) file.copy("dev/Makevars.save","~/.R/Makevars")
  missing.imports <- setdiff(missing.imports, pkg.eigen)
} else {
  print("> No ultra-verbose packages!")
}

print(">>> Installing missing CRAN/BioConductor packages...")
message( ">>> Installing ",length(missing.imports), " missing ",
        "CRAN/BioConductor packages:\n", paste(missing.imports,collapse="\n"))
if(length(missing.imports>0)) {
  BiocManager::install( missing.imports, ask=FALSE, force=FALSE )
}

message(">>> Installing ",length(missing.remotes)," missing remote packages")
if(length(missing.remotes>0)) {
  missing.remotes <- sub("^url::","",missing.remotes)
  for(url in missing.remotes) {
    message("Installing remote: ", url)
    remotes::install_url( url, ask=FALSE, force=FALSE )
  }
}


## ---------------------------------------------------------------------
## remove unneccessary big packages??
## ---------------------------------------------------------------------
print(">>> removing not needed packages...")
BIG.NOTUSED <- c(
    "reactome.db", ## >2GB!!!
    "terra",
    "RNAseqData.HNRNPC.bam.chr14",
    "tximportData"
    ## "EnsDb.Hsapiens.v86",
    ## "EnsDb.Mmusculus.v79",
    ## "TxDb.Hsapiens.UCSC.hg19.knownGene",  ## need for import
    ## "TxDb.Mmusculus.UCSC.mm10.knownGene"  ## need for import
)
remove.pkgs <- function(pkgs) {
  for(i in 1:length(pkgs)) {
    if(require(pkgs[i])) try(remove.packages(pkgs[i]))
  }
}
BIG.NOTUSED <- intersect(BIG.NOTUSED, installed.packages())
if(length(BIG.NOTUSED)>0) {
    remove.pkgs(BIG.NOTUSED)
}


message("**** FINISHED INSTALLING REQUIREMENTS *****")
