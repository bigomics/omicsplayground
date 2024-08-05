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
source("dev/write_description.R")

## some are perhaps already installed by playbase or other packages
P <- installed.packages()
installed.pkg <- sort(as.character(P[,"Package"]))
missing.imports <- setdiff( pkg.imports, installed.pkg )
missing.remotes <- pkg.remotes[!names(pkg.remotes) %in% installed.pkg]

print(">>> Checking for ultra-verbose packages...")
pkg.eigen <- rownames(P)[grep("RcppEigen",as.character(P[,"LinkingTo"]))]

if(any(missing.imports %in% c("RccpEigen",pkg.eigen))) {
  pkg.sel <- intersect(missing.imports,c("RccpEigen",pkg.eigen))
  message("> Pre-installing ultra-verbose packages: ", paste(pkg.sel,collapse=" "))
  if(!dir.exists("~/.R")) dir.create("~/.R")
  if(!file.exists("~/.R/Makevars")) file.copy("dev/Makevars","~/.R/Makevars.save")
  file.copy("dev/Makevars.no-error","~/.R/Makevars")
  BiocManager::install(pkg.sel,ask=FALSE)
  file.remove("~/.R/Makevars")
  if(!file.exists("~/.R/Makevars.save")) file.copy("dev/Makevars.save","~/.R/Makevars")  
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
## Install Kaleido for plotly
## ---------------------------------------------------------------------

if(0) {
  print(">>> Installing Kaleido/plotly packages...")
  ## Install a clean reticulate and miniconda
  # install.packages('reticulate', force=TRUE) # remove reticulate install since its already done.. and we get checksum error for some reason at this step
  unlink("~/.local/share/r-miniconda", recursive = TRUE)
  reticulate::install_miniconda()
  reticulate::conda_install("r-reticulate", "python-kaleido")
  reticulate::conda_install("r-reticulate", "plotly", channel = "plotly")
  reticulate::use_miniconda("r-reticulate")
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
