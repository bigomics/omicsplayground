# This file is supposed to run from the root Playground folder
if (basename(getwd()) != "omicsplayground") {
    stop("Please run from the OmicsPlayground root folder")
}

if(!require("renv")) install.packages("renv")
if(!require("BiocManager")) install.packages("BiocManager")
if(!require("remotes")) install.packages("remotes")
if(!require("devtools")) install.packages("devtools")
if(!require("reticulate")) install.packages("reticulate")

## ---------------------------------------------------------------------
## Automatically scan all used packages and install
## ---------------------------------------------------------------------
## We use renv to detect dependencies. Renv is looking for library and
## require statements in the r/R source files.
source("dev/write_description.R")

## some are perhaps already installed by playbase or other packages
installed.pkg <- sort(as.character(installed.packages()[,"Package"]))
missing.imports <- setdiff( pkg.imports, installed.pkg )
missing.remotes <- pkg.remotes[!names(pkg.remotes) %in% installed.pkg]

print(">>> Installing missing CRAN/BioConductor packages...")
message( ">>> Installing ",length(missing.imports),
        " missing CRAN/BioConductor packages:\n",
        paste(missing.imports,collapse="\n"))
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

print(">>> Installing Kaleido/plotly packages...")
## Install a clean reticulate and miniconda
# install.packages('reticulate', force=TRUE) # remove reticulate install since its already done.. and we get checksum error for some reason at this step
unlink("~/.local/share/r-miniconda", recursive = TRUE)
reticulate::install_miniconda()
reticulate::conda_install("r-reticulate", "python-kaleido")
reticulate::conda_install("r-reticulate", "plotly", channel = "plotly")
reticulate::use_miniconda("r-reticulate")

## ---------------------------------------------------------------------
## remove unneccessary big packages??
## ---------------------------------------------------------------------
print(">>> removing not needed packages...")

BIG.NOTUSED <- c(
    "reactome.db", ## >2GB!!!
    "terra",
    ## "DeMixT", ## purify
    "RNAseqData.HNRNPC.bam.chr14",
    ## "RSpectra",  ## ???
    "tximportData"
    ## "EnsDb.Hsapiens.v86",
    ## "EnsDb.Mmusculus.v79",
    ## "TxDb.Hsapiens.UCSC.hg19.knownGene",  ## need for import
    ## "TxDb.Mmusculus.UCSC.mm10.knownGene"  ## need for import
)
remove.pkgs <- function(pkgs) {
  for(i in 1:length(pkgs)) try(remove.packages(pkgs[i])
}
if(length(BIG.NOTUSED)>0) remove.pkgs(BIG.NOTUSED)


message("**** FINISHED INSTALLING REQUIREMENTS *****")
