## 
## Write license file of the used/installed packages
## 

# This file is supposed to run from the root Playground folder
if (basename(getwd()) != "omicsplayground") {
    stop("Please run from the OmicsPlayground root folder")
}

## packages used in omicplayground 
renv1 <- renv::dependencies(path = "components", root = getwd(), errors = "ignored")
opg.pkg <- unique(renv1$Package)

## packages used in playbase
renv2 <- renv::dependencies(path = "../playbase/R", root = getwd(), errors = "ignored")
playbase.pkg <- unique(renv2$Package)

bigomics.pkg <- c("playbase","playdata","bigdash","bigLoaders","wizardR")
INSTALLED.PKGS <- sort(setdiff(c(opg.pkg,playbase.pkg), bigomics.pkg))

lisc <- installed.packages(fields = "License")
sel <- which(lisc[, "Package"] %in% INSTALLED.PKGS)
lisc1 <- lisc
lisc1 <- lisc[sel, ]
lisc1 <- lisc1[order(lisc1[, "Package"]), ]
lisc1 <- lisc1[!duplicated(lisc1[, "Package"]), ]
lisc2 <- lisc1[, c("Package", "Version", "License")]
write.table(lisc2, "RPackageLicenses.txt", sep='\t', row.names=FALSE)


