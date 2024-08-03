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

## create fixed length sentences
fixstr <- function(s, n = 30) {
  substring(paste0(s, paste(rep(" ", n), collapse = "")), 1, n)
}
lisc2 <- rbind(c("PACKAGE", "VERSION", "LICENSE"), lisc2)
lisc2 <- cbind(fixstr(lisc2[, 1], 36), fixstr(lisc2[, 2], 15), fixstr(lisc2[, 3], 30))
lisc.text <- apply(lisc2, 1, function(s) paste0(s, collapse = ""))
head(lisc.text)
write(lisc.text, "RPackageLicenses.txt")
##  lisc2[grep("LGPL|AGPL|GPL-3",lisc2[,"License"]),]

