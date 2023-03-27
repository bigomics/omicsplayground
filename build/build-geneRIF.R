##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

require(org.Hs.eg.db)
source("../R/pgx-functions.R")

cat("************************************************************************\n")
cat("******************* BUILDING GENERIF DATABASE **************************\n")
cat("************************************************************************\n")

##
## http://www.sthda.com/english/wiki/text-mining-and-word-cloud-
##   fundamentals-in-r-5-simple-steps-you-should-know
## 
##install.packages("tm")
library(tm)
library(data.table)

GENERIFS <- fread("../libx/generifs_basic.txt", quote=FALSE, sep="\t")
dim(GENERIFS)

rif.text <- GENERIFS$"GeneRIF text"
docs <- Corpus(VectorSource(rif.text))
##inspect(docs)

toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
docs <- tm_map(docs, toSpace, "/")
docs <- tm_map(docs, toSpace, "@")
docs <- tm_map(docs, toSpace, "\\|")

## Convert the text to lower case
docs <- tm_map(docs, content_transformer(tolower))
## Remove numbers
##docs <- tm_map(docs, removeNumbers)
## Remove english common stopwords
docs <- tm_map(docs, removeWords, stopwords("english"))

## Remove your own stop word specify your stopwords as a character vector
docs <- tm_map(docs, removeWords, c("expression","cell")) 
## Remove punctuations
docs <- tm_map(docs, removePunctuation, preserve_intra_word_dashes=TRUE)
## Eliminate extra white spaces
docs <- tm_map(docs, stripWhitespace)
## Text stemming
## docs <- tm_map(docs, stemDocument)

dtm <- TermDocumentMatrix(docs)
m <- sparseMatrix( dtm$j, dtm$i, x=dtm$v,
                  dims=c(dtm$ncol, dtm$nrow),
                  dimnames=rev(dtm$dimnames))
dim(m)
rownames(m) <- paste0(rif.text," (PMID:",GENERIFS$"PubMed ID (PMID) list",")")

sum(duplicated(rownames(m)))
is.dup <- (duplicated(rownames(m)))
table(is.dup)
m <- m[!is.dup,]
dim(m)

v <- sort(Matrix::colSums(m),decreasing=TRUE)
d <- data.frame(word=names(v),freq=v)
head(d, 200)
tail(d, 10)

stop.words <- c("study","cells","association","huge","navigator","role","observational",
                "may","gene","results","associated","data","patients","suggest",
                "levels","factor","function","via","novel","development","show",
                "important","indicate","response","genes","risk","findings","involved",
                "can","significantly","complex","high","mechanism","found","required",
                "plays","potential","regulated","domain","play","type","effect","identified",
                "induced","gene-enviroment","expressed","also","two","significant","might",
                "effects","essential","critical","higher","demonstrate","loss",
                "showed","analysis","early","functions","contribute","vivo","regulate",
                "model","new","different","compared","provide","population","first",
                "roles","specific","tissue","region","key","susceptibility",
                "factors","targeting","induces","studies","shows","axis","active","identify",
                "leads","part","sites","likely","acts","using","single","process",
                "independently","involves","responses","thus","tils","bcc","cows",
                "exacerbate","argue","operate","lower","leading","chinese","form","mediates",
                "feeback","anti","illuminate","meal","decreased","conditions","bring","critically",
                "limit","perpetuate","demonstrated","contributes","suggesting",
                "primary","induces","well","related","multiple","tissues","independent",
                "observed","studies","present","cause","due","caused","suggested",
                "poor","binds","three","one","two","","","","","","","","","","","","","","","",
                "","","","","","","","","","","","","","","","","","","","")

stop.words <- c(stop.words, names(which(v <= 3)))
m <- m[,!colnames(m) %in% stop.words]
dim(m)

cat("writing to geneRIF-matrix.rds")
saveRDS(m, file="../lib/geneRIF-matrix.rds")

cat("************************************************************************\n")
cat("************************ BUILD ANNOTATIONS DONE ************************\n")
cat("************************************************************************\n")
