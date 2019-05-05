require(glmnet)
require(sva)
source("../R/gx-heatmap.r")
source("../R/gx-limma.r")

##BiocManager::install("survcomp")
##install.packages("plsRcox")
library(survival)
library(plsRcox)

methods=c("glmnet","randomforest","boruta","xgboost")
methods=c("glmnet","randomforest","xgboost")

source("../R/pgx-predict.R")

##----------------------------------------------------------------------
## combine all importance values
##----------------------------------------------------------------------

load("../pgx/GSE10846-dlbcl-mRNA-8k.pgx",verbose=1)
y <- factor(ngs$samples$dlbcl.type)
time <- as.numeric(ngs$samples$OS.years)
status <- (ngs$samples$OS.status=="DEAD")
##y <- Surv( data=ngs$samples)

load("../pgx/sallusto2019-th1star-cf-12k-LT.pgx",verbose=1)
y <- factor(ngs$samples$cell.type)
y <- factor(ngs$samples$treatment)
table(y)

X <- as.matrix(ngs$X)
res <- ngs$gx.meta$meta[[1]]

##X <- as.matrix(ngs$gsetX)
##res <- ngs$gset.meta$meta[[1]]

##jj <- which(res[,"adj.P.Val"] < 0.05 & abs(res[,"logFC"]) > 1)
jj <- which(res[,"meta.q"] < 0.05 & abs(res[,"meta.fx"]) > 1)
length(jj)
X <- X[rownames(res)[jj],]
X <- head(X[order(-apply(X,1,sd)),],200)
dim(X)


##----------------------------------------------------------------------
## tests (miXomics,pls,...)
##----------------------------------------------------------------------
library(mixOmics)
y <- as.character(y)
y[sample(150,40)] <- "class3"
res <- splsda( t(X), y, keepX = c(25, 25))  
plotIndiv(res)                                      
plotVar(res)                                        
cim(res, comp = 1)
network(res, comp = 1)


library(plsRcox)
vars <- list(x=t(X),time=time,status=status)
(res <- cv.coxpls(vars,nt=10))
res <- plsRcox(t(X),time=time,event=status,nt=10)
summary(res)
cf <- res$Coeffs[,1] 

##----------------------------------------------------------------------
## compute importance values
##----------------------------------------------------------------------
table(y)

methods=c("glmnet","randomforest","boruta","xgboost")
methods=c("glmnet","randomforest","xgboost","splsda")
P <- pgx.multiclassVariableImportance(X, y, methods=methods)
P <- pgx.variableImportance(X, y, methods=methods)

P <- pgx.survivalVariableImportance(
    X, time, status, methods=c("glmnet","randomforest","xgboost","plsRcox"))

P <- t( t(P) / apply(P,2,max,na.rm=TRUE))
P <- pmax(P,0.1)
P <- P[order(-rowSums(P,na.rm=TRUE)),]
head(P)

par(mfrow=c(2,2), mar=c(8,4,2,2))
frame()
barplot( t(head(P,30)), las=1, horiz=TRUE)
klr <- grey.colors(ncol(P))
legend("topright", legend=colnames(P), fill=klr)

R <- (apply(P,2,rank)/nrow(P))**4
R <- R[order(-rowSums(R)),]
frame()
barplot( t(head(R,30)), las=1, horiz=TRUE)
legend("topright", legend=colnames(R), fill=klr)

par(mfrow=c(4,5), mar=c(4,3,2,1))
for(i in 1:20) {
    g <- rownames(R)[i]
    boxplot( X[g,] ~ y, col="grey90", main=g, las=3)
}
head(res[rownames(R),])

##install.packages("rpart.plot")
require(rpart)
require(rpart.plot)
sel <- rownames(X)
sel <- head(rownames(R),20)
##sel <- head(rownames(R),10)

jj <- 1:length(y)
## jj <- which(y %in% c("Th1","Th1star"))
jj <- c(jj,jj)
##jj <- c(jj,jj,jj)
rX <- X[sel,jj]
rX <- rX + 1e-1*matrix(rnorm(length(rX)),nrow(rX),ncol(rX))
df <- data.frame( y=as.character(y[jj]), t(rX[,] ))

rf <- rpart( y ~ ., data=df)
par(mfrow=c(1,1), mar=c(4,4,2,2)*0)
rpart.plot(rf)


library(party)
ct <- ctree( y ~ ., data=df)
plot(ct)
plot(ct, type="simple")

cf <- cforest( y ~ ., data=df)
pt <- prettytree(cf@ensemble[[1]], names(cf@data@get("input"))) 
nt <- new("BinaryTree")
nt@tree <- pt 
nt@data <- cf@data 
nt@responses <- cf@responses 
nt@where <- cf@where[[1]]
nt@weights <- cf@weights[[1]]
##nt@get_where <- cf@get_where

plot(nt)
plot(nt, type="simple")

