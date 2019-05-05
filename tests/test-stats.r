
library(genefilter)
A = matrix(rnorm(1000),100,10)
B = matrix(rnorm(800),100,8)
y <- factor(c(rep("A",ncol(A)), rep("B",ncol(B))))
out <- rowttests( cbind(A,B), y)
head(out)



A = mn[i,6:86]
B = ctrl[i,6:45]
y <- factor(c(rep("mn",ncol(A)), rep("ctrl",ncol(B))))
out <- rowttests( cbind(A,B), y)
head(out)
p_val <- out$p.value
log_e <- out$dm
plot( log_e, -log10(p_val) )

