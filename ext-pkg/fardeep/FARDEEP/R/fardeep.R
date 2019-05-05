#' Using the idea of least trimmed square to detect and remove outliers before estimating the coefficients.
#' A robust method for gene-expression deconvolution.
#' 
#' @param x input matrix of predictors with n rows and p columns.
#' @param y input vector of dependent variable.
#' @param alpha1  parameter used to adjust the upper bound of outliers. Take value from 0 to 1, default 0.1.
#' @param alpha2  parameter used to adjust the lower bound of outliers. Take value larger than 1, default 1.5.
#' @param up  upper bound of parameter k in function alts, default 10.
#' @param low lower bound of parameter k in function alts, default 1.
#' @param nn  whether coefficients are non-negative,default TRUE.
#' @param intercept whether intercept is included in model, default TRUE.
#' @param lognorm  whether noise is log-normal distributed, default TRUE.
#' @return beta:  estimation of coefficients.
#' @return number_outlier:  number of outliers.
#' @return outlier_detect:  index of detected outliers.
#' @return X.new:  good observed points for independent variables.
#' @return Y.new:  good observed points for dependent variables.
#' @return k:  modified k (if the input value is not appropriate).
#' @author YuningHao, Ming Yan, Yu Lei and YuyingXie
#' @references YuningHao, Ming Yan, Yu Lei and YuyingXie. Fast and Robust Deconvolution of Tumor Infiltrating Lymphocyte from Expression Profiles using Least Trimmed Squares.
#' @examples
#' library(FARDEEP)
#' data(LM22)
#' data(mixture)
#' result = fardeep(LM22, mixture, alpha1 = 0.1, alpha2 = 1.5, up = 10, low = 1, nn = TRUE, intercept = TRUE, lognorm = TRUE, permn = 100, QN = TRUE)
#' coef = result$abs.beta
#' @export
fardeep = function(X, Y, alpha1 = 0.1, alpha2 = 1.5, up = 10, low = 1, nn = TRUE, intercept = TRUE, lognorm = TRUE, permn = 100, QN = TRUE){
    ##source("https://bioconductor.org/biocLite.R", echo = FALSE)
    ##biocLite("preprocessCore", suppressUpdates=TRUE, verbose = FALSE, quiet = TRUE)
    library(preprocessCore)
    X = data.matrix(X, rownames.force = TRUE)
    Y = data.matrix(Y, rownames.force = TRUE)
    X = X[order(rownames(X)),]
    Y = Y[order(rownames(Y)),]
    temp.c = colnames(Y)
    temp.r = rownames(Y)
    sig.gene = rownames(X)
    if(QN){
        Y = normalize.quantiles(Y)
        colnames(Y) = temp.c
        rownames(Y) = temp.r
    }
    YinX = temp.r %in% sig.gene
    Y = Y[YinX,]
    XinY = sig.gene %in% rownames(Y)
    X = X[XinY,]
    n = dim(X)[1]
    p = dim(X)[2]
    if(permn > 0) {
        itor = 1
        Ylist = as.list(data.matrix(Y))
        dist = NULL
        while(itor <= permn){
            y.samp = as.numeric(Ylist[sample(length(Ylist), dim(X)[1])])
            ktwo = tuningBIC(x = X, y = y.samp, alpha1 = alpha1, alpha2 = alpha2, up = up, low = low, nn = nn, intercept = intercept, lognorm = lognorm)
            mod = alts(x = X, y = y.samp, alpha1 = alpha1, alpha2 = alpha2, k = ktwo, nn = nn, intercept = intercept)
            if(intercept){
                yhat = cbind(1, X) %*% mod$beta
            }else{
                yhat = X %*% mod$beta
            }
            corr = cor(yhat, y.samp)
            dist = rbind(dist, corr)
            itor = itor + 1
        }
        nulldist = sort(dist)
    }
    nsamp = dim(Y)[2]
    para.k = c()
    pval = c()
    abs.beta = NULL
    for (i in 1:nsamp){
        ys = data.matrix(Y[, i])
        xs = data.matrix(X)
        k  = tuningBIC(x = xs, y = ys, alpha1 = alpha1, alpha2 = alpha2, up = up, low = low, nn = nn, intercept = intercept)
        reg   = alts(x = xs, y = ys, alpha1 = alpha1, alpha2 = alpha2, k = k, nn = nn, intercept = intercept)
        para.k = c(para.k, reg$k)
        if(intercept){
            abs.beta = rbind(abs.beta, reg$beta[-1])
            yh = cbind(1, xs) %*% reg$beta
        }else{
            abs.beta = rbind(abs.beta, reg$beta)
            yh = xs %*% reg$beta
        }
        if(permn > 0){
            pcor = cor(yh, ys)
            pv = 1 - (which.min(abs(nulldist - pcor)) / length(nulldist))
            pval = c(pval, pv)
        }
    }
    colnames(abs.beta) = colnames(X)
    rownames(abs.beta) = colnames(Y)
    if(permn > 0){
        pval = t(as.matrix(pval))
        colnames(pval) = colnames(Y)
    }
    para.k = t(as.matrix(para.k))
    colnames(para.k) = colnames(Y)
    result = list(abs.beta = abs.beta, pval = pval, k.value = para.k, Y.used = Y)
    return(result)
}
