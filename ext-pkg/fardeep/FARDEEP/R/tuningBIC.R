#' Tuning parameter k in function alts using Bayesian Information Criterion (BIC) with some adjustment.
#'
#' @param x input matrix of predictors with n rows and p columns.
#' @param y input vector of dependent variable with length n.
#' @param alpha1  parameter used to adjust the upper bound of outliers. Take value from 0 to 1, default 0.1.
#' @param alpha2  parameter used to adjust the lower bound of outliers. Take value larger than 1, default 1.5.
#' @param up  upper bound of parameter k in function alts, default 10.
#' @param low lower bound of parameter k in function alts, default 1.
#' @param nn  whether coefficients are non-negative, default TRUE.
#' @param intercept whether intercept is included in model, default TRUE.
#' @param lognorm  whether noise is log-normal distributed, default TRUE.
#' @return k:   tuning result of parameter k for function alts.
#' @author YuningHao, Ming Yan, Yu Lei and YuyingXie
#' @references YuningHao, Ming Yan, Yu Lei and YuyingXie. Fast and Robust Deconvolution of Tumor Infiltrating Lymphocyte from Expression Profiles using Least Trimmed Squares.
#'
#' @examples
#' library(FARDEEP)
#' samp = sample.sim(n = 500, p = 20, sig = 1, a1 = 0.1, a2 = 0.2, nn = TRUE, intercept = TRUE)
#' k = tuningBIC(samp$x, samp$y, alpha1 = 0.1, alpha2 = 1.5, up = 10, low = 1, nn = TRUE, intercept = TRUE, lognorm = FALSE)
#' @export

tuningBIC = function(x, y, alpha1 = 0.1, alpha2 = 1.5, up = 10, low = 1, nn = TRUE,
                      intercept = TRUE, lognorm = TRUE){
  n = dim(x)[1]
  p = dim(x)[2]
  para  = NULL
  BIC.alts = NULL
  for (j in seq (low, up, 0.1)) {
    reg1  = alts(x = x, y = y, alpha1 = alpha1, alpha2 = alpha2, k = j, nn = nn, intercept = intercept)
    if (intercept){
      if (lognorm){
        res = log(abs(reg1$Y.new - cbind (1, reg1$X.new) %*% reg1$beta), 2)
      }else{
        res = reg1$Y.new - cbind (1, reg1$X.new) %*% reg1$beta
      }
    }else{
      if (lognorm){
        res = log(abs(reg1$Y.new - reg1$X.new %*% reg1$beta), 2)
      }else{
        res = reg1$Y.new - reg1$X.new %*% reg1$beta
      }
    }
    sse   = t (res) %*% res
    t     = reg1$number_outlier + p + 1
    no    = reg1$number_outlier
    BIC2  = (n - no) * log (sse / (n - no)) + t * (log (n - no) + 1)
    BIC.alts = rbind (BIC.alts, BIC2)
  }
  ind2    = which.min (BIC.alts)
  seq     = seq (low, up, 0.1)
  k       = seq[ind2]
  return (k)
}
