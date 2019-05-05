#' Generate random sample with different proportion of outliers and leverage points
#'
#' @param n number of observations.
#' @param p number of independent variables (predictors).
#' @param sig variance of dependent variable.
#' @param a1  proportion of outliers.
#' @param a2  proportion of leverage points in outliers.
#' @param nn  whether coefficients are non-negative, default TRUE.
#' @param intercept whether intercept is included in model, default TRUE.
#' @return y:  vector of dependent variable.
#' @return x:  matrix of predictors with n rows and p columns.
#' @return loc:  index of added outliers.
#' @return beta: vector of coefficients.
#' @author YuningHao, Ming Yan, Yu Lei and YuyingXie
#' @references YuningHao, Ming Yan, Yu Lei and YuyingXie. Fast and Robust Deconvolution of Tumor Infiltrating Lymphocyte from Expression Profiles using Least Trimmed Squares.
#' @examples
#' library(FARDEEP)
#' samp = sample.sim(n = 500, p = 20, sig = 1, a1 = 0.1, a2 = 0.2, nn = TRUE, intercept = TRUE)
#' @export
sample.sim = function (n, p, sig, a1, a2, nn = TRUE, intercept = FALSE){
  u     = matrix (stats::runif (n * p, 0, 30), n, p)
  sigma = matrix (0.5, p, p)
  diag (sigma) = 1
  x     = u %*% (sigma ^ (1/2))
  e     = stats::rnorm (n, 0, sig)
  tau   = matrix(0, n, 1)
  loc.a1  = sample (n, floor(a1 * n))
  tau[loc.a1] = stats::rnorm (length(loc.a1), 20, 5)
  loc.a2  = sample (loc.a1, a2 * length(loc.a1))
  lev = 2 * max(x)
  x[loc.a2, ] = stats::rnorm(length(loc.a2) * p, lev, 1)
  if (nn){
    if (intercept){
      beta  = stats::runif(p + 1, 0, 1)
      y     = cbind(1, x) %*% beta + tau + e
    }else{
      beta  = stats::runif(p, 0, 1)
      y     = x %*% beta + tau + e
      }
  }else{
    if(intercept){
      beta  = stats::rnorm(p + 1, 0, 1)
      y     = cbind(1, x) %*% beta + tau + e
    }else{
      beta  = stats::rnorm(p, 0, 1)
      y     = x %*% beta + tau + e
    }
  }
  loc   = loc.a1
  result  = list (y = y, x = x, loc = loc, beta = beta)
  return(result)
}




